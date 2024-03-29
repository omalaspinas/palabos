/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 *
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * The most recent release of Palabos can be downloaded at
 * <https://palabos.unige.ch/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \file
 * Helper classes for parallel 2D multiblock lattice -- generic implementation.
 */

#include "parallelism/parallelBlockCommunicator2D.h"

#include <algorithm>

#include "atomicBlock/atomicBlock2D.h"
#include "core/plbDebug.h"
#include "core/plbProfiler.h"
#include "multiBlock/multiBlock2D.h"
#include "multiBlock/multiBlockManagement2D.h"

namespace plb {

#ifdef PLB_MPI_PARALLEL

CommunicationStructure2D::CommunicationStructure2D(
    std::vector<Overlap2D> const &overlaps, MultiBlockManagement2D const &originManagement,
    MultiBlockManagement2D const &destinationManagement, plint sizeOfCell)
{
    plint fromEnvelopeWidth = originManagement.getEnvelopeWidth();
    plint toEnvelopeWidth = destinationManagement.getEnvelopeWidth();
    SparseBlockStructure2D const &fromSparseBlock = originManagement.getSparseBlockStructure();
    SparseBlockStructure2D const &toSparseBlock = destinationManagement.getSparseBlockStructure();

    SendRecvPool sendPool, recvPool;
    for (pluint iOverlap = 0; iOverlap < overlaps.size(); ++iOverlap) {
        Overlap2D const &overlap = overlaps[iOverlap];
        CommunicationInfo2D info;

        info.fromBlockId = overlap.getOriginalId();
        info.toBlockId = overlap.getOverlapId();

        SmartBulk2D originalBulk(fromSparseBlock, fromEnvelopeWidth, info.fromBlockId);
        SmartBulk2D overlapBulk(toSparseBlock, toEnvelopeWidth, info.toBlockId);

        Box2D originalCoordinates(overlap.getOriginalCoordinates());
        Box2D overlapCoordinates(overlap.getOverlapCoordinates());
        info.fromDomain = originalBulk.toLocal(originalCoordinates);
        info.toDomain = overlapBulk.toLocal(overlapCoordinates);
        info.absoluteOffset = Dot2D(
            overlapCoordinates.x0 - originalCoordinates.x0,
            overlapCoordinates.y0 - originalCoordinates.y0);

        plint lx = info.fromDomain.x1 - info.fromDomain.x0 + 1;
        plint ly = info.fromDomain.y1 - info.fromDomain.y0 + 1;
        PLB_PRECONDITION(lx == info.toDomain.x1 - info.toDomain.x0 + 1);
        PLB_PRECONDITION(ly == info.toDomain.y1 - info.toDomain.y0 + 1);

        plint numberOfCells = lx * ly;

        ThreadAttribution const &fromAttribution = originManagement.getThreadAttribution();
        ThreadAttribution const &toAttribution = destinationManagement.getThreadAttribution();
        info.fromProcessId = fromAttribution.getMpiProcess(info.fromBlockId);
        info.toProcessId = toAttribution.getMpiProcess(info.toBlockId);

        if (fromAttribution.isLocal(info.fromBlockId) && toAttribution.isLocal(info.toBlockId)) {
            sendRecvPackage.push_back(info);
        } else if (fromAttribution.isLocal(info.fromBlockId)) {
            sendPackage.push_back(info);
            sendPool.subscribeMessage(info.toProcessId, numberOfCells * sizeOfCell);
        } else if (toAttribution.isLocal(info.toBlockId)) {
            recvPackage.push_back(info);
            recvPool.subscribeMessage(info.fromProcessId, numberOfCells * sizeOfCell);
        }
    }

    sendComm = SendPoolCommunicator(sendPool);
    recvComm = RecvPoolCommunicator(recvPool);
}

////////////////////// Class ParallelBlockCommunicator2D /////////////////////

ParallelBlockCommunicator2D::ParallelBlockCommunicator2D() :
    overlapsModified(true), communication(0)
{ }

ParallelBlockCommunicator2D::ParallelBlockCommunicator2D(ParallelBlockCommunicator2D const &) :
    overlapsModified(true), communication(0)
{ }

ParallelBlockCommunicator2D::~ParallelBlockCommunicator2D()
{
    delete communication;
}

ParallelBlockCommunicator2D &ParallelBlockCommunicator2D::operator=(
    ParallelBlockCommunicator2D const &rhs)
{
    ParallelBlockCommunicator2D(rhs).swap(*this);
    return *this;
}

void ParallelBlockCommunicator2D::swap(ParallelBlockCommunicator2D &rhs)
{
    std::swap(overlapsModified, rhs.overlapsModified);
    std::swap(communication, rhs.communication);
}

ParallelBlockCommunicator2D *ParallelBlockCommunicator2D::clone() const
{
    return new ParallelBlockCommunicator2D(*this);
}

void ParallelBlockCommunicator2D::duplicateOverlaps(
    MultiBlock2D &multiBlock, modif::ModifT whichData) const
{
    MultiBlockManagement2D const &multiBlockManagement = multiBlock.getMultiBlockManagement();
    PeriodicitySwitch2D const &periodicity = multiBlock.periodicity();

    // Implement a caching mechanism for the communication structure.
    if (overlapsModified) {
        overlapsModified = false;
        LocalMultiBlockInfo2D const &localInfo = multiBlockManagement.getLocalInfo();
        std::vector<Overlap2D> overlaps(multiBlockManagement.getLocalInfo().getNormalOverlaps());
        for (pluint iOverlap = 0; iOverlap < localInfo.getPeriodicOverlaps().size(); ++iOverlap) {
            PeriodicOverlap2D const &pOverlap = localInfo.getPeriodicOverlaps()[iOverlap];
            if (periodicity.get(pOverlap.normalX, pOverlap.normalY)) {
                overlaps.push_back(pOverlap.overlap);
            }
        }
        delete communication;
        communication = new CommunicationStructure2D(
            overlaps, multiBlockManagement, multiBlockManagement, multiBlock.sizeOfCell());
    }

    communicate(*communication, multiBlock, multiBlock, whichData);
}

void ParallelBlockCommunicator2D::communicate(
    std::vector<Overlap2D> const &overlaps, MultiBlock2D const &originMultiBlock,
    MultiBlock2D &destinationMultiBlock, modif::ModifT whichData) const
{
    PLB_PRECONDITION(originMultiBlock.sizeOfCell() == destinationMultiBlock.sizeOfCell());

    CommunicationStructure2D communication(
        overlaps, originMultiBlock.getMultiBlockManagement(),
        destinationMultiBlock.getMultiBlockManagement(), originMultiBlock.sizeOfCell());
    global::profiler().start("mpiCommunication");
    communicate(communication, originMultiBlock, destinationMultiBlock, whichData);
    global::profiler().stop("mpiCommunication");
}

void ParallelBlockCommunicator2D::communicate(
    CommunicationStructure2D &communication, MultiBlock2D const &originMultiBlock,
    MultiBlock2D &destinationMultiBlock, modif::ModifT whichData) const
{
    bool staticMessage = whichData == modif::staticVariables;
    // 1. Non-blocking receives.
    communication.recvComm.startBeingReceptive(staticMessage);

    // 2. Non-blocking sends.
    for (unsigned iSend = 0; iSend < communication.sendPackage.size(); ++iSend) {
        CommunicationInfo2D const &info = communication.sendPackage[iSend];
        AtomicBlock2D const &fromBlock = originMultiBlock.getComponent(info.fromBlockId);
        fromBlock.getDataTransfer().send(
            info.fromDomain, communication.sendComm.getSendBuffer(info.toProcessId), whichData);
        communication.sendComm.acceptMessage(info.toProcessId, staticMessage);
    }

    // 3. Local copies which require no communication.
    for (unsigned iSendRecv = 0; iSendRecv < communication.sendRecvPackage.size(); ++iSendRecv) {
        CommunicationInfo2D const &info = communication.sendRecvPackage[iSendRecv];
        AtomicBlock2D const &fromBlock = originMultiBlock.getComponent(info.fromBlockId);
        AtomicBlock2D &toBlock = destinationMultiBlock.getComponent(info.toBlockId);
        plint deltaX = info.fromDomain.x0 - info.toDomain.x0;
        plint deltaY = info.fromDomain.y0 - info.toDomain.y0;
        toBlock.getDataTransfer().attribute(
            info.toDomain, deltaX, deltaY, fromBlock, whichData, info.absoluteOffset);
    }

    // 4. Finalize the receives.
    for (unsigned iRecv = 0; iRecv < communication.recvPackage.size(); ++iRecv) {
        CommunicationInfo2D const &info = communication.recvPackage[iRecv];
        AtomicBlock2D &toBlock = destinationMultiBlock.getComponent(info.toBlockId);
        toBlock.getDataTransfer().receive(
            info.toDomain, communication.recvComm.receiveMessage(info.fromProcessId, staticMessage),
            whichData, info.absoluteOffset);
    }

    // 5. Finalize the sends.
    communication.sendComm.finalize(staticMessage);
}

void ParallelBlockCommunicator2D::signalPeriodicity() const
{
    overlapsModified = true;
}

#endif  // PLB_MPI_PARALLEL

}  // namespace plb
