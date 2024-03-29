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
 * 2D couplings -- generic implementation
 */
#ifndef COUPLING_2D_HH
#define COUPLING_2D_HH

#include "core/globalDefs.h"
#include "multiBlock/coupling2D.h"

namespace plb {

/* ********************  FullDomainCollideAndStreamAction2D ************************** */

template <typename T, template <typename U> class Descriptor>
FullDomainCollideAndStreamAction2D<T, Descriptor>::FullDomainCollideAndStreamAction2D(
    plint blockId_) :
    blockId(blockId_)
{ }

template <typename T, template <typename U> class Descriptor>
FullDomainCollideAndStreamAction2D<T, Descriptor>
    *FullDomainCollideAndStreamAction2D<T, Descriptor>::clone() const
{
    return new FullDomainCollideAndStreamAction2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void FullDomainCollideAndStreamAction2D<T, Descriptor>::execute(std::vector<id_t> &allMultiBlocks)
{
    PLB_ASSERT(blockId < (plint)allMultiBlocks.size());
    MultiBlock2D *block = multiBlockRegistration2D().find(allMultiBlocks[blockId]);
    PLB_ASSERT(block);
    MultiBlockLattice2D<T, Descriptor> *lattice =
        dynamic_cast<MultiBlockLattice2D<T, Descriptor> *>(block);
    PLB_ASSERT(lattice);
    lattice->externalCollideAndStream();
}

// QUESTION: Does northing. Normal?
template <typename T, template <typename U> class Descriptor>
void FullDomainCollideAndStreamAction2D<T, Descriptor>::regenerate(std::vector<id_t> &)
{ }

/* ********************  CollideAndStreamAction2D ************************** */

template <typename T, template <typename U> class Descriptor>
CollideAndStreamAction2D<T, Descriptor>::CollideAndStreamAction2D(plint blockId_, Box2D domain_) :
    blockId(blockId_), domain(domain_)
{ }

template <typename T, template <typename U> class Descriptor>
CollideAndStreamAction2D<T, Descriptor> *CollideAndStreamAction2D<T, Descriptor>::clone() const
{
    return new CollideAndStreamAction2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void CollideAndStreamAction2D<T, Descriptor>::execute(std::vector<id_t> &allMultiBlocks)
{
    PLB_ASSERT(blockId < (plint)allMultiBlocks.size());
    MultiBlock2D *block = multiBlockRegistration2D().find(allMultiBlocks[blockId]);
    PLB_ASSERT(block);
    MultiBlockLattice2D<T, Descriptor> *lattice =
        dynamic_cast<MultiBlockLattice2D<T, Descriptor> *>(block);
    PLB_ASSERT(lattice);
    lattice->collideAndStream(domain);
}

// QUESTION: Does nothing. Normal?
template <typename T, template <typename U> class Descriptor>
void CollideAndStreamAction2D<T, Descriptor>::regenerate(std::vector<id_t> &)
{ }

/* ********************  FullDomainStreamAction2D ************************** */

template <typename T, template <typename U> class Descriptor>
FullDomainStreamAction2D<T, Descriptor>::FullDomainStreamAction2D(plint blockId_) :
    blockId(blockId_)
{ }

template <typename T, template <typename U> class Descriptor>
FullDomainStreamAction2D<T, Descriptor> *FullDomainStreamAction2D<T, Descriptor>::clone() const
{
    return new FullDomainStreamAction2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void FullDomainStreamAction2D<T, Descriptor>::execute(std::vector<id_t> &allMultiBlocks)
{
    PLB_ASSERT(blockId < (plint)allMultiBlocks.size());
    MultiBlock2D *block = multiBlockRegistration2D().find(allMultiBlocks[blockId]);
    PLB_ASSERT(block);
    MultiBlockLattice2D<T, Descriptor> *lattice =
        dynamic_cast<MultiBlockLattice2D<T, Descriptor> *>(block);
    PLB_ASSERT(lattice);
    lattice->externalStream();
}

// QUESTION: Does nothing. Normal?
template <typename T, template <typename U> class Descriptor>
void FullDomainStreamAction2D<T, Descriptor>::regenerate(std::vector<id_t> &)
{ }

/* ********************  StreamAction2D ************************** */

template <typename T, template <typename U> class Descriptor>
StreamAction2D<T, Descriptor>::StreamAction2D(plint blockId_, Box2D domain_) :
    blockId(blockId_), domain(domain_)
{ }

template <typename T, template <typename U> class Descriptor>
StreamAction2D<T, Descriptor> *StreamAction2D<T, Descriptor>::clone() const
{
    return new StreamAction2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void StreamAction2D<T, Descriptor>::execute(std::vector<id_t> &allMultiBlocks)
{
    PLB_ASSERT(blockId < (plint)allMultiBlocks.size());
    MultiBlock2D *block = multiBlockRegistration2D().find(allMultiBlocks[blockId]);
    PLB_ASSERT(block);
    MultiBlockLattice2D<T, Descriptor> *lattice =
        dynamic_cast<MultiBlockLattice2D<T, Descriptor> *>(block);
    PLB_ASSERT(lattice);
    lattice->stream(domain);
}

// QUESTION: DOes nothing. Normal?
template <typename T, template <typename U> class Descriptor>
void StreamAction2D<T, Descriptor>::regenerate(std::vector<id_t> &)
{ }

/* ********************  IncrementTimeAction2D ************************** */

template <typename T, template <typename U> class Descriptor>
IncrementTimeAction2D<T, Descriptor>::IncrementTimeAction2D(plint blockId_) : blockId(blockId_)
{ }

template <typename T, template <typename U> class Descriptor>
IncrementTimeAction2D<T, Descriptor> *IncrementTimeAction2D<T, Descriptor>::clone() const
{
    return new IncrementTimeAction2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void IncrementTimeAction2D<T, Descriptor>::execute(std::vector<id_t> &allMultiBlocks)
{
    PLB_ASSERT(blockId < (plint)allMultiBlocks.size());
    MultiBlock2D *block = multiBlockRegistration2D().find(allMultiBlocks[blockId]);
    PLB_ASSERT(block);
    MultiBlockLattice2D<T, Descriptor> *lattice =
        dynamic_cast<MultiBlockLattice2D<T, Descriptor> *>(block);
    PLB_ASSERT(lattice);
    lattice->incrementTime();
}

// QUESTION: Empty. Is this normal?
template <typename T, template <typename U> class Descriptor>
void IncrementTimeAction2D<T, Descriptor>::regenerate(std::vector<id_t> &)
{ }

/* ********************  Actions2D ************************** */

template <typename T, template <typename U> class Descriptor>
plint Actions2D::addCollideAndStream(plint blockNum, Box2D domain)
{
    actions.push_back(new CollideAndStreamAction2D<T, Descriptor>(blockNum, domain));
    return (plint)actions.size() - 1;
}

template <typename T, template <typename U> class Descriptor>
plint Actions2D::addCollideAndStream(plint blockNum)
{
    actions.push_back(new FullDomainCollideAndStreamAction2D<T, Descriptor>(blockNum));
    return (plint)actions.size() - 1;
}

template <typename T, template <typename U> class Descriptor>
plint Actions2D::addStream(plint blockNum, Box2D domain)
{
    actions.push_back(new StreamAction2D<T, Descriptor>(blockNum, domain));
    return (plint)actions.size() - 1;
}

template <typename T, template <typename U> class Descriptor>
plint Actions2D::addStream(plint blockNum)
{
    actions.push_back(new FullDomainStreamAction2D<T, Descriptor>(blockNum));
    return (plint)actions.size() - 1;
}

template <typename T, template <typename U> class Descriptor>
plint Actions2D::addIncrementTime(plint blockNum)
{
    actions.push_back(new IncrementTimeAction2D<T, Descriptor>(blockNum));
    return (plint)actions.size() - 1;
}

}  // namespace plb

#endif  // COUPLING_2D_HH
