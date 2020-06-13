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
  * Flow in a lid-driven 3D cavity. The cavity is square and has no-slip walls,
  * except for the top wall which is diagonally driven with a constant
  * velocity. The benchmark is challenging because of the velocity
  * discontinuities on corner nodes.
  **/

#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor

template<typename T>
class CylinderShapeDomain3D : public plb::DomainFunctional3D {
public:
    CylinderShapeDomain3D(plb::plint cx_, plb::plint cz_, plb::plint radius)
        : cx(cx_),
          cz(cz_),
          radiusSqr(plb::util::sqr(radius))
    { }
    virtual bool operator() (plb::plint iX, plb::plint iY, plb::plint iZ) const {
        return plb::util::sqr(iX-cx) + plb::util::sqr(iZ-cz) <= radiusSqr;
    }
    virtual CylinderShapeDomain3D<T>* clone() const {
        return new CylinderShapeDomain3D<T>(*this);
    }
private:
    plb::plint cx;
    plb::plint cz;
    plb::plint radiusSqr;
};

void cylinderSetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                  IncomprFlowParam<T> const& parameters,
                  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition,
                  Array<plint,3> &forceIds )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();

    Box3D inlet = Box3D(0, 0, 0, ny-1, 1, nz-2);
    boundaryCondition.addVelocityBoundary0N(inlet, lattice, boundary::dirichlet);
    Box3D outlet = Box3D(nx-1, nx-1, 0, ny-1, 1, nz-2);
    boundaryCondition.addVelocityBoundary0P(outlet, lattice, boundary::neumann);
    Box3D top = Box3D(1, nx-2, 0, ny-1, nz-1, nz-1);
    boundaryCondition.addVelocityBoundary2P(top, lattice, boundary::freeslip);
    Box3D bottom = Box3D(1, nx-2, 0, ny-1, 0, 0);
    boundaryCondition.addVelocityBoundary2N(bottom, lattice, boundary::freeslip);

    Box3D edge_bot_left = Box3D(0, 0, 0, ny-1, 0, 0);
    boundaryCondition.addExternalVelocityEdge1NN (
            edge_bot_left, lattice, boundary::dirichlet );
    Box3D edge_top_left = Box3D(0, 0, 0, ny-1, nz-1, nz-1);
    boundaryCondition.addExternalVelocityEdge1PN (
            edge_top_left, lattice, boundary::dirichlet );
    Box3D edge_bot_right = Box3D(nx-1, nx-1, 0, ny-1, 0, 0);
    boundaryCondition.addExternalVelocityEdge1NP (
            edge_bot_right, lattice, boundary::neumann );
    Box3D edge_top_right = Box3D(nx-1, nx-1, 0, ny-1, nz-1, nz-1);
    boundaryCondition.addExternalVelocityEdge1PP (
            edge_top_right, lattice, boundary::neumann );

    setBoundaryVelocity(lattice, lattice.getBoundingBox(), Array<T,3>(parameters.getLatticeU(), 0.0, 0.0) );
    // Box3D topLid = Box3D(0, nx-1, ny-1, ny-1, 0, nz-1);
    // Box3D everythingButTopLid = Box3D(0, nx-1, 0, ny-2, 0, nz-1);

    // boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, lattice.getBoundingBox(), boundary::dirichlet);

    lattice.periodicity().toggle(1, true);
    
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), 
        (T) 1., Array<T,3>(parameters.getLatticeU(),(T)0.,(T)0.) );

    plint cx     = nx/2;
    plint cz     = nz/2+2; // cz is slightly offset to avoid full symmetry,
                          //   and to get a Von Karman Vortex street.
    plint radius = parameters.getResolution() / 2;

    
    lattice.toggleInternalStatistics(true);
    forceIds[0] = lattice.internalStatSubscription().subscribeSum();
    forceIds[1] = lattice.internalStatSubscription().subscribeSum();
    forceIds[2] = lattice.internalStatSubscription().subscribeSum();

    defineDynamics(lattice, lattice.getBoundingBox(),
                   new CylinderShapeDomain3D<T>(cx, cz, radius),
                   new plb::MomentumExchangeBounceBack<T,DESCRIPTOR>(forceIds));
    initializeMomentumExchange (lattice, lattice.getBoundingBox() );

    lattice.initialize();
}

template<class BlockLatticeT>
void writeGifs(BlockLatticeT& lattice,
               IncomprFlowParam<T> const& parameters, plint iter)
{
    const plint imSize = 600;
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    const plint xComponent = 0;

    Box3D slice(0, nx-1, ny/2, ny/2, 0, nz-1);
    ImageWriter<T> imageWriter("leeloo");

    imageWriter.writeScaledGif( createFileName("ux", iter, 6),
                                *computeVelocityComponent (lattice, slice, xComponent),
                                imSize, imSize );

    imageWriter.writeScaledGif( createFileName("uNorm", iter, 6),
                                *computeVelocityNorm (lattice, slice),
                                imSize, imSize );
    imageWriter.writeScaledGif( createFileName("omega", iter, 6),
                                *computeNorm(*computeVorticity (
                                        *computeVelocity(lattice) ), slice ),
                                imSize, imSize );
}

template<class BlockLatticeT>
void writeVTK(BlockLatticeT& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();

    //VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    ParallelVtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), 3, dx);

    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", dx/dt);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1./dt);
}


int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    IncomprFlowParam<T> parameters(
            (T) 1e-2,  // uMax
            (T) 10.,   // Re
            10,        // N
            50.,        // lx
            0.5,        // ly
            50.         // lz
    );
    const T logT     = (T)1/(T)10;
    const T imSave   = (T)1/(T)10;
    const T vtkSave  = (T)1/(T)100;
    const T maxT     = (T)10.1;

    pcout << "omega= " << parameters.getOmega() << std::endl;
    writeLogFile(parameters, "3D cylinder");

    T omega = parameters.getOmega();

    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
            parameters.getNx(), parameters.getNy(), parameters.getNz(),
            new BGKdynamics<T,DESCRIPTOR>(omega) );


    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createInterpBoundaryCondition3D<T,DESCRIPTOR>();

    Array<plint,3> forceIds;
    cylinderSetup(lattice, parameters, *boundaryCondition, forceIds);

    T previousIterationTime = T();
    // Loop over main time iteration.
    for (plint iT=0; iT<parameters.nStep(maxT); ++iT) {
        global::timer("mainLoop").restart();

        if (iT%parameters.nStep(imSave)==0) {
            pcout << "Writing Gif ..." << endl;
            writeGifs(lattice, parameters, iT);
        }

        if (iT%parameters.nStep(vtkSave)==0 && iT>0) {
            pcout << "Saving VTK file ..." << endl;
            writeVTK(lattice, parameters, iT);
        }

        if (iT%parameters.nStep(logT)==0) {
            pcout << "step " << iT
                  << "; t=" << iT*parameters.getDeltaT();
        }

        // Execute a time iteration.
        lattice.collideAndStream();

        // Access averages from internal statistics ( their value is defined
        //   only after the call to lattice.collideAndStream() )
        if (iT%parameters.nStep(logT)==0) {
            pcout << "; av energy="
                  << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                  << "; av rho="
                  << setprecision(10) << getStoredAverageDensity<T>(lattice) << endl;
            pcout << "Time spent during previous iteration: "
                  << previousIterationTime << endl;

            T drag = lattice.getInternalStatistics().getSum(forceIds[0]);
            T lift = lattice.getInternalStatistics().getSum(forceIds[2]);

            T length = parameters.getNy(); // Because of periodicity!
            T diameter = parameters.getResolution() + 1;

            T vel = parameters.getLatticeU();

            T avgRho = computeAverage(*computeDensity(lattice));

            T drag_coef_factor = 1.0 / (length * util::sqr(vel) * diameter * avgRho);
            pcout << "Drag coefficient: " << drag * drag_coef_factor << std::endl;
            pcout << "Lift coefficient: " << lift * drag_coef_factor << std::endl;
            pcout << "Average rho: " << avgRho << std::endl;
            pcout << std::endl;
        }

        previousIterationTime = global::timer("mainLoop").stop();
    }

    delete boundaryCondition;
}

