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
 * Flow in a lid-driven 2D cavity. The cavity is square and has no-slip walls,
 * except for the top wall which is driven to the right with a constant
 * velocity. The benchmark is challenging because of the velocity
 * discontinuities on corner nodes. The code on the other hand is very simple.
 * It could for example be used as a first example, to get familiar with Palabos.
 **/

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "palabos2D.h"
#include "palabos2D.hh"  // include full template code

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9Descriptor

void cavitySetup(
    MultiBlockLattice2D<T, DESCRIPTOR> &lattice, IncomprFlowParam<T> const &parameters,
    OnLatticeBoundaryCondition2D<T, DESCRIPTOR> &boundaryCondition)
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();

    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

    setBoundaryVelocity(lattice, lattice.getBoundingBox(), Array<T, 2>((T)0., (T)0.));
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), (T)1., Array<T, 2>((T)0., (T)0.));

    T u = parameters.getLatticeU();
    setBoundaryVelocity(lattice, Box2D(1, nx - 2, ny - 1, ny - 1), Array<T, 2>(u, (T)0.));
    initializeAtEquilibrium(
        lattice, Box2D(1, nx - 2, ny - 1, ny - 1), (T)1., Array<T, 2>(u, (T)0.));

    lattice.initialize();
}

template <class BlockLatticeT>
void writeGif(BlockLatticeT &lattice, plint iter)
{
    const plint imSize = 600;

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(
        createFileName("uNorm", iter, 6), *computeVelocityNorm(lattice), imSize, imSize);
    imageWriter.writeScaledGif(
        createFileName("logUnorm", iter, 6),
        *computeLog(*add((T)1.e-8, *computeVelocityNorm(lattice))), imSize, imSize);
}

template <class BlockLatticeT>
void writeVTK(BlockLatticeT &lattice, IncomprFlowParam<T> const &parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx / dt);
    vtkOut.writeData<2, float>(*computeVelocity(lattice), "velocity", dx / dt);
}

enum class CollisionType { BGK, TRT, BGKma1, TRTma1 };

/**
 * This function returns a pointer to the IsoThermalBulkDynamics base class. Depending on the
 * collision_type parameter it assign to the returned pointer a pointer to a different type of
 * dynamics, either BGK or TRT. The TRT class has a variant defined using a linearized equilibrium,
 * its name is Ma1TRTdynamics. This functions set the omega- parameter of the TRT class starting
 * from the magic parameter Lambda. The omega- is automaticcally recomputed when setting the magic
 * parameter.
 * @param basicString
 * @param tau_plus
 * @param magic
 * @param collision_type
 * @return
 */
IsoThermalBulkDynamics<T, DESCRIPTOR> *selectDynamics(
    const CollisionType collision_type, double tau_plus, double magic)
{
    IsoThermalBulkDynamics<T, DESCRIPTOR> *dynamics = nullptr;
    switch (collision_type) {
    case CollisionType::BGK:
        dynamics = new BGKdynamics<T, DESCRIPTOR>(1. / tau_plus);
        break;
    case CollisionType::BGKma1:  // in this we use a special Ma1TRT with tau-=tau+
        dynamics = new Ma1TRTdynamics<T, DESCRIPTOR>(1. / tau_plus, 1. / tau_plus, false);
        break;
    case CollisionType::TRT:
        dynamics = new TRTdynamics<T, DESCRIPTOR>(1. / tau_plus);
        dynamics->setParameter(dynamicParams::magicParameter, magic);
        break;
    case CollisionType::TRTma1:
        dynamics = new Ma1TRTdynamics<T, DESCRIPTOR>(1. / tau_plus);
        dynamics->setParameter(dynamicParams::magicParameter, magic);
        break;
    default:
        pcout << "ERROR! This collision model has not set up.\n";
        abort();
    }
    return dynamics;
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");

    IncomprFlowParam<T> parameters(
        (T)4e-3,  // uMax
        (T)1,     // Re
        64,       // N
        1.,       // lx
        1.        // ly
    );
    const T logT = (T)0.1;
#ifndef PLB_REGRESSION
    const T imSave = (T)0.01;
    const T vtkSave = (T)1.;
    const T maxT = (T)1.1;
#else
    const T maxT = (T)0.51;
#endif

    writeLogFile(parameters, "2D cavity");
    T tau_plus = parameters.getTau();
    if (tau_plus > 2.0) {
        pcout << "ERROR! The computed value of relaxation parameter tau_plus+= " << tau_plus
              << " is too high.\n";
        abort();
    }

    // Try different types of collision_models!
    CollisionType collision_model = CollisionType::TRTma1;
    // Check out how to use TRTdynamics inside the selectDynamics function.
    auto dynamics = selectDynamics(collision_model, tau_plus, 1. / 8.);

    MultiBlockLattice2D<T, DESCRIPTOR> lattice(parameters.getNx(), parameters.getNy(), dynamics);

    OnLatticeBoundaryCondition2D<T, DESCRIPTOR> *boundaryCondition =
        createInterpBoundaryCondition2D<T, DESCRIPTOR>();

    cavitySetup(lattice, parameters, *boundaryCondition);

#ifndef PLB_REGRESSION
    T previousIterationTime = T();
#endif

    // Main loop over time iterations.
    for (plint iT = 0; iT * parameters.getDeltaT() < maxT; ++iT) {
#ifndef PLB_REGRESSION
        global::timer("mainLoop").restart();

        if (iT % parameters.nStep(imSave) == 0 && iT > 0) {
            pcout << "Saving Gif ..." << endl;
            writeGif(lattice, iT);
            pcout << endl;
        }

        if (iT % parameters.nStep(vtkSave) == 0 && iT > 0) {
            pcout << "Saving VTK file ..." << endl;
            writeVTK(lattice, parameters, iT);
        }
#endif

        if (iT % parameters.nStep(logT) == 0) {
            pcout << "step " << iT << "; t=" << iT * parameters.getDeltaT();
        }

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();

        if (iT % parameters.nStep(logT) == 0) {
            pcout << "; av energy=" << setprecision(10) << getStoredAverageEnergy(lattice)
                  << "; av rho=" << getStoredAverageDensity(lattice) << endl;
#ifndef PLB_REGRESSION
            pcout << "Time spent during previous iteration: " << previousIterationTime << endl;
#endif
        }

#ifndef PLB_REGRESSION
        previousIterationTime = global::timer("mainLoop").stop();
#endif
    }

    delete boundaryCondition;
}
