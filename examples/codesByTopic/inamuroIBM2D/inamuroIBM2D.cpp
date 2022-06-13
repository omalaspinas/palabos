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
 * Flow around a 2D cylinder inside a channel, with the creation of a von
 * Karman vortex street. This example makes use of bounce-back nodes to
 * describe the shape of the cylinder. The outlet is modeled through a
 * Neumann (zero velocity-gradient) condition.
 */

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "palabos2D.h"
#include "palabos2D.hh"

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9Descriptor

plint largeEnvelopeWidth = 4;  // Because of immersed walls.

/// Velocity on the parabolic Poiseuille profile
T poiseuilleVelocity(plint iY, IncomprFlowParam<T> const &parameters)
{
    T y = (T)iY / parameters.getResolution();
    return 4. * parameters.getLatticeU() * (y - y * y);
}

/// Linearly decreasing pressure profile
T poiseuillePressure(plint iX, IncomprFlowParam<T> const &parameters)
{
    T Lx = parameters.getNx() - 1;
    T Ly = parameters.getNy() - 1;
    return 8. * parameters.getLatticeNu() * parameters.getLatticeU() / (Ly * Ly)
           * (Lx / (T)2 - (T)iX);
}

/// Convert pressure to density according to ideal gas law
T poiseuilleDensity(plint iX, IncomprFlowParam<T> const &parameters)
{
    return poiseuillePressure(iX, parameters) * DESCRIPTOR<T>::invCs2 + (T)1;
}

/// A functional, used to initialize the velocity for the boundary conditions
template <typename T>
class PoiseuilleVelocity {
public:
    PoiseuilleVelocity(IncomprFlowParam<T> parameters_) : parameters(parameters_) { }
    void operator()(plint iX, plint iY, Array<T, 2> &u) const
    {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }

private:
    IncomprFlowParam<T> parameters;
};

/// A functional, used to initialize a pressure boundary to constant density
template <typename T>
class ConstantDensity {
public:
    ConstantDensity(T density_) : density(density_) { }
    T operator()(plint iX, plint iY) const
    {
        return density;
    }

private:
    T density;
};

/// A functional, used to create an initial condition for the density and velocity
template <typename T>
class PoiseuilleVelocityAndDensity {
public:
    PoiseuilleVelocityAndDensity(IncomprFlowParam<T> parameters_) : parameters(parameters_) { }
    void operator()(plint iX, plint iY, T &rho, Array<T, 2> &u) const
    {
        rho = poiseuilleDensity(iX, parameters);
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }

private:
    IncomprFlowParam<T> parameters;
};

class SurfaceVelocity {
public:
    SurfaceVelocity(T t_) : t(t_) { }
    Array<T, 2> operator()(Array<T, 2> const &pos)
    {
        return Array<T, 2>(T(), T());
    }

private:
    T t;
};

void boundarySetup(
    MultiBlockLattice2D<T, DESCRIPTOR> &lattice, std::vector<MultiBlock2D *> &rhoBarJarg,
    IncomprFlowParam<T> const &parameters,
    OnLatticeBoundaryCondition2D<T, DESCRIPTOR> &boundaryCondition)
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    Box2D outlet(nx - 1, nx - 1, 1, ny - 2);

    // Create Velocity boundary conditions everywhere
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, Box2D(0, 0, 1, ny - 2));
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, Box2D(0, nx - 1, 0, 0));
    boundaryCondition.setVelocityConditionOnBlockBoundaries(
        lattice, Box2D(0, nx - 1, ny - 1, ny - 1));
    // .. except on right boundary, where we prefer an outflow condition
    //    (zero velocity-gradient).
    boundaryCondition.setVelocityConditionOnBlockBoundaries(
        lattice, Box2D(nx - 1, nx - 1, 1, ny - 2), boundary::outflow);

    setBoundaryVelocity(lattice, lattice.getBoundingBox(), PoiseuilleVelocity<T>(parameters));
    setBoundaryDensity(lattice, outlet, ConstantDensity<T>(1.));
    initializeAtEquilibrium(
        lattice, lattice.getBoundingBox(), PoiseuilleVelocityAndDensity<T>(parameters));

    // Initialization
    applyProcessingFunctional(
        new BoxRhoBarJfunctional2D<T, DESCRIPTOR>(), lattice.getBoundingBox(), rhoBarJarg);
}

void cylinderSetup(
    IncomprFlowParam<T> const &parameters, vector<Array<T, 2>> &vertices, vector<T> &areas, T dV)
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    plint cx = nx / 4;
    plint cy = ny / 2 + 2;  // cy is slightly offset to avoid full symmetry,
                            //   and to get a Von Karman Vortex street.
    plint radius = cy / 4;
    // S/N*dx = dV \in [0.5, 1], so N = S/dx/dV
    // S the surface area of the cylinder, N the number of nodes
    plint nodeNum = util::roundToInt(2 * M_PI * radius / dV) + 1;
    for (int i = 0; i < nodeNum; ++i) {
        T angle = (T)i * 2. * M_PI / (T)nodeNum;
        vertices.push_back(
            Array<T, 2>(cx + radius * std::cos(angle), cy + radius * std::sin(angle)));
        areas.push_back(dV);
    }
}

void writeGif(MultiBlockLattice2D<T, DESCRIPTOR> &lattice, plint iter)
{
    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("u", iter, 6), *computeVelocityNorm(lattice));
}

void writeVTK(
    MultiBlockLattice2D<T, DESCRIPTOR> &lattice, IncomprFlowParam<T> const &parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx / dt);
    vtkOut.writeData<2, float>(*computeVelocity(lattice), "velocity", dx / dt);
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");

    IncomprFlowParam<T> parameters(
        (T)1e-2,  // uMax
        (T)600.,  // Re
        100,      // N
        6.,       // lx
        1.        // ly
    );
    const T logT = (T)0.02;
#ifdef PLB_REGRESSION
    const T maxT = (T)0.1;
#else
    const T imSave = (T)0.06;
    const T vtkSave = (T)1.;
    const T maxT = (T)20.1;
#endif

    writeLogFile(parameters, "Poiseuille flow");

    // The immersed boundary condition works with the external macroscopic variable
    MultiBlockLattice2D<T, DESCRIPTOR> *lattice = new MultiBlockLattice2D<T, DESCRIPTOR>(
        parameters.getNx(), parameters.getNy(),
        new BGKdynamics<T, DESCRIPTOR>(parameters.getOmega()));
    lattice->toggleInternalStatistics(false);

    MultiScalarField2D<T> *rhoBar =
        generateMultiScalarField<T>((MultiBlock2D &)*lattice, largeEnvelopeWidth).release();
    rhoBar->toggleInternalStatistics(false);
    MultiTensorField2D<T, 2> *j =
        generateMultiTensorField<T, 2>((MultiBlock2D &)*lattice, largeEnvelopeWidth).release();
    j->toggleInternalStatistics(false);

    std::vector<MultiBlock2D *> rhoBarJarg;
    rhoBarJarg.push_back(lattice);
    rhoBarJarg.push_back(rhoBar);
    rhoBarJarg.push_back(j);

    integrateProcessingFunctional(
        new ExternalRhoJcollideAndStream2D<T, DESCRIPTOR>(), lattice->getBoundingBox(), rhoBarJarg,
        0);
    integrateProcessingFunctional(
        new BoxRhoBarJfunctional2D<T, DESCRIPTOR>(), lattice->getBoundingBox(), rhoBarJarg, 1);

    OnLatticeBoundaryCondition2D<T, DESCRIPTOR> *boundaryCondition =
        createLocalBoundaryCondition2D<T, DESCRIPTOR>();

    vector<Array<T, 2>> vertices;
    vector<T> areas;
    T dV = 0.45;
    MultiContainerBlock2D container(*lattice);

    boundarySetup(*lattice, rhoBarJarg, parameters, *boundaryCondition);
    cylinderSetup(parameters, vertices, areas, dV);  // create boundary nodes

    // Instantiate the immersed wall data
    instantiateImmersedWallData(vertices, areas, container);

    // Main loop over time iterations.
    for (plint iT = 0; iT * parameters.getDeltaT() < maxT; ++iT) {
        // At this point, the state of the lattice corresponds to the
        //   discrete time iT. However, the stored averages (getStoredAverageEnergy
        //   and getStoredAverageDensity) correspond to the previous time iT-1.

#ifndef PLB_REGRESSION
        if (iT % parameters.nStep(imSave) == 0) {
            pcout << "Saving Gif ..." << endl;
            writeGif(*lattice, iT);
        }

        if (iT % parameters.nStep(vtkSave) == 0 && iT > 0) {
            pcout << "Saving VTK file ..." << endl;
            writeVTK(*lattice, parameters, iT);
        }
#endif

        if (iT % parameters.nStep(logT) == 0) {
            pcout << "step " << iT << "; t=" << iT * parameters.getDeltaT();
        }

        // Lattice Boltzmann iteration step.
        lattice
            ->executeInternalProcessors();  // Execute all processors and communicate appropriately.

        // Perform the immersed boundary iterations.
        plint ibIter = 5;
        for (int i = 0; i < ibIter; ++i) {
            inamuroIteration(
                SurfaceVelocity(iT + 1.0), *rhoBar, *j, container, (T)parameters.getTau(), false);
        }

        // At this point, the state of the lattice corresponds to the
        //   discrete time iT+1, and the stored averages are upgraded to time iT.
        if (iT % parameters.nStep(logT) == 0) {
            pcout << "; av energy =" << setprecision(10) << computeAverageEnergy<T>(*lattice)
                  << "; av rho =" << computeAverageDensity<T>(*lattice) << endl;
        }
    }

    delete boundaryCondition;
    delete j;
    delete rhoBar;
    delete lattice;
}
