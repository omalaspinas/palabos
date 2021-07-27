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

/* This file was written by Christophe Coreixas and Jonas Latt for the
 * Palabos Summer School 2021.
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "palabos2D.h"
#include "palabos2D.hh"           // include full template code
#include "utility_dsl2d_param.h"  // code for parameters and log file

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9Descriptor

//////// Select the collision model from data in the .xml file.
template <typename T, template <typename U> class Descriptor>
Dynamics<T, Descriptor> *getDynamics(Param<T> &param)
{
    Dynamics<T, Descriptor> *dyn;
    if (param.dynName == "BGK_Ma2") {  // BGK with second-order equilibrium
        dyn = new BGKdynamics<T, Descriptor>(param.omega);
    } else if (param.dynName == "RM")
    {  // Collision based on raw moment space (equivalent to Complete_BGK if hoOmega=SRT)
        dyn = new RMdynamics<T, Descriptor>(param.omega);
    } else if (param.dynName == "HM")
    {  // Collision based on Hermite moment space (equivalent to Complete_BGK if hoOmega=SRT)
        dyn = new HMdynamics<T, Descriptor>(param.omega);
    } else if (param.dynName == "CM")
    {  // Collision based on central moment space (equivalent to Complete_BGK if hoOmega=SRT)
        dyn = new CMdynamics<T, Descriptor>(param.omega);
    } else if (param.dynName == "CHM") {  // Collision based on central Hermite moment space
                                          // (equivalent to Complete_BGK if hoOmega=SRT)
        dyn = new CHMdynamics<T, Descriptor>(param.omega);
    } else if (param.dynName == "K") {  // Collision based on cumulant space
        dyn = new Kdynamics<T, Descriptor>(param.omega);
    } else if (param.dynName == "GH")
    {  // Collision based on Gauss-Hermite quadrature (HM with weighted scalar product, equivalent
       // to Complete_BGK if hoOmega=SRT)
        dyn = new GHdynamics<T, Descriptor>(param.omega);
    } else if (param.dynName == "RR") {  // Recursive regularization of populations (equivalent to
                                         // Complete_Regularized_BGK if hoOmega=SRT)
        dyn = new RRdynamics<T, Descriptor>(param.omega);
    } else {
        pcout << "Error: Dynamics name does not exist, please choose among BGK_Ma2, RM, HM, CM, "
                 "CHM, K, GH and RR."
              << std::endl;
        exit(-1);
    }

    return dyn;
}

//////// A functional, used to initialize the simulation
template <typename T>
class DoubleShearLayerInitialVelocityField {
public:
    DoubleShearLayerInitialVelocityField(Param<T> const &param_) : param(param_) { }

    void operator()(plint iX, plint iY, T &rho, Array<T, 2> &u) const
    {
        rho = (T)1;

        const plint nx = param.nx;
        const plint ny = param.ny;

        T ux = 0.;
        T uy = 0.;

        T kappa = 80.;
        T delta = 0.05;
        T u0    = param.u0;

        T x = (T)iX / T(nx);
        T y = (T)iY / T(ny);
        if (y <= 0.5) {
            ux = u0 * tanh(kappa * (y - 0.25));
            uy = u0 * delta * sin(2. * M_PI * (x + 0.25));
        } else {
            ux = u0 * tanh(kappa * (0.75 - y));
            uy = u0 * delta * sin(2. * M_PI * (x + 0.25));
        }

        u[0] = ux;
        u[1] = uy;
    }

private:
    Param<T> param;
};

//////// Initialize the simulation
void simulationSetup(MultiBlockLattice2D<T, DESCRIPTOR> &lattice, const Param<T> &param)
{
    // Set periodic boundaries.
    lattice.periodicity().toggleAll(true);

    // Initialize the simulation domain.
    initializeAtEquilibrium(
        lattice, lattice.getBoundingBox(), DoubleShearLayerInitialVelocityField<T>(param));

    // Call initialize to get the lattice ready for the simulation.
    lattice.initialize();
}

//////// Post processing
/// Produce a GIF snapshot of the vorticity.
void writeGif(MultiBlockLattice2D<T, DESCRIPTOR> &lattice, plint iter)
{
    const plint imSize = 600;

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(
        createFileName("vorticity", iter, 6), *computeVorticity(*computeVelocity(lattice)), imSize,
        imSize);
}
/// Output data into a VTK file.
void writeVTK(MultiBlockLattice2D<T, DESCRIPTOR> &lattice, Param<T> &param, plint iter)
{
    T dx = param.dx;
    T dt = param.dt;
    VtkStructuredImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeDensity(lattice), "density", 1.0);
    vtkOut.writeData<2, float>(*computeVelocity(lattice), "velocity", dx / dt);
    vtkOut.writeData<float>(
        *computeVelocityNorm(lattice), "mach", 1.0 / std::sqrt(DESCRIPTOR<T>::cs2));
    vtkOut.writeData<float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", (T)1 / dt);
}

/////////////////////
/// MAIN PROGRAMM ///
/////////////////////
int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);

    if (argc != 5) {
        pcout << argc << std::endl;
        pcout << "Error! Wrong number of parameters." << std::endl;
        pcout << "Syntax: " << (std::string)global::argv(0) << " config.xml Re N Ma" << std::endl;
        pcout << "Example: " << (std::string)global::argv(0) << " config.xml 10000 128 0.1"
              << std::endl;
        pcout << "Example (mpi with 4 cores): mpirun -np 4 " << (std::string)global::argv(0)
              << " config.xml 10000 128 0.1" << std::endl;
        exit(1);
    }

    ///// Read and print the parameters of the simulation
    Param<T> param(argv[1]);
    param.writeLogFile();
    ///// Output directory
    global::directories().setOutputDir(param.outDirName);

    ///// Initialize the relaxation matrix from the .xml file.
    Array<T, DESCRIPTOR<T>::numRelaxationTimes> allOmega;
    if (param.lbm == "D2Q9") {
        allOmega[0] = param.omega;      // relaxation of M20 and M02
        allOmega[1] = param.omega;      // relaxation of M11
        allOmega[2] = param.omega3;     // relaxation of M21 and M12
        allOmega[3] = param.omega4;     // relaxation of M22
        allOmega[4] = param.omegaBulk;  // relaxation of bulk moment (M20 + M02)
    } else {
        pcout << "Error: lattice does not exist." << std::endl;
        return EXIT_FAILURE;
    }

    ///// Generate the dynamics and the corresponding lattice from the .xml file.
    Dynamics<T, DESCRIPTOR> *dyn = getDynamics<T, DESCRIPTOR>(param);
    MultiBlockLattice2D<T, DESCRIPTOR> lattice(param.nx, param.ny, dyn);

    ///// Initialize relexation frequencies from those in "param"
    if (param.dynName == "RM") {
        RMdynamics<T, DESCRIPTOR>::allOmega = allOmega;
    } else if (param.dynName == "HM") {
        HMdynamics<T, DESCRIPTOR>::allOmega = allOmega;
    } else if (param.dynName == "CM") {
        CMdynamics<T, DESCRIPTOR>::allOmega = allOmega;
    } else if (param.dynName == "CHM") {
        CHMdynamics<T, DESCRIPTOR>::allOmega = allOmega;
    } else if (param.dynName == "K") {
        Kdynamics<T, DESCRIPTOR>::allOmega = allOmega;
    } else if (param.dynName == "GH") {
        GHdynamics<T, DESCRIPTOR>::allOmega = allOmega;
    } else if (param.dynName == "RR") {
        RRdynamics<T, DESCRIPTOR>::allOmega = allOmega;
    }

    ///// Initialization from analytical profiles
    simulationSetup(lattice, param);

#ifndef PLB_REGRESSION
    ///// Initial state is outputed
    pcout << "Writing VTK file at iteration = 0" << std::endl;
    writeVTK(lattice, param, 0);
    writeGif(lattice, 0);
#endif

    ///// Simulation maximal time, and output frequency (in terms of convective times).
#ifndef PLB_REGRESSION
    plint vtkTout = param.vtkT * param.tc;
#endif
    plint tMax    = param.tAdim * param.tc;

    ///// Output the evolution of the kinetic energy
    plb_ofstream statsOut((param.outDirName + "/stats_" + param.fnameBase + ".dat").c_str());
    T initial_energy = std::numeric_limits<T>::max();

    ///// Main loop.
    for (plint iT = 0; iT <= tMax; ++iT) {
        ///// Lattice Boltzmann iteration step.
        lattice.collideAndStream();

#ifndef PLB_REGRESSION
        ///// Output.
        if (iT % vtkTout == 0 and iT != 0) {
            pcout << "Writing VTK file at iteration = " << iT << std::endl;
            writeGif(lattice, iT);
            writeVTK(lattice, param, iT);
        }
#endif

        ///// Compute kinetic energy for stats and stability criterion
        T kinEnergy     = computeAverageEnergy(lattice);
        T kinEnergyAdim = kinEnergy / (0.5 * param.u0 * param.u0);
#ifndef PLB_REGRESSION
        statsOut << iT / param.tc << " " << kinEnergyAdim << std::endl;
#else
        pcout << iT / param.tc << " " << kinEnergyAdim << std::endl;
#endif

        ///// Stability test based on the kinetic energy.
        // Here, the kinetic energy should be smaller than 2 times its initial value
        if (iT == 0)
            initial_energy = kinEnergyAdim;
        if (kinEnergyAdim > 2. * initial_energy) {
            pcout << "Catastrophic error: energy has increased or is NaN!" << std::endl;
#ifndef PLB_REGRESSION
            writeGif(lattice, iT);
            writeVTK(lattice, param, iT);
#endif
            return EXIT_FAILURE;
        }
    }
    /// Close stats file
    statsOut.close();

    return EXIT_SUCCESS;
}
