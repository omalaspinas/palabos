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

#include "palabos3D.h"
#include "palabos3D.hh"          // include full template code
#include "utility_dsl3d_param.h" // code for parameters and log file
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
//////// Lattice (Make sure it matches "lbm" in the .xml file!)
#define DESCRIPTOR D3Q27Descriptor
// #define DESCRIPTOR D3Q19Descriptor

//////// Select the collision model from data in the .xml file.
template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor> *getDynamics(Param<T> &param)
{
    Dynamics<T,Descriptor> *dyn;
    if (param.dynName == "BGK_Ma2") {   // BGK with second-order equilibrium
        dyn = new BGKdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "RM") { // Collision based on raw moment space (equivalent to Complete_BGK if hoOmega=SRT)
        dyn = new RMdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "HM") { // Collision based on Hermite moment space (equivalent to Complete_BGK if hoOmega=SRT)
        dyn = new HMdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "CM") { // Collision based on central moment space (equivalent to Complete_BGK if hoOmega=SRT)
        dyn = new CMdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "CHM") {// Collision based on central Hermite moment space (equivalent to Complete_BGK if hoOmega=SRT)
        dyn = new CHMdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "K") {  // Collision based on cumulant space
        dyn = new Kdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "GH") { // Collision based on Gauss-Hermite quadrature (HM with weighted scalar product, equivalent to Complete_BGK if hoOmega=SRT)
        dyn = new GHdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "RR") { // Recursive regularization of populations (equivalent to Complete_Regularized_BGK if hoOmega=SRT)
        dyn = new RRdynamics<T,Descriptor>(param.omega);
    } else {
        pcout << "Error: Dynamics name does not exist, please choose among BGK_Ma2, RM, HM, CM, CHM, K, GH and RR." << std::endl;
        exit(-1);
    }

    return dyn;
}


//////// A functional, used to initialize the simulation
template<typename T>
class DoubleShearLayerInitialVelocityField {
public:
    DoubleShearLayerInitialVelocityField(Param<T> const& param_)
        : param(param_)
    { }
    
    void operator()(plint iX, plint iY, plint iZ, T &rho, Array<T,3>& u) const {
        rho = (T)1;

        const plint nx = param.nx;
        const plint ny = param.ny;
        
        T ux = 0.;
        T uy = 0.;
        T uz = 0.;

        T kappa =80.;
        T delta = 0.05;
        T u0 = param.u0;

        T x = (T)iX/T(nx);
        T y = (T)iY/T(ny);
        if (y <= 0.5){
            ux   = u0*tanh(kappa*(y-0.25));
            uy   = u0*delta*sin(2.*M_PI*(x+0.25));
        }
        else{
            ux   = u0*tanh(kappa*(0.75-y));
            uy   = u0*delta*sin(2.*M_PI*(x+0.25));
        }

        u[0] = ux;
        u[1] = uy;
        u[2] = uz;
    }
    
private:
    Param<T> param;
};


//////// Initialize the simulation
void simulationSetup(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                     const Param<T> &param)
{
    // Set periodic boundaries.
    lattice.periodicity().toggleAll(true); 

    // Initialize the simulation domain.
    initializeAtEquilibrium (
        lattice, lattice.getBoundingBox(),
        DoubleShearLayerInitialVelocityField<T>(param) );

    // Call initialize to get the lattice ready for the simulation.
    lattice.initialize();
}


//////// Post processing
/// Produce a GIF snapshot of the vorticity-norm.
void writeGif(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
              Param<T> &param, plint iter)
{
    const plint imSize = 600;
    const plint nx = param.nx;
    const plint ny = param.ny;

    Box3D slice(0, nx-1, 0, ny-1, 0, 0);
    ImageWriter<T> imageWriter("leeloo");

    imageWriter.writeScaledGif( createFileName("vorticity", iter, 6),
                                *computeNorm(*computeVorticity (
                                        *computeVelocity(lattice) ), slice ),
                                imSize, imSize );
}
/// Output data into a VTK file.
void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
              Param<T> &param, plint iter)
{
    T dx = param.dx;
    T dt = param.dt;
    VtkStructuredImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeDensity(lattice), "density", 1.0);
    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", dx/dt);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "mach", 1.0/std::sqrt(DESCRIPTOR<T>::cs2));
    vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", (T)1/dt);
}


/////////////////////
/// MAIN PROGRAMM ///
/////////////////////
int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    if (argc != 5)
    {
        pcout << argc << std::endl;
        pcout << "Error! Wrong number of parameters." << std::endl;
        pcout << "Syntax: " << (std::string)global::argv(0) << " config.xml Re N Ma" << std::endl;
        pcout << "Example: " << (std::string)global::argv(0) << " config.xml 10000 128 0.1" << std::endl;
        pcout << "Example (mpi with 4 cores): mpirun -np 4 " << (std::string)global::argv(0) << " config.xml 10000 128 0.1" << std::endl;
        exit(1);
    }

    ///// Read and print the parameters of the simulation
    Param<T> param(argv[1]);
    param.writeLogFile();
    ///// Output directory
    global::directories().setOutputDir(param.outDirName);

    ///// Initialize the diagonal relaxation matrix from the .xml file.
    // Q19 and Q27 do not have the same number of relaxation parameters!
    Array<T, DESCRIPTOR<T>::numRelaxationTimes> allOmega;
    if (param.lbm == "D3Q19"){
        allOmega[0] = param.omega;     // relaxation of M200 and cyclic permutations
        allOmega[1] = param.omega;     // relaxation of M110 and cyclic permutations
        allOmega[2] = param.omega3;    // relaxation of M210 and cyclic permutations
        allOmega[3] = param.omega4;    // relaxation of M220 and cyclic permutations
        allOmega[4] = param.omegaBulk; // relaxation of bulk moment (M200 + M020 + M002)
    } else if (param.lbm == "D3Q27"){
        allOmega[0] = param.omega;     // relaxation of M200 and cyclic permutations
        allOmega[1] = param.omega;     // relaxation of M110 and cyclic permutations
        allOmega[2] = param.omega3;    // relaxation of M210 and cyclic permutations
        allOmega[3] = param.omega3;    // relaxation of M111 and cyclic permutations
        allOmega[4] = param.omega4;    // relaxation of M220 and cyclic permutations
        allOmega[5] = param.omega4;    // relaxation of M211 and cyclic permutations
        allOmega[6] = param.omega5;    // relaxation of M221 and cyclic permutations
        allOmega[7] = param.omega6;    // relaxation of M222 and cyclic permutations
        allOmega[8] = param.omegaBulk; // relaxation of bulk moment (M200 + M020 + M002)
    } else {
        pcout << "Error: lbm name does not exist." << std::endl;
        return EXIT_FAILURE;
    }
    
    ///// Generate the dynamics and the corresponding lattice from the .xml file.
    Dynamics<T,DESCRIPTOR> *dyn = getDynamics<T,DESCRIPTOR>(param);
    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
              param.nx, param.ny, param.nz, dyn);

    if (param.dynName == "RM") { 
        RMdynamics<T,DESCRIPTOR>::allOmega = allOmega;
    } else if (param.dynName == "HM") { 
        HMdynamics<T,DESCRIPTOR>::allOmega = allOmega;
    } else if (param.dynName == "CM") { 
        CMdynamics<T,DESCRIPTOR>::allOmega = allOmega;
    } else if (param.dynName == "CHM") { 
        CHMdynamics<T,DESCRIPTOR>::allOmega = allOmega;
    } else if (param.dynName == "K") { 
        Kdynamics<T,DESCRIPTOR>::allOmega = allOmega;
    } else if (param.dynName == "GH") { 
        GHdynamics<T,DESCRIPTOR>::allOmega = allOmega;
    } else if (param.dynName == "RR") { 
        RRdynamics<T,DESCRIPTOR>::allOmega = allOmega;
    }

    ///// Initialization from analytical profiles
    simulationSetup(lattice, param);

#ifndef PLB_REGRESSION
    ///// Initial state is outputed
    pcout << "Writing VTK file at iteration = 0" << std::endl;
    writeVTK(lattice, param, 0);
    writeGif(lattice, param, 0);
#endif

    ///// Simulation maximal time, and output frequency (in terms of convective times).
#ifndef PLB_REGRESSION
    plint vtkTout = param.vtkT * param.tc;
#endif
    plint tMax = param.tAdim * param.tc;
    
    ///// Output the evolution of the kinetic energy      
    plb_ofstream statsOut((param.outDirName+"/stats_" + param.fnameBase + ".dat").c_str());
    T initial_energy = std::numeric_limits<T>::max();

    ///// Main loop.
    for (plint iT=0; iT<=tMax; ++iT) {

        ///// Lattice Boltzmann iteration step.
        lattice.collideAndStream();
        
#ifndef PLB_REGRESSION
        ///// Output.
        if (iT % vtkTout == 0 and iT != 0) {
            pcout << "Writing VTK file at iteration = " << iT << std::endl;
            writeVTK(lattice, param, iT);
            writeGif(lattice, param, iT);
        }
#endif

        ///// Compute kinetic energy for stats and stability criterion
        T kinEnergy = computeAverageEnergy(lattice);
        T kinEnergyAdim = kinEnergy/(0.5*param.u0*param.u0);
#ifndef PLB_REGRESSION
        statsOut << iT/param.tc << " " << kinEnergyAdim << std::endl;    
#else
        pcout << iT/param.tc << " " << kinEnergyAdim << std::endl;    
#endif

        ///// Stability test based on the kinetic energy.
        // Here, the kinetic energy should be smaller than 2 times its initial value
        if (iT == 0) initial_energy = kinEnergyAdim;
        if (kinEnergyAdim > 2.*initial_energy) {
            pcout << "Catastrophic error: energy has increased or is NaN!" << std::endl;
#ifndef PLB_REGRESSION
            writeVTK(lattice, param, iT);
            writeGif(lattice, param, iT);
#endif
            return EXIT_FAILURE;
        }
    }
    /// Close stats file
    statsOut.close();

    return EXIT_SUCCESS;
}
