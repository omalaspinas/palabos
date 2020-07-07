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
 * Convection of a 2D isothermal vortex by a uniform mean flow in a periodic domain. 
 * Used as a validation case for the comparison of collision models in terms of stability 
 * and accuracy (Coreixas et al. 'Impact of collision models on the physical properties 
 * and the stability of lattice Boltzmann methods', 2020, Phil. Trans. R. Soc. A., 378).
 */

#include "palabos2D.h"
#include "palabos2D.hh"   // include full template code
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
#define DESCRIPTOR D2Q9Descriptor

// This structure holds all the user-defined parameters, and some
// derived values needed for the simulation.
template<typename T>
class Param {
public:
    plint N;                    // Resolution
    plint nx, ny;               // Domain size (LB units)
    T lx, ly;                   // Domain size (Physical units)
    T dx, dt;                   // Space and time steps
    T Re, Ma, tc;               // Dimensionless parameters
    T u0;                       // Velocity (LB units)
    T soundSpeed;               // Speed of sound (Physical units)
    T cSmago;                   // Smagorinsky constant (turbulent flows)
    T omega, tau, nu;           // Collision parameters (LB units)
    T omega3, omega4;           // Relaxation frequencies of 3rd and 4th order moments
    bool iniFneq;               // Init parameters (include nonequilibrium part in the init)
    plint intIniFneq;

    T tAdim, vtsT;              // Simulation time and frequency of .vti outputs
    std::string hoOmega;

    std::string lbm;
    std::string dynName, outDirName, logfileName;
    std::string fnameBase;

    Param() { }
    Param(std::string xmlFname) {

        //////// Parameters from the xml file
        XMLreader document(xmlFname);
        document["lattice"]["lbm"].read(lbm);
        document["lattice"]["dynName"].read(dynName);
        document["lattice"]["hoOmega"].read(hoOmega);
        document["simuParam"]["soundSpeed"].read(soundSpeed);
        document["simuParam"]["dx"].read(dx);
        document["simuParam"]["cSmago"].read(cSmago);
        document["initialization"]["intIniFneq"].read(intIniFneq);
        document["io"]["output"].read(outDirName);
        document["io"]["tAdim"].read(tAdim);
        document["io"]["vtsT"].read(vtsT);

        //////// Numerical discretization 
        // dx = 0.01;                                // Space step in physical units [m] (now read from .xml file  to avoid har coded value)
        T cs = ::sqrt(DESCRIPTOR<T>::cs2);           // Sound speed in LB units
        // T gamma = 1.4;                            // Specific heat ratio of air (diatomic gas)
        // T rGas = 287;                             // Gas constant                       
        // T Tref = 273.15 + 20.;                    // Ambient temperature (20°C)
        // soundSpeed = ::sqrt(gamma*rGas*Tref;      // Sound speed in physical units (m/s)
        //            = 343.20208332701009;
        // soundSpeed = 340.;                        // Simplified value (now read from .xml file to avoid hard coded value)       
        dt = (cs/soundSpeed)*dx;                     // Time step in physical units [s]
        
        //////// Simulation domain parameters
        global::argv(3).read(N);
        nx = N; ny = N; 
        lx = nx * dx; ly = ny * dx;

        //////// Dimensionless parameters
        global::argv(4).read(Ma); 
        u0 = (T)(Ma * cs);
        tc = (T)(N/u0);
        global::argv(2).read(nu);
        Re = (u0*N)/nu;
        tau = nu/DESCRIPTOR<T>::cs2;
        omega = 1./(tau + 0.5);
        if (hoOmega == "SRT") { 
            omega3 = omega;
            omega4 = omega;
        } else if (hoOmega == "REG") { 
            omega3 = 1.;
            omega4 = 1.;
        } else {
            pcout << "Error: Relaxation of high-order moments not correct." << std::endl;
            exit(-1);
        }

        //////// Improved initialization step
        iniFneq = (intIniFneq == 0 ? false : true);

        //////// File name for log and stats
        std::stringstream fnameBaseStr;
        fnameBaseStr << std::setprecision(7) << floor(1e6*nu); //Keeps only the first 6 decimals as an interger
        fnameBaseStr << "_";
        fnameBaseStr << N;
        fnameBaseStr << "_0_";
        fnameBaseStr << std::setprecision(7) << (int)(Ma*100.);
        fnameBaseStr << "_";
        fnameBaseStr << (iniFneq ? 1 : 0);
        fnameBaseStr >> fnameBase;
    }
    
    void writeLogFile() {
        plb_ofstream fout((outDirName+"/log_"+fnameBase+".dat").c_str());//========== Numerical Parameters ===========//

        fout << " //======== LBM Parameters ===============// " << std::endl;
        fout << "Lattice  -->           "<< lbm << std::endl;
        fout << "Dynamics -->           "<< dynName << std::endl;
        fout << "HO Relaxation Rype --> "<< hoOmega << std::endl;
        fout << std::endl;

        fout << " //======== Physical Parameters ==========// " << std::endl;
        fout << "Flow properties (dimensionless):    " << std::endl;
        fout << "Re = " << Re << std::endl;
        fout << "Ma = " << Ma << std::endl;
        fout << "Flow properties (physical units):    " << std::endl;
        fout << "nu = " << nu*dx*dx/dt << " [m2/s]" << std::endl;
        fout << "c  = " << soundSpeed << " [m/s]" << std::endl;
        fout << "u0 = " << Ma * ::sqrt(DESCRIPTOR<T>::cs2) * dx/dt << " [m/s]" << std::endl;
        fout << "tc = " << tc * dt << " [s]" << std::endl;
        fout << "Geometry (physical units):    " << std::endl;
        fout << "lx = " << lx << " [m]" << std::endl;
        fout << "ly = " << ly << " [m]" << std::endl;
        fout << std::endl;

        fout << " //======== Numerical Parameters =========// " << std::endl;
        fout << "Numerical discretization (physical units):    " << std::endl;        
        fout << "dx = " << dx << " [m]" << std::endl;
        fout << "dt = " << dt << " [s]" << std::endl;
        fout << "Geometry (LB units):    " << std::endl;
        fout << "N  = " << N << " (resolution)" << std::endl;
        fout << "nx = " << nx << std::endl;
        fout << "ny = " << ny << std::endl;
        fout << "Flow properties (LB units):    " << std::endl;
        fout << "nuLB = " << nu << std::endl;
        fout << "u0LB = " << Ma * ::sqrt(DESCRIPTOR<T>::cs2) << std::endl;
        fout << "tcLB = " << round(tc) << " (" << tc << ")" << std::endl;
        fout << "Collision parameters (LB units):    " << std::endl;
        fout << "tau = " << tau << std::endl;
        fout << "omega = " << omega << std::endl;
        if (lbm == "D2Q9"){
            fout << "omega3 = " << omega3 << std::endl;
            fout << "omega4 = " << omega4 << std::endl;
        }
        fout << "Large Eddy Simulation parameters:    " << std::endl;
        fout << "cSmago = " << cSmago << std::endl;
        fout << std::endl;

        fout << " //======== Improved Initialization ======// " << std::endl;
        fout << "Initializes fNeq (using FD gradients): ";
        if (iniFneq) fout << "True" << std::endl;
        else fout << "False" << std::endl;
        fout << std::endl;

        fout << " //======== Simulation parameters ======// " << std::endl;
        fout << "output= " << outDirName << std::endl;
        fout << "tAdim = " << tAdim << " * tc" << std::endl;
        fout << "      = " << (int)(tAdim * tc) * dt << " [s]" << std::endl;
        fout << "      = " << (int)(tAdim * tc) << " [iterations]" << std::endl;
        fout << "vtsT  = " << vtsT << " * tc" << std::endl;
        fout << "      = " << (int)(vtsT * tc) * dt << " [s]" << std::endl;
        fout << "      = " << (int)(vtsT * tc) << " [iterations]" << std::endl;
    }
};

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor> *getDynamics(Param<T> &param)
{
    Dynamics<T,Descriptor> *dyn;
    if (param.dynName == "BGK_Ma2") { 
        dyn = new BGKdynamics<T,Descriptor>(param.omega); // Second-order equilibrium
    } else if (param.dynName == "RM") { 
        dyn = new RMdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "HM") { 
        dyn = new HMdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "CM") { 
        dyn = new CMdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "CHM") { 
        dyn = new CHMdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "K") { 
        dyn = new Kdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "GH") { 
        dyn = new GHdynamics<T,Descriptor>(param.omega);
    } else if (param.dynName == "RR") { 
        dyn = new RRdynamics<T,Descriptor>(param.omega);
    } else {
        pcout << "Error: dynamics name does not exist." << std::endl;
        exit(-1);
    }

    return dyn;
}


//////// A functional, used to initialize the simulation
template<typename T>
class IsothermalCoVoInitialDensityAndVelocityField {
public:
    IsothermalCoVoInitialDensityAndVelocityField(Param<T> const& param_)
        : param(param_)
    { }
    
    void operator()(plint iX, plint iY, T &rho, Array<T,2>& u) const {

        const plint nx = param.nx;
        const plint ny = param.ny;
        
        T beta =0.5;
        T u0 = param.u0;
        T bU = beta*u0;
        T cs = ::sqrt(DESCRIPTOR<T>::cs2);

        T xCentered = (T)(iX-(nx-1.)/2.)/(T)(nx);
        T yCentered = (T)(iY-(ny-1.)/2.)/(T)(ny);
        T R = 1./20.;
        T r = (xCentered*xCentered+yCentered*yCentered)/(R*R);

        rho = (T)( (1.-0.5*(bU)*(bU)*exp(-r)/(cs*cs)) );
        u[0] = u0*(1.-beta*(yCentered)/R * exp(-r/2.));
        u[1] = bU*(xCentered)/R * exp(-r/2.);
    }
    
private:
    Param<T> param;
};


//////// Initialize the simulation
void simulationSetup(MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                     const Param<T> &param)
{
    // Set periodic boundaries.
    lattice.periodicity().toggleAll(true); 

    // Initialize the simulation domain.
    initializeAtEquilibrium (
        lattice, lattice.getBoundingBox(),
        IsothermalCoVoInitialDensityAndVelocityField<T>(param) );

    // Call initialize to get the lattice ready for the simulation.
    lattice.initialize();
}
void accurateInitialCondition(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, const Param<T> &param, T epsilon)
{
    std::unique_ptr<MultiScalarField2D<T> > density;
    density = computeDensity(lattice);
    if (param.iniFneq) {
        std::unique_ptr<MultiTensorField2D<T,2> > velocity = computeVelocity(lattice);
        std::unique_ptr<MultiTensorField2D<T,3> > S = computeStrainRate(*velocity);
        recomposeFromFlowVariables(lattice, *density, *velocity, *S);
    }
}


//////// Post processing
/// Produce a GIF snapshot of the velocity-norm.
void writeGif(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, plint iter)
{
    const plint imSize = 600;

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("u", iter, 6),
                            *computeVelocityNorm(lattice), imSize, imSize);
}
/// Write the full velocity and the velocity-norm into a VTK file.
void writeVTS(MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
              Param<T> &param, plint iter)
{
    T dx = param.dx;
    T dt = param.dt;
    VtkStructuredImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "machNorm", 1.0/std::sqrt(DESCRIPTOR<T>::cs2));
    vtkOut.writeData<float>(*computeDensity(lattice), "density", 1.0);
    vtkOut.writeData<2,float>(*computeGradient(*computeDensity(lattice)), "gradRho", 1.0/dx);
    vtkOut.writeData<2,float>(*computeVelocity(lattice), "velocity", dx/dt);
    vtkOut.writeData<float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", (T)1/dt);
    vtkOut.writeData<3,float>(*computeStrainRate(*computeVelocity(lattice)), "strain", (T)1/dt);
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    if (argc != 5)
    {
        pcout << argc << std::endl;
        pcout << "Error! Wrong number of parameters." << std::endl;
        pcout << "Syntax: " << (std::string)global::argv(0) << " config.xml nu N Ma" << std::endl;
        pcout << "Example: " << (std::string)global::argv(0) << " config.xml 1e-5 256 0.1" << std::endl;
        pcout << "Example (mpi): mpirun -np 2 " << (std::string)global::argv(0) << " config.xml 1e-5 256 0.1" << std::endl;
        exit(1);
    }

    ///// Read and print the parameters of the simulation
    Param<T> param(argv[1]);
    param.writeLogFile();
    ///// Output directory
    global::directories().setOutputDir(param.outDirName);

    ///// Initialize the diagonal relaxation matrix from the .xml file.
    Array<T, DESCRIPTOR<T>::numRelaxationTimes> allOmega;
    if (param.lbm == "D2Q9"){
        allOmega[0] = param.omega;  // relaxation of M20 and M02
        allOmega[1] = param.omega;  // relaxation of M11 
        allOmega[2] = param.omega3; // relaxation of M21 and M12
        allOmega[3] = param.omega4; // relaxation of M22
    } else {
        pcout << "Error: lbm name does not exist." << std::endl;
        exit(-1);
    }
    
    ///// Generate the dynamics and the corresponding lattice from the .xml file.
    Dynamics<T,DESCRIPTOR> *dyn = getDynamics<T,DESCRIPTOR>(param);
    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
              param.nx, param.ny, dyn);

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

    lattice.toggleInternalStatistics(false);

    ///// Initialization from analytical profiles
    simulationSetup(lattice, param);

    //// Either add fNeqor do nothing
    accurateInitialCondition(lattice, param,1.e-4);

    ///// Initial state is outputed
    writeVTS(lattice, param, 0);    

    ///// Simulation maximal time, and output frequency (in terms of iterations).
    plint vtsTout = param.vtsT * param.tc;
    plint tmax = param.tAdim * param.tc;

    ///// Main loop over time iterations.
    for (plint iT=0; iT<=tmax; ++iT) {

        ///// Lattice Boltzmann iteration step.
        lattice.collideAndStream();
        ///// Output.
        if (iT % vtsTout == 0) {
            // pcout << "Writing VTS file at iteration = " << iT << std::endl;
            pcout << "Writing VTS file at t = " << iT/param.tc << " (iteration " << iT << ")" << std::endl;
            // pcout << "maximal velocity = " << computeMax(*computeVelocityNorm(lattice)) << std::endl;
            writeVTS(lattice, param, iT);
        }
    }

    return 0;
}
