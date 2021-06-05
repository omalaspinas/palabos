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

#ifndef UTILITY_DSL3D_PARAM_H
#define UTILITY_DSL3D_PARAM_H

namespace plb {
    
// This structure holds all the user-defined parameters, and some
// derived values needed for the simulation.
template<typename T>
class Param {
public:
    plint N;                    // Resolution
    plint nx, ny, nz;           // Domain size (LB units)
    T lx, ly, lz;               // Domain size (Physical units)
    T dx, dt;                   // Space and time steps
    T Re, Ma, tc;               // Dimensionless parameters
    T u0;                       // Velocity (LB units)
    T soundSpeed;               // Speed of sound (Physical units)
    T cs;                       // Speed of sound (LB units)
    T omega, tau, nu;           // Collision parameters (LB units)
    T omegaBulk;                // Relaxation frequency for the bulk viscosity
    T omega3, omega4;           // Relaxation frequencies of 3rd and 4th order moments
    T omega5, omega6;           // Relaxation frequencies of 5th and 6th order moments

    T tAdim, vtkT;              // Simulation time and frequency of .vti outputs
    std::string hoOmega;
    bool bulk;
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
        document["lattice"]["bulk"].read(bulk);
        document["simuParam"]["soundSpeed"].read(soundSpeed);
        document["simuParam"]["dx"].read(dx);
        document["simuParam"]["nz"].read(nz);
        document["io"]["output"].read(outDirName);
        document["io"]["tAdim"].read(tAdim);
        document["io"]["vtkT"].read(vtkT);

        //////// Numerical discretization 
        // dx = 0.01;                             // Space step in physical units [m] (now read from .xml file  to avoid hard coded value)
        T cs  = ::sqrt(1./3.);                    // Sound speed in LB units (aka, lattice constant) 
        // T gamma = 1.4;                         // Specific heat ratio of air (diatomic gas)
        // T rGas = 287;                          // Gas constant                       
        // T Tref = 273.15 + 20.;                 // Ambient temperature (20°C)
        // soundSpeed = ::sqrt(gamma*rGas*Tref);  // Sound speed in physical units (m/s)
        //            = 343.20208332701009;
        // soundSpeed = 343.;                     // Simplified value (now read from .xml file to avoid hard coded value)       
        dt = (cs/soundSpeed)*dx;                  // Time step in physical units [s] computed from the speed of sound (acoustic scaling)
        
        //////// Simulation domain parameters
        global::argv(3).read(N);
        nx = N; ny = N; //nz is given in the xml. file to avoid any hard coded value 
        lx = nx * dx; ly = ny * dx; lz = nz * dx;

        //////// Dimensionless parameters
        global::argv(4).read(Ma); 
        u0 = (T)(Ma * cs);                        // Ma = u0/cs
        tc = (T)(N/u0);
        global::argv(2).read(Re);
        nu = (u0*N)/Re;                           // Re = (u0*N)/nu
        tau = nu/(cs*cs);
        omega = 1./(tau + 0.5);

        //////// Bulk and high-order relaxation parameters
        omegaBulk = bulk? (T)1.: omega;   // Relaxation of bulk viscosity related moments (M20 + M02)
        if (hoOmega == "SRT") {           // Single relaxation time formulation
            omega3 = omega;
            omega4 = omega;
            omega5 = omega;
            omega6 = omega;
        } else if (hoOmega == "REG") {    // Regularization of high-order moments
            omega3 = 1.;
            omega4 = 1.;
            omega5 = 1.;
            omega6 = 1.;
        } else {
            pcout << "Error: Relaxation of high-order moments not correct, please choose either SRT or REG." << std::endl;
            exit(-1);
        }

        //////// File name for log and stats
        std::stringstream fnameBaseStr;
        fnameBaseStr << std::setprecision(7) << Re;
        fnameBaseStr << "_";
        fnameBaseStr << N;
        fnameBaseStr << "_0_";
        fnameBaseStr << std::setprecision(7) << (int)(Ma*100.);
        fnameBaseStr >> fnameBase;
    }
    
    void writeLogFile() {
        plb_ofstream fout((outDirName+"/log_"+fnameBase+".dat").c_str());//========== Numerical Parameters ===========//

        fout << " //======== LBM Parameters ===============// " << std::endl;
        fout << "Lattice  -->           "<< lbm << std::endl;
        fout << "Dynamics -->           "<< dynName << std::endl;
        fout << "HO Relaxation Type --> "<< hoOmega << std::endl;
        fout << std::endl;

        fout << " //======== Physical Parameters ==========// " << std::endl;
        fout << "Flow properties (dimensionless):    " << std::endl;
        fout << "Re = " << Re << std::endl;
        fout << "Ma = " << Ma << std::endl;
        fout << "Flow properties (physical units):    " << std::endl;
        fout << "nu = " << nu*dx*dx/dt << " [m2/s]" << std::endl;
        fout << "c  = " << soundSpeed << " [m/s]" << std::endl;
        fout << "u0 = " << Ma * cs * dx/dt << " [m/s]" << std::endl;
        fout << "tc = " << tc * dt << " [s]" << std::endl;
        fout << "Geometry (physical units):    " << std::endl;
        fout << "lx = " << lx << " [m]" << std::endl;
        fout << "ly = " << ly << " [m]" << std::endl;
        fout << "lz = " << lz << " [m]" << std::endl;
        fout << std::endl;

        fout << " //======== Numerical Parameters =========// " << std::endl;
        fout << "Numerical discretization (physical units):    " << std::endl;        
        fout << "dx = " << dx << " [m]" << std::endl;
        fout << "dt = " << dt << " [s]" << std::endl;
        fout << "Geometry (LB units):    " << std::endl;
        fout << "N  = " << N << " (resolution)" << std::endl;
        fout << "nx = " << nx << std::endl;
        fout << "ny = " << ny << std::endl;
        fout << "nz = " << nz << std::endl;
        fout << "Flow properties (LB units):    " << std::endl;
        fout << "nuLB = " << nu << std::endl;
        fout << "u0LB = " << Ma * cs << std::endl;
        fout << "tcLB = " << round(tc) << " (" << tc << ")" << std::endl;
        fout << "Collision parameters (LB units):    " << std::endl;
        fout << "tau = " << tau << std::endl;
        fout << "omega = " << omega << std::endl;
        fout << "omegaBulk = " << omegaBulk << std::endl;
        if (lbm == "D3Q19"){
            fout << "omega3 = " << omega3 << std::endl;
            fout << "omega4 = " << omega4 << std::endl;
        } else if (lbm == "D3Q27"){
            fout << "omega3 = " << omega3 << std::endl;
            fout << "omega4 = " << omega4 << std::endl;
            fout << "omega5 = " << omega5 << std::endl;
            fout << "omega6 = " << omega6 << std::endl;
        }

        fout << " //======== Simulation parameters ======// " << std::endl;
        fout << "output= " << outDirName << std::endl;
        fout << "tAdim = " << tAdim << " * tc" << std::endl;
        fout << "      = " << (int)(tAdim * tc) * dt << " [s]" << std::endl;
        fout << "      = " << (int)(tAdim * tc) << " [iterations]" << std::endl;
        fout << "vtkT  = " << vtkT << " * tc" << std::endl;
        fout << "      = " << (int)(vtkT * tc) * dt << " [s]" << std::endl;
        fout << "      = " << (int)(vtkT * tc) << " [iterations]" << std::endl;
    }
};

}  // namespace plb

#endif  // UTILITY_DSL3D_PARAM_H
