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

#include "palabos3D.h"
#include "palabos3D.hh"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor

#define NMAX 150

const T pi = (T)4.*std::atan((T)1.);

struct Param
{
    plint numOutletSpongeCells;  // Number of the lattice nodes contained in the outlet sponge zone.
    int outletSpongeZoneType;    // Type of the outlet sponge zone (Viscosity or Smagorinsky).
    T targetSpongeCSmago;        // Target Smagorinsky parameter at the end of the Smagorinsky sponge Zone.
};

Param param;

template<typename T>
class CylinderShapeDomain3D : public plb::DomainFunctional3D {
public:
    CylinderShapeDomain3D(plb::plint cx_, plb::plint cz_, plb::plint radius)
        : cx(cx_),
        cz(cz_),
        radiusSqr(plb::util::sqr(radius))
    { }
    virtual bool operator() (plb::plint iX, plb::plint iY, plb::plint iZ) const {
        return plb::util::sqr(iX - cx) + plb::util::sqr(iZ - cz) <= radiusSqr;
    }
    virtual CylinderShapeDomain3D<T>* clone() const {
        return new CylinderShapeDomain3D<T>(*this);
    }
private:
    plb::plint cx;
    plb::plint cz;
    plb::plint radiusSqr;
};

static T poiseuillePressure(IncomprFlowParam<T> const &parameters, plint maxN)
{
    const T a = parameters.getNx()-1;
    const T b = parameters.getNy()-1;

    const T nu = parameters.getLatticeNu();
    const T uMax = parameters.getLatticeU();

    T sum = T();
    for (plint iN = 0; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum += ((T)1 / (std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a))));
    }
    for (plint iN = 1; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum -= ((T)1 / (std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a))));
    }

    T alpha = -(T)8 * uMax * pi * pi * pi / (a*a*(pi*pi*pi-(T)32*sum)); // alpha = -dp/dz / mu

    T deltaP = - (alpha * nu);

    return deltaP;
}

T poiseuilleVelocity(plint iX, plint iY, IncomprFlowParam<T> const& parameters, plint maxN)
{
    const T a = parameters.getNx()-1;
    const T b = parameters.getNy()-1;

    const T x = (T)iX - a / (T)2;
    const T y = (T)iY - b / (T)2;

    const T alpha = - poiseuillePressure(parameters,maxN) / parameters.getLatticeNu();

    T sum = T();

    for (plint iN = 0; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;

        sum += (std::cos(twoNplusOne*pi*x/a)*std::cosh(twoNplusOne*pi*y/a)
             / ( std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a)) ));
    }
    for (plint iN = 1; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;

        sum -= (std::cos(twoNplusOne*pi*x/a)*std::cosh(twoNplusOne*pi*y/a)
             / ( std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a)) ));
    }

    sum *= ((T)4 * alpha * a *a /std::pow(pi,(T)3));
    sum += (alpha / (T)2 * (x * x - a*a / (T)4));
    
    return sum;
}

template <typename T>
class SquarePoiseuilleDensityAndVelocity {
public:
    SquarePoiseuilleDensityAndVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
        : parameters(parameters_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, T &rho, Array<T,3>& u) const {
        rho = (T)1;
        u[0] = T();
        u[1] = T();
        u[2] = poiseuilleVelocity(iX, iY, parameters, maxN);
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

template <typename T>
class SquarePoiseuilleVelocity {
public:
    SquarePoiseuilleVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
        : parameters(parameters_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const  {
        u[0] = T();
        u[1] = T();
        u[2] = poiseuilleVelocity(iX, iY, parameters, maxN);
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

void simulationSetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                      IncomprFlowParam<T> const& parameters,
                      OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition,
                      Array<plint,3> &forceIds )
{
    // No periodic boundaries
    lattice.periodicity().toggleAll(false);

    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();

    // Flow direction: z
    // Top/Bottom direction: y
    Box3D top    = Box3D(0, nx-1, ny-1, ny-1, 0, nz-1); // Full Area
    Box3D bottom = Box3D(0, nx-1, 0, 0, 0, nz-1); // Full Area
    
    Box3D inlet  = Box3D(0, nx-1, 0, ny-1, 0, 0); // Full Area
    Box3D outlet = Box3D(1, nx-2, 1, ny-2, nz-1, nz-1); // Offset from wall boundaries by 1 lattice unit
    
    Box3D right  = Box3D(0, 0, 0, ny-1, 0, nz-1); // Full Area
    Box3D left   = Box3D(nx-1, nx-1, 0, ny-1, 0, nz-1); // Full Area
    
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, inlet);
    boundaryCondition.addVelocityBoundary2P(outlet, lattice, boundary::neumann);

    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, top);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, bottom);
    
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, left);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, right);
    
    setBoundaryVelocity(lattice, inlet, SquarePoiseuilleVelocity<T>(parameters, NMAX));

    setBoundaryVelocity(lattice, top, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, bottom, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, left, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, right, Array<T,3>((T)0.0,(T)0.0,(T)0.0));

    /*
    * Implement the outlet sponge zone.
    */
    if (param.numOutletSpongeCells > 0)
    {
        T bulkValue;
        Array<plint, 6> numSpongeCells;

        if (param.outletSpongeZoneType == 0) {
            pcout << "Generating an outlet viscosity sponge zone." << std::endl;
            bulkValue = parameters.getOmega();
        }
        else if (param.outletSpongeZoneType == 1) {
            pcout << "Generating an outlet Smagorinsky sponge zone." << std::endl;
            bulkValue = 0.14; // Parameter for the Smagorinsky LES model
        }
        else {
            pcout << "Error: unknown type of sponge zone." << std::endl;
            exit(-1);
        }

        // Number of sponge zone lattice nodes at all the outer domain boundaries.
        // So: 0 means the boundary at x = 0
        //     1 means the boundary at x = nx-1
        //     2 means the boundary at y = 0
        //     and so on...
        numSpongeCells[0] = 0;
        numSpongeCells[1] = 0;
        numSpongeCells[2] = 0;
        numSpongeCells[3] = 0;
        numSpongeCells[4] = 0;
        numSpongeCells[5] = param.numOutletSpongeCells;

        std::vector<MultiBlock3D*> args;
        args.push_back(&lattice);

        if (param.outletSpongeZoneType == 0) {
            applyProcessingFunctional(new ViscositySpongeZone3D<T, DESCRIPTOR>(
                nx, ny, nz, bulkValue, numSpongeCells),
                lattice.getBoundingBox(), args);
        }
        else {
            applyProcessingFunctional(new SmagorinskySpongeZone3D<T, DESCRIPTOR>(
                nx, ny, nz, bulkValue, param.targetSpongeCSmago, numSpongeCells),
                lattice.getBoundingBox(), args);
        }
    }

    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), SquarePoiseuilleDensityAndVelocity<T>(parameters, NMAX));

    // Add the obstacle: cylinder 3d
    plint cx     = nx/2 + 2; // cx is slightly offset to avoid full symmetry,
                             // and to get a Von Karman Vortex street.
    plint cz     = nz/6;
    plint radius = parameters.getResolution() / 2; // the diameter is the reference length

    lattice.toggleInternalStatistics(true);
    forceIds[0] = lattice.internalStatSubscription().subscribeSum();
    forceIds[1] = lattice.internalStatSubscription().subscribeSum();
    forceIds[2] = lattice.internalStatSubscription().subscribeSum();

    defineDynamics(lattice, lattice.getBoundingBox(),
                   new CylinderShapeDomain3D<T>(cx, cz, radius),
                   new plb::MomentumExchangeBounceBack<T,DESCRIPTOR>(forceIds));
    initializeMomentumExchange(lattice, lattice.getBoundingBox());

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

    Box3D slice(0, nx - 1, ny / 2, ny / 2, 0, nz - 1);
    ImageWriter<T> imageWriter("leeloo");

    imageWriter.writeScaledGif(createFileName("uNorm", iter, 6),
        *computeVelocityNorm(lattice, slice),
        imSize, imSize);
}

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    if (argc != 10)
    {
        pcout << "Error: the parameters are wrong. \n";
        pcout << "Give Re Resolution CylinderDiameter ChannelWidth ChannelHeight ChannhelLength outletSpongeZoneType numOutletSpongeCells targetSpongeCSmago \n";
        pcout << "Example: ./executable 100 20 6 50 50 75 1 1 0.6 \n";
        exit(1);
    }

    T Re_    = atof(argv[1]);
    plint N_ = atoi(argv[2]);
    // Physical units, reference length is the cylinder diameter (D_)
    T D_     = atof(argv[3]);
    T W_     = atof(argv[4]);
    T h_     = atof(argv[5]);
    T L_     = atof(argv[6]);
    // Sponge Zone related
    param.outletSpongeZoneType = atoi(argv[7]); // Outlet sponge zone type ("Viscosity":0 or "Smagorinsky":1)
    param.numOutletSpongeCells = atoi(argv[8]); // Number of the lattice nodes contained in the outlet sponge zone.
    param.targetSpongeCSmago   = atof(argv[9]); // This parameter is relevant only if the chosen sponge zone is of "Smagorinsky" type.

    // Use the class IncomprFlowParam to convert from
    // dimensionless variables to lattice units, in the
    // context of incompressible flows.
    IncomprFlowParam<T> parameters(
        0.1,    // Reference velocity (the maximum velocity in the Poiseuille profile) in lattice units.
        Re_,    // Reynolds number
        N_,     // Resolution of the reference length (cylinder diameter)
        W_/D_,  // dimensionless: channel lateral length
        h_/D_,  // dimensionless: channel height
        L_/D_   // dimensionless: channel length
    );

    const T vtkSave = (T) 0.1; // Time intervals at which to save GIF VTKs, in dimensionless time units
    const T maxT    = (T)200.0; // Total simulation time, in dimensionless time units

    pcout << "omega= " << parameters.getOmega() << std::endl;
    writeLogFile(parameters, "3D square Poiseuille with Cylinder as an obstacle");

    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
        parameters.getNx(), parameters.getNy(), parameters.getNz(), 
        new BGKdynamics<T, DESCRIPTOR>(parameters.getOmega()) );

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    Array<plint,3> forceIds;
    simulationSetup(lattice, parameters, *boundaryCondition, forceIds);

    // Loop over main time iteration.
    for (plint iT=0; iT<parameters.nStep(maxT); ++iT)
    {
        if (iT%parameters.nStep(vtkSave)==0 && iT>0)
        {
            pcout << "step " << iT << "; t=" << iT*parameters.getDeltaT() << std::endl;
            
            writeGifs(lattice, parameters, iT);
        }

        // Execute a time iteration.
        lattice.collideAndStream();
    }
    
    delete boundaryCondition;
}
