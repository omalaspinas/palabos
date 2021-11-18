#ifndef SIMPARAMS // or some other unique name
#define SIMPARAMS

#include "palabos3D.h"
#include "palabos3D.hh"
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>      /* srand, rand */
#include <ctime>        /* time */

#define DESCRIPTOR descriptors::RhoBarJD3Q19Descriptor

typedef double T;
typedef plb::Array<T,3> Velocity;

using namespace plb;


// actually counts nx based on clot file
double countNXY(std::string clotGeometryFile)
{
    std::fstream clotFile;
    clotFile.open(clotGeometryFile.c_str(), std::fstream::in);
    std::string line;
    getline(clotFile,line);
    std::stringstream ss(line);

    std::string val;

    size_t count=0;
    while(ss>>val){ // ss is used more like cin
        count++;
    }
    clotFile.close();
    return count;
}

// functional that describes the pipe domain
template<typename T>
class Pipe : public DomainFunctional3D {
public:
    Pipe(Array<T,3> center_={0,0,0}, plint radius_=0, plint length_=0, bool isSquarePipe_=false)
        : center(center_),
          radius(radius_),
          length(length_),
          isSquarePipe(isSquarePipe_)
    { }

    virtual bool operator() (Array<T,3>const & pos) const
    {
        return !this->operator()(pos[0], pos[1], pos[2]);
    }

    virtual bool operator() (plint iX, plint iY, plint iZ) const
    {
        if (!isSquarePipe)
            return(util::sqr((T) iX - center[0]) + util::sqr((T) iY - center[1]) > util::sqr(radius) || (iZ < center[2]) || (iZ>center[2]+length));
        else
            return (util::sqr((T) iX - center[0]) > util::sqr(radius) || util::sqr((T) iY - center[1]) > util::sqr(radius) || (iZ < center[2]) || (iZ>center[2]+length)); // square pipe
    }

    virtual Pipe<T>* clone() const
    {
        return new Pipe<T>(*this);
    }

    const bool inPipeSection(plint iX, plint iY) const
    {
        if (!isSquarePipe)
            return(util::sqr((T) iX - center[0]) + util::sqr((T) iY - center[1]) < util::sqr(radius));
        else
            return (util::sqr((T) iX - center[0]) < util::sqr(radius) && util::sqr((T) iY - center[1]) < util::sqr(radius));
    }

    const plint getRadius() const {return radius;}
    const Array<T,3> getCenter() const{return center;}
    const plint getLength() const {return length;}
private:
    Array<T,3> center;
    plint radius;
    plint length;
    bool isSquarePipe;
};

/// Linearly decreasing pressure profile
T poiseuillePressure(plint iZ, const Pipe<T>& pipe, T DPAdim) {
    T Lz = pipe.getLength();
    T z = iZ; 
    return DPAdim * (Lz-z+pipe.getCenter()[2])/Lz;
}

/// Convert pressure to density according to ideal gas law
T poiseuilleDensity(plint iZ, const Pipe<T>& pipe, T DPAdim) {
    return poiseuillePressure(iZ,pipe, DPAdim)*DESCRIPTOR<T>::invCs2 + (T)1;
}

/// A functional, used to create an initial condition for with zero velocity,
///   and linearly decreasing pressure.
template<typename T>
class PoiseuilleDensityAndZeroVelocity {
public:
    PoiseuilleDensityAndZeroVelocity(const Pipe<T>& pipe_, T DPAdim_)
        : pipe(pipe_), DPAdim(DPAdim_)
    { }
    void operator()(plint iX, plint iY, plint iZ, T& rho, Array<T,3>& u) const {
        rho = poiseuilleDensity(iZ, pipe, DPAdim);
        u[0] = T();
        u[1] = T();
        u[2] = T();
    }
private:
    Pipe<T> pipe;
    T DPAdim;
};

class SimulationParameters {
    /*
     * Parameters set by the user.
     */
public:
    Cuboid<T> fullDomain;   // Full simulation domain.in physical units
    T DeltaP = 1;//200, k1 5e3 OK         // Pressure drop along pipe, in Pa
    T nu = 3e-6;            // Fluid kinematic viscosity in physical units.
    T rho = 1e3;            // Fluid density in physical units
    
    std::string clotGeometryFile = "clot_0.0003.txt";   // file from which to read the solidFraction of the clot   
    double k1 = 5e3;    // degradation rate of fibrin
    bool interactNeighbor = true;   // particles interact also with diagonal clot voxels if true
    double clotBeginZ = 0.;     // ratio of z length where clot starts
    double clotEndZ;    // determined by the clot geometry file
    // quantities for computation of Rf(x,t), for relation between qty of Fn and porosity locally (Diamond and Anand)
    double Rf0 = 250; // fibrin fiber radius nanometers
    double p0 = 0.01116; // fibrin fiber density Fb/nm^2
    std::string permeModel = "Davies"; // gives relation k(ns*) (which uses value of Rf0). Options : Davies, Clague, JJ (Jackson-James), or None (gamma = ns)
    
    bool usesParticles = true;     // no particles = no fibrinolysis
    plint particlesNum = 100.;     // how many particles to inject initially (NB: doesn't affect concentration of anti-fibrin)
    double antiFnConcentration = 10.e-4;     // initial concentration of anti-fibrin, in mol/m³ (Franck's particle model). 1nM = 1e-6 mol/m³
    double k2 = 1.;     // rate of degradation of anti-fibrin
    bool constantConc = true;
    
    T mach = 0.1;       // mach number
    T maxTime = 50;   // Max simulation time, in s.
    double startClotTime = 0.;      // when to put the clot, in s
    double startParticleTime = .2;  // when to start injecting particles, in s
    double stopParticleTime = 1000.;    // when to stop injecting particles, in s

    plint outIter = 1000;                                  // Number of iterations for .vtk output.
    plint statIter = 1000;                                 // Number of iterations for screen and .dat output.
    bool writeVtk = true;                                  // write vtk output or not
    std::string outDir = "tmp/";
    std::string clotStatFilename = "clotStat.dat";
    std::string partStatFilename = "partStat.dat";
    bool overwriteVtk = false;  // if true, simulation will only produce one set of vtk files at start of reaction, and one at the end (to reduce disk usage)

    double avogadro = 6.022141e23;
    double pi = (double)4.*std::atan((double)1.);

    /*
     * Parameters NOT set by the user.
     */

    // for the particles
    plint commEnvelope;

    // domain
    plint nx, ny, nz;
    T dx, dt; // physical units, the corresponding lb are units
    T u;    // Fluid velocity in physical units.
    T rho_LB, nu_LB, u_LB;  // quantities in lattice units
    T reynolds; // reynolds number, indicative value
    Array<T,3> physicalLocation; // The origin of the system
    T omega;    // relaxation rate

    plint antiFnQtyPerCell;
    double V0antiFn;    // volume before clot where particles are injected
    double rhoWK_t_in;    // inlet density

    double C_P; // conversion factor for pressure : P = C_P*P_adim ; and rho_adim = P/(C_P*cs2) + rho_LB = P_adim/cs2 + rho_LB

    // params for checkpoint
    plint iniIter = 0;

    SimulationParameters(){}
    
    void printSimulationParameters()
    {
        pcout << "u (m/s) = " << u << std::endl;
        pcout << "u_LB = " << u_LB << std::endl;
        pcout << "maxSimTime (s) = " << maxTime << std::endl;
        pcout << "rho (kg/m-3) = " << rho << std::endl;
        pcout << "nu = " << nu << std::endl;
        pcout << "nu_LB = " << nu_LB << std::endl;
        pcout << "statIter = " << statIter << std::endl;
        pcout << "outIter = " << outIter << std::endl;
        pcout << "rho_LB = " << rho_LB << std::endl;
        pcout << "omega = " << omega << std::endl;
        pcout << "tau = " << (T) 1 / omega << std::endl;
        pcout << "dx (m) = " << dx<< std::endl;
        pcout << "dt (s) = " << dt << std::endl;
        pcout << "nx = " << nx << std::endl;
        pcout << "ny = " << ny << std::endl;
        pcout << "nz = " << nz << std::endl;
        pcout << "clotBeginZ = " << clotBeginZ << std::endl;
        pcout << "DeltaP (Pa) = " << DeltaP << std::endl;
        pcout << "Reynolds = " << reynolds << std::endl;
        pcout << "Mach = " << mach << std::endl;
        pcout << "antiFibrinQtyPerVoxel = " << antiFnQtyPerCell << std::endl;
        pcout << "Rf0 (nm) = " << Rf0 << std::endl;
        pcout << "p0 = " << p0 << std::endl;
        pcout << "k1 = " << k1 << std::endl;
        pcout << "clotFile = " << clotGeometryFile << std::endl;
        pcout << std::endl;
    }

    std::unique_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > setupProblem(
        MultiParticleField3D<DenseParticleField3D<T,DESCRIPTOR>>& particles, 
        MultiScalarField3D<plint> & clotFlags, MultiScalarField3D<T> & clotSolidFraction, MultiScalarField3D<T> & clotSolidFractionPhys)
    {

        // initialization of the domain
        fullDomain.lowerLeftCorner[0] = 0;          // lower left x bound, in m
        fullDomain.lowerLeftCorner[1] = 0;          // lower left y bound, in m
        fullDomain.lowerLeftCorner[2] = 0;          // lower left z bound, in m
        fullDomain.upperRightCorner[0] = 0.0058;    // upper right x bound, in m
        fullDomain.upperRightCorner[1] = 0.0058;    // upper right y bound, in m
        fullDomain.upperRightCorner[2] = 0.006495652;      // upper right z bound, in m
        physicalLocation = Array<T,3>(fullDomain.x0(), fullDomain.y0(), fullDomain.z0());
        dx = (fullDomain.x1() - fullDomain.x0()) / (countNXY(clotGeometryFile));
        pcout << "domain : " << fullDomain.x1() - fullDomain.x0() << "; countNXY : " << (countNXY(clotGeometryFile)) <<
        "; nx : " << T(fullDomain.x1() - fullDomain.x0()) / T(dx) << std::endl;
        commEnvelope = (plint) 1.5 + 1;
        nx = T(fullDomain.x1() - fullDomain.x0()) / T(dx) -2;
        ny = T(fullDomain.y1() - fullDomain.y0()) / T(dx) -2;
        nz = T(fullDomain.z1() - fullDomain.z0()) / T(dx) -2;

        reynolds = sqrt(DeltaP*pow(nx*dx,3.)*rho/(nz*dx))/(nu*rho);
        u_LB = mach*sqrt(DESCRIPTOR<T>::cs2);
        dt = sqrt(dx*dx/(3.*DESCRIPTOR<T>::cs2));
        nu_LB = nu * dt / (dx*dx);
        u=u_LB*dx/dt;
        rho_LB = 1.; 
        omega = (T) 1 / (DESCRIPTOR<T>::invCs2 * nu_LB + (T) 0.5);
        C_P = rho * (dx * dx) / (dt * dt);  // conversion factor numerical -> physical pressure

        // definition of the lattice
        std::unique_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > lattice = 
        generateMultiBlockLattice<T,DESCRIPTOR>(Box3D(0,nx-1, 0,ny-1, 0,nz-1), 
            new PartialBBdynamics<T,DESCRIPTOR>(omega), commEnvelope);
        lattice->toggleInternalStatistics(false);
        lattice->periodicity().toggleAll(false);
        
        // Pipe
        plint minN = std::min(nx, ny);
        T pipeRadius = floor((T) 0.5 * (T) minN);
        T pipeLength = (T) 1.*(nz);
        Array<T,3> pipeCenter(pipeRadius, pipeRadius, (T) 0.0 * (nz - 1));
        Pipe<T> poiseuillePipe(pipeCenter, pipeRadius, pipeLength,false);

        // declare the clot flags and SF for checkpointing load
        clotFlags = MultiScalarField3D<plint>(nx, ny, nz);
        clotSolidFraction = MultiScalarField3D<T>(nx, ny, nz);
        clotSolidFractionPhys = MultiScalarField3D<T>(nx, ny, nz);

        //initialize the flags and SF to 0
        setToConstant<plint>(clotFlags,clotFlags.getBoundingBox(),0);
        setToConstant<T>(clotSolidFraction,clotSolidFraction.getBoundingBox(),0);
        setToConstant<T>(clotSolidFractionPhys,clotSolidFractionPhys.getBoundingBox(),0);


        // Initialization.
        initializeAtEquilibrium(*(lattice.get()), lattice->getBoundingBox(), PoiseuilleDensityAndZeroVelocity<T>(poiseuillePipe, DeltaP/C_P)/*1., Array<T,3>(0,0,simParam.u)*/);

        // Bounce-Back walls for the pipe.
        defineDynamics(*(lattice.get()), lattice->getBoundingBox(), new Pipe<T>(poiseuillePipe), new BounceBack<T, DESCRIPTOR>());


        pcout << "z range : " << fullDomain.z0() <<"; " <<fullDomain.z1() <<std::endl;
        // Concentration  of Fn and tPA based on the scale of the domain
        plint numVoxelsBeforeClot = floor(pow(pipeRadius, 2)*3.1416*floor(clotBeginZ*nz));
        V0antiFn = (numVoxelsBeforeClot+1) * dx*dx*dx;
        pcout << "Volume before clot (in m³) = " << V0antiFn << std::endl;
        // antiFnQtyPerCell = nb parts in 1 super part of anti-fibrin
        antiFnQtyPerCell = (double)antiFnConcentration*V0antiFn*avogadro/(double)particlesNum;      // but consider the whole volume before clot for initial antiFn

        // initialize the pointer that will contain the inlet density
        rhoWK_t_in = double(DeltaP/C_P * DESCRIPTOR<T>::invCs2 + 1);

        lattice.get()->initialize();
        
        // create and initialize the particle field
        particles = MultiParticleField3D<DenseParticleField3D<T,DESCRIPTOR>> (nx,ny,nz);

        return lattice;
    } // setupProblem
};  // SimulationParameters

// just a tag for voxels with fibrin totally consumed
enum flagNums{fibrinDestroyed=-1};


#endif
