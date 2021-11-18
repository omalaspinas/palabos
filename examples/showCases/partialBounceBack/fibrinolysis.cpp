/* This file is part of the Palabos library.
 *
 * Author: Remy Petkantchin
 * remy.petkantchin@unige.ch
 */

#include <fstream>
#include <cstdio>
#include <cstdlib>     /* srand, rand */
#include <ctime>       /* time */
#include "palabos3D.h"
#include "palabos3D.hh"
#include "simParams.h"
#include "clot.h"
#include "particles.h"

#define DESCRIPTOR descriptors::RhoBarJD3Q19Descriptor

typedef double T;
typedef DenseParticleField3D<T,DESCRIPTOR> ParticleFieldT;

using namespace plb;

SimulationParameters simParam;  // class that contains the parameters of the simulation to be ran

// function that writes .vti files readable with ParaView
void writeVtk(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, MultiScalarField3D<plint>& clotFlags
    , MultiScalarField3D<double>& clotSolidFraction, MultiScalarField3D<double>& clotSolidFractionPhys
    , plint iter)
{
    plint startPartIter = simParam.startParticleTime/simParam.dt;

	std::string fileName = createFileName("tmp/vtk/voldata", 
        (simParam.overwriteVtk ? ((iter>=startPartIter && iter<=startPartIter+simParam.outIter+1) ? 0 : 1) : iter), 8);
    VtkImageOutput3D<T> vtkOut(fileName, simParam.dx, simParam.physicalLocation);

    // z velocity in clot
    std::unique_ptr<MultiTensorField3D<double,3> > vel(computeVelocity(lattice));
    // porosity (=(1-n_s)) field in clot
    std::unique_ptr<MultiScalarField3D<double> > porosityField(subtract(1., clotSolidFraction));
    // macroscopic velocity in clot, based on Walsh ((1-n_s) * u)
    std::unique_ptr<MultiTensorField3D<double,3> > velMacroWalsh(multiply(*(porosityField.get()),*(vel.get())));

    // velocity
    vtkOut.writeData<3,float>(*(velMacroWalsh.get()), "velocity", simParam.dx / simParam.dt);
    vtkOut.writeData<3,float>(*(velMacroWalsh.get()), "velocityLB", 1.);

    // pressure
    T pressureScale = simParam.rho * (simParam.dx * simParam.dx) / (simParam.dt * simParam.dt) * DESCRIPTOR<T>::cs2;
    T pressureOffset = - pressureScale * simParam.rho_LB;

    vtkOut.writeData<float>(*computeDensity(lattice), "pressure", pressureScale, pressureOffset);
    vtkOut.writeData<float>(*computeDensity(lattice), "densityLB", 1., -simParam.rho_LB);
   

    std::string fileNameClot = createFileName("tmp/vtk/voldataClot", 
        (simParam.overwriteVtk ? ((iter>=startPartIter && iter<=startPartIter+simParam.outIter+1) ? 0 : 1) : iter), 8);
    VtkImageOutput3D<plint> clotOut(fileNameClot, simParam.dx, simParam.physicalLocation);
    clotOut.writeData<plint>(clotFlags, "clotFlags");

    std::string fileNameClotSF = createFileName("tmp/vtk/voldataSF", 
        (simParam.overwriteVtk ? ((iter>=startPartIter && iter<=startPartIter+simParam.outIter+1) ? 0 : 1) : iter), 8);
    VtkImageOutput3D<double> clotOutSF(fileNameClotSF, simParam.dx, simParam.physicalLocation);
    clotOutSF.writeData<double>(clotSolidFraction, "clotSolidFraction (gamma)");
    clotOutSF.writeData<double>(clotSolidFractionPhys, "clotSolidFractionPhys (ns*)");

}

double computeVSeepage(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, MultiScalarField3D<T>& sFField)
{
    Box3D clotDomain(lattice.getBoundingBox());
    clotDomain.z0 = simParam.clotBeginZ * (simParam.nz - 1)+2;
    clotDomain.z1 = simParam.clotEndZ * (simParam.nz - 1)-2;
    // z velocity in clot
    std::unique_ptr<MultiScalarField3D<T> > velZ(computeVelocityComponent(lattice, clotDomain, 2));
    // porosity (=(1-n_s)) field in clot
    std::unique_ptr<MultiScalarField3D<T> > porosityField(subtract(1., sFField, clotDomain));
    // macroscopic velocity in clot, based on Walsh ((1-n_s) * u)
    std::unique_ptr<MultiScalarField3D<T> > velMacroWalsh(multiply(*(porosityField.get()),*(velZ.get()), clotDomain));
    // average of Walsh macro velocity
    return computeAverage(*(velMacroWalsh.get()), clotDomain);
}

// returns true if flag value is fibrin
bool clotBoolMask(plint val)
{
    return (val!=flagNums::fibrinDestroyed && val!=0);
}

// TODO changed solid fraction here
/// Copy Unknown Populations Plus Impose Constant Pressure
template<typename T, template<typename U> class Descriptor, int direction, int orientation>
class FluidPressureInlet3D : public BoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    FluidPressureInlet3D(double rho_in_t_) : rho_in_t(rho_in_t_) {}

    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice) {
        std::vector<plint> const& unknownIndices = indexTemplates::subIndex<Descriptor<T>, direction, -orientation>();
    enum {
        normalX = direction==0 ? orientation : 0,
        normalY = direction==1 ? orientation : 0,
        normalZ = direction==2 ? orientation : 0
    };

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                Cell<T,Descriptor>& neighbor = lattice.get(iX-normalX, iY-normalY, iZ-normalZ);
                for (pluint fIndex=0; fIndex<unknownIndices.size(); ++fIndex) {
                    plint iPop = unknownIndices[fIndex];
                    cell[iPop] = neighbor[iPop];
                }
                T rhoBar;
                Array<T,3> j;
                cell.getDynamics().computeRhoBarJ(cell, rhoBar, j);
                Array<T,Descriptor<T>::q> oldFeq, newFeq;
                T jSqr = normSqr(j);
                cell.getDynamics().computeEquilibria(oldFeq, rhoBar, j, jSqr);
                cell.getDynamics().computeEquilibria(newFeq, rho_in_t-1, j, jSqr);
                for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                    cell[iPop] += newFeq[iPop]-oldFeq[iPop];
                }
            }
        }
    }
    }
    virtual FluidPressureInlet3D<T,Descriptor,direction,orientation>* clone() const {
        return new FluidPressureInlet3D<T,Descriptor,direction,orientation>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }

private:
    double rho_in_t;
};

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
class MyFluidPressureOutlet3D : public BoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice) {
      
        std::vector<plint> const& unknownIndices = indexTemplates::subIndex<Descriptor<T>, direction, -orientation>();
        enum {
            normalX = direction==0 ? orientation : 0,
            normalY = direction==1 ? orientation : 0,
            normalZ = direction==2 ? orientation : 0
        };

        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                    Cell<T,Descriptor>& neighbor = lattice.get(iX-normalX, iY-normalY, iZ-normalZ);
                    for (pluint fIndex=0; fIndex<unknownIndices.size(); ++fIndex) {
                        plint iPop = unknownIndices[fIndex];
                        cell[iPop] = neighbor[iPop];
                    }
                    T rhoBar;
                    Array<T,3> j;
                    cell.getDynamics().computeRhoBarJ(cell, rhoBar, j);
                    Array<T,Descriptor<T>::q> oldFeq, newFeq;
                    T jSqr = normSqr(j);
                    cell.getDynamics().computeEquilibria(oldFeq, rhoBar, j, jSqr);
                    cell.getDynamics().computeEquilibria(newFeq, T(), j, jSqr);
                    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                        cell[iPop] += newFeq[iPop]-oldFeq[iPop];
                    }
                }
            }
        }
    }
    virtual MyFluidPressureOutlet3D<T,Descriptor,direction,orientation>* clone() const {
        return new MyFluidPressureOutlet3D<T,Descriptor,direction,orientation>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
};


int main(int argc, char* argv[])
{    
    plbInit(&argc, &argv);

    pcout << std::endl;
    pcout << "Initialization." << std::endl;
    pcout << std::endl;

    global::directories().setOutputDir(simParam.outDir);

    global::IOpolicy().activateParallelIO(true);


    // Container for the lattice of particles
    MultiParticleField3D<ParticleFieldT> parts(1,1,1);
    // declare the flags for the bouce-back pixels of the clot
    MultiScalarField3D<plint> clotFlags(1,1,1);
    // Clot solid fraction initialization (for PBB dynamics) 0 : fluid node, 1 : solid node
    MultiScalarField3D<double> clotSolidFraction(1,1,1);    // corresponds to ns in Walsh : model sF (gamma)
    MultiScalarField3D<double> clotSolidFractionPhys(1,1,1);    // corresponds to ns* in Walsh : physical sF
    
    MultiBlockLattice3D<T,DESCRIPTOR>* lattice;
    lattice = simParam.setupProblem(parts, clotFlags, clotSolidFraction, clotSolidFractionPhys).release();
    MultiParticleField3D<ParticleFieldT>* particles = new MultiParticleField3D<ParticleFieldT>(parts);

    simParam.printSimulationParameters();


    // ------------------------ Clot zone -------------------------------- 
    // convert ratio to voxel position
    plint clotBeginZ = (T) simParam.clotBeginZ * (simParam.nz);
    plint clotEndZ = clotBeginZ + countClotZLength()-1;
    simParam.clotEndZ = double(clotEndZ)/double(simParam.nz);

    Clot clot(*lattice, clotBeginZ, clotEndZ);

    if(simParam.usesParticles)
    {
        setupCouplings(lattice, particles, &clotSolidFraction);
    }

    plb_ofstream statFile;	// file that will contain data statistics on the clot
    statFile.open(std::string(simParam.outDir+simParam.clotStatFilename).c_str(), fstream::out);

    plb_ofstream partStatFile;   // file that will contain data statistics on the particles    
    partStatFile.open(std::string(simParam.outDir+simParam.partStatFilename).c_str(), fstream::out);

    plb_ofstream flowFile;    // file that will contain flow information
    flowFile.open(std::string(simParam.outDir+"flow.dat").c_str(), fstream::out);

    // pre-computation for inlet/outlet
    Array<T,3> velInlet, velOutlet;
    T rhoInlet = 1.;
    T rhoOutlet = 1.;
    Array<T,3> inletCenter, outletCenter;
    T radiusIn;



    // declare outside to be able to use it in Zou-He pulsating pressure
    Pipe<T> poiseuillePipe;

    // Pipe
    plint minN = std::min(simParam.nx, simParam.ny);
    T pipeRadius = floor((T) 0.5 * (T) minN);
    radiusIn = pipeRadius;
    T pipeLength = (T) 1.*(simParam.nz);
    Array<T,3> pipeCenter(pipeRadius, pipeRadius, (T) 0.0 * (simParam.nz - 1));
    poiseuillePipe = Pipe<T>(pipeCenter, pipeRadius, pipeLength);
    inletCenter = pipeCenter;
    outletCenter = pipeCenter;
    outletCenter[2] = pipeLength-1;

    T convRhoToPressure = simParam.rho * (simParam.dx * simParam.dx) / (simParam.dt * simParam.dt) * DESCRIPTOR<T>::cs2;
    T convRhoToPressureAdim = simParam.rho_LB * DESCRIPTOR<T>::cs2;
    T offsetRhoToPressure = - convRhoToPressure * simParam.rho_LB;

    plint FnInitial = 0;
    MultiScalarField3D<plint> clotFlagsInitial(simParam.nx, simParam.ny, simParam.nz);
    MultiScalarField3D<double> L_x_t(simParam.nx, simParam.ny, simParam.nz);     // from Diamond and Anand "Inner clot diffusion and permeation during fibrinolysis, in µM
    MultiScalarField3D<double> R_x_t(simParam.nx, simParam.ny, simParam.nz);     // from Diamond and Anand "Inner clot diffusion and permeation during fibrinolysis, in nm
    
    // time at which clot should be done
    double tfClot = simParam.startClotTime+0.5*simParam.dt;

    // ------------ Starting iterations ------------------------

    pcout << "Starting the iterations" << std::endl;
    pcout << std::endl;

#ifdef PLB_REGRESSION
    simParam.maxTime = 3.0;
    simParam.startParticleTime = 1.0;
    simParam.startClotTime = 2.0;
    simParam.outIter = 1000;
#endif

    for (plint iter = simParam.iniIter; iter < simParam.maxTime/simParam.dt; iter++) {

		// generate clot
    	if (double(iter*simParam.dt) >= simParam.startClotTime && iter*simParam.dt < tfClot){
            clot.generateFibrin(clotFlags, clotSolidFraction, clotSolidFractionPhys, 1.);
                
            FnInitial += computeSum(clotFlags);
            clotFlagsInitial = clotFlags;
        }

    /******** SETTING BCs **********/

        T t = iter*simParam.dt;
        
        // constant inlet density
        if (t==0) {
            simParam.rhoWK_t_in = simParam.DeltaP/simParam.C_P * DESCRIPTOR<T>::invCs2 + 1;
        }

    /******* BCs set ******/

        if (simParam.usesParticles && iter*simParam.dt >= simParam.startParticleTime)
        {
            // inject particles
            if (iter*simParam.dt < simParam.stopParticleTime)
            {
        	    MultiScalarField3D<plint> tPAField(simParam.nx, simParam.ny, simParam.nz);
                std::vector<MultiBlock3D*> particlesTPAArg;
                particlesTPAArg.push_back(particles);
                particlesTPAArg.push_back(&tPAField);
                Box3D domainTPA(lattice->getBoundingBox());
                domainTPA.z1 = floor(simParam.clotBeginZ*simParam.nz)+1;     // stop the computation of tPA before the clot
                // first compute the tPA scalar field based on the particles
                applyProcessingFunctional(new ComputeTPAScalarField<T, DESCRIPTOR>(), domainTPA, particlesTPAArg);
                plint totalTPABeforeClot = computeIntSum(tPAField, domainTPA);
                
                while (double(totalTPABeforeClot)/double(simParam.V0antiFn*simParam.avogadro) < simParam.antiFnConcentration)
                {
                    tPAField = MultiScalarField3D<plint>(simParam.nx, simParam.ny, simParam.nz);
                    injectParticles(iter, 1, particles);
                    applyProcessingFunctional(new ComputeTPAScalarField<T, DESCRIPTOR>(), domainTPA, particlesTPAArg);
                    totalTPABeforeClot = computeIntSum(tPAField, domainTPA);
                }
            }

            // Handle interaction between particles and clot
            if (iter*simParam.dt >= simParam.startClotTime)
            {
                std::vector<MultiBlock3D*> particleFlagSFArg;
            	particleFlagSFArg.push_back(particles);
            	particleFlagSFArg.push_back(&clotFlags);
                particleFlagSFArg.push_back(&clotSolidFraction);
                particleFlagSFArg.push_back(&clotSolidFractionPhys);
                //// ------------------------- Interact ---------------------
                clot.interact(particleFlagSFArg,iter);  // apply functional to handle the interaction with the clot

            	if (iter % simParam.statIter == 0){		// count and write the clot qty of matter
            		statFile << count(clotFlags, clotBoolMask) << " " << computeSum(clotFlags) << " ";
                }
            }

        }

        // compute porosity(x,t) based on Diamond & Anand, with R_f(x,t), when there is lysis
        if (simParam.usesParticles)
        {

            std::vector<MultiBlock3D*> sFFlags;
            sFFlags.push_back(&clotSolidFraction);
            sFFlags.push_back(&R_x_t);
            sFFlags.push_back(&L_x_t);
            sFFlags.push_back(&clotFlagsInitial);
            sFFlags.push_back(&clotFlags);
            sFFlags.push_back(&clotSolidFractionPhys);
            applyProcessingFunctional(new updatesF_R_L_x_t<double, DESCRIPTOR>(), clotSolidFraction.getBoundingBox(), sFFlags);
            if (iter % simParam.statIter == 0 && iter*simParam.dt >= simParam.startParticleTime && iter*simParam.dt >= simParam.startClotTime){
                Box3D clotDomain(lattice->getBoundingBox());
                clotDomain.z0 = simParam.clotBeginZ * (simParam.nz);
                clotDomain.z1 = simParam.clotEndZ * (simParam.nz)-1;
                clotDomain.x0 += radiusIn/sqrt(2.);
                clotDomain.x1 -= radiusIn/sqrt(2.);
                clotDomain.y0 += radiusIn/sqrt(2.);
                clotDomain.y1 -= radiusIn/sqrt(2.);
                // output average solid fraction remaining
                statFile << computeAverage(clotSolidFractionPhys, clotDomain) << " " << iter*simParam.dt<<endl;
            }

            std::vector<MultiBlock3D*> latticeClotArg;
            latticeClotArg.push_back(lattice);
            latticeClotArg.push_back(&clotSolidFraction);
            applyProcessingFunctional(new ApplySolidFractionToRhoBar<double, DESCRIPTOR>(), clotSolidFraction.getBoundingBox(), latticeClotArg);

        }

        // ----------------------------------- OUTPUT -----------------------------------------

        // energy check
        T energy = computeAverageEnergy(*lattice) * simParam.rho * (simParam.dx * simParam.dx) / (simParam.dt * simParam.dt);
        if (!util::isFiniteNumber(energy)) {
            pcerr <<  "The simulation is unstable. Aborting ..." << std::endl;
            global::mpi().barrier();
            exit(1);
        }

#ifndef PLB_REGRESSION
        // write vtk
        if (simParam.writeVtk && iter % simParam.outIter == 0) {
            pcout << "Output (vtk) to disk at iteration: " << iter << ", t = " << iter * simParam.dt << std::endl;
            writeVtk(*lattice, clotFlags, clotSolidFraction, clotSolidFractionPhys, iter);
            pcout << std::endl;
        }
#endif
        // count particles and write in output file
        if (simParam.writeVtk && simParam.usesParticles && iter % simParam.outIter == 0){
            findNumParticlesPerCell(particles,iter);
        }

        lattice->executeInternalProcessors();

        // output flow.dat data
        if (iter % simParam.outIter == 0)
        {
            Cell<T,DESCRIPTOR>& cellInlet = lattice->get(round(inletCenter[0]), round(inletCenter[1]), round(inletCenter[2]));
            cellInlet.computeVelocity(velInlet); 
            rhoInlet = cellInlet.computeDensity();   
            Cell<T,DESCRIPTOR>& cellOutlet = lattice->get(round(outletCenter[0]), round(outletCenter[1]), round(outletCenter[2]));
            cellOutlet.computeVelocity(velOutlet);    
            rhoOutlet = cellOutlet.computeDensity(); 
            Cell<T,DESCRIPTOR>& cellInlet_1 = lattice->get(round(inletCenter[0]), round(inletCenter[1])+1, round(inletCenter[2]));
            Array<T,3> velInlet_1;
            cellInlet_1.computeVelocity(velInlet_1); 
            T rhoInlet_1 = cellInlet_1.computeDensity();
            Cell<T,DESCRIPTOR>& cellOutlet_1 = lattice->get(round(outletCenter[0]), round(outletCenter[1]), round(outletCenter[2])-1);
            Array<T,3> velOutlet_1;
            cellOutlet_1.computeVelocity(velOutlet_1);
            T rhoOutlet_1 = cellOutlet_1.computeDensity();
            Cell<T,DESCRIPTOR>& cellClotStart = lattice->get(round(outletCenter[0]), round(outletCenter[1]), round(clotBeginZ+2));
            Array<T,3> velClotStart;
            cellClotStart.computeVelocity(velClotStart);
            Cell<T,DESCRIPTOR>& cellClotEnd = lattice->get(round(outletCenter[0]), round(outletCenter[1]), round(clotEndZ-2));
            Array<T,3> velClotEnd;
            cellClotEnd.computeVelocity(velClotEnd);

            T gradPClot = double((lattice->get(round(outletCenter[0]), round(outletCenter[1]), round(clotBeginZ)).computeDensity()*convRhoToPressure+offsetRhoToPressure)
                                - (lattice->get(round(outletCenter[0]), round(outletCenter[1]), round(clotEndZ)).computeDensity()*convRhoToPressure+offsetRhoToPressure)) / double((clotEndZ-clotBeginZ)*simParam.dx);

            T gradPClotAdim = double((lattice->get(round(outletCenter[0]), round(outletCenter[1]), round(clotBeginZ)).computeDensity())
                                - (lattice->get(round(outletCenter[0]), round(outletCenter[1]), round(clotEndZ)).computeDensity()))
                                *convRhoToPressureAdim / double(clotEndZ-clotBeginZ);

            T permeability = computeVSeepage(*lattice, clotSolidFraction)*(simParam.dx/simParam.dt)*simParam.rho*simParam.nu/gradPClot;
            T permeAdim = computeVSeepage(*lattice, clotSolidFraction)*simParam.rho_LB*simParam.nu_LB/gradPClotAdim;

#ifndef PLB_REGRESSION
            flowFile << velInlet[0]*simParam.dx/simParam.dt << " " << velInlet[1]*simParam.dx/simParam.dt<< " " 
#else
            pcout << velInlet[0]*simParam.dx/simParam.dt << " " << velInlet[1]*simParam.dx/simParam.dt<< " " 
#endif
            << velInlet[2]*simParam.dx/simParam.dt<< " "<< rhoInlet*convRhoToPressure+offsetRhoToPressure << " "
            << velInlet_1[0]*simParam.dx/simParam.dt << " " << velInlet_1[1]*simParam.dx/simParam.dt<< " " 
            << velInlet_1[2]*simParam.dx/simParam.dt<< " "<< rhoInlet_1*convRhoToPressure+offsetRhoToPressure<<" "
            << velOutlet[0]*simParam.dx/simParam.dt << " " << velOutlet[1]*simParam.dx/simParam.dt<< " "
            << velOutlet[2]*simParam.dx/simParam.dt << " " << rhoOutlet*convRhoToPressure+offsetRhoToPressure << " "
            << velOutlet_1[0]*simParam.dx/simParam.dt << " " << velOutlet_1[1]*simParam.dx/simParam.dt<< " "
            << velOutlet_1[2]*simParam.dx/simParam.dt << " " << rhoOutlet_1*convRhoToPressure+offsetRhoToPressure << " "
            << computeVSeepage(*lattice, clotSolidFraction)*simParam.dx/simParam.dt << " " << permeability << " " << gradPClot << " "
            << computeVSeepage(*lattice, clotSolidFraction) << " " << permeAdim << " " << gradPClotAdim << " "
            << simParam.dx << " " << simParam.dt << " " << iter*simParam.dt << endl;
        }
        // ---------- COLLIDE AND STREAM ------------------
	    lattice->collideAndStream();
        
        // here, one could adapt the BCs based on the evolution of the lysis
        Box3D outlet(0,simParam.nx-1, 0,simParam.ny-1, simParam.nz-1,simParam.nz-1);
        applyProcessingFunctional(new MyFluidPressureOutlet3D<T,DESCRIPTOR,2,1>(), outlet, *lattice);
        Box3D inlet(0,simParam.nx-1, 0,simParam.ny-1, 0,0);
        applyProcessingFunctional(new FluidPressureInlet3D<T,DESCRIPTOR,2,-1>(simParam.rhoWK_t_in), inlet, *lattice);
        
        
    }

    statFile.close();
    partStatFile.close();
    flowFile.close();

    delete particles;
    delete lattice;


    return 0;
}

