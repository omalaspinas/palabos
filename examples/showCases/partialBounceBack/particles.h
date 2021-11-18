#ifndef PARTICLES // or some other unique name
#define PARTICLES

#include "palabos3D.h"
#include "palabos3D.hh"
#include "simParams.h"
#include <fstream>
#include <cstdio>
#include <cstdlib>      /* srand, rand */
#include <ctime>        /* time */
#include <algorithm>    // std::min

#define DESCRIPTOR descriptors::RhoBarJD3Q19Descriptor
typedef double T;
typedef DenseParticleField3D<T,DESCRIPTOR> ParticleFieldT;
extern SimulationParameters simParam;
#define PADDING 8

/// Execute the particle-fluid interaction step, and take into account the porous medium
/// to compute the particles' velocity
template<typename T, template<typename U> class Descriptor>
class PBBFluidToParticleCoupling3D : public BoxProcessingFunctional3D
{
public:
    /// Particle speed = scaling*fluid speed.
    PBBFluidToParticleCoupling3D(T scaling_) : scaling(scaling_) { }
    /// Arguments: [0] Particle-field; [1] Fluid.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks) {
        PLB_PRECONDITION( blocks.size()==3 );
        DenseParticleField3D<T,Descriptor>& particleField =
            *dynamic_cast<DenseParticleField3D<T,Descriptor>*>(blocks[0]);
        BlockLattice3D<T,DESCRIPTOR>& fluid =
            *dynamic_cast<BlockLattice3D<T,DESCRIPTOR>*>(blocks[1]);
        ScalarField3D<T>& sF = *dynamic_cast<ScalarField3D<T>*>(blocks[2]);
        particleField.fluidToParticleCoupling(domain, fluid, sF, scaling);
    }
    virtual PBBFluidToParticleCoupling3D<T,Descriptor>* clone() const {
        return new PBBFluidToParticleCoupling3D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::dynamicVariables; // Particle field.
        modified[1] = modif::nothing;  // Fluid.
        modified[2] = modif::nothing;  // Solid fraction.
    }
private:
    T scaling;
};

template<typename T, template<typename U> class Descriptor>
class ComputeTPAScalarField : public BoxProcessingFunctional3D
{   
public:
    ComputeTPAScalarField() {}

    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks)
    {
        PLB_PRECONDITION( blocks.size()==2 );
        ParticleField3D<T,Descriptor>& particleField = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
        ScalarField3D<plint>& numParticlefield = *dynamic_cast<ScalarField3D<plint>*>(blocks[1]);   // ScalarField that will contain the qty of tPA in the domain
        Dot3D offsetNumPart = computeRelativeDisplacement(particleField, numParticlefield);

        std::vector<Particle3D<T,Descriptor>*> particles;
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    particleField.findParticles(Box3D(iX,iX,iY,iY,iZ,iZ), particles);
                    for (plint iPart=0 ; iPart<particles.size() ; ++iPart)
                        numParticlefield.get(iX+offsetNumPart.x,iY+offsetNumPart.y,iZ+offsetNumPart.z) += particles[iPart]->getTag();
                }
            }
        }
    }

    virtual ComputeTPAScalarField<T,Descriptor>* clone() const
    {
        return new ComputeTPAScalarField<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;

    }
};


// Put in place the fluid-particle interactions,
// and the algorithms for the Verlet integration of the particle motion.
void setupCouplings(MultiBlockLattice3D<T,DESCRIPTOR>* lattice, MultiParticleField3D<ParticleFieldT>* particles, 
    MultiScalarField3D<T>* sF)
{
    // In the following the data processors for the equations of motion and the
    // interaction terms are manually added to the particle field. 

    // Prepare the arguments to be provided to the data processors: the particle field
    // for Verlet integration, and particles+fluid for the coupling terms.
    std::vector<MultiBlock3D*> particleFluidArg;
    particleFluidArg.push_back(particles);
    particleFluidArg.push_back(lattice);
    particleFluidArg.push_back(sF);

    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(particles);

    // Compute fluid->particle force (which includes friction), taking into account the porosity
    integrateProcessingFunctional (
            new PBBFluidToParticleCoupling3D<T,DESCRIPTOR>(1.),
            lattice->getBoundingBox(), particleFluidArg, 0 );

    // Calls the function advance() on each particle which, for Verlet particles,
    // updates the position according to v(t+0.5)
    integrateProcessingFunctional (
            new AdvanceParticlesFunctional3D<T,DESCRIPTOR>(-1.0),
            lattice->getBoundingBox(), particleArg, 0);

    // absorb particles exiting domain (= deleting them)
    Box3D absorbtionDomain(sF->getBoundingBox());
    absorbtionDomain.z0 = simParam.nz;
    absorbtionDomain.z1 = simParam.nz;

    integrateProcessingFunctional (
            new AbsorbParticlesFunctional3D<T,DESCRIPTOR>(),
            absorbtionDomain, particleArg, 0 );


    // Execute all data processors once to start the simulation off with well-defined initial values.
    particles->executeInternalProcessors();
}

void findNumParticlesPerCell(MultiParticleField3D<ParticleFieldT>* particles,plint iT)
{

    pcout << "# total particles = " << countParticles(*particles, particles->getBoundingBox()) << ", writing number of particles per cell VTK...\n";
    plint startPartIter = simParam.startParticleTime/simParam.dt;

    writeParticleVtk(*particles, std::string(simParam.outDir+"vtk/particles"+std::to_string((simParam.overwriteVtk ? ((iT>=startPartIter && iT<=startPartIter+simParam.outIter+1) ? 0 : 1) : iT))+".vtk").c_str(), simParam.dx);
}


// data processor to create a given number of particles randomly in a circular area, and setting their age
template<typename T, template<typename U> class Descriptor, class DomainFunctional>
class InjectParticlesAtInletSlice : public BoxProcessingFunctional3D
{   
public:
    InjectParticlesAtInletSlice(plint iter_, plint numParticles_, DomainFunctional functional_) : iter(iter_), numParticles(numParticles_), functional(functional_) {}
    InjectParticlesAtInletSlice(InjectParticlesAtInletSlice<T, Descriptor, DomainFunctional> const& rhs) : iter(rhs.iter), numParticles(rhs.numParticles), functional(rhs.functional) {}
    InjectParticlesAtInletSlice<T, Descriptor, DomainFunctional>& operator=(InjectParticlesAtInletSlice<T, Descriptor, DomainFunctional>const& rhs){
        InjectParticlesAtInletSlice<T, Descriptor, DomainFunctional>(rhs).swap(*this); return this; }
    void swap(InjectParticlesAtInletSlice<T, Descriptor, DomainFunctional>& rhs) {
        std::swap(iter, rhs.iter);
        std::swap(numParticles, rhs.numParticles);
        std::swap(functional, rhs.functional);
    }

    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks)
    {
        ParticleField3D<T,Descriptor>& particleField = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
        Dot3D offset = particleField.getLocation();
        PointParticle3D<T,DESCRIPTOR> *particleTemplate=0;
        random_device rd;       
        mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)

        uniform_real_distribution<double> uni_x(domain.x0, domain.x1); 
        uniform_real_distribution<double> uni_y(domain.y0,domain.y1);  
        uniform_real_distribution<double> uni_z(domain.z0,domain.z0+1.e-6);  
        double cx, cy, cz;
        plint radius = floor(0.5*(simParam.nx-2));
        plint numVoxelsBeforeClot = floor(pow(radius, 2)*3.1416*(floor(simParam.clotBeginZ*simParam.nz)+1)); //number of voxels inside the artery, before the clot. Useful to compute how many particles to inject in a given subdomain (parallelization).
        // plint numVoxelsHere = floor(pow(radius, 2)*3.1416*(std::min(abs(double(domain.z1+offset.z)),abs(floor(simParam.clotBeginZ*simParam.nz)-(domain.z0+offset.z)))));
        plint numVoxelsHere = floor(pow(radius, 2)*3.1416*double(domain.z1+offset.z));
        
        plint partsToInjectInThisSubdomain = ceil(numParticles*(double)numVoxelsHere/(double)numVoxelsBeforeClot);
        // pcout << "Parts to inject : " << double(numParticles)*(double)numVoxelsHere/(double)numVoxelsBeforeClot << " " << numParticles << " " << numVoxelsHere << " " << numVoxelsBeforeClot << endl;
        for (size_t i = 0 ; i < partsToInjectInThisSubdomain ; ++i)
        {
             
            // find a random position within tube
            
            do{
                cx = uni_x(rng);
                cy = uni_y(rng);
                cz = uni_z(rng);
            }while(!functional(Array<T,3>(cx+offset.x,cy+offset.y,cz+offset.z)));   // do this until we are in the domain functional


            particleTemplate = new PointParticle3D<T,DESCRIPTOR>(simParam.antiFnQtyPerCell, Array<T,3>(cx+offset.x,cy+offset.y,cz+offset.z), Array<T,3>(0.,0.,0.));
            particleField.addParticle(domain, particleTemplate->clone());
            delete particleTemplate; 
        }

    }

    virtual InjectParticlesAtInletSlice<T,Descriptor, DomainFunctional>* clone() const
    {
        return new InjectParticlesAtInletSlice<T,Descriptor, DomainFunctional>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::dynamicVariables;

    }
private:
    plint iter, numParticles;
    DomainFunctional functional;
};


void injectParticles(plint iter, plint numPart, MultiParticleField3D<ParticleFieldT>* particles)
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(particles);
    // pcout << "Inject particles" << std::endl;

    Box3D injectionDomain(2, simParam.nx-3, 2, simParam.ny-3, 0, floor(simParam.clotBeginZ*(simParam.nz))+1);
    plint minN = std::min(simParam.nx, simParam.ny);
    T pipeRadius = floor((T) 0.5 * (T) minN)-1;
    T pipeLength = (T) simParam.clotBeginZ*(simParam.nz-1);
    Array<T,3> pipeCenter(pipeRadius, pipeRadius, 0);
    Pipe<T> poiseuillePipe(pipeCenter, pipeRadius, pipeLength,false);

    applyProcessingFunctional (new InjectParticlesAtInletSlice<T,DESCRIPTOR,Pipe<T>> (iter, numPart, poiseuillePipe),
                     injectionDomain, particleArg);     

}


#endif

