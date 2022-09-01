#include <cstdlib>
#include <iostream>
#include <random>

#include "SDCfunctions.h"
#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;
using namespace std;

typedef double T;

#define NSDESCRIPTOR descriptors::ForcedD3Q19Descriptor
#define NSDYNAMICS   GuoExternalForceBGKdynamics

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialization of the particle volume fraction field.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename T, template <typename NSU> class nsDescriptor>
struct IniVolFracProcessor3D : public BoxProcessingFunctional3D_S<T> {
    IniVolFracProcessor3D(RayleighTaylorFlowParam<T, nsDescriptor> parameters_) :
        parameters(parameters_)
    { }
    virtual void process(Box3D domain, ScalarField3D<T> &volfracField)
    {
        Dot3D absoluteOffset = volfracField.getLocation();

        T nz = parameters.getNz();
        T up = 0.135;  // upper layer thickness
        T low = 0.25;  // lower layer thickness

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint absoluteZ = absoluteOffset.z + iZ;

                    T VolFrac;
                    if (absoluteZ < (nz - 1) - util::roundToInt(((up / (up + low)) * nz)))
                        VolFrac = 0.0;
                    else
                        VolFrac = 0.02 * volfracField.get(iX, iY, iZ) + 0.99;
                    if (absoluteZ > 0.95 * (nz - 1))
                        VolFrac = (0.02 * volfracField.get(iX, iY, iZ) + 0.99) * 0.5
                                  * (1 - erf((absoluteZ - 0.98 * (nz - 1)) / 6));
                    volfracField.get(iX, iY, iZ) = VolFrac;
                }
            }
        }
    }

    virtual IniVolFracProcessor3D<T, nsDescriptor> *clone() const
    {
        return new IniVolFracProcessor3D<T, nsDescriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }
    virtual BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulkAndEnvelope;
    }

private:
    RayleighTaylorFlowParam<T, nsDescriptor> parameters;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialization of the density field.
///////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename T, template <typename NSU> class nsDescriptor>
struct IniDensityProcessor3D : public BoxProcessingFunctional3D_S<T> {
    IniDensityProcessor3D(RayleighTaylorFlowParam<T, nsDescriptor> parameters_) :
        parameters(parameters_)
    { }
    virtual void process(Box3D domain, ScalarField3D<T> &densityField)
    {
        Dot3D absoluteOffset = densityField.getLocation();

        T nz = parameters.getNz();

        T up = 0.135;  // upper layer thickness
        T low = 0.25;  // lower layer thickness

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint absoluteZ = absoluteOffset.z + iZ;
                    T dens;
                    if (absoluteZ < (nz - 1) - util::roundToInt(((up / (up + low)) * nz)))
                        dens = 1008.4;
                    else
                        dens = 1000.;
                    densityField.get(iX, iY, iZ) = dens;
                }
            }
        }
    }
    virtual IniDensityProcessor3D<T, nsDescriptor> *clone() const
    {
        return new IniDensityProcessor3D<T, nsDescriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }
    virtual BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulkAndEnvelope;
    }

private:
    RayleighTaylorFlowParam<T, nsDescriptor> parameters;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fields and lattice initialization
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void ExpSetup(
    MultiBlockLattice3D<T, NSDESCRIPTOR> &nsLattice, MultiScalarField3D<T> &volfracField,
    MultiScalarField3D<T> &densityField, RayleighTaylorFlowParam<T, NSDESCRIPTOR> &parameters)
{
    initializeAtEquilibrium(
        nsLattice, nsLattice.getBoundingBox(), (T)1., Array<T, 3>((T)0., (T)0., (T)0.));

    T up = 0.135;  // upper layer thickness
    T low = 0.25;  // lower layer thickness
    T nx = parameters.getNx();
    T ny = parameters.getNy();
    T nz = parameters.getNz();
    T nupper = (nz - 1) - util::roundToInt(((up / (up + low)) * nz));
    Box3D upper(0, nx - 1, 0, ny - 1, nupper, nz - 1);
    uint32_t seed = 1;
    setToRandom(volfracField, upper, seed);
    applyProcessingFunctional(
        new IniVolFracProcessor3D<T, NSDESCRIPTOR>(parameters), volfracField.getBoundingBox(),
        volfracField);

    applyProcessingFunctional(
        new IniDensityProcessor3D<T, NSDESCRIPTOR>(parameters), densityField.getBoundingBox(),
        densityField);

    nsLattice.initialize();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// VTK outputs
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void writeVTK(
    MultiBlockLattice3D<T, NSDESCRIPTOR> &nsLattice, MultiScalarField3D<T> &volfracField,
    MultiScalarField3D<T> &densityField, RayleighTaylorFlowParam<T, NSDESCRIPTOR> const &parameters,
    plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();

    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<3, float>(*computeVelocity(nsLattice), "NSvelocity", dx / dt);
    vtkOut.writeData<float>(volfracField, "VolFrac", (T)1);
    vtkOut.writeData<float>(densityField, "dens", (T)1);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Image outputs
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void writeGif(
    MultiBlockLattice3D<T, NSDESCRIPTOR> &nsLattice, MultiScalarField3D<T> &volfracField,
    MultiScalarField3D<T> &densityField, int iT)
{
    const plint imSize = 600;
    const plint nx = nsLattice.getNx();
    const plint ny = nsLattice.getNy();
    const plint nz = nsLattice.getNz();
    Box3D slice((nx - 1) / 2, (nx - 1) / 2, 0, ny - 1, 0, nz - 1);
    ImageWriter<T> imageWriter("leeloo.map");
    imageWriter.writeScaledGif(
        createFileName("u", iT, 6), *computeVelocityNorm(nsLattice, slice), imSize, imSize);

    imageWriter.writeScaledGif(
        createFileName("VolFrac", iT, 6), *extractSubDomain(volfracField, slice), imSize, imSize);

    imageWriter.writeScaledGif(
        createFileName("Density", iT, 6), *extractSubDomain(densityField, slice), imSize, imSize);

    imageWriter.writeScaledGif(
        createFileName("Vorticity", iT, 6),
        *computeNorm(*computeVorticity(*computeVelocity(nsLattice)), slice), imSize, imSize);
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::timer("simTime").start();

    const T lx = 0.075;      // domain extent in the x-direction
    const T ly = 0.303;      // domain extent in the y-direction
    const T lz = 0.385;      // domain extent in the z-direction
    const T uCar = 0.01;     // physical characteristic velocity
    const T uMax = 8e-4;     // characteristic velocity in lattice units
    const T Gr = 1.48225e7;  // Grashof number
    const T Ri = 1.0;        // Richardson number
    const T Di = 3e-9;       // Physical diffusion coefficient

    const T Dp = 40e-6;                 // Equivalent particle diameter
    const T TotalVolFrac = 0.00198472;  // Particle volume fraction

    const plint resolution = 250;  // Mesh resolution

    global::directories().setOutputDir("./tmp/");

    RayleighTaylorFlowParam<T, NSDESCRIPTOR> parameters(
        Ri, Gr, uMax, uCar, resolution, lx, ly, lz, Di);

    T rho0 = 1000;                                // Fresh water density
    T rhoP = 2519.24;                             // Particle density
    T g = 9.81;                                   // Gravity acceleration
    T g_lb = g * parameters.getLatticeGravity();  // Gravity acceleration in lattice units
    T rho0f = 0.0;
    T eps = 1e-6;

    writeLogFile(parameters, "palabos.log");

    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    plint nz = parameters.getNz();

    T nsOmega = parameters.getSolventOmega();
    T convers = parameters.getDeltaT() / parameters.getDeltaX();
    T mu = 1e-3;  // Dynamic viscosity

    plint envelopeWidth = 2;
    SparseBlockStructure3D blockStructure(createRegularDistribution3D(nx, ny, nz));

    Dynamics<T, NSDESCRIPTOR> *nsdynamics = new NSDYNAMICS<T, NSDESCRIPTOR>(nsOmega);
    MultiBlockLattice3D<T, NSDESCRIPTOR> nsLattice(nx, ny, nz, nsdynamics->clone());
    defineDynamics(nsLattice, nsLattice.getBoundingBox(), nsdynamics->clone());
    delete nsdynamics;
    nsdynamics = 0;

    MultiScalarField3D<T> volfracField(
        MultiBlockManagement3D(
            blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(), envelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<T>(), T());

    MultiScalarField3D<T> densityField(volfracField);
    MultiScalarField3D<T> T_t(volfracField), T_tp1(volfracField), Q(volfracField),
        v_sedimentation(volfracField);
    MultiScalarField3D<T> volfracField_RK(volfracField), phi_1(volfracField), phi_2(volfracField),
        phi_n_adv(volfracField), phi_1_adv(volfracField), phi_2_adv(volfracField);
    MultiScalarField3D<T> D_t(densityField), D_tp1(densityField), Q_d(densityField);

    MultiTensorField3D<T, 3> velocity(
        MultiBlockManagement3D(
            blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(), envelopeWidth),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiTensorAccess<T, 3>());

    Box3D bottom(0, nx - 1, 0, ny - 1, 0, 0);
    Box3D top(0, nx - 1, 0, ny - 1, nz - 1, nz - 1);

    Box3D front(nx - 1, nx - 1, 0, ny - 1, 1, nz - 2);
    Box3D back(0, 0, 0, ny - 1, 1, nz - 2);

    Box3D left(0, nx - 1, 0, 0, 1, nz - 2);
    Box3D right(0, nx - 1, ny - 1, ny - 1, 1, nz - 2);

    nsLattice.toggleInternalStatistics(false);

    ExpSetup(nsLattice, volfracField, densityField, parameters);

    Dynamics<T, NSDESCRIPTOR> *nsbbDynamics = new BounceBack<T, NSDESCRIPTOR>(rho0f);
    defineDynamics(nsLattice, bottom, nsbbDynamics->clone());
    defineDynamics(nsLattice, top, nsbbDynamics->clone());
    defineDynamics(nsLattice, left, nsbbDynamics->clone());
    defineDynamics(nsLattice, right, nsbbDynamics->clone());
    defineDynamics(nsLattice, front, nsbbDynamics->clone());
    defineDynamics(nsLattice, back, nsbbDynamics->clone());
    delete nsbbDynamics;
    nsbbDynamics = 0;

    Array<T, NSDESCRIPTOR<T>::d> forceOrientation(T(), T(), (T)1);

    std::vector<MultiBlock3D *> args_f;
    args_f.push_back(&nsLattice);
    args_f.push_back(&volfracField);
    args_f.push_back(&densityField);

    integrateProcessingFunctional(
        new ScalarBuoyancyTermProcessor3D<T, NSDESCRIPTOR>(
            g_lb, rho0, rhoP, TotalVolFrac, parameters.getDeltaT(), forceOrientation),
        nsLattice.getBoundingBox(), args_f, 1);

    T tIni = global::timer("simTime").stop();
#ifndef PLB_REGRESSION
    pcout << "time elapsed for ExpSetup:" << tIni << endl;
#endif
    global::timer("simTime").start();

    plint evalTime = 2000;
    plint iT = 0;
    plint saveIter = 0.05 / parameters.getDeltaT();
    plint saveIterVtk = 0.05 / parameters.getDeltaT();
#ifndef PLB_REGRESSION
    plint maxT = 20 / parameters.getDeltaT();
#else
    plint maxT = 0.5 / parameters.getDeltaT();
#endif
    util::ValueTracer<T> converge((T)1, (T)100, 1.0e-3);
    pcout << "Max Number of iterations: " << maxT << endl;
    pcout << "Number of saving iterations: " << saveIter << endl;
    for (iT = 0; iT <= maxT; ++iT) {
        if (iT == (evalTime)) {
            T tEval = global::timer("simTime").stop();
            T remainTime = (tEval - tIni) / (T)evalTime * (T)maxT / (T)3600;
            global::timer("simTime").start();
            pcout << "Remaining " << (plint)remainTime << " hours, and ";
            pcout << (plint)((T)60 * (remainTime - (T)((plint)remainTime)) + 0.5) << " minutes."
                  << endl;
        }

        if (iT % saveIterVtk == 0) {
            pcout << iT * parameters.getDeltaT() << " : Writing VTK." << endl;
            writeVTK(nsLattice, volfracField, densityField, parameters, iT);
            pcout << "Average energy = "
                  << computeAverageEnergy(nsLattice, nsLattice.getBoundingBox());
            pcout << " ; Average volume fraction = " << TotalVolFrac * computeAverage(volfracField)
                  << endl;
        }

        if (iT % saveIter == 0) {
            pcout << iT * parameters.getDeltaT() << " : Writing gif." << endl;
            writeGif(nsLattice, volfracField, densityField, iT);
        }

        bool upwind = true;
        bool neumann = true;

        // Lattice Boltzmann iteration step.
        nsLattice.collideAndStream();

        computeVelocity(nsLattice, velocity, nsLattice.getBoundingBox());

        Actions3D cycle;
        plint D_t_ID = cycle.addBlock(D_t);
        plint D_tp1_ID = cycle.addBlock(D_tp1);
        plint densityField_ID = cycle.addBlock(densityField);
        plint velocity_ID = cycle.addBlock(velocity);
        plint Q_d_ID = cycle.addBlock(Q_d);
        plint T_t_ID = cycle.addBlock(T_t);
        plint T_tp1_ID = cycle.addBlock(T_tp1);
        plint volfracField_ID = cycle.addBlock(volfracField);
        plint volfracField_RK_ID = cycle.addBlock(volfracField_RK);
        plint v_sedimentation_ID = cycle.addBlock(v_sedimentation);
        plint Q_ID = cycle.addBlock(Q);
        plint phi_n_adv_ID = cycle.addBlock(phi_n_adv);
        plint phi_1_ID = cycle.addBlock(phi_1);
        plint phi_2_ID = cycle.addBlock(phi_2);
        plint phi_1_adv_ID = cycle.addBlock(phi_1_adv);
        plint phi_2_adv_ID = cycle.addBlock(phi_2_adv);

        ///////////////////////////////////////////////////////////////////////////////////////////////
        // Solve the advection diffusion for the density field with 1st order upwind finite
        // difference
        //////////////////////////////////////////////////////////////////////////////////////////////
        cycle.addProcessor(
            new CopyConvertScalarFunctional3D<T, T>(), densityField_ID, D_t_ID,
            densityField.getBoundingBox());
        cycle.addCommunication(D_t_ID, modif::staticVariables);
        cycle.addProcessor(
            new CopyConvertScalarFunctional3D<T, T>(), densityField_ID, D_tp1_ID,
            densityField.getBoundingBox());
        cycle.addCommunication(D_tp1_ID, modif::staticVariables);
        cycle.addProcessor(
            new AdvectionDiffusionNeumannFd3D<T>(
                Di * parameters.getLatticeKappa(), upwind, neumann, nx, ny, nz),
            D_t_ID, D_tp1_ID, densityField_ID, velocity_ID, Q_d_ID, densityField.getBoundingBox());
        cycle.addCommunication(densityField_ID, modif::staticVariables);

        ///////////////////////////////////////////////////////////////////////////////////////////////
        // Solve the advection-diffusion-sedimentation for the particle field with 3rd order WENO
        //////////////////////////////////////////////////////////////////////////////////////////////
        cycle.addProcessor(
            new ComputeSedimentationVelocity3D<T>(rhoP, Dp, convers, mu, g), densityField_ID,
            volfracField_ID, v_sedimentation_ID, v_sedimentation.getBoundingBox());
        cycle.addCommunication(v_sedimentation_ID, modif::staticVariables);
        cycle.addProcessor(
            new CopyConvertScalarFunctional3D<T, T>(), volfracField_ID, T_t_ID,
            volfracField.getBoundingBox());
        cycle.addCommunication(T_t_ID, modif::staticVariables);
        cycle.addProcessor(
            new CopyConvertScalarFunctional3D<T, T>(), volfracField_ID, T_tp1_ID,
            volfracField.getBoundingBox());
        cycle.addCommunication(T_tp1_ID, modif::staticVariables);
        std::vector<plint> arg1;
        arg1.push_back(T_t_ID);
        arg1.push_back(T_tp1_ID);
        arg1.push_back(phi_n_adv_ID);
        arg1.push_back(velocity_ID);
        arg1.push_back(Q_ID);
        arg1.push_back(v_sedimentation_ID);
        cycle.addProcessor(
            new WENO3<T>(Di * parameters.getLatticeKappa(), eps, neumann, nx, ny, nz), arg1,
            volfracField.getBoundingBox());
        cycle.addCommunication(phi_n_adv_ID, modif::staticVariables);
        cycle.addProcessor(
            new RK3Step1Functional3D<T>(), volfracField_ID, phi_n_adv_ID, phi_1_ID,
            volfracField.getBoundingBox());
        cycle.addCommunication(phi_1_ID, modif::staticVariables);
        cycle.addProcessor(
            new CopyConvertScalarFunctional3D<T, T>(), phi_1_ID, T_t_ID,
            volfracField.getBoundingBox());
        cycle.addCommunication(T_t_ID, modif::staticVariables);
        cycle.addProcessor(
            new CopyConvertScalarFunctional3D<T, T>(), phi_1_ID, T_tp1_ID,
            volfracField.getBoundingBox());
        cycle.addCommunication(T_tp1_ID, modif::staticVariables);
        std::vector<plint> arg2;
        arg2.push_back(T_t_ID);
        arg2.push_back(T_tp1_ID);
        arg2.push_back(phi_1_adv_ID);
        arg2.push_back(velocity_ID);
        arg2.push_back(Q_ID);
        arg2.push_back(v_sedimentation_ID);
        cycle.addProcessor(
            new WENO3<T>(Di * parameters.getLatticeKappa(), eps, neumann, nx, ny, nz), arg2,
            volfracField.getBoundingBox());
        cycle.addCommunication(phi_1_adv_ID, modif::staticVariables);
        cycle.addProcessor(
            new RK3Step2Functional3D<T>(), volfracField_ID, phi_1_ID, phi_1_adv_ID, phi_2_ID,
            volfracField.getBoundingBox());
        cycle.addCommunication(phi_1_ID, modif::staticVariables);
        cycle.addProcessor(
            new CopyConvertScalarFunctional3D<T, T>(), phi_2_ID, T_t_ID,
            volfracField.getBoundingBox());
        cycle.addCommunication(T_t_ID, modif::staticVariables);
        cycle.addProcessor(
            new CopyConvertScalarFunctional3D<T, T>(), phi_2_ID, T_tp1_ID,
            volfracField.getBoundingBox());
        cycle.addCommunication(T_tp1_ID, modif::staticVariables);
        std::vector<plint> arg3;
        arg3.push_back(T_t_ID);
        arg3.push_back(T_tp1_ID);
        arg3.push_back(phi_2_adv_ID);
        arg3.push_back(velocity_ID);
        arg3.push_back(Q_ID);
        arg3.push_back(v_sedimentation_ID);
        cycle.addProcessor(
            new WENO3<T>(Di * parameters.getLatticeKappa(), eps, neumann, nx, ny, nz), arg3,
            volfracField.getBoundingBox());
        cycle.addCommunication(phi_2_adv_ID, modif::staticVariables);
        cycle.addProcessor(
            new RK3Step3Functional3D<T>(), volfracField_ID, phi_2_ID, phi_2_adv_ID,
            volfracField_RK_ID, volfracField.getBoundingBox());
        cycle.addCommunication(volfracField_RK_ID, modif::staticVariables);
        cycle.addProcessor(
            new CopyConvertScalarFunctional3D<T, T>(), volfracField_RK_ID, volfracField_ID,
            volfracField.getBoundingBox());
        cycle.addCommunication(volfracField_ID, modif::staticVariables);

        cycle.execute();
    }
#ifndef PLB_REGRESSION
    saveBinaryBlock(nsLattice, "checkpoint_fluid.dat");
    saveBinaryBlock(volfracField, "checkpoint_vf.dat");
    saveBinaryBlock(densityField, "checkpoint_dens.dat");

    writeGif(nsLattice, volfracField, densityField, iT);

    T tEnd = global::timer("simTime").stop();

    T totalTime = tEnd - tIni;
    T nx100 = nsLattice.getNx() / (T)100;
    T ny100 = nsLattice.getNy() / (T)100;
    T nz100 = nsLattice.getNz() / (T)100;
    pcout << "N=" << resolution << endl;
    pcout << "number of processors: " << global::mpi().getSize() << endl;
    pcout << "simulation time: " << totalTime << endl;
    pcout << "total time: " << tEnd << endl;
    pcout << "total iterations: " << iT << endl;
    pcout << "Msus: " << nx100 * ny100 * nz100 * (T)iT / totalTime << endl;
#endif

    return 0;
}
