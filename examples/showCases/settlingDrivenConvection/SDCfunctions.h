#ifndef SDCFUNCTIONS_H
#define SDCFUNCTIONS_H

#include <cmath>
#include <vector>

#include "palabos3D.h"
#include "palabos3D.hh"

#ifndef M_PI
#define M_PI 3.14159265358979323846L
#endif

namespace plb {

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// 1st order Upwind finite differences for advection-diffuion including Neumann boundary conditions
///////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
class AdvectionDiffusionNeumannFd3D : public BoxProcessingFunctional3D {
public:
    AdvectionDiffusionNeumannFd3D(
        T d_, bool upwind_, bool neumann_, plint nx_, plint ny_, plint nz_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual AdvectionDiffusionNeumannFd3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T d;
    bool upwind, neumann;
    plint nx, ny, nz;
};

template <typename T>
AdvectionDiffusionNeumannFd3D<T>::AdvectionDiffusionNeumannFd3D(
    T d_, bool upwind_, bool neumann_, plint nx_, plint ny_, plint nz_) :
    d(d_), upwind(upwind_), neumann(neumann_), nx(nx_), ny(ny_), nz(nz_)
{ }

template <typename T>
void AdvectionDiffusionNeumannFd3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_PRECONDITION(fields.size() == 5);
    ScalarField3D<T> *phi_t = dynamic_cast<ScalarField3D<T> *>(fields[0]);
    ScalarField3D<T> *phi_tp1 = dynamic_cast<ScalarField3D<T> *>(fields[1]);
    ScalarField3D<T> *result = dynamic_cast<ScalarField3D<T> *>(fields[2]);
    TensorField3D<T, 3> *uField = dynamic_cast<TensorField3D<T, 3> *>(fields[3]);
    ScalarField3D<T> *Q = dynamic_cast<ScalarField3D<T> *>(fields[4]);

    Dot3D ofs1 = computeRelativeDisplacement(*phi_t, *phi_tp1);
    Dot3D ofs2 = computeRelativeDisplacement(*phi_t, *result);
    Dot3D ofs3 = computeRelativeDisplacement(*phi_t, *uField);
    Dot3D ofs4 = computeRelativeDisplacement(*phi_t, *Q);

    Dot3D absoluteOffset = phi_tp1->getLocation();

    if (upwind) {
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint absoluteX = absoluteOffset.x + iX + ofs1.x;
                    plint absoluteY = absoluteOffset.y + iY + ofs1.y;
                    plint absoluteZ = absoluteOffset.z + iZ + ofs1.z;

                    T phiC = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    T phiE = phi_tp1->get(iX + 1 + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    T phiW = phi_tp1->get(iX - 1 + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    T phiN = phi_tp1->get(iX + ofs1.x, iY + 1 + ofs1.y, iZ + ofs1.z);
                    T phiS = phi_tp1->get(iX + ofs1.x, iY - 1 + ofs1.y, iZ + ofs1.z);
                    T phiT = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ + 1 + ofs1.z);
                    T phiB = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ - 1 + ofs1.z);

                    Array<T, 3> const &u = uField->get(iX + ofs3.x, iY + ofs3.y, iZ + ofs3.z);

                    Array<T, 3> adv;
                    T diffX, diffY, diffZ;

                    adv[0] =
                        (util::greaterThan(u[0], (T)0)
                             ? (phiC - phiW)
                             : (util::lessThan(u[0], (T)0) ? (phiE - phiC)
                                                           : (T)0.5 * (phiE - phiW)));
                    adv[1] =
                        (util::greaterThan(u[1], (T)0)
                             ? (phiC - phiS)
                             : (util::lessThan(u[1], (T)0) ? (phiN - phiC)
                                                           : (T)0.5 * (phiN - phiS)));
                    adv[2] =
                        (util::greaterThan(u[2], (T)0)
                             ? (phiC - phiB)
                             : (util::lessThan(u[2], (T)0) ? (phiT - phiC)
                                                           : (T)0.5 * (phiT - phiB)));

                    diffX = phiW + phiE - (T)2 * phiC;
                    diffY = phiS + phiN - (T)2 * phiC;
                    diffZ = phiT + phiB - (T)2 * phiC;

                    if (neumann) {
                        if (absoluteX == 0) {
                            adv[0] = 0;
                            diffX = (T)2 * phiE - (T)2 * phiC;
                        }

                        if (absoluteX == nx - 1) {
                            adv[0] = 0;
                            diffX = (T)2 * phiW - (T)2 * phiC;
                        }

                        if (absoluteY == 0) {
                            adv[1] = 0;
                            diffY = (T)2 * phiN - (T)2 * phiC;
                        }

                        if (absoluteY == ny - 1) {
                            adv[1] = 0;
                            diffY = (T)2 * phiS - (T)2 * phiC;
                        }

                        if (absoluteZ == 0) {
                            adv[2] = 0;
                            diffZ = (T)2 * phiT - (T)2 * phiC;
                        }

                        if (absoluteZ == nz - 1) {
                            adv[2] = 0;
                            diffZ = (T)2 * phiB - (T)2 * phiC;
                        }
                    }

                    T advection = u[0] * adv[0] + u[1] * adv[1] + u[2] * adv[2];

                    T diffusion = d * (diffX + diffY + diffZ);

                    result->get(iX + ofs2.x, iY + ofs2.y, iZ + ofs2.z) =
                        phi_t->get(iX, iY, iZ) + diffusion - advection
                        + Q->get(iX + ofs4.x, iY + ofs4.y, iZ + ofs4.z);
                }
            }
        }
    } else {
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint absoluteX = absoluteOffset.x + iX + ofs1.x;
                    plint absoluteY = absoluteOffset.y + iY + ofs1.y;
                    plint absoluteZ = absoluteOffset.z + iZ + ofs1.z;

                    T phiC = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    T phiE = phi_tp1->get(iX + 1 + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    T phiW = phi_tp1->get(iX - 1 + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    T phiN = phi_tp1->get(iX + ofs1.x, iY + 1 + ofs1.y, iZ + ofs1.z);
                    T phiS = phi_tp1->get(iX + ofs1.x, iY - 1 + ofs1.y, iZ + ofs1.z);
                    T phiT = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ + 1 + ofs1.z);
                    T phiB = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ - 1 + ofs1.z);

                    Array<T, 3> const &u = uField->get(iX + ofs3.x, iY + ofs3.y, iZ + ofs3.z);
                    Array<T, 3> adv;
                    T diffX, diffY, diffZ;

                    adv[0] = (T)0.5 * (phiE - phiW);
                    adv[1] = (T)0.5 * (phiN - phiS);
                    adv[2] = (T)0.5 * (phiT - phiB);

                    diffX = phiW + phiE - (T)2 * phiC;
                    diffY = phiS + phiN - (T)2 * phiC;
                    diffZ = phiT + phiB - (T)2 * phiC;

                    if (neumann) {
                        if (absoluteX == 0) {
                            adv[0] = 0;
                            diffX = (T)2 * phiE - (T)2 * phiC;
                        }

                        if (absoluteX == nx - 1) {
                            adv[0] = 0;
                            diffX = (T)2 * phiW - (T)2 * phiC;
                        }

                        if (absoluteY == 0) {
                            adv[1] = 0;
                            diffY = (T)2 * phiN - (T)2 * phiC;
                        }

                        if (absoluteY == ny - 1) {
                            adv[1] = 0;
                            diffY = (T)2 * phiS - (T)2 * phiC;
                        }

                        if (absoluteZ == 0) {
                            adv[2] = 0;
                            diffZ = (T)2 * phiT - (T)2 * phiC;
                        }

                        if (absoluteZ == nz - 1) {
                            adv[2] = 0;
                            diffZ = (T)2 * phiB - (T)2 * phiC;
                        }
                    }

                    T advection = u[0] * adv[0] + u[1] * adv[1] + u[2] * adv[2];

                    T diffusion = d * (diffX + diffY + diffZ);

                    result->get(iX + ofs2.x, iY + ofs2.y, iZ + ofs2.z) =
                        phi_t->get(iX, iY, iZ) + diffusion - advection
                        + Q->get(iX + ofs4.x, iY + ofs4.y, iZ + ofs4.z);
                }
            }
        }
    }
}

template <typename T>
AdvectionDiffusionNeumannFd3D<T> *AdvectionDiffusionNeumannFd3D<T>::clone() const
{
    return new AdvectionDiffusionNeumannFd3D<T>(*this);
}

template <typename T>
void AdvectionDiffusionNeumannFd3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // phi_t
    modified[1] = modif::nothing;          // phi_tp1
    modified[2] = modif::staticVariables;  // result
    modified[3] = modif::nothing;          // u
    modified[4] = modif::nothing;          // Q
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Buoyant force term
///////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T, template <typename U> class FluidDescriptor>
class ScalarBuoyancyTermProcessor3D : public BoxProcessingFunctional3D {
public:
    ScalarBuoyancyTermProcessor3D(
        T gravity_, T rho0_, T rhoP_, T TotalVolFrac_, T dt_, Array<T, FluidDescriptor<T>::d> dir_);

    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual ScalarBuoyancyTermProcessor3D<T, FluidDescriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T gravity, rho0, rhoP, TotalVolFrac, dt;
    Array<T, FluidDescriptor<T>::d> dir;
};

template <typename T, template <typename U> class FluidDescriptor>
ScalarBuoyancyTermProcessor3D<T, FluidDescriptor>::ScalarBuoyancyTermProcessor3D(
    T gravity_, T rho0_, T rhoP_, T TotalVolFrac_, T dt_, Array<T, FluidDescriptor<T>::d> dir_) :
    gravity(gravity_), rho0(rho0_), rhoP(rhoP_), TotalVolFrac(TotalVolFrac_), dt(dt_), dir(dir_)
{
    // We normalize the direction of the force vector.
    T normDir = std::sqrt(VectorTemplate<T, FluidDescriptor>::normSqr(dir));
    for (pluint iD = 0; iD < FluidDescriptor<T>::d; ++iD) {
        dir[iD] /= normDir;
    }
}

template <typename T, template <typename U> class FluidDescriptor>
void ScalarBuoyancyTermProcessor3D<T, FluidDescriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    typedef FluidDescriptor<T> D;
    enum { forceOffset = FluidDescriptor<T>::ExternalField::forceBeginsAt };

    PLB_PRECONDITION(fields.size() == 3);
    BlockLattice3D<T, FluidDescriptor> *fluid =
        dynamic_cast<BlockLattice3D<T, FluidDescriptor> *>(fields[0]);
    ScalarField3D<T> *volfracfield = dynamic_cast<ScalarField3D<T> *>(fields[1]);
    ScalarField3D<T> *densityfield = dynamic_cast<ScalarField3D<T> *>(fields[2]);

    Dot3D offset1 = computeRelativeDisplacement(*fluid, *volfracfield);
    Dot3D offset2 = computeRelativeDisplacement(*fluid, *densityfield);

    Array<T, D::d> gravOverrho0(
        gravity * dir[0] / rho0, gravity * dir[1] / rho0, gravity * dir[2] / rho0);

    T maxiT = 1.0 / dt;
    T iT = fluid->getTimeCounter().getTime();
    T gain = util::sinIncreasingFunction(iT, maxiT);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T localVolfrac = volfracfield->get(iX + offset1.x, iY + offset1.y, iZ + offset1.z);
                // Computation of the Boussinesq force
                T *force = fluid->get(iX, iY, iZ).getExternal(forceOffset);
                T dens = densityfield->get(iX + offset2.x, iY + offset2.y, iZ + offset2.z);
                // volfracfield is the order-0 moment of the advection-diffusion lattice.
                const T diffT = rhoP - rho0;
                for (pluint iD = 0; iD < D::d; ++iD) {
                    force[iD] = -gain * gravOverrho0[iD]
                                * (diffT * localVolfrac * TotalVolFrac
                                   + (dens - rho0) * (1 - localVolfrac * TotalVolFrac));
                }
            }
        }
    }
}

template <typename T, template <typename U> class FluidDescriptor>
ScalarBuoyancyTermProcessor3D<T, FluidDescriptor>
    *ScalarBuoyancyTermProcessor3D<T, FluidDescriptor>::clone() const
{
    return new ScalarBuoyancyTermProcessor3D<T, FluidDescriptor>(*this);
}

template <typename T, template <typename U> class FluidDescriptor>
void ScalarBuoyancyTermProcessor3D<T, FluidDescriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Settling velocity field
///////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
class ComputeSedimentationVelocity3D : public BoxProcessingFunctional3D {
public:
    ComputeSedimentationVelocity3D(T rhoP_, T Dp_, T convers_, T mu_, T g_);

    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual ComputeSedimentationVelocity3D<T> *clone() const;

private:
    T rhoP;
    T Dp;
    T convers;
    T mu;
    T g;
};

template <typename T>
ComputeSedimentationVelocity3D<T>::ComputeSedimentationVelocity3D(
    T rhoP_, T Dp_, T convers_, T mu_, T g_) :
    rhoP(rhoP_), Dp(Dp_), convers(convers_), mu(mu_), g(g_)
{ }

template <typename T>
void ComputeSedimentationVelocity3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    // typedef DensityDescriptor<T> D;
    PLB_PRECONDITION(atomicBlocks.size() == 3);
    ScalarField3D<T> *densityField = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[0]);
    ScalarField3D<T> *volfracField = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[1]);
    ScalarField3D<T> *v_sed = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[2]);

    Dot3D offset1 = computeRelativeDisplacement(*densityField, *volfracField);
    Dot3D offset2 = computeRelativeDisplacement(*densityField, *v_sed);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T dens = densityField->get(iX, iY, iZ);
                T localvolfrac = volfracField->get(iX + offset1.x, iY + offset1.y, iZ + offset1.z);
                T vel_sed;

                if (localvolfrac > 0)
                    vel_sed = -convers * (0.5 * Dp * Dp * g * (rhoP - dens)) / (9 * mu);
                else
                    vel_sed = 0;

                v_sed->get(iX + offset2.x, iY + offset2.y, iZ + offset2.z) = vel_sed;
            }
        }
    }
}

template <typename T>
ComputeSedimentationVelocity3D<T> *ComputeSedimentationVelocity3D<T>::clone() const
{
    return new ComputeSedimentationVelocity3D<T>(*this);
}

template <typename T>
void ComputeSedimentationVelocity3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Units conversion class
///////////////////////////////////////////////////////////////////////////////////////////////////////////

/// A useful class for the conversion between dimensionless and lattice units.
template <typename T, template <typename NSU> class nsDescriptor>
class RayleighTaylorFlowParam {
public:
    RayleighTaylorFlowParam(
        T Ri_, T Gr_, T uMax_, T uCar_, T resolution_, T lx_, T ly_, T lz_, T Di_ = T()) :
        Ri(Ri_),
        Gr(Gr_),
        uMax(uMax_),
        uCar(uCar_),
        resolution(resolution_),
        lx(lx_),
        ly(ly_),
        lz(lz_),
        Di(Di_)
    { }
    /// Rayleigh number
    T getRi() const
    {
        return Ri;
    }
    /// Grasshoff number
    T getGr() const
    {
        return Gr;
    }

    T getUcar() const
    {
        return uCar;
    }

    T getResolution() const
    {
        return resolution;
    }
    /// x-length in dimensionless units
    T getLx() const
    {
        return lx;
    }
    /// y-length in dimensionless units
    T getLy() const
    {
        return ly;
    }
    /// z-length in dimensionless units
    T getLz() const
    {
        return lz;
    }
    /// lattice spacing in dimensionless units
    T getDeltaX() const
    {
        return 1 / (T)resolution;
    }
    T getDeltaT() const
    {
        return getLatticeU() * getDeltaX() / getUcar();
    }
    /// conversion from dimensionless to lattice units for space coordinate
    plint nCell(T l) const
    {
        return (plint)(l / getDeltaX() + (T)0.5);
    }
    /// conversion from dimensionless to lattice units for time coordinuLbate
    plint nStep(T t) const
    {
        return (plint)(t / getDeltaT() + (T)0.5);
    }
    /// number of lattice cells in x-direction
    plint getNx() const
    {
        return nCell(lx) + 1;
    }
    /// number of lattice cells in y-direction
    plint getNy() const
    {
        return nCell(ly) + 1;
    }
    /// number of lattice cells in z-direction
    plint getNz() const
    {
        return nCell(lz) + 1;
    }
    /// velocity in lattice units (proportional to Mach number)
    T getLatticeU() const
    {
        return uMax;
    }
    /// Reynolds number
    T getRe() const
    {
        return std::sqrt(getGr() / getRi());
    }
    /// viscosity in lattice units
    T getLatticeNu() const
    {
        return getLatticeU() * (getNz() - 1) / getRe();
    }
    /// Diffusivity in lattice units
    T getLatticeKappa() const
    {
        return (getDeltaT() / (getDeltaX() * getDeltaX()));
    }
    /// viscosity in lattice units
    T getLatticeGravity() const
    {
        return getDeltaT() * getDeltaT() / getDeltaX();
    }
    /// relaxation time
    T getSolventTau() const
    {
        return nsDescriptor<T>::invCs2 * getLatticeNu() + (T)0.5;
    }
    /// relaxation frequency
    T getSolventOmega() const
    {
        return (T)1 / getSolventTau();
    }

private:
    T Ri, Gr, uMax, uCar, resolution, lx, ly, lz, Di;
};

template <typename T, template <typename NSU> class nsDescriptor>
void writeLogFile(
    RayleighTaylorFlowParam<T, nsDescriptor> const &parameters, std::string const &title)
{
    std::string fullName = global::directories().getLogOutDir() + "plbLog.dat";
    std::ofstream ofile(fullName.c_str());
    ofile << title << "\n\n";
    ofile << "Reynolds number:           Re=" << parameters.getRe() << "\n";
    ofile << "Richardson number:          Ri=" << parameters.getRi() << "\n";
    ofile << "Grasshoff number:            Gr=" << parameters.getGr() << "\n";
    ofile << "Kinematic viscosity:       Nu=" << parameters.getLatticeNu() << "\n";
    ofile << "Diffusivity:               Kappa=" << parameters.getLatticeKappa() << "\n";
    ofile << "Lattice resolution:         N=" << parameters.getResolution() << "\n";
    ofile << "Extent of the system:      lx=" << parameters.getLx() << "\n";
    ofile << "Extent of the system:      ly=" << parameters.getLy() << "\n";
    ofile << "Extent of the system:      lz=" << parameters.getLz() << "\n";
    ofile << "Grid spacing deltaX:       dx=" << parameters.getDeltaX() << "\n";
    ofile << "Time step deltaT:          dt=" << parameters.getDeltaT() << "\n";
    ofile << "Solvent omega:        omega_S=" << parameters.getSolventOmega() << "\n";
    ofile << "Caracteristic vel:       uLb=" << parameters.getLatticeU() << "\n";
    ofile << "Number of cells x:       Nx=" << parameters.getNx() << "\n";
    ofile << "Number of cells y:       Ny=" << parameters.getNy() << "\n";
    ofile << "Number of cells z:       Nz=" << parameters.getNz() << "\n";
}

}  // namespace plb

#endif
