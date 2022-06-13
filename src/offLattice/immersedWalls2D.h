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

#ifndef IMMERSED_WALLS_2D_H
#define IMMERSED_WALLS_2D_H

#include <memory>

#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "core/array.h"
#include "core/globalDefs.h"
#include "dataProcessors/dataAnalysisWrapper2D.h"
#include "multiBlock/multiBlockGenerator2D.h"

namespace plb {

// TODO: Several functionals could be added based on the 3D implementations
// from immersedWalls3H.h/hh

/* ******** FavierDeltaFunction ************************************ */

/* ******** InamuroDeltaFunction2D ************************************ */

template <typename T>
class InamuroDeltaFunction2D {
public:
    InamuroDeltaFunction2D(int N_) : N(N_), dx(4. / (T)N), invDx(1. / dx)
    {
        sampleFunction();
    }
    T rawValue(T r) const
    {
        T rabs = std::fabs(r);
        T rsqr = r * r;
        if (rabs < 1.) {
            return 0.125 * (3. - 2. * rabs + std::sqrt(1. + 4. * rabs - 4. * rsqr));
        } else if (rabs < 2.) {
            return 0.125 * (5. - 2. * rabs - std::sqrt(-7. + 12. * rabs - 4. * rsqr));
        } else {
            return 0.;
        }
    }
    T w(T r) const
    {
        int position = (int)((r + 2.0) * invDx + 0.5);
        if (position <= 0) {
            return 0.;
        }
        if (position >= N) {
            return 0.;
        }
        return samples[position];
    }
    T W(Array<T, 2> const &r) const
    {
        return w(r[0]) * w(r[1]);
    }

private:
    void sampleFunction()
    {
        samples.resize(N + 1);
        for (int i = 0; i <= N; ++i) {
            samples[i] = rawValue(-2. + dx * i);
        }
    }

private:
    int N;
    T dx, invDx;
    std::vector<T> samples;
};

template <typename T>
inline InamuroDeltaFunction2D<T> const &inamuroDeltaFunction2D()
{
    static InamuroDeltaFunction2D<T> deltaFunction(1000);
    return deltaFunction;
}

/* ******** ImmersedWallData2D ************************************ */

template <typename T>
struct ImmersedWallData2D : public ContainerBlockData {
    std::vector<Array<T, 2> > vertices;
    std::vector<T> areas;
    std::vector<Array<T, 2> > normals;
    std::vector<Array<T, 2> > g;
    std::vector<int> flags;  // Flag for each vertex used to distinguish between vertices for
                             // conditional reduction operations.
    std::vector<pluint> globalVertexIds;
    virtual ImmersedWallData2D<T> *clone() const
    {
        return new ImmersedWallData2D<T>(*this);
    }
};

/* ******** SurfaceBlockData2D ************************************ */

/* ******** Utility functions ************************************ */

template <typename T>
inline bool closedOpenContained(Array<T, 2> const &x, Box2D const &box)
{
    return x[0] >= (box.x0 - 0.5) && x[0] < (box.x1 + 0.5) && x[1] >= (box.y0 - 0.5)
           && x[1] < (box.y1 + 0.5);
    // in order to count correctly the particles, a 0.5 must be added
}

/* ******** ReduceImmersedTorque2D ************************************ */
/* ******** ReduceImmersedForce3D ************************************ */
/* ******** ReduceImmersedArea3D ************************************ */

/* ******** InamuroIteration2D ************************************ */

template <typename T, class VelFunction>
class InamuroIteration2D : public BoxProcessingFunctional2D {
public:
    InamuroIteration2D(VelFunction velFunction_, T tau_, bool incompressibleModel_);
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual InamuroIteration2D<T, VelFunction> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    VelFunction velFunction;
    T tau;
    bool incompressibleModel;
};

template <typename T, class VelFunction>
void inamuroIteration(
    VelFunction velFunction, MultiScalarField2D<T> &rhoBar, MultiTensorField2D<T, 2> &j,
    MultiContainerBlock2D &container, T tau, bool incompressibleModel)
{
    std::vector<MultiBlock2D *> args;
    args.push_back(&rhoBar);
    args.push_back(&j);
    args.push_back(&container);
    applyProcessingFunctional(
        new InamuroIteration2D<T, VelFunction>(velFunction, tau, incompressibleModel),
        rhoBar.getBoundingBox(), args);
}

/* ******** IndexedInamuroIteration2D ************************************ */

// This is the same as InamuroIteration2D, with the difference that
// the VelFunction accepts as argument a global vertex index instead of
// a 2D position in space.
template <typename T, class VelFunction>
class IndexedInamuroIteration2D : public BoxProcessingFunctional2D {
public:
    IndexedInamuroIteration2D(VelFunction velFunction_, T tau_, bool incompressibleModel_);
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual IndexedInamuroIteration2D<T, VelFunction> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    VelFunction velFunction;
    T tau;
    bool incompressibleModel;
};

template <typename T, class VelFunction>
void indexedInamuroIteration(
    VelFunction velFunction, MultiScalarField2D<T> &rhoBar, MultiTensorField2D<T, 2> &j,
    MultiContainerBlock2D &container, T tau, bool incompressibleModel)
{
    std::vector<MultiBlock2D *> args;
    args.push_back(&rhoBar);
    args.push_back(&j);
    args.push_back(&container);
    applyProcessingFunctional(
        new IndexedInamuroIteration2D<T, VelFunction>(velFunction, tau, incompressibleModel),
        rhoBar.getBoundingBox(), args);
}

/* ******** IndexedImmersedBoundaryIteration2D ************************************ */
/* ******** ConstVelInamuroIteration2D ************************************ */
/* ******** ComputeImmersedBoundaryForce3D ************************************ */

/* ******** InstantiateImmersedWallData2D ************************************ */

template <typename T>
class InstantiateImmersedWallData2D : public BoxProcessingFunctional2D {
public:
    InstantiateImmersedWallData2D(
        std::vector<Array<T, 2> > const &vertices_, std::vector<T> const &areas_,
        std::vector<Array<T, 2> > const &normals_);
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual InstantiateImmersedWallData2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    std::vector<Array<T, 2> > const &vertices;
    std::vector<T> const &areas;
    std::vector<Array<T, 2> > const &normals;
};

template <typename T>
void instantiateImmersedWallData(
    std::vector<Array<T, 2> > const &vertices, std::vector<T> const &areas,
    MultiContainerBlock2D &container)
{
    static std::vector<Array<T, 2> > dummyNormals;
    std::vector<MultiBlock2D *> args;
    args.push_back(&container);
    applyProcessingFunctional(
        new InstantiateImmersedWallData2D<T>(vertices, areas, dummyNormals),
        container.getBoundingBox(), args);
}

template <typename T>
void instantiateImmersedWallData(
    std::vector<Array<T, 2> > const &vertices, std::vector<T> const &areas,
    std::vector<Array<T, 2> > const &normals, MultiContainerBlock2D &container)
{
    std::vector<MultiBlock2D *> args;
    args.push_back(&container);
    applyProcessingFunctional(
        new InstantiateImmersedWallData2D<T>(vertices, areas, normals), container.getBoundingBox(),
        args);
}

/* ******** InstantiateImmersedWallDataWithTagging2D ************************************ */
/* ******** InstantiateImmersedWallDataWithIndexedTagging2D ************************************ */
/* ******** InstantiateSurfaceBlockData2D ****************************** */
/* ******** SurfaceOnLattice2D **************************************** */
/* ******** SurfaceOnLattice2D_N **************************************** */
/* ******** ResetForceStatistics2D ************************************ */
/* ******** RecomputeImmersedForce2D ************************************ */
/* ******** OpenSurfaceImmersedForce2D ************************************ */

}  // namespace plb

#endif  // IMMERSED_WALLS_2D_H
