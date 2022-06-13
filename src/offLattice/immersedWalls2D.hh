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

#ifndef IMMERSED_WALLS_2D_HH
#define IMMERSED_WALLS_2D_HH

#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "core/array.h"
#include "core/globalDefs.h"
#include "immersedWalls2D.h"

namespace plb {

/* ******** InamuroIteration2D ************************************ */

template <typename T, class VelFunction>
InamuroIteration2D<T, VelFunction>::InamuroIteration2D(
    VelFunction velFunction_, T tau_, bool incompressibleModel_) :
    velFunction(velFunction_), tau(tau_), incompressibleModel(incompressibleModel_)
{ }

template <typename T, class VelFunction>
void InamuroIteration2D<T, VelFunction>::processGenericBlocks(
    [[maybe_unused]] Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 3);
    ScalarField2D<T> *rhoBar = dynamic_cast<ScalarField2D<T> *>(blocks[0]);
    TensorField2D<T, 2> *j = dynamic_cast<TensorField2D<T, 2> *>(blocks[1]);
    AtomicContainerBlock2D *container = dynamic_cast<AtomicContainerBlock2D *>(blocks[2]);
    PLB_ASSERT(rhoBar);
    PLB_ASSERT(j);
    PLB_ASSERT(container);
    Dot2D location = rhoBar->getLocation();
    Dot2D ofsJ = computeRelativeDisplacement(*rhoBar, *j);
    ImmersedWallData2D<T> *wallData = dynamic_cast<ImmersedWallData2D<T> *>(container->getData());
    PLB_ASSERT(wallData);

    std::vector<Array<T, 2> > const &vertices = wallData->vertices;
    std::vector<T> const &areas = wallData->areas;
    PLB_ASSERT(vertices.size() == areas.size());
    std::vector<Array<T, 2> > deltaG(vertices.size());
    std::vector<Array<T, 2> > &g = wallData->g;
    PLB_ASSERT(vertices.size() == g.size());

    // In this iteration, the force is computed for every vertex.
    if (incompressibleModel) {
        for (pluint i = 0; i < vertices.size(); ++i) {
            Array<T, 2> const &vertex = vertices[i];
            Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
            const Array<plint, 2> xLim(
                (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<plint, 2> yLim(
                (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<T, 2> fracPos(util::frac(vertex[0]), util::frac(vertex[1]));
            Array<T, 2> averageJ;
            averageJ.resetToZero();
            // Use the weighting function to compute the average momentum
            // and the average density on the surface vertex.
            // x   x . x   x
            for (plint dx = xLim[0]; dx <= xLim[1]; ++dx) {
                for (plint dy = yLim[0]; dy <= yLim[1]; ++dy) {
                    Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                    Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
                    Array<T, 2> r((T)dx - fracPos[0], (T)dy - fracPos[1]);
                    T W = inamuroDeltaFunction2D<T>().W(r);
                    averageJ += W * nextJ;
                }
            }
            // averageJ += (T)0.5*g[i];
            Array<T, 2> wallVelocity = velFunction(vertex);
            deltaG[i] = areas[i] * (wallVelocity - averageJ);
            g[i] += deltaG[i];
        }
    } else {  // Compressible model.
        for (pluint i = 0; i < vertices.size(); ++i) {
            Array<T, 2> const &vertex = vertices[i];
            Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
            const Array<plint, 2> xLim(
                (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<plint, 2> yLim(
                (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<T, 2> fracPos(util::frac(vertex[0]), util::frac(vertex[1]));
            Array<T, 2> averageJ;
            averageJ.resetToZero();
            T averageRhoBar = T();
            // Use the weighting function to compute the average momentum
            // and the average density on the surface vertex.
            // x   x . x   x
            for (plint dx = xLim[0]; dx <= xLim[1]; ++dx) {
                for (plint dy = yLim[0]; dy <= yLim[1]; ++dy) {
                    Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                    T nextRhoBar = rhoBar->get(pos[0], pos[1]);
                    Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
                    Array<T, 2> r((T)dx - fracPos[0], (T)dy - fracPos[1]);
                    T W = inamuroDeltaFunction2D<T>().W(r);
                    averageJ += W * nextJ;
                    averageRhoBar += W * nextRhoBar;
                }
            }
            // averageJ += (T)0.5*g[i];
            Array<T, 2> wallVelocity = velFunction(vertex);
            deltaG[i] = areas[i] * ((averageRhoBar + (T)1.) * wallVelocity - averageJ);
            // g[i] += deltaG[i];
            g[i] += deltaG[i] / ((T)1.0 + averageRhoBar);
        }
    }

    // In this iteration, the force is applied from every vertex to the grid nodes.
    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 2> const &vertex = vertices[i];
        Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
        const Array<plint, 2> xLim(
            (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> yLim(
            (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<T, 2> fracPos(util::frac(vertex[0]), util::frac(vertex[1]));
        for (plint dx = xLim[0]; dx <= xLim[1]; ++dx) {
            for (plint dy = yLim[0]; dy <= yLim[1]; ++dy) {
                Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
                Array<T, 2> r((T)dx - fracPos[0], (T)dy - fracPos[1]);
                T W = inamuroDeltaFunction2D<T>().W(r);
                nextJ += tau * W * deltaG[i];
                j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y) = nextJ;
            }
        }
    }
}

template <typename T, class VelFunction>
InamuroIteration2D<T, VelFunction> *InamuroIteration2D<T, VelFunction>::clone() const
{
    return new InamuroIteration2D<T, VelFunction>(*this);
}

template <typename T, class VelFunction>
void InamuroIteration2D<T, VelFunction>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // RhoBar
    modified[1] = modif::staticVariables;  // J
    modified[2] = modif::nothing;          // Container Block with triangle data.
}

template <typename T, class VelFunction>
BlockDomain::DomainT InamuroIteration2D<T, VelFunction>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** IndexedInamuroIteration2D ************************************ */

template <typename T, class VelFunction>
IndexedInamuroIteration2D<T, VelFunction>::IndexedInamuroIteration2D(
    VelFunction velFunction_, T tau_, bool incompressibleModel_) :
    velFunction(velFunction_), tau(tau_), incompressibleModel(incompressibleModel_)
{ }

template <typename T, class VelFunction>
void IndexedInamuroIteration2D<T, VelFunction>::processGenericBlocks(
    [[maybe_unused]] Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 3);
    ScalarField2D<T> *rhoBar = dynamic_cast<ScalarField2D<T> *>(blocks[0]);
    TensorField2D<T, 2> *j = dynamic_cast<TensorField2D<T, 2> *>(blocks[1]);
    AtomicContainerBlock2D *container = dynamic_cast<AtomicContainerBlock2D *>(blocks[2]);
    PLB_ASSERT(rhoBar);
    PLB_ASSERT(j);
    PLB_ASSERT(container);
    Dot2D location = rhoBar->getLocation();
    Dot2D ofsJ = computeRelativeDisplacement(*rhoBar, *j);
    ImmersedWallData2D<T> *wallData = dynamic_cast<ImmersedWallData2D<T> *>(container->getData());
    PLB_ASSERT(wallData);

    std::vector<Array<T, 2> > const &vertices = wallData->vertices;
    std::vector<T> const &areas = wallData->areas;
    PLB_ASSERT(vertices.size() == areas.size());
    std::vector<Array<T, 2> > deltaG(vertices.size());
    std::vector<Array<T, 2> > &g = wallData->g;
    PLB_ASSERT(vertices.size() == g.size());
    std::vector<pluint> const &globalVertexIds = wallData->globalVertexIds;
    PLB_ASSERT(vertices.size() == globalVertexIds.size());

    if (incompressibleModel) {
        for (pluint i = 0; i < vertices.size(); ++i) {
            Array<T, 2> const &vertex = vertices[i];
            Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
            const Array<plint, 2> xLim(
                (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<plint, 2> yLim(
                (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<T, 2> fracPos(util::frac(vertex[0]), util::frac(vertex[1]));
            Array<T, 2> averageJ;
            averageJ.resetToZero();
            // x   x . x   x
            for (plint dx = xLim[0]; dx <= xLim[1]; ++dx) {
                for (plint dy = yLim[0]; dy <= yLim[1]; ++dy) {
                    Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                    Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
                    Array<T, 2> r((T)dx - fracPos[0], (T)dy - fracPos[1]);
                    T W = inamuroDeltaFunction2D<T>().W(r);
                    averageJ += W * nextJ;
                }
            }
            // averageJ += (T)0.5*g[i];
            Array<T, 2> wallVelocity = velFunction(globalVertexIds[i]);
            deltaG[i] = areas[i] * (wallVelocity - averageJ);
            g[i] += deltaG[i];
        }
    } else {  // Compressible model.
        for (pluint i = 0; i < vertices.size(); ++i) {
            Array<T, 2> const &vertex = vertices[i];
            Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
            const Array<plint, 2> xLim(
                (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<plint, 2> yLim(
                (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<T, 2> fracPos(util::frac(vertex[0]), util::frac(vertex[1]));
            Array<T, 2> averageJ;
            averageJ.resetToZero();
            T averageRhoBar = T();
            // x   x . x   x
            for (plint dx = xLim[0]; dx <= xLim[1]; ++dx) {
                for (plint dy = yLim[0]; dy <= yLim[1]; ++dy) {
                    Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                    T nextRhoBar = rhoBar->get(pos[0], pos[1]);
                    if (nextRhoBar == -1)
                        continue;
                    Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
                    Array<T, 2> r((T)dx - fracPos[0], (T)dy - fracPos[1]);
                    T W = inamuroDeltaFunction2D<T>().W(r);
                    averageJ += W * nextJ;
                    averageRhoBar += W * nextRhoBar;
                }
            }
            // averageJ += (T)0.5*g[i];
            Array<T, 2> wallVelocity = velFunction(globalVertexIds[i]);
            deltaG[i] = areas[i] * ((averageRhoBar + (T)1.) * wallVelocity - averageJ);
            // g[i] += deltaG[i];
            g[i] += deltaG[i] / ((T)1.0 + averageRhoBar);
        }
    }

    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 2> const &vertex = vertices[i];
        Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
        const Array<plint, 2> xLim(
            (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> yLim(
            (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<T, 2> fracPos(util::frac(vertex[0]), util::frac(vertex[1]));
        for (plint dx = xLim[0]; dx <= xLim[1]; ++dx) {
            for (plint dy = yLim[0]; dy <= yLim[1]; ++dy) {
                Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));

                T nextRhoBar = rhoBar->get(pos[0], pos[1]);
                if (nextRhoBar == -1)
                    continue;

                Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
                Array<T, 2> r((T)dx - fracPos[0], (T)dy - fracPos[1]);
                T W = inamuroDeltaFunction2D<T>().W(r);
                nextJ += tau * W * deltaG[i];
                j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y) = nextJ;
            }
        }
    }
}

template <typename T, class VelFunction>
IndexedInamuroIteration2D<T, VelFunction> *IndexedInamuroIteration2D<T, VelFunction>::clone() const
{
    return new IndexedInamuroIteration2D<T, VelFunction>(*this);
}

template <typename T, class VelFunction>
void IndexedInamuroIteration2D<T, VelFunction>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // RhoBar
    modified[1] = modif::staticVariables;  // J
    modified[2] = modif::nothing;          // Container Block with triangle data.
}

template <typename T, class VelFunction>
BlockDomain::DomainT IndexedInamuroIteration2D<T, VelFunction>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** InstantiateImmersedWallData2D ************************************ */

template <typename T>
InstantiateImmersedWallData2D<T>::InstantiateImmersedWallData2D(
    std::vector<Array<T, 2> > const &vertices_, std::vector<T> const &areas_,
    std::vector<Array<T, 2> > const &normals_) :
    vertices(vertices_), areas(areas_), normals(normals_)
{
    PLB_ASSERT(vertices.size() == areas.size());
    PLB_ASSERT(normals.size() == 0 || normals.size() == areas.size());
}

template <typename T>
void InstantiateImmersedWallData2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    AtomicContainerBlock2D *container = dynamic_cast<AtomicContainerBlock2D *>(blocks[0]);
    PLB_ASSERT(container);
    bool useNormals = normals.size() > 0;
    Dot2D location = container->getLocation();
    Array<T, 2> offset(location.x, location.y);
    ImmersedWallData2D<T> *wallData = new ImmersedWallData2D<T>;
    Box2D extendedEnvelope(domain.enlarge(2).shift(location.x, location.y));

    for (pluint i = 0; i < vertices.size(); ++i) {
        // Vertices which are close to the boundaries of the extendedEnvelope
        // are irrelevant, because they will act upon the bulk of the computational
        // domain through an Inamuro kernel, which at this distance is close to zero.
        // It is therefore OK, numerically speaking to exclude an epsilon-margin close
        // to these boundaries. Plus, it is required for technical reasons, because if
        // later on we pass across the boundaries of the extendedEnvelope because
        // of roundoff errors, the code will crash.
        static const T epsilon = 1.e-4;
        if (contained(vertices[i], extendedEnvelope, epsilon)) {
            wallData->vertices.push_back(vertices[i]);
            wallData->areas.push_back(areas[i]);
            if (useNormals) {
                wallData->normals.push_back(normals[i]);
            }
            wallData->g.push_back(Array<T, 2>((T)0., (T)0.));
            wallData->globalVertexIds.push_back(i);
        }
    }
    wallData->flags = std::vector<int>(wallData->vertices.size(), 0);
    container->setData(wallData);
}

template <typename T>
InstantiateImmersedWallData2D<T> *InstantiateImmersedWallData2D<T>::clone() const
{
    return new InstantiateImmersedWallData2D<T>(*this);
}

template <typename T>
void InstantiateImmersedWallData2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Container Block with triangle data.
}

template <typename T>
BlockDomain::DomainT InstantiateImmersedWallData2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

}  // namespace plb

#endif  // IMMERSED_WALLS_2D_HH