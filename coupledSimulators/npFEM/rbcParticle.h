///////////////////////////////////////////////////////////////////////////////
/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 * 
 * Contact for Palabos:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 * 
 * Contact for npFEM:
 * Christos Kotsalos
 * kotsaloscv@gmail.com
 * Computer Science Department
 * University of Geneva
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
///////////////////////////////////////////////////////////////////////////////
#pragma once

#include "palabos3D.h"
#include "palabos3D.hh"
#include "rbcGlobal.h"
#include "npfemConstants.h"

namespace plb {
namespace npfem {
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template<typename T, template<typename U> class Descriptor>
class SurfaceParticleWithAreaBase3D : public PointParticle3D<T,Descriptor> {
public:
    SurfaceParticleWithAreaBase3D();
    SurfaceParticleWithAreaBase3D(plint tag_, Array<T,3> const& position_, Array<T,3> const& velocity_, plint surfaceId_);

    virtual void serialize(HierarchicSerializer& serializer) const;
    virtual void unserialize(HierarchicUnserializer& unserializer);

    plint getSurfaceId() const { return surfaceId; }
    plint& getSurfaceId() { return surfaceId; }

    virtual T getArea(RawConnectedTriangleMesh<T>* wallMesh) const = 0;
private:
    plint surfaceId;
};

template<typename T, template<typename U> class Descriptor>
SurfaceParticleWithAreaBase3D<T,Descriptor>::SurfaceParticleWithAreaBase3D()
    : PointParticle3D<T,Descriptor>(),
      surfaceId(-1)
{ }

template<typename T, template<typename U> class Descriptor>
SurfaceParticleWithAreaBase3D<T,Descriptor>::SurfaceParticleWithAreaBase3D(plint tag_, Array<T,3> const& position_,
        Array<T,3> const& velocity_, plint surfaceId_)
    : PointParticle3D<T,Descriptor>(tag_, position_, velocity_),
      surfaceId(surfaceId_)
{ }

template<typename T, template<typename U> class Descriptor>
void SurfaceParticleWithAreaBase3D<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    PointParticle3D<T,Descriptor>::serialize(serializer);
    serializer.addValue(surfaceId);
}

template<typename T, template<typename U> class Descriptor>
void SurfaceParticleWithAreaBase3D<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    PointParticle3D<T,Descriptor>::unserialize(unserializer);
    unserializer.readValue(surfaceId);
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template <typename T, template <typename U> class Descriptor>
class bodyVertexParticle : public SurfaceParticleWithAreaBase3D<T, Descriptor>
{
public:
    bodyVertexParticle();
    // vertexID aka "tag", bodyID aka "surfaceID"
    // position in lattice units
    bodyVertexParticle(pluint vertexID_, Array<T, 3> const& position_, pluint bodyID_, Array<T, 3> const& velocity_);
    virtual void serialize(HierarchicSerializer& serializer) const;
    virtual void unserialize(HierarchicUnserializer& unserializer);
    virtual int getId() const { return id; }
    bodyVertexParticle<T, Descriptor>* clone() const
    {
        return new bodyVertexParticle<T, Descriptor>(*this);
    }
    virtual T getArea(RawConnectedTriangleMesh<T>* wallMesh) const;
private:
    static plint id;
};

template <typename T, template <typename U> class Descriptor>
plint bodyVertexParticle<T, Descriptor>::id = meta::registerGenericParticle3D<T,
    Descriptor, bodyVertexParticle<T, Descriptor>>("bodyVertex");

template <typename T, template <typename U> class Descriptor>
void bodyVertexParticle<T, Descriptor>::serialize(
    HierarchicSerializer& serializer) const
{
    SurfaceParticleWithAreaBase3D<T, Descriptor>::serialize(serializer);
}

template <typename T, template <typename U> class Descriptor>
void bodyVertexParticle<T, Descriptor>::unserialize(
    HierarchicUnserializer& unserializer)
{
    SurfaceParticleWithAreaBase3D<T, Descriptor>::unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
bodyVertexParticle<T, Descriptor>::bodyVertexParticle()
{
}

template <typename T, template <typename U> class Descriptor>
bodyVertexParticle<T, Descriptor>::bodyVertexParticle(pluint vertexID_,
    Array<T, 3> const& position_, pluint bodyID_, Array<T, 3> const& velocity_)
    : SurfaceParticleWithAreaBase3D<T, Descriptor>(vertexID_, position_, velocity_, bodyID_)
{
}

template <typename T, template <typename U> class Descriptor>
T bodyVertexParticle<T, Descriptor>::getArea(RawConnectedTriangleMesh<T>* wallMesh) const
{
    pluint bodyID = (pluint)this->getSurfaceId();

    if (bodyID == IDforWall)
    {
        return wallMesh->vertex(this->getTag())->area();
    }
    else
    {
        typename std::map<pluint, LocalMesh<T>*>::iterator it = LocalMeshes<T>().find(bodyID);
        return it->second->mesh.vertex(this->getTag())->area();
    }
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Construct local meshes for all bodies which have at least one vertex
// in the current processor or its envelope, and update all vertices
// for which a particle is found.
template <typename T, template <typename U> class Descriptor>
class ConstructLocalMeshesFromParticles : public BoxProcessingFunctional3D
{
public:
    ConstructLocalMeshesFromParticles(RawConnectedTriangleMesh<T>* rbcTemplate_, RawConnectedTriangleMesh<T>* pltTemplate_,
            std::map<pluint, pluint> bodyToType_,
            plint particleEnvelopeWidth_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks);
    virtual ConstructLocalMeshesFromParticles<T, Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    RawConnectedTriangleMesh<T> *rbcTemplate, *pltTemplate;
    std::map<pluint, pluint> bodyToType;
    plint particleEnvelopeWidth;
};

template <typename T, template <typename U> class Descriptor>
ConstructLocalMeshesFromParticles<T, Descriptor>::
    ConstructLocalMeshesFromParticles(RawConnectedTriangleMesh<T>* rbcTemplate_, RawConnectedTriangleMesh<T>* pltTemplate_,
            std::map<pluint, pluint> bodyToType_,
            plint particleEnvelopeWidth_)
    : rbcTemplate(rbcTemplate_)
    , pltTemplate(pltTemplate_)
    , bodyToType(bodyToType_)
    , particleEnvelopeWidth(particleEnvelopeWidth_)
{
}

template <typename T, template <typename U> class Descriptor>
void ConstructLocalMeshesFromParticles<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D*> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    ParticleField3D<T, Descriptor>& particleField = *dynamic_cast<ParticleField3D<T, Descriptor>*>(blocks[0]);

    // Bulk & Envelope
    std::vector<Particle3D<T, Descriptor>*> found;
    particleField.findParticles(domain.enlarge(particleEnvelopeWidth), found);

    // Loop over all local particles
    for (pluint iParticle = 0; iParticle < found.size(); ++iParticle)
    {
        bodyVertexParticle<T, Descriptor>* particle = dynamic_cast<bodyVertexParticle<T, Descriptor>*>(found[iParticle]);
        PLB_ASSERT(particle);

        pluint bodyID = particle->getSurfaceId();

        if (bodyID == IDforWall)
            continue;

        pluint vertexID = particle->getTag();

        typename std::map<pluint, LocalMesh<T>*>::iterator it = LocalMeshes<T>().find(bodyID);
        if (it == LocalMeshes<T>().end())
        {
            // if RBC
            if (bodyToType[bodyID] == 0)
                LocalMeshes<T>()[bodyID] = new LocalMesh<T>(*rbcTemplate, bodyID);
            // if PLT
            if (bodyToType[bodyID] == 1)
                LocalMeshes<T>()[bodyID] = new LocalMesh<T>(*pltTemplate, bodyID);

            it = LocalMeshes<T>().find(bodyID);
        }

        it->second->mesh.vertex(vertexID)->get() = particle->getPosition();
    }
}

template <typename T, template <typename U> class Descriptor>
ConstructLocalMeshesFromParticles<T, Descriptor>*
    ConstructLocalMeshesFromParticles<T, Descriptor>::clone() const
{
    return new ConstructLocalMeshesFromParticles<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ConstructLocalMeshesFromParticles<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template <typename T, template <typename U> class Descriptor>
class CollisionsForcesCombo : public BoxProcessingFunctional3D
{
public:
    CollisionsForcesCombo(T dx_, T omega_, T densityOffset_,
                          T Cf_, T Cp_, T Ca_, T collisions_threshold_rep_, T collisions_threshold_nonRep_,
                          std::vector<Array<T, 3>> wallVertexNormals_, std::map<pluint, pluint> bodyToType_,
                          bool CellPacking_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks);
    virtual CollisionsForcesCombo<T, Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T dx_p;
    T omega;
    T rho0;
    T Cf, Cp, Ca;
    T collisions_threshold_rep, collisions_threshold_nonRep;
    std::vector<Array<T, 3>> wallVertexNormals;
    std::map<pluint, pluint> bodyToType;
    bool CellPacking;
};

template <typename T, template <typename U> class Descriptor>
CollisionsForcesCombo<T, Descriptor>::CollisionsForcesCombo(T dx_, T omega_, T densityOffset_,
                                                            T Cf_, T Cp_, T Ca_, T collisions_threshold_rep_, T collisions_threshold_nonRep_,
                                                            std::vector<Array<T, 3>> wallVertexNormals_, std::map<pluint, pluint> bodyToType_,
                                                            bool CellPacking_)
    : dx_p(dx_)
    , omega(omega_)
    , rho0(densityOffset_)
    , Cf(Cf_)
    , Cp(Cp_)
    , Ca(Ca_)
    , collisions_threshold_rep(collisions_threshold_rep_)
    , collisions_threshold_nonRep(collisions_threshold_nonRep_)
    , wallVertexNormals(wallVertexNormals_)
    , bodyToType(bodyToType_)
    , CellPacking(CellPacking_)
{
}

template <typename T, template <typename U> class Descriptor>
void CollisionsForcesCombo<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D*> blocks)
{
    PLB_PRECONDITION(blocks.size() == 3);
    ParticleField3D<T, Descriptor>* particleField = dynamic_cast<ParticleField3D<T, Descriptor>*>(blocks[0]);
    ScalarField3D<T>* rhoBar = dynamic_cast<ScalarField3D<T>*>(blocks[1]);
    TensorField3D<T, SymmetricTensorImpl<T, 3>::n>* PiNeq = dynamic_cast<TensorField3D<T, SymmetricTensorImpl<T, 3>::n>*>(blocks[2]);

    PLB_ASSERT(particleField);
    PLB_ASSERT(rhoBar);
    PLB_ASSERT(PiNeq);

    Dot3D location = particleField->getLocation();

    Dot3D ofsRB = computeRelativeDisplacement(*particleField, *rhoBar);
    Dot3D ofsPN = computeRelativeDisplacement(*particleField, *PiNeq);

    // Get all particles in the bulk of the domain (no envelope)
    std::vector<Particle3D<T, Descriptor>*> found;
    // domain is in local coords
    particleField->findParticles(domain, found);

    // For mesh info
    plint isWrittenTag;

    // For collisions and periodicity
    std::map<pluint, Array<T, 3>> localCOMs; // COM: Center of mass
    std::map<pluint, pluint> localcnt;

    // Loop over all particles in the current bulk.
    for (pluint iParticle = 0; iParticle < found.size(); ++iParticle)
    {
        bodyVertexParticle<T, Descriptor>* particle = dynamic_cast<bodyVertexParticle<T, Descriptor>*>(found[iParticle]);
        PLB_ASSERT(particle);

        pluint bodyID = particle->getSurfaceId();

        if (bodyID == IDforWall)
            continue;

        // global coordinates
        Array<T, 3> vertex = particle->getPosition();

        if (localCOMs.find(bodyID) == localCOMs.end())
        {
            localCOMs[bodyID] = Array<T, 3>(0., 0., 0.);
            localcnt[bodyID] = 0;
        }

        // Periodicity & Collisions Handling
        localCOMs[bodyID] += vertex;
        localcnt[bodyID] += 1;
    }

    for (auto const& value : localCOMs)
    {
        pluint bodyID = value.first;

        // Center of Mass calculation (global coordinates - lattice units)
        localCOMs[bodyID] /= (T)localcnt[bodyID];

        typename std::map<pluint, LocalMesh<T>*>::iterator it = LocalMeshes<T>().find(bodyID);

        // Clean containers
        it->second->vertexIDs.clear();
        it->second->shearForces.clear();
        it->second->normals.clear();
        it->second->pressure.clear();
        it->second->area.clear();
        it->second->collisionNeighbors.clear();
        it->second->collisionNeighborsNormals.clear();
    }

    bool areaWeighted = true;
    for (pluint iParticle = 0; iParticle < found.size(); ++iParticle)
    {
        bodyVertexParticle<T, Descriptor>* particle = dynamic_cast<bodyVertexParticle<T, Descriptor>*>(found[iParticle]);
        PLB_ASSERT(particle);

        pluint bodyID = particle->getSurfaceId();

        // The Wall particles are represented in the same particle field as the body
        // particles, so we now exclude them based on their tag.
        if (bodyID == IDforWall)
            continue;

        pluint vertexId = particle->getTag();

        typename std::map<pluint, LocalMesh<T>*>::iterator it = LocalMeshes<T>().find(bodyID);

///////////////////////////////////////////////////////////////////////////////
        // Collision Handling First
///////////////////////////////////////////////////////////////////////////////

        T collisions_threshold = std::max(collisions_threshold_rep, collisions_threshold_nonRep);
        // The 2 comes from the support of phi_4 kernel
        plint offset_ = (plint)collisions_threshold - (plint)2;

        // If there are no particles of other bodies in the IBM kernel
        // (same kernel as the one that we compute the forces on the bodies)
        // then the fluid by itself can resolve the lubrication zone. Thus there 
        // is NO need for a collision framework, which resolves essentially 
        // this lubrication zone. This is why, the search for nearest neighbors
        // is done inside this volume. If on the contary we detect particles
        // of other bodies, then the collision framework has to take over 
        // because the fluid resolution is not enough any more.

        // global coordinates - lattice units
        Array<T, 3> vertex = particle->getPosition();

        // intPos in local coordinates
        Array<plint, 3> intPos((plint) vertex[0] - location.x,
                               (plint) vertex[1] - location.y,
                               (plint) vertex[2] - location.z);
        const Array<plint, 2> xLim((vertex[0] < (T) 0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> yLim((vertex[1] < (T) 0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> zLim((vertex[2] < (T) 0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));

        // Same kernel as in IBM (~)
        Array<plint, 3> corner1(intPos + Array<plint, 3>(xLim[0], yLim[0], zLim[0]) - offset_);
        Array<plint, 3> corner2(intPos + Array<plint, 3>(xLim[1], yLim[1], zLim[1]) + offset_);
        Box3D IBM_kernel(corner1[0], corner2[0], corner1[1], corner2[1], corner1[2], corner2[2]);

        std::vector<Particle3D<T, Descriptor>*> particlesInKernel;
        particleField->findParticles(IBM_kernel, particlesInKernel);

        // Calculate the nearest neighbor and then proceed to PIK
        // PIK: Particle In Kernel

        T minDist = std::numeric_limits<T>::max();
        T Dist;
        plint index_ = -1;
        bool rep = true;
        bool rep_ = true;
        for (pluint ipart = 0; ipart < particlesInKernel.size(); ++ipart)
        {
            bodyVertexParticle<T, Descriptor>* part = dynamic_cast<bodyVertexParticle<T, Descriptor>*>(particlesInKernel[ipart]);

            pluint candidateBodyId = part->getSurfaceId();

            // 1st filter - kernel "contaminated" by another body
            // "contamination" in the sense that we are going to interpolate
            // quantities from the interior of another body.
            if (candidateBodyId != bodyID)
            {
                Array<T, 3> candidateVertex = part->getPosition();

                Dist = sqrt((vertex[0] - candidateVertex[0])*(vertex[0] - candidateVertex[0]) + 
                            (vertex[1] - candidateVertex[1])*(vertex[1] - candidateVertex[1]) + 
                            (vertex[2] - candidateVertex[2])*(vertex[2] - candidateVertex[2]));

                if (candidateBodyId == IDforWall)
                {
                    // if RBC
                    if (bodyToType[bodyID] == 0)
                    {
                        rep = false;
                        collisions_threshold = collisions_threshold_nonRep;
                    }
                    // if PLT
                    else
                    {
                        rep = false;
                        collisions_threshold = collisions_threshold_nonRep;
                    }
                }
                else
                {
                    if (bodyToType[bodyID] == bodyToType[candidateBodyId])
                    {
                        // if RBC
                        if (bodyToType[bodyID] == 0)
                        {
                            rep = false;
                            collisions_threshold = collisions_threshold_nonRep;
                        }
                        // if PLT
                        else
                        {
                            rep = false;
                            collisions_threshold = collisions_threshold_nonRep;
                        }
                    }
                    else
                    {
                        rep = true;
                        collisions_threshold = collisions_threshold_rep;
                    }
                }
                
                if (Dist < minDist && Dist <= collisions_threshold)
                {
                    rep_ = rep;
                    minDist = Dist;
                    index_ = (plint)ipart;
                }
            }
        }

        bool particleInKernel = false;
        if (index_ > 0)
        {
            bodyVertexParticle<T, Descriptor>* part = dynamic_cast<bodyVertexParticle<T, Descriptor>*>(particlesInKernel[index_]);

            pluint candidateBodyId = part->getSurfaceId();
            Array<T, 3> candidateVertex = part->getPosition();
            pluint verId = part->getTag();

            Array<T, 3> normal;
            typename std::map<pluint, LocalMesh<T>*>::iterator it_tmp;
            if (candidateBodyId == IDforWall)
            {
                normal = wallVertexNormals[verId];
            }
            else
            {
                it_tmp = LocalMeshes<T>().find(candidateBodyId);
                normal = it_tmp->second->mesh.vertex(verId)->normal(areaWeighted);
            }

            // Encode repulsion in the normal
            if (rep_)
                normal *= 2.0;

            it->second->collisionNeighbors.push_back((candidateVertex - localCOMs[bodyID]) * dx_p);
            it->second->collisionNeighborsNormals.push_back(normal);

            collisions_threshold = rep_ ? collisions_threshold_rep : collisions_threshold_nonRep;
            offset_ = (plint)collisions_threshold - (plint)2;
            
            if (minDist <= 1.)
                particleInKernel = true;
        }

///////////////////////////////////////////////////////////////////////////////
        // ForceToLocalMesh
///////////////////////////////////////////////////////////////////////////////

        isWrittenTag = it->second->mesh.getVertexTag("isWritten");

        // Interpolate rhoBar and PiNeq on the particle position.
        T area = it->second->mesh.vertex(vertexId)->area();
        // "normal" points towards the fluid
        Array<T, 3> normal = it->second->mesh.vertex(vertexId)->normal(areaWeighted);

        const Array<T, 3> fracPos(util::frac(vertex[0]), util::frac(vertex[1]), util::frac(vertex[2]));

        // 1 : outside, 2 : inside
        T singleValRhoBar = 0.0; //averageRhoBar1 = 0.0, averageRhoBar2 = 0.0;
        Array<T, SymmetricTensorImpl<T, 3>::n> singleValPiNeq; //averagePiNeq1, averagePiNeq2;
        singleValPiNeq.resetToZero();
        //averagePiNeq1.resetToZero();
        //averagePiNeq2.resetToZero();
        //plint n1 = 0, n2 = 0;
        //T w1 = 0.0, w2 = 0.0;

        Array<T, 3> shearForce;
        T pressure;
        if (!particleInKernel && !CellPacking)
        {
            T maxDotProduct = -1.;
            for (plint dx = xLim[0]; dx <= xLim[1]; dx++)
            {
                for (plint dy = yLim[0]; dy <= yLim[1]; dy++)
                {
                    for (plint dz = zLim[0]; dz <= zLim[1]; dz++)
                    {
                        Array<plint, 3> pos(intPos + Array<plint, 3>(dx, dy, dz));
                        Array<T, 3> r((T)dx - fracPos[0], (T)dy - fracPos[1], (T)dz - fracPos[2]);
                        //T W = inamuroDeltaFunction<T>().W(r);

                        T dotProd = dot(r, normal);

                        // Simplest "interpolation"
                        // Keep the values of the most distant site
                        if (dotProd > maxDotProduct)
                        {
                            T nextRhoBar = rhoBar->get(pos[0] + ofsRB.x, pos[1] + ofsRB.y, pos[2] + ofsRB.z);
                            Array<T, SymmetricTensorImpl<T, 3>::n>& nextPiNeq = PiNeq->get(pos[0] + ofsPN.x, pos[1] + ofsPN.y, pos[2] + ofsPN.z);

                            maxDotProduct = dotProd;
                            singleValRhoBar = nextRhoBar;
                            singleValPiNeq = nextPiNeq;
                        }

                        // Inside/ Outside interpolation
                        /*
                        if (dotProd > 0) {
                        averageRhoBar1 += W * nextRhoBar;
                        averagePiNeq1 += W * nextPiNeq;
                        w1 += W;
                        ++n1;
                        } else {
                        averageRhoBar2 += W * nextRhoBar;
                        averagePiNeq2 += W * nextPiNeq;
                        w2 += W;
                        ++n2;
                        }
                        */

                        // For rigid bodies the whole kernel interpolation
                        // works fine. However for deformable bodies, the
                        // density variations cause huge instabilities.
                    }
                }
            }

            T singleValRho = Descriptor<T>::fullRho(singleValRhoBar);
            /*
            T averageRho1 = rho0;
            T averageRho2 = rho0;
            if (n1 > 0) {
            averageRhoBar1 /= w1;
            averageRho1 = Descriptor<T>::fullRho(averageRhoBar1);
            averagePiNeq1 /= w1;
            }
            if (n2 > 0) {
            averageRhoBar2 /= w2;
            averageRho2 = Descriptor<T>::fullRho(averageRhoBar2);
            averagePiNeq2 /= w2;
            }
            */

            Array<T, 3> singleValPi_n;
            SymmetricTensorImpl<T, 3>::matVectMult(singleValPiNeq, normal, singleValPi_n);
            /*
            Array<T, 3> averagePi_n1;
            SymmetricTensorImpl<T, 3>::matVectMult(averagePiNeq1, normal, averagePi_n1);
            // Important: on the reverse side, take the negative normal.
            Array<T, 3> averagePi_n2;
            SymmetricTensorImpl<T, 3>::matVectMult(averagePiNeq2, -normal, averagePi_n2);
            */

            // Pressure (Not yet the force since we will correct it later)
            pressure = -(singleValRho - rho0)*Descriptor<T>::cs2;
            
            // Shear Force
            shearForce = area*(Descriptor<T>::invRho(singleValRhoBar)*(omega / (T)2. - (T)1.)*singleValPi_n);
        }
        else
        {
            pressure = 0.;
            shearForce.resetToZero();
        }

        // Palabos: lb units
        // ShapeOp: physical units
        it->second->vertexIDs.push_back(vertexId);
        it->second->shearForces.push_back(Cf * shearForce);
        it->second->normals.push_back(normal);
        it->second->pressure.push_back(Cp * pressure);
        it->second->area.push_back(Ca * area);

        // We need to know which vertices belong to the bulk of the local
        // domain, for which the force has been computed, so we can, later
        // on, send them to the ShapeOp solver.
        it->second->mesh.vertex(vertexId)->setTag(isWrittenTag, 1);
    }
}

template <typename T, template <typename U> class Descriptor>
CollisionsForcesCombo<T, Descriptor>*
CollisionsForcesCombo<T, Descriptor>::clone() const
{
    return new CollisionsForcesCombo<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void CollisionsForcesCombo<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing; // Particle Field
    modified[1] = modif::nothing; // RhoBar
    modified[2] = modif::nothing; // PiNeq
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT CollisionsForcesCombo<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template <typename T, template <typename U> class Descriptor>
class LocalMeshToParticleVelocity3D : public BoxProcessingFunctional3D
{
public:
    LocalMeshToParticleVelocity3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks);
    virtual LocalMeshToParticleVelocity3D<T, Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
LocalMeshToParticleVelocity3D<T, Descriptor>::LocalMeshToParticleVelocity3D()
{
}

template <typename T, template <typename U> class Descriptor>
void LocalMeshToParticleVelocity3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D*> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    
    ParticleField3D<T, Descriptor>* particleField = dynamic_cast<ParticleField3D<T, Descriptor>*>(blocks[0]);
    
    PLB_ASSERT(particleField);

    std::vector<Particle3D<T, Descriptor>*> particles;
    particleField->findParticles(domain, particles);

    // Loop over all particles in the bulk.
    for (pluint iParticle = 0; iParticle < particles.size(); ++iParticle)
    {
        bodyVertexParticle<T, Descriptor>* particle = dynamic_cast<bodyVertexParticle<T, Descriptor>*>(particles[iParticle]);
        PLB_ASSERT(particle);

        pluint bodyID = particle->getSurfaceId();
        
        if (bodyID == IDforWall)
            continue;

        pluint vertexId = particle->getTag();

        typename std::map<pluint, LocalMesh<T>*>::const_iterator it = LocalMeshes<T>().find(bodyID);

        // The particle "velocity" is simply a position offset which will be
        // applied the next time "AdvanceParticles" is applied, to guarantee
        // that all particles move to the position "vertexPosition" computed by
        // ShapeOp in a way that is compatible with parallelism.
        particle->getVelocity() = it->second->velocities[vertexId];
    }
}

template <typename T, template <typename U> class Descriptor>
LocalMeshToParticleVelocity3D<T, Descriptor>*
    LocalMeshToParticleVelocity3D<T, Descriptor>::clone() const
{
    return new LocalMeshToParticleVelocity3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void LocalMeshToParticleVelocity3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::dynamicVariables; // Particle Field
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT
    LocalMeshToParticleVelocity3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template<typename T, template<typename U> class Descriptor>
class MultiDirectForcingImmersedBoundaryIteration3D : public BoxProcessingFunctional3D {
public:
    MultiDirectForcingImmersedBoundaryIteration3D(T tau_, bool incompressibleModel_,
                                                  RawConnectedTriangleMesh<T>* wallMesh_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks);
    virtual MultiDirectForcingImmersedBoundaryIteration3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T tau;
    bool incompressibleModel;
    RawConnectedTriangleMesh<T>* wallMesh;
};

template<typename T, template<typename U> class Descriptor>
MultiDirectForcingImmersedBoundaryIteration3D<T,Descriptor>::MultiDirectForcingImmersedBoundaryIteration3D(
        T tau_, bool incompressibleModel_,
        RawConnectedTriangleMesh<T>* wallMesh_)
    : tau(tau_),
      incompressibleModel(incompressibleModel_),
      wallMesh(wallMesh_)
{ }

template<typename T, template<typename U> class Descriptor>
void MultiDirectForcingImmersedBoundaryIteration3D<T,Descriptor>::processGenericBlocks(Box3D domain,
        std::vector<AtomicBlock3D*> blocks)
{
    PLB_PRECONDITION(blocks.size() == 3);

    ParticleField3D<T,Descriptor>* particleField = dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    ScalarField3D<T>* rhoBar = dynamic_cast<ScalarField3D<T>*>(blocks[1]);
    TensorField3D<T,3>* j = dynamic_cast<TensorField3D<T,3>*>(blocks[2]);
    PLB_ASSERT(particleField);
    PLB_ASSERT(rhoBar);
    PLB_ASSERT(j);

    Dot3D ofsRB = computeRelativeDisplacement(*particleField, *rhoBar);
    Dot3D ofsJ  = computeRelativeDisplacement(*particleField, *j);

    Dot3D location = particleField->getLocation();

    Box3D extendedDomain(domain.enlarge(2));

    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField->findParticles(extendedDomain, particles);
    std::vector<Array<T,3>> deltaG(particles.size());

    if (incompressibleModel) {
        for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
            SurfaceParticleWithAreaBase3D<T,Descriptor>* particle =
                dynamic_cast<SurfaceParticleWithAreaBase3D<T,Descriptor>*>(particles[iParticle]);
            PLB_ASSERT(particle);

            Array<T,3> particlePos = particle->getPosition();
            Array<plint,3> intPos((plint) particlePos[0] - location.x,
                                  (plint) particlePos[1] - location.y,
                                  (plint) particlePos[2] - location.z);
            const Array<plint,2> xLim((particlePos[0] < (T) 0 ? Array<plint,2>(-2, 1) : Array<plint,2>(-1, 2)));
            const Array<plint,2> yLim((particlePos[1] < (T) 0 ? Array<plint,2>(-2, 1) : Array<plint,2>(-1, 2)));
            const Array<plint,2> zLim((particlePos[2] < (T) 0 ? Array<plint,2>(-2, 1) : Array<plint,2>(-1, 2)));
            const Array<T,3> fracPos(util::frac(particlePos[0]), util::frac(particlePos[1]), util::frac(particlePos[2]));
            Array<T,3> averageJ; averageJ.resetToZero();
            // x   x . x   x
            for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
                for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                    for (plint dz = zLim[0]; dz <= zLim[1]; dz++) {
                        Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                        Array<T,3> nextJ = j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y, pos[2]+ofsJ.z);
                        Array<T,3> r((T)dx-fracPos[0],(T)dy-fracPos[1],(T)dz-fracPos[2]);
                        T W = inamuroDeltaFunction<T>().W(r);
                        averageJ += W*nextJ;
                    }
                }
            }

            Array<T,3> const& wallVelocity = particle->getVelocity();
            T area = particle->getArea(wallMesh);

            deltaG[iParticle] = area * (wallVelocity - averageJ);
        }
    } else { // Compressible model.
        for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
            SurfaceParticleWithAreaBase3D<T,Descriptor>* particle =
                dynamic_cast<SurfaceParticleWithAreaBase3D<T,Descriptor>*>(particles[iParticle]);
            PLB_ASSERT(particle);

            Array<T,3> particlePos = particle->getPosition();
            Array<plint,3> intPos((plint) particlePos[0] - location.x,
                                  (plint) particlePos[1] - location.y,
                                  (plint) particlePos[2] - location.z);
            const Array<plint,2> xLim((particlePos[0] < (T) 0 ? Array<plint,2>(-2, 1) : Array<plint,2>(-1, 2)));
            const Array<plint,2> yLim((particlePos[1] < (T) 0 ? Array<plint,2>(-2, 1) : Array<plint,2>(-1, 2)));
            const Array<plint,2> zLim((particlePos[2] < (T) 0 ? Array<plint,2>(-2, 1) : Array<plint,2>(-1, 2)));
            const Array<T,3> fracPos(util::frac(particlePos[0]), util::frac(particlePos[1]), util::frac(particlePos[2]));
            Array<T,3> averageJ; averageJ.resetToZero();
            T averageRhoBar = T();
            // x   x . x   x
            for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
                for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                    for (plint dz = zLim[0]; dz <= zLim[1]; dz++) {
                        Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                        T nextRhoBar = rhoBar->get(pos[0]+ofsRB.x, pos[1]+ofsRB.y, pos[2]+ofsRB.z);
                        Array<T,3> nextJ = j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y, pos[2]+ofsJ.z);
                        Array<T,3> r((T)dx-fracPos[0],(T)dy-fracPos[1],(T)dz-fracPos[2]);
                        T W = inamuroDeltaFunction<T>().W(r);
                        averageJ += W*nextJ;
                        averageRhoBar += W*nextRhoBar;
                    }
                }
            }

            Array<T,3> const& wallVelocity = particle->getVelocity();
            T area = particle->getArea(wallMesh);

            deltaG[iParticle] = area * (Descriptor<T>::fullRho(averageRhoBar) * wallVelocity - averageJ);
        }
    }
    
    for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
        SurfaceParticleWithAreaBase3D<T,Descriptor>* particle =
            dynamic_cast<SurfaceParticleWithAreaBase3D<T,Descriptor>*>(particles[iParticle]);
        PLB_ASSERT(particle);

        Array<T,3> particlePos = particle->getPosition();
        Array<plint,3> intPos((plint) particlePos[0] - location.x,
                              (plint) particlePos[1] - location.y,
                              (plint) particlePos[2] - location.z);
        const Array<plint,2> xLim((particlePos[0] < (T) 0 ? Array<plint,2>(-2, 1) : Array<plint,2>(-1, 2)));
        const Array<plint,2> yLim((particlePos[1] < (T) 0 ? Array<plint,2>(-2, 1) : Array<plint,2>(-1, 2)));
        const Array<plint,2> zLim((particlePos[2] < (T) 0 ? Array<plint,2>(-2, 1) : Array<plint,2>(-1, 2)));
        const Array<T,3> fracPos(util::frac(particlePos[0]), util::frac(particlePos[1]), util::frac(particlePos[2]));

        for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
            for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                for (plint dz = zLim[0]; dz <= zLim[1]; dz++) {
                    Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                    Array<T,3> nextJ = j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y, pos[2]+ofsJ.z);
                    Array<T,3> r((T)dx-fracPos[0],(T)dy-fracPos[1],(T)dz-fracPos[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    nextJ += tau*W*deltaG[iParticle];
                    j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y, pos[2]+ofsJ.z) = nextJ;
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
MultiDirectForcingImmersedBoundaryIteration3D<T,Descriptor>*
MultiDirectForcingImmersedBoundaryIteration3D<T,Descriptor>::clone() const
{
    return new MultiDirectForcingImmersedBoundaryIteration3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void MultiDirectForcingImmersedBoundaryIteration3D<T,Descriptor>::getTypeOfModification(
        std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;           // Particle Field
    modified[1] = modif::nothing;           // RhoBar
    modified[2] = modif::staticVariables;   // J
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT MultiDirectForcingImmersedBoundaryIteration3D<T,Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
} // namespace plb
}