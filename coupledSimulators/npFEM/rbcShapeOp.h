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

#include <numeric>

#include "shapeOpWrapper.h"
#include "rbcFiltering.h"
#include "npfemConstants.h"

#define FILTERING

namespace plb {
namespace npfem {


template <typename T>
inline T getPeriodicPosition1D(const T pos, const pluint domainSize, bool& corrected)
{
    if (pos > domainSize)
    {
        corrected = true;
        return pos - ((plint)(std::abs(pos) / (T)domainSize)) * (T)domainSize;
    }
    else if (pos < 0.)
    {
        corrected = true;
        return pos + ((plint)(std::abs(pos) / (T)domainSize) + 1) * (T)domainSize;
    }
    
    corrected = false;
    return pos;
}

template <typename T>
inline void controlPositionPeriodic(Array<T, 3>& pos, const pluint lx, const pluint ly, const pluint lz, bool& corrected_x, bool& corrected_y, bool& corrected_z)
{
    pos[0] = getPeriodicPosition1D<T>(pos[0], lx, corrected_x);
    pos[1] = getPeriodicPosition1D<T>(pos[1], ly, corrected_y);
    pos[2] = getPeriodicPosition1D<T>(pos[2], lz, corrected_z);
}

template <typename T>
class ShapeOpBody
{
public:
	ShapeOpBody(){};

	void init(pluint bodyID_, RawConnectedTriangleMesh<T>* templateMesh_, std::vector<int> *mesh_graph = NULL)
    {
		int nb = solver.getPoints().cols();
		commForces = plb::npfem::Matrix3X::Zero(3, nb);
		shearForces = plb::npfem::Matrix3X::Zero(3, nb);
		normals = plb::npfem::Matrix3X::Zero(3, nb);
		numReceived = 0;
		pressure = plb::npfem::VectorX::Zero(nb);
		area = plb::npfem::VectorX::Zero(nb);
		bodyID = bodyID_ ;
		solver.bodyID_ = bodyID;
		templateMesh = templateMesh_;
		this->mesh_graph = mesh_graph;
		numVertices_onSurface_ = numVertices_onSurface();
	}

    ShapeOpBody(pluint bodyID_, ShapeOp_Solver &solver_, RawConnectedTriangleMesh<T>* templateMesh_)
        : bodyID(bodyID_)
        , solver(solver_)
        , templateMesh(templateMesh_)
        , numVertices_onSurface_(numVertices_onSurface())
        , numReceived(0)
        , commForces(3, numVertices())
        , shearForces(3, numVertices())
        , normals(3, numVertices())
        , pressure(numVertices())
        , area(numVertices())
    { solver.bodyID_ = bodyID; }

	void setId(pluint id)
    {
		bodyID = id;
		solver.bodyID_ = id;
	};

    pluint numVertices() const
    {
        return (pluint)solver.getPoints().cols();
    }

    pluint numVertices_onSurface() const
    {
        const std::vector<bool>& onSurfaceParticle = solver.get_onSurfaceParticle();

        pluint cnt = 0;
        for (pluint i = 0; i < numVertices(); ++i)
            if (onSurfaceParticle[i])
                ++cnt;

        PLB_ASSERT(cnt <= numVertices());
        return cnt;
    }

    pluint numTasksToSendVelocities() const
    {
        return (pluint)taskIDs.size();
    }

    pluint getID() const
    {
        return bodyID;
    }

    ShapeOp_Solver& getSolver()
    {
        return solver;
    }

    pluint getTask(pluint taskIndex) const
    {
        PLB_ASSERT(taskIndex < taskIDs.size());
        return taskIDs[taskIndex];
    }

    void resetBeforeReceiving()
    {
        numReceived = 0;
        
        commForces.setZero();
        shearForces.setZero();
        normals.setZero();
        pressure.setZero();
        area.setZero();
        
        vertexIDsPerProcessor.clear();
        taskIDs.clear();
        
        velocities.clear();

        // Collisions Handling
        collisionNeighbors.clear();
        collisionNeighborsNormals.clear();
    }

    bool communicationComplete()
    {
#ifdef ENABLE_LOGS
        plb::global::logfile_nonparallel("localProgress.log").flushEntry("bodyID: " + util::val2str(bodyID) +
        ", Vertices on Surface: " + util::val2str(numVertices_onSurface_) +
        ", Received Vertices: " + util::val2str(numReceived));
#endif // ENABLE_LOGS

        return (numReceived == numVertices_onSurface_);
    }

    // - Vertex-IDs and forces for a group of vertices
    // - The taskID of the processor from which this data came (not needed
    //   immediately, but is being stored to know which velocities
    //   to send to which processor afterwards).
    void receiveContainers(
        // Force Containers
        std::vector<pluint> const& vertexIDs,
        std::vector<Array<T, 3>> const& plbShearForces,
        std::vector<Array<T, 3>> const& plbNormals,
        std::vector<T> const& plbPressure, std::vector<T> const& plbArea,
        pluint taskID,
        // Collision Containers
        std::vector<Array<T, 3>> newNeighbors, std::vector<Array<T, 3>> newNeighborsNormals)
    {
        // Force Containers

        // ShapeOp side Center of Mass
        SOLocalMeshCenterOfMass = Array<T, 3>(0., 0., 0.);
        const plb::npfem::Matrix3X& points = solver.getPoints();

        for (pluint i = 0; i < (pluint)vertexIDs.size(); ++i) 
        {
            plb::npfem::Vector3 shearForce;
            shearForce[0] = plbShearForces[i][0];
            shearForce[1] = plbShearForces[i][1];
            shearForce[2] = plbShearForces[i][2];
            shearForces.col(vertexIDs[i]) = shearForce;

            plb::npfem::Vector3 normal;
            normal[0] = plbNormals[i][0];
            normal[1] = plbNormals[i][1];
            normal[2] = plbNormals[i][2];
            normals.col(vertexIDs[i]) = normal;

            pressure[vertexIDs[i]] = plbPressure[i];
            area[vertexIDs[i]] = plbArea[i];

            SOLocalMeshCenterOfMass += Array<T, 3>(
                                        points.col(vertexIDs[i])[0]
                                      , points.col(vertexIDs[i])[1]
                                      , points.col(vertexIDs[i])[2] );
        }
        
        // For collisions & periodicity
        SOLocalMeshCenterOfMass /= (T)vertexIDs.size();

        // Count the number of vertices on which the force has been
        // received, to make sure that at the end all vertices are up to date
        // (onSurface ones).
        numReceived += vertexIDs.size();

        // Store vertex-IDs and processor ID for later, to communicate back the
        // positions & velocities.
        vertexIDsPerProcessor.push_back(vertexIDs);
        taskIDs.push_back(taskID);

        ///////////////////////////////////////////////////////////////////////

        // Collision Containers

        if (newNeighbors.size() > 0)
        {
            // Transformation first (periodicity handling)
            for (pluint i = 0; i < (pluint)newNeighbors.size(); i++)
                newNeighbors[i] += SOLocalMeshCenterOfMass;

            collisionNeighbors.insert(collisionNeighbors.end(), newNeighbors.begin(), newNeighbors.end());
            collisionNeighborsNormals.insert(collisionNeighborsNormals.end(), newNeighborsNormals.begin(), newNeighborsNormals.end());
        }
    }

    void applyForces(bool CellPacking)
    {
        const std::vector<bool>& onSurfaceParticle = solver.get_onSurfaceParticle();

        // Pressure normalization
        T pressureIntegral = 0.;
        T totalArea = 0.;
        for (pluint i = 0; i < (pluint)pressure.size(); ++i)
        {
            pressureIntegral += pressure[i] * area[i];
            totalArea += area[i];
        }

        T p0 = pressureIntegral / totalArea;

        for (pluint i = 0; i < numVertices(); ++i)
        {
            if (onSurfaceParticle[i])
            {
                // Force Capping
                // For the specific IBM that we use, a force capping is rather useless.
                // Instead the velocity capping moderates the instabilites (if any and needed).
                commForces.col(i) = shearForces.col(i) + (pressure[i] - p0) * area[i] * normals.col(i).normalized();
            }
        }

#ifdef FILTERING
        forceFiltered = commForces;
        
        // On the GPU side it is done directly in the SO solver
#ifndef NPFEM_CUDA
        // Spatial Filtering
        plb::npfem::Matrix3X tmpForces = forceFiltered;

        for (pluint i = 0; i < (pluint)forceFiltered.cols(); ++i)
        {
            if (!onSurfaceParticle[i])
                continue;
				
            plb::npfem::VectorX vx(mesh_graph[i].size());
            plb::npfem::VectorX vy(mesh_graph[i].size());
            plb::npfem::VectorX vz(mesh_graph[i].size());

            int j = 0;
            for (auto v : mesh_graph[i])
            {
                vx[j] = forceFiltered.col(v)[0];
                vy[j] = forceFiltered.col(v)[1];
                vz[j] = forceFiltered.col(v)[2];
                ++j;
            }

            T mx, my, mz;
            median(vx, mx);
            median(vy, my);
            median(vz, mz);

            tmpForces.col(i) = plb::npfem::Vector3(mx, my, mz);
        }

        forceFiltered = tmpForces;
#endif // not GPU

        if (CellPacking)
            forceFiltered.setZero();

        editVertexForce(solver, forceFiltered);
#else
        if (CellPacking)
            commForces.setZero();

        editVertexForce(solver, commForces);
#endif // FILTERING
    }

    void applyCollidingNeighbors()
    {
        plb::npfem::Matrix3X neighbors(3, collisionNeighbors.size());
        plb::npfem::Matrix3X neighborsNormals(3, collisionNeighbors.size());
        for (pluint i = 0; i < (pluint)neighbors.cols(); ++i) {
            neighbors.col(i) = plb::npfem::Vector3(collisionNeighbors[i][0],
                                                collisionNeighbors[i][1],
                                                collisionNeighbors[i][2]);

            neighborsNormals.col(i) = plb::npfem::Vector3(collisionNeighborsNormals[i][0],
                                                       collisionNeighborsNormals[i][1],
                                                       collisionNeighborsNormals[i][2]);
        }

        // Delete previous (t-1) neighbors
        solver.verticesOfCollidingPieces_.clear();
        solver.normalsOfCollidingPieces_.clear();

        // Apply new ones (t)
        solver.verticesOfCollidingPieces_.push_back(neighbors);
        solver.normalsOfCollidingPieces_.push_back(neighborsNormals);
    }

    void sendVelocities()
    {
        // Since ShapeOp received something from Palabos,
        // it knows where to send back.
        for (pluint i = 0; i < numTasksToSendVelocities(); ++i)
        {
            // Prepare the container
            pluint taskID = taskIDs[i];
            std::vector<Array<T, 3>> newVelocities;
            fillVelocityContainers(i, newVelocities);
            velocities.push_back(newVelocities);

            // Send only if the data is not going to the current processor.
            if (newVelocities.size() > 0 && taskID != (pluint)global::mpi().getRank())
            {
                mpiRequests.push_back(new MPI_Request);
                MPI_Isend(&velocities.back()[0][0],
                    (int)(3 * velocities.back().size()), MPI_DOUBLE, (int)taskID,
                    shapeOpMPItagVelocities, MPI_COMM_WORLD, mpiRequests.back());
            }
        }
    }

    void fillVelocityContainers(pluint i, std::vector<Array<T, 3>>& velocities_)
    {
        PLB_ASSERT(i < numTasksToSendVelocities());
        velocities_.clear();

        plb::npfem::Matrix3X& vels = solver.getVelocities();

        for (pluint j = 0; j < vertexIDsPerProcessor[i].size(); ++j)
        {
            pluint vertexID = vertexIDsPerProcessor[i][j];
            Array<T, 3> vel(vels.col(vertexID)[0], vels.col(vertexID)[1], vels.col(vertexID)[2]);
            velocities_.push_back(vel);
        }
    }

    void completeRequests()
    {
#ifdef ENABLE_LOGS
        plb::global::logfile_nonparallel("localProgress.log").flushEntry("bodyID=" + util::val2str(bodyID) +
            ", Start MPI Completion (Velocities)");
#endif // ENABLE_LOGS

        MPI_Status status;
        for (pluint i = 0; i < mpiRequests.size(); ++i) {
            MPI_Wait(mpiRequests[i], &status);
            delete mpiRequests[i];
        }
        mpiRequests.clear();

#ifdef ENABLE_LOGS
        plb::global::logfile_nonparallel("localProgress.log").flushEntry("bodyID=" + util::val2str(bodyID) +
            ", Finish MPI Completion (Velocities)");
#endif // ENABLE_LOGS
    }

    /*
    // Create a Palabos Mesh from a GPU_solver.
    template <typename T>
	RawConnectedTriangleMesh<T> generateMeshGPU(const plb::npfem::MatrixXXCuda *points_t, 
												const std::vector<std::vector<int>> &triangles, const int nb_cells, double dx) {
		typedef Array<Array<T, 3>, 3> RawTriangle;

		// All points (OnSurface + interior)
		int npoint = points_t->cols();

		//std::cout << points_t->col(0) << std::endl;
		//std::cout << points_t->col(1) << std::endl;

		std::vector<Array<T, 3>> vertices(npoint*nb_cells);
		std::vector<Array<plint, 3>> plb_triangles(nb_cells*triangles.size());

		// Connectivity of the Surface only
		int k = 0;
		for (int c = 0; c < nb_cells; c++) {
			for (const auto &i : triangles) {
				plb_triangles[k] = Array<plint, 3>(i[0] + c*npoint, i[1] + c*npoint, i[2] + c*npoint);
				//std::cout << triangles[k][0]<< " " << triangles[k][1]<< " " << triangles[k][2] <<" |" << k <<"\n";
				k++;
			}
		}
	
		k = 0;
		for (int j = 0; j < nb_cells; j++) {
			for (int i = 0; i < points_t->cols(); ++i) {
				Array<T, 3> vertex(points_t->col(i)[3 * j] / dx, points_t->col(i)[3 * j + 1] / dx,
					points_t->col(i)[3 * j + 2] / dx);
				vertices[k++] = vertex;
			}
		}
		
        RawConnectedTriangleMesh<T> mesh(vertices, plb_triangles);
		return mesh;
	}
    */

	RawConnectedTriangleMesh<T> mergeMultipleMeshes(std::vector<plb::npfem::Solver*> &rbc, std::vector<plb::npfem::Solver*> &plt, T dx,
        T dump_RBCs_ratio = 1.0, T dump_PLTs_ratio = 1.0)
    {
		std::vector<Array<T, 3>> vertices;
		std::vector<Array<plint, 3>> triangles;


		for (int j = 0; j < (int)((T)rbc.size()*dump_RBCs_ratio); j++)
        {
			const plb::npfem::Matrix3X& points = rbc[j]->getPoints();

			for (pluint i = 0; i < (pluint)points.cols(); ++i)
            {
				Array<T, 3> vertex(points.col(i)[0], points.col(i)[1], points.col(i)[2]);
				vertices.push_back(vertex / dx);
			}
        }

		for (int j = 0; j < (int)((T)plt.size()*dump_PLTs_ratio); j++)
        {
			const plb::npfem::Matrix3X& points = plt[j]->getPoints();

			for (pluint i = 0; i < (pluint)points.cols(); ++i)
            {
				Array<T, 3> vertex(points.col(i)[0], points.col(i)[1], points.col(i)[2]);
				vertices.push_back(vertex / dx);
			}
		}

		int last_p = 0;
		// Connectivity of the Surface only
		for (int j = 0; j < (int)((T)rbc.size()*dump_RBCs_ratio); j++)
        {
			for (auto i = rbc[j]->getConnectivityList().begin(); i < rbc[j]->getConnectivityList().end(); ++i)
				triangles.push_back(Array<plint, 3>((*i)[0] + last_p, (*i)[1] + last_p, (*i)[2] + last_p));

			last_p += rbc[j]->getPoints().cols();
		}

		for (int j = 0; j < (int)((T)plt.size()*dump_PLTs_ratio); j++)
        {
			for (auto i = plt[j]->getConnectivityList().begin(); i < plt[j]->getConnectivityList().end(); ++i)
				triangles.push_back(Array<plint, 3>((*i)[0] + last_p, (*i)[1] + last_p, (*i)[2] + last_p));
			
			last_p += plt[j]->getPoints().cols();
		}

		RawConnectedTriangleMesh<T> mesh(vertices, triangles);
        /*
		mesh.registerVertexProperty("shearForce_0");
		mesh.registerVertexProperty("shearForce_1");
		mesh.registerVertexProperty("shearForce_2");
		mesh.registerVertexProperty("normal_0");
		mesh.registerVertexProperty("normal_1");
		mesh.registerVertexProperty("normal_2");
		mesh.registerVertexProperty("pressure");
		mesh.registerVertexProperty("area");
        */
		mesh.registerVertexTag("isWritten");
		
        return mesh;
	}

    RawConnectedTriangleMesh<T> generatePalabosMesh(T dx, T dt, T rho, pluint nx,
        pluint ny, pluint nz, pluint iT, pluint dt_ShapeOp)
    {
        typedef typename RawConnectedTriangleMesh<T>::PVertexIterator PVertexIterator;
        typedef typename RawConnectedTriangleMesh<T>::PVertex PVertex;

        // Whatever we get from ShapeOp is in physical units
        plb::npfem::Matrix3X points = solver.getPoints(); // Do a copy
        //const plb::npfem::Matrix3X& Palabos_Forces = solver.get_Palabos_Forces();
        const plb::npfem::Matrix3X& vels = solver.getVelocities();

        // ShapeOp has already moved the points at the next
        // ShapeOp timestep (which is in principal different
        // from Palabos timestep). However in the visualization,
        // we do not want to see steps but something continuous
        // Below is a small correction for the visualization only
        if (dt_ShapeOp != 1) {
            T delta = solver.getTimeStep();
            plb::npfem::Matrix3X oldPoints = points - delta * vels;

            T delta_Palabos = delta / ((T)dt_ShapeOp);
            pluint proceedBy = iT % dt_ShapeOp;

            // TODO: Maybe proceedBy+1 instead of proceedBy
            points = oldPoints + vels * delta_Palabos * (T)(proceedBy+1);
        }

        RawConnectedTriangleMesh<T> mesh(*templateMesh);

        pluint i = 0;
        PVertexIterator vertexIterator(mesh.vertexIterator());
        while (!vertexIterator->end())
        {
            PVertex vertex(vertexIterator->next());
            Array<T, 3> v(points.col(i)[0], points.col(i)[1], points.col(i)[2]);
            v /= dx;

            /*
            The periodicity is fully handled by the Palabos DenseParticles
            interface. However, the periodicity on the visualization of the 
            bodies is handled by this simple function, since ShapeOp does not 
            care about the domain boundaries
            */
            bool x = false, y = false, z = false;
            controlPositionPeriodic<T>(v, nx - 1, ny - 1, nz - 1, x, y, z);

            (*vertex)[0] = v[0];
            (*vertex)[1] = v[1];
            (*vertex)[2] = v[2];
            
            ++i;
        }

        return mesh;
    }

    // Move ShapeOp solvers inside the domain
    void imposePlbPeriodicityToSO(T dx, pluint nx, pluint ny, pluint nz)
    {
        const plb::npfem::Matrix3X& points = solver.getPoints();
        plb::npfem::Vector3 centroid = points.rowwise().mean();
        Array<T, 3> c_lb(centroid[0] / dx, centroid[1] / dx, centroid[2] / dx);
        Array<T, 3> c_lb_init = c_lb;

        bool x = false, y = false, z = false;
        controlPositionPeriodic<T>(c_lb, nx - 1, ny - 1, nz - 1, x, y, z);

        if (x || y || z)
        {
            Array<T, 3> displace = c_lb - c_lb_init;
            plb::npfem::Vector3 displaceSO(displace[0] * dx, displace[1] * dx, displace[2] * dx);
            solver.shiftPoints(displaceSO);
        }
    }

	std::vector<int> *mesh_graph;

private:
    pluint bodyID;
    ShapeOp_Solver solver;
    RawConnectedTriangleMesh<T>* templateMesh;
    pluint numVertices_onSurface_;
    
    // A body may belong in more than one processors in Palabos side.
    // Consequently, it is going to receive the forces
    // from the corresponding processors.
    pluint numReceived;
    plb::npfem::Matrix3X commForces;
    plb::npfem::Matrix3X shearForces, normals;
	plb::npfem::VectorX pressure, area;

    std::vector<pluint> taskIDs; // from where it receives data
    std::vector<std::vector<pluint>> vertexIDsPerProcessor;
    

    // ShapeOp sends back velocities
    // We use these containers, because we are using Isend
    // and we need to keep the buffers safe
    std::vector<std::vector<Array<T, 3>>> velocities;

    // MPI Completion
    std::vector<MPI_Request*> mpiRequests;

    // Collisions Handling
    Array<T, 3> SOLocalMeshCenterOfMass;
    std::vector<Array<T, 3>> collisionNeighbors, collisionNeighborsNormals;

    // Force Filtering
    plb::npfem::Matrix3X forceFiltered;
};

// For a given body ID, returns the index inside the vector of ShapeOp solvers
// that has the solver of that body.
template <typename T>
pluint getSolverID(pluint bodyID, std::vector<ShapeOpBody<T>> const& shapeOpBodies)
{
    plint solverID = -1;
    for (plint i = 0; i < (plint)shapeOpBodies.size(); ++i) {
        if (shapeOpBodies[i].getID() == bodyID)
            solverID = i;
    }
    PLB_ASSERT(solverID != -1);
    return (pluint)solverID;
}

// Checks if all local body solvers have received the newest force data.
template <typename T>
bool communicationComplete(std::vector<ShapeOpBody<T>>& shapeOpBodies)
{
    pluint cnt = 0;
    for (pluint i = 0; i < shapeOpBodies.size(); ++i)
    {
        if (shapeOpBodies[i].communicationComplete())
            cnt += 1;
    }
    
    if (cnt == (pluint)shapeOpBodies.size())
        return true;
    else
        return false;
}

// Create a Palabos Mesh from a ShapeOp solver.
template <typename T>
RawConnectedTriangleMesh<T> generateMesh(ShapeOp_Solver& solver, T dx)
{
    std::vector<Array<T, 3>> vertices;
    std::vector<Array<plint, 3>> triangles;

    // ShapeOp in physical units
    // All points (OnSurface + interior)
    const plb::npfem::Matrix3X& points = solver.getPoints();
    for (pluint i = 0; i < (pluint)points.cols(); ++i)
    {
        Array<T, 3> vertex(points.col(i)[0], points.col(i)[1], points.col(i)[2]);
        // Convert to lb units
        vertices.push_back(vertex / dx);
    }

    // Connectivity of the Surface only
    for (auto i = solver.getConnectivityList().begin(); i < solver.getConnectivityList().end(); ++i)
        triangles.push_back(Array<plint, 3>((*i)[0], (*i)[1], (*i)[2]));

    RawConnectedTriangleMesh<T> mesh(vertices, triangles);
    
    /*
    mesh.registerVertexProperty("shearForce_0");
    mesh.registerVertexProperty("shearForce_1");
    mesh.registerVertexProperty("shearForce_2");
    mesh.registerVertexProperty("normal_0");
    mesh.registerVertexProperty("normal_1");
    mesh.registerVertexProperty("normal_2");
    mesh.registerVertexProperty("pressure");
    mesh.registerVertexProperty("area");
    */
    mesh.registerVertexTag("isWritten");
    
    return mesh;
}


/*
template <typename T>
RawConnectedTriangleMesh<T> generateNeighbors(ShapeOp_Solver& solver, T dx,
    pluint nx, pluint ny, pluint nz)
{
    std::vector<Array<T, 3>> vertices;
    std::vector<Array<plint, 3>> triangles;

    const plb::npfem::Matrix3X& points = solver.verticesOfCollidingPieces_.back();

    for (pluint i = 0; i < (pluint)points.cols(); ++i) {
        Array<T, 3> vertex(points.col(i)[0], points.col(i)[1], points.col(i)[2]);
        vertex /= dx;
        
        bool x = false, y = false, z = false;
        controlPositionPeriodic<T>(vertex, nx - 1, ny - 1, nz - 1, x, y, z);

        vertices.push_back(vertex);
    }

    if (points.cols() >= 3) {
        for (pluint i = 2; i < (pluint)points.cols(); ++i) {
            triangles.push_back(Array<plint, 3>(i-2, i-1, i));
        }
    }

    RawConnectedTriangleMesh<T> mesh(vertices, triangles);
    mesh.registerVertexProperty("shearForce_0");
    mesh.registerVertexProperty("shearForce_1");
    mesh.registerVertexProperty("shearForce_2");
    mesh.registerVertexProperty("normal_0");
    mesh.registerVertexProperty("normal_1");
    mesh.registerVertexProperty("normal_2");
    mesh.registerVertexProperty("pressure");
    mesh.registerVertexProperty("area");
    mesh.registerVertexTag("isWritten");
    
    return mesh;
}
*/

}
}