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

namespace plb {
namespace npfem {

// Do not use MPI_TAG = 0 because it is reserved for Palabos
const int MPItagIniPalabosParticles = 2709;

const int shapeOpMPItagForcesAndCollisionContainersBodyID = 2710;
const int shapeOpMPItagForcesAndCollisionContainersNumVertices = 2711;
const int shapeOpMPItagForcesAndCollisionContainersVertexIDs = 2712;
const int shapeOpMPItagForcesAndCollisionContainersShearForces = 2713;
const int shapeOpMPItagForcesAndCollisionContainersNormals = 2714;
const int shapeOpMPItagForcesAndCollisionContainersPressure = 2715;
const int shapeOpMPItagForcesAndCollisionContainersArea = 2716;
const int shapeOpMPItagForcesAndCollisionContainersNumCollidingNeighbors = 2717;
const int shapeOpMPItagForcesAndCollisionContainersCollisionNeighbors = 2718;
const int shapeOpMPItagForcesAndCollisionContainerscollisionNeighborsNormals = 2719;

const int shapeOpMPItagVelocities = 2730;

const pluint IDforWall = std::numeric_limits<pluint>::max();

}
}