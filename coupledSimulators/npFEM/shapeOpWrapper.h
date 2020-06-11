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
#ifndef SHAPEOPWRAPPER_H
#define SHAPEOPWRAPPER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "Solver.h"
#ifndef NPFEM_SA
#ifdef NPFEM_CUDA
#include "Solver_GPU.h"
#endif
#endif // !NPFEM_SA
#include "Constraint.h"
#include "Force.h"

namespace plb {
namespace npfem {

typedef plb::npfem::Solver ShapeOp_Solver;

void setPointsFromCSV(ShapeOp_Solver& s, std::string filename, bool best_fit_T = false);
void setVelsFromCSV(ShapeOp_Solver& s, std::string filename);

void savePointsToCSV(ShapeOp_Solver& s, std::string filename, size_t iT = 1, size_t dt_ShapeOp = 1);
void saveVelsToCSV(ShapeOp_Solver& s, std::string filename);

void setConstraintsFromCSV(ShapeOp_Solver& s, std::string filename);

void setConnectivityListFromCSV(ShapeOp_Solver& s, std::string filename);

void setForcesFromCSV(ShapeOp_Solver& s, std::string filename);

void saveForcesToCSV(ShapeOp_Solver& s, std::string filename);

void setOnSurfaceParticle(ShapeOp_Solver& s, std::string filename);

// This runs only once as an initialization. Afterwards, modify forces with
// editVertexForce
void addVertexForce(ShapeOp_Solver& s, const plb::npfem::Matrix3X& forces, const int cell_id = 0);

void editVertexForce(ShapeOp_Solver& s, const plb::npfem::Matrix3X& forces, const int cell_id = 0);

}
}

///////////////////////////////////////////////////////////////////////////////
#ifdef SHAPEOP_HEADER_ONLY
#include "shapeOpWrapper.cpp"
#endif
///////////////////////////////////////////////////////////////////////////////
#endif
