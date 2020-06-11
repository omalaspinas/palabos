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
#ifndef CUDASTRCT
#define CUDASTRCT

#include "sparse_matrix.h"
#include "common.h"

namespace plb {
namespace npfem {

struct matrix3{
	float data[9];
};

//MESH
struct Mesh_info {
	float rho = 1125.0;
	float miu = 35.;
	float lambda = 5;
	float kappa = 0.;
	float dx_p = 0.5;
    float threshold_nonRep = 1.;
    float weight_col_nonRep = 20.;
    float threshold_rep = 1.;
    float weight_col_rep = 20.;
    float beta_morse = 1.;
	float volume_weight = 250.;
	float calpha = 0.7;
	float vel_cap = 0.04;
	float vel_cap_fin = 0.004;
	cuda_scalar volume = 100;
	ShapeOpScalar mass_total = 0;
	int nb_cells = 1;
	int n_points = 0;
	int idO = 0;
	int n_constraints = 0;
	int n_triangles = 0;
	int sum_points = 256;
	int sum_triangles = 512;
};

struct graph_data
{
	int *indexs = NULL;
	int *edges = NULL;
	int nb_edges = 0;
};

//MESH AND CONSTRAINTS data; this data is the same for every cell and doesnt change over time.
struct Mesh_data
{
	//ntriangle
	short       *triangles;
	cuda_scalar *A;
	//nconstraint
	int *ConstraintType;
	int *idO;
	ShapeOpScalar *rangeMin;
	ShapeOpScalar *rangeMax;
	ShapeOpScalar *Scalar1;
	ShapeOpScalar *weight;
	//nconstraint*4;
	int *idI;
	ShapeOpScalar *vectorx;
	ShapeOpScalar *matrix22;
	//nconstraint*8
	ShapeOpScalar *matrix33;
	//npoint
	ShapeOpScalar *mass;
	//npoint*degree
	sparse_matrix_cuda L;
	sparse_matrix_cuda J;
	sparse_matrix_cuda M_star;
	//npoint^2
	ShapeOpScalar  *N_inv_dense;
	ShapeOpScalar  *M_tild_inv;
};

struct Collision_data
{
	int total_neighbours = 0;
	int *nb_neighbours;
    ShapeOpScalar *colid_points;
    ShapeOpScalar *colid_normals;
};

//SIMULATION data; that is to say data different for each cells, wich gonna change over time.
struct Simulation_input
{
	//3*npoint*ncell
	ShapeOpScalar *forces_ex;
	ShapeOpScalar *points;
	ShapeOpScalar *velocities;
	//ncell
	int *nb_colid_search;
};
//Simulation data that are intermediate value which doesnt need to be initialize on the CPU.
struct Simulation_data
{
	//3*npoint*ncell
	ShapeOpScalar *oldPoints;
	ShapeOpScalar *points_last_iter;
	ShapeOpScalar *projections;
	ShapeOpScalar *momentum;
	cuda_scalar   *force_npd;

	ShapeOpScalar *gradient;
	ShapeOpScalar *gradient_prev;
	ShapeOpScalar *descent_direction;
	ShapeOpScalar *q;
	ShapeOpScalar *energies;
	ShapeOpScalar *energies_prev;
	ShapeOpScalar *nearest_normals;
	ShapeOpScalar *nearest_points;
  
	//nb constraint
	ShapeOpScalar *E_nonePD;
	//nbcell
	ShapeOpScalar *descent_direction_dot_gradient;

	//memesize*3*npoint*ncell 
	ShapeOpScalar *s;
	ShapeOpScalar *t;
	//memesize*ncell
	ShapeOpScalar *rho;
	ShapeOpScalar *alpha;
	ShapeOpScalar *center;
	//ncell
	int *has_converged;
};

}
}

#endif
