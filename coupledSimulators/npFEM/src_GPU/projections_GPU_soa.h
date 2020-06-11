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
#ifndef PROJECTIONS_GPU_H
#define PROJECTIONS_GPU_H
///////////////////////////////////////////////////////////////////////////////
#include <cuda_runtime.h>

#include "common.h"
#include "sum_cuda.h"
#include "projections_GPU_MATH.h"

namespace plb {
namespace npfem {
////////////////////////////////////////////////////////////////////////////////
__device__ void project_volume_d(int n_points, int n_tri, int n_constraints, double *points_d, short *triangles_d, cuda_scalar *force, 
								 double *E_nonePD, cuda_scalar volume0, const int nb_sum_thread, double *buffer_double, cuda_scalar *grad_c, const int id, const float volume_weight) {

	grad_c[threadIdx.x             ] = 0;
	grad_c[threadIdx.x +   n_points] = 0;
	grad_c[threadIdx.x + 2*n_points] = 0;

	const int cell_shift = blockIdx.x*3*n_points;
	int tp_id;

	//float *buffer = (float*)buffer_double;

	buffer_double[threadIdx.x] = 0;
	__syncthreads();

	for(int tid = threadIdx.x; tid < n_tri; tid += blockDim.x){
		//printf("tid %d \n", tid);
						  
		short id0 = triangles_d[IDX(tid, 0, n_tri)];
		short id1 = triangles_d[IDX(tid, 1, n_tri)];
		short id2 = triangles_d[IDX(tid, 2, n_tri)];

		//printf("ido %d %d %d\n", id0, id1, id2);
		tp_id = cell_shift + id0;

		cuda_scalar p0x = points_d[tp_id             ];
		cuda_scalar p0y = points_d[tp_id +   n_points];
		cuda_scalar p0z = points_d[tp_id + 2*n_points];

		tp_id = cell_shift + id1;

		cuda_scalar p1x = points_d[tp_id             ];
		cuda_scalar p1y = points_d[tp_id +   n_points];
		cuda_scalar p1z = points_d[tp_id + 2*n_points];

		tp_id = cell_shift + id2;

		cuda_scalar p2x = points_d[tp_id             ];
		cuda_scalar p2y = points_d[tp_id +   n_points];
		cuda_scalar p2z = points_d[tp_id + 2*n_points];

		cuda_scalar edge0_x = p1x - p0x;
		cuda_scalar edge0_y = p1y - p0y;
		cuda_scalar edge0_z = p1z - p0z;

		cuda_scalar edge1_x = p2x - p0x;
		cuda_scalar edge1_y = p2y - p0y;
		cuda_scalar edge1_z = p2z - p0z;
	
		cuda_scalar n1, n0, n2;
		cross(edge0_x, edge0_y, edge0_z, edge1_x, edge1_y, edge1_z, n0, n1, n2);
		cuda_scalar area = Normalize_3(n0, n1, n2)/2;
	
		buffer_double[threadIdx.x]  += area*(n0*(p0x + p1x + p2x) + n1*(p0y + p1y + p2y) + n2*(p0z + p1z + p2z))/9;

		atomicAdd(grad_c + id0             , (area*n0)/3);
		atomicAdd(grad_c + id0 +   n_points, (area*n1)/3);
		atomicAdd(grad_c + id0 + 2*n_points, (area*n2)/3);

		atomicAdd(grad_c + id1             , (area*n0)/3);
		atomicAdd(grad_c + id1 +   n_points, (area*n1)/3);
		atomicAdd(grad_c + id1 + 2*n_points, (area*n2)/3);

		atomicAdd(grad_c + id2             , (area*n0)/3);
		atomicAdd(grad_c + id2 +   n_points, (area*n1)/3);
		atomicAdd(grad_c + id2 + 2*n_points, (area*n2)/3);
	}
	__syncthreads();
	//if (threadIdx.x == 224)printf("nb_sum_thread %d \n", nb_sum_thread);
	//if (threadIdx.x == 0)printf("v0 %f v1 %f \n", volume0, buffer[0]);

	cuda_scalar dx = grad_c[threadIdx.x             ];
	cuda_scalar dy = grad_c[threadIdx.x +   n_points];
	cuda_scalar dz = grad_c[threadIdx.x + 2*n_points];

	__syncthreads();
	sum(buffer_double, nb_sum_thread, n_points, threadIdx.x);
	__syncthreads();

	//if (threadIdx.x == 0)printf("norm pourri %d %d\n", nb_sum_thread, blockDim.x);
	cuda_scalar C = (buffer_double[0] - volume0);
	__syncthreads();
	buffer_double[threadIdx.x] = dx*dx + dy*dy + dz*dz;
	sum(buffer_double, nb_sum_thread, n_points, threadIdx.x);
	__syncthreads();

	dx = (C/buffer_double[0])*dx;
	dy = (C/buffer_double[0])*dy;
	dz = (C/buffer_double[0])*dz;

	force[id             ] -= volume_weight*dx;
	force[id +   n_points] -= volume_weight*dy;
	force[id + 2*n_points] -= volume_weight*dz;
	
	E_nonePD[blockIdx.x*n_constraints + threadIdx.x] += 0.5*volume_weight*(dx*dx + dy*dy + dz*dz);
	
	//if(threadIdx.x == 224)printf("force volume %f %f %f | %d\n", force[id], force[id + n_points], force[id + 2 * n_points], threadIdx.x);
}
///////////////////////////////////////////////////////////////////////////////
__device__
void project_Bending(int tid, int n_points, int n_constraints, int n_projected_points, int nb_cells,
	ShapeOpScalar *points_d, ShapeOpScalar *projections_d, cuda_scalar *f_int_nonePD_d,
	int *ConstraintType_d, int *idO_d, ShapeOpScalar *rangeMin_d, ShapeOpScalar *rangeMax_d, ShapeOpScalar *Scalar1_d, ShapeOpScalar *weight_d, ShapeOpScalar *E_nonePD_d,
	int *idI_d, ShapeOpScalar *vectorx_d, ShapeOpScalar *matrix22_d, ShapeOpScalar *matrix33_d)
{
	cuda_scalar rangeMin_ = rangeMin_d[tid];
	cuda_scalar rangeMax_ = rangeMax_d[tid];
	cuda_scalar weight_ = weight_d[tid];
	int idO_ = idO_d[tid] + blockIdx.x*3*n_projected_points;
	int cell_shift;
	cuda_scalar n_ = Scalar1_d[tid];

	cuda_scalar w0_, w1_, w2_, w3_;
	int idI0_, idI1_, idI2_, idI3_;

	w0_ = vectorx_d[IDX(tid, 0, n_constraints)];
	w1_ = vectorx_d[IDX(tid, 1, n_constraints)];
	w2_ = vectorx_d[IDX(tid, 2, n_constraints)];
	w3_ = vectorx_d[IDX(tid, 3, n_constraints)];

	cell_shift = blockIdx.x*3*n_points;
	// Involved points (indices)
	idI0_ = idI_d[tid                  ] + cell_shift;
	idI1_ = idI_d[tid +   n_constraints] + cell_shift;
	idI2_ = idI_d[tid + 2*n_constraints] + cell_shift;
	idI3_ = idI_d[tid + 3*n_constraints] + cell_shift;

	cuda_scalar e0=0, e1=0, e2=0;
	if (n_ > 1e-6){
		e0 = w0_ * points_d[IDX(idI0_, 0, n_points)] + w1_ * points_d[IDX(idI1_, 0, n_points)] + w2_ * points_d[IDX(idI2_, 0, n_points)] + w3_ * points_d[IDX(idI3_, 0, n_points)];
		e1 = w0_ * points_d[IDX(idI0_, 1, n_points)] + w1_ * points_d[IDX(idI1_, 1, n_points)] + w2_ * points_d[IDX(idI2_, 1, n_points)] + w3_ * points_d[IDX(idI3_, 1, n_points)];
		e2 = w0_ * points_d[IDX(idI0_, 2, n_points)] + w1_ * points_d[IDX(idI1_, 2, n_points)] + w2_ * points_d[IDX(idI2_, 2, n_points)] + w3_ * points_d[IDX(idI3_, 2, n_points)];

		cuda_scalar l = Norm_3(e0, e1, e2);

		if (l > 1e-6){
			e0 /= l;
			e1 /= l;
			e2 /= l;
			l = n_ * CLAMP(l / n_, rangeMin_, rangeMax_);
			e0 *= l;
			e1 *= l;
			e2 *= l;
		}
	}

	projections_d[IDX(idO_, 0, n_projected_points)] = weight_ * e0;
	projections_d[IDX(idO_, 1, n_projected_points)] = weight_ * e1;
	projections_d[IDX(idO_, 2, n_projected_points)] = weight_ * e2;
	/*
	if (tid == 107){
		printf("projection %d %d %d %d \n", idI0_, idI1_, idI2_, idI3_);
		printf("weight %f %f %f %f n %f\n", w0_, w1_, w2_, w3_, n_);
		printf("e %f %f %f \n", e0, e1, e2);
		printf("projection %f %f %f\n", points_d[IDX(idI0_, 0, n_points)], points_d[IDX(idI0_, 1, n_points)], points_d[IDX(idI0_, 2, n_points)]);
		printf("projection %f %f %f\n", points_d[IDX(idI1_, 0, n_points)], points_d[IDX(idI1_, 1, n_points)], points_d[IDX(idI1_, 2, n_points)]);
		printf("projection %f %f %f\n", points_d[IDX(idI2_, 0, n_points)], points_d[IDX(idI2_, 1, n_points)], points_d[IDX(idI2_, 2, n_points)]);
		printf("projection %f %f %f\n", points_d[IDX(idI3_, 0, n_points)], points_d[IDX(idI3_, 1, n_points)], points_d[IDX(idI3_, 2, n_points)]);
		printf("projection %f %f %f \n", projections_d[IDX(idO_, 0, n_projected_points)], projections_d[IDX(idO_, 1, n_projected_points)], projections_d[IDX(idO_, 2, n_projected_points)]);
	}
	*/
}
///////////////////////////////////////////////////////////////////////////////

__device__ 
void project_Area(int tid, int n_points, int n_constraints, int n_projected_points, int nb_cells,
	ShapeOpScalar *points_d, ShapeOpScalar *projections_d, cuda_scalar *f_int_nonePD_d,
	int *ConstraintType_d, int *idO_d, ShapeOpScalar *rangeMin_d, ShapeOpScalar *rangeMax_d, ShapeOpScalar *Scalar1_d, ShapeOpScalar *weight_d, ShapeOpScalar *E_nonePD_d,
	int *idI_d, ShapeOpScalar *vectorx_d, ShapeOpScalar *matrix22_d, ShapeOpScalar *matrix33_d)
{
	cuda_scalar rangeMin_ = rangeMin_d[tid];
	cuda_scalar rangeMax_ = rangeMax_d[tid];
	cuda_scalar weight_ = weight_d[tid];
	int idO_ = idO_d[tid] + blockIdx.x*3*n_projected_points;

	cuda_scalar rest00_, rest01_,
		rest10_, rest11_;
	int idI0_, idI1_, idI2_, cell_shift;

	// Involved points (indices)
	cell_shift = blockIdx.x*3*n_points;
	idI0_ = idI_d[IDX(tid, 0, n_constraints)] + cell_shift;
	idI1_ = idI_d[IDX(tid, 1, n_constraints)] + cell_shift;
	idI2_ = idI_d[IDX(tid, 2, n_constraints)] + cell_shift;
	//if (tid == 340)printf("ids %d %d %d tid %d\n", idI0_, idI1_, idI2_, tid);

	rest00_ = matrix22_d[IDX(tid, 0, n_constraints)];
	rest10_ = matrix22_d[IDX(tid, 1, n_constraints)];
	rest01_ = matrix22_d[IDX(tid, 2, n_constraints)];
	rest11_ = matrix22_d[IDX(tid, 3, n_constraints)];
	//if (tid == 768)printf("matrix2 %f %f %f %f tid %d\n", rest00_, rest10_, rest01_, rest11_, tid);

	cuda_scalar edges00, edges01,
		edges10, edges11,
		edges20, edges21;
	cuda_scalar P00, P01,
		P10, P11,
		P20, P21;

	edges00 = points_d[IDX(idI1_, 0, n_points)] - points_d[IDX(idI0_, 0, n_points)];
	edges10 = points_d[IDX(idI1_, 1, n_points)] - points_d[IDX(idI0_, 1, n_points)];
	edges20 = points_d[IDX(idI1_, 2, n_points)] - points_d[IDX(idI0_, 2, n_points)];
	// edges.col(1)
	edges01 = points_d[IDX(idI2_, 0, n_points)] - points_d[IDX(idI0_, 0, n_points)];
	edges11 = points_d[IDX(idI2_, 1, n_points)] - points_d[IDX(idI0_, 1, n_points)];
	edges21 = points_d[IDX(idI2_, 2, n_points)] - points_d[IDX(idI0_, 2, n_points)];

	P00 = edges00;
	P10 = edges10;
	P20 = edges20;
	Normalize_3(P00, P10, P20);

	cuda_scalar dot = DOT_3(edges01, edges11, edges21, P00, P10, P20);
	P01 = edges01 - dot * P00;
	P11 = edges11 - dot * P10;
	P21 = edges21 - dot * P20;
	Normalize_3(P01, P11, P21);
	// P.col(1)

	cuda_scalar tmp00, tmp01,
				tmp10, tmp11;
	cuda_scalar F00, F01,
				F10, F11;

	Matrix_Product_23_32(P00, P10, P20,
						 P01, P11, P21,
		//
						edges00, edges01,
						edges10, edges11,
						edges20, edges21,
		//
						tmp00, tmp01,
						tmp10, tmp11);

	Matrix_Product_22_22(tmp00, tmp01,
						 tmp10, tmp11,
		//
						rest00_, rest01_,
						rest10_, rest11_,
		//
						F00, F01,
						F10, F11);

	// SVD
	cuda_scalar U00, U01,
		U10, U11;
	cuda_scalar SIG00, SIG11;
	cuda_scalar V00, V01,
		V10, V11;

	svd_joel22(F00, F01,
			F10, F11,
			//
			U00, U01,
			U10, U11,
			//
			SIG00, SIG11,
			//
			V00, V01,
			V10, V11);

	cuda_scalar S0 = SIG00, S1 = SIG11;
	cuda_scalar d0 = 0., d1 = 0.;
	for (int i = 0; i < SHAPEOP_INNER_ITERATIONS; ++i)
	{
		cuda_scalar v = S0 * S1;
		cuda_scalar f = v - CLAMP(v, rangeMin_, rangeMax_);
		cuda_scalar g0 = S1, g1 = S0;
		cuda_scalar dot_gd = DOT_2(g0, g1, d0, d1);
		cuda_scalar dot_gg = DOT_2(g0, g1, g0, g1);
		cuda_scalar frac = -((f - dot_gd) / dot_gg);
		d0 = frac * g0; d1 = frac * g1;
		S0 = SIG00 + d0; S1 = SIG11 + d1;
	}
	SIG00 = S0; 
    SIG11 = S1;

	Matrix_Product_22_22(U00, U01,
		U10, U11,
		//
		SIG00, (cuda_scalar)0,
		(cuda_scalar)0, SIG11,
		//
		tmp00, tmp01,
		tmp10, tmp11);

	Matrix_Product_22_22(tmp00, tmp01,
						 tmp10, tmp11,

						V00, V10,
						V01, V11,

						F00, F01,
						F10, F11);

	cuda_scalar PF00, PF01,
		PF10, PF11,
		PF20, PF21;
	Matrix_Product_32_22(P00, P01,
		P10, P11,
		P20, P21,
		//
		F00, F01,
		F10, F11,
		//
		PF00, PF01,
		PF10, PF11,
		PF20, PF21);

	projections_d[IDX(idO_, 0, n_projected_points)] = weight_ * PF00; projections_d[IDX(idO_ + 1, 0, n_projected_points)] = weight_ * PF01;
	projections_d[IDX(idO_, 1, n_projected_points)] = weight_ * PF10; projections_d[IDX(idO_ + 1, 1, n_projected_points)] = weight_ * PF11;
	projections_d[IDX(idO_, 2, n_projected_points)] = weight_ * PF20; projections_d[IDX(idO_ + 1, 2, n_projected_points)] = weight_ * PF21;

	//if(tid == 340)printf("proj %f %f %f | %d %d\n", projections_d[IDX(idO_, 0, n_projected_points)] , projections_d[IDX(idO_, 1, n_projected_points)], projections_d[IDX(idO_, 2, n_projected_points)], idO_, n_projected_points);

	//#ifdef DEBUG
	//if (threadIdx.x == 109)printf("proj  %f %f %f | %d \n", projections_d[IDX(idO_, 0, n_projected_points)], projections_d[IDX(idO_, 1, n_projected_points)], projections_d[IDX(idO_, 2, n_projected_points)], blockIdx.x);
	//#endif
	
}

///////////////////////////////////////////////////////////////////////////////
__device__ 
void project_surface_material_and_area(int tid, int n_points, int n_constraints, int n_projected_points, int nb_cells,
	ShapeOpScalar *points_d, ShapeOpScalar *projections_d, cuda_scalar *f_int_nonePD_d,
	int *ConstraintType_d, int *idO_d, ShapeOpScalar *rangeMin_d, ShapeOpScalar *rangeMax_d, ShapeOpScalar *Scalar1_d, ShapeOpScalar *weight_d, ShapeOpScalar *E_nonePD_d,
	int *idI_d, ShapeOpScalar *vectorx_d, ShapeOpScalar *matrix22_d, ShapeOpScalar *matrix33_d, cuda_scalar *A_d, const float miu, const float lambda, const float kappa)
{
  //printf("called material \n");

  cuda_scalar rangeMin_ = rangeMin_d[tid];
  cuda_scalar rangeMax_ = rangeMax_d[tid];
  cuda_scalar weight_   = weight_d[tid];
  int idO_              = idO_d[tid + 512];

  cuda_scalar rest00_, rest01_,
              rest10_, rest11_;
  int idI0_, idI1_, idI2_;

  // Involved points (indices)
  idI0_ = idI_d[IDX(tid, 0, n_constraints)];
  idI1_ = idI_d[IDX(tid, 1, n_constraints)];
  idI2_ = idI_d[IDX(tid, 2, n_constraints)];

  // ColMajor Order
  rest00_ = matrix22_d[IDX(tid, 0, n_constraints)];
  rest10_ = matrix22_d[IDX(tid, 1, n_constraints)];
  rest01_ = matrix22_d[IDX(tid, 2, n_constraints)];
  rest11_ = matrix22_d[IDX(tid, 3, n_constraints)];

  cuda_scalar edges00, edges01,
              edges10, edges11,
              edges20, edges21;

  cuda_scalar P00, P01,
              P10, P11,
              P20, P21;

  edges00 = points_d[IDX(idI1_, 0, n_points)] - points_d[IDX(idI0_, 0, n_points)];
  edges10 = points_d[IDX(idI1_, 1, n_points)] - points_d[IDX(idI0_, 1, n_points)];
  edges20 = points_d[IDX(idI1_, 2, n_points)] - points_d[IDX(idI0_, 2, n_points)];
  // edges.col(1)
  edges01 = points_d[IDX(idI2_, 0, n_points)] - points_d[IDX(idI0_, 0, n_points)];
  edges11 = points_d[IDX(idI2_, 1, n_points)] - points_d[IDX(idI0_, 1, n_points)];
  edges21 = points_d[IDX(idI2_, 2, n_points)] - points_d[IDX(idI0_, 2, n_points)];

  //debug to be deleled!!
  //edges00 += 0.01;

  P00 = edges00;
  P10 = edges10;
  P20 = edges20;
  Normalize_3(P00, P10, P20);

  cuda_scalar dot = DOT_3(edges01, edges11, edges21, P00, P10, P20);
  P01 = edges01 - dot * P00;
  P11 = edges11 - dot * P10;
  P21 = edges21 - dot * P20;
  Normalize_3(P01, P11, P21);
  // P.col(1)

  cuda_scalar tmp00, tmp01,
			  tmp10, tmp11;
  cuda_scalar F00, F01,
			  F10, F11;

  Matrix_Product_23_32(P00, P10, P20,
					   P01, P11, P21,
	//
					edges00, edges01,
					edges10, edges11,
					edges20, edges21,
	//
					tmp00, tmp01,
					tmp10, tmp11);

  Matrix_Product_22_22(tmp00, tmp01,
					   tmp10, tmp11,
	//
					   rest00_, rest01_,
					   rest10_, rest11_,
	//
					   F00, F01,
					   F10, F11);


  // [U SIG V] = SVD(F)
  // SVD
  cuda_scalar U00, U01,
              U10, U11;

  cuda_scalar SIG00, SIG11; 
 
  cuda_scalar V00, V01,
              V10, V11;

     svd_joel22(F00, F01,
             F10, F11,
          //
             U00, U01,
             U10, U11,
          //
             SIG00, SIG11,
          //
             V00, V01,
             V10, V11);
  //save those value for later
  cuda_scalar S0 = SIG00, S1 = SIG11;
/*
if( idI0_ == 0 && idI1_ == 1 && idI2_ == 2){
	printf("rest %f %f \nrest %f %f \n", rest00_, rest01_, rest10_, rest11_);
	printf("edge %f %f \nedge %f %f  \nedge %f %f\n", edges00, edges01, edges10, edges11, edges20, edges21);
	printf("P %f %f \nP %f %f \nP %f %f  \n", P00, P01, P10, P11, P20, P21);
	printf("F %f %f \nF %f %f \n", F00, F01, F10, F11);
    printf("U %f %f \nU %f %f \n\nV %f %f \nV %f %f \n---------------------------------\n", U00, U01, U10, U11, V00, V01, V10, V11);
	printf("singular %f %f \n", SIG00, SIG11);
}
*/

  cuda_scalar A = -0.5*A_d[tid];
 //none linear business
  E_nonePD_d[tid] = A_d[tid]*(f_tr(SIG00, miu, lambda, kappa) + f_tr(SIG11, miu, lambda, kappa))/2;
	//printf("e_non_pd %f %f\n", E_nonePD_d[tid], A_d[tid]);
  /*
  if (idI0_ == 0 && idI1_ == 1 && idI2_ == 2) {
	  printf("singular %f %f \n", SIG00, SIG11);
  }
  */
 //F = U*S*VT = PIOLA
  Matrix_Product_22_22(U00, U01,
                       U10, U11,
                       //
						f_prime_tr(SIG00, miu, lambda, kappa), (cuda_scalar)0.,
                       (cuda_scalar)0., f_prime_tr(SIG11, miu, lambda, kappa),
                       //
                       tmp00, tmp01,
                       tmp10, tmp11);

  Matrix_Product_22_22(tmp00, tmp01,
                       tmp10, tmp11,
                       //
                       V00, V10,
                       V01, V11,
                       //
                       F00, F01,
                       F10, F11);
  /*
  if (idI0_ == 0 && idI1_ == 1 && idI2_ == 2) {
	  printf("Piola %f %f \nPiola %f %f A %f \n", F00, F01, F10, F11, A_d[tid]);
  }
  */
  // tmp =  Piola*rest_T
  Matrix_Product_22_22(F00, F01,
					   F10, F11,
					   rest00_, rest10_,
					   rest01_, rest11_,
					   tmp00, tmp01,
					   tmp10, tmp11);
	/*
  if (idI0_ == 0 && idI1_ == 1 && idI2_ == 2) {
	  printf("Piola*rest %f %f \nPiola*rest %f %f\n", tmp00, tmp01, tmp10, tmp11);
  }
  */
  cuda_scalar F20, F21;
  //F = P*Piola*rest_T;
  Matrix_Product_32_22(P00, P01,
					   P10, P11,
	                   P20, P21,
					   tmp00, tmp01,
					   tmp10, tmp11,
	                   F00, F01,
	                   F10, F11,
	                   F20, F21);
 
 /*
  if (idI0_ == 0 && idI1_ == 1 && idI2_ == 2) {
	  printf("P*Piola*rest %f %f \nP*Piola*rest %f %f\n P*Piola*rest %f %f\n", A*F00, A*F01, A*F10, A*F11, A*F20, A*F21);
  }
 */
  atomicAdd(f_int_nonePD_d + IDX(idI0_, 0, n_points), A*(-F00 - F01));
  atomicAdd(f_int_nonePD_d + IDX(idI0_, 1, n_points), A*(-F10 - F11));
  atomicAdd(f_int_nonePD_d + IDX(idI0_, 2, n_points), A*(-F20 - F21));

  atomicAdd(f_int_nonePD_d + IDX(idI1_, 0, n_points), A*F00);
  atomicAdd(f_int_nonePD_d + IDX(idI1_, 1, n_points), A*F10);
  atomicAdd(f_int_nonePD_d + IDX(idI1_, 2, n_points), A*F20);

  atomicAdd(f_int_nonePD_d + IDX(idI2_, 0, n_points), A*F01);
  atomicAdd(f_int_nonePD_d + IDX(idI2_, 1, n_points), A*F11);
  atomicAdd(f_int_nonePD_d + IDX(idI2_, 2, n_points), A*F21);
 
  cuda_scalar d0 = 0., d1 = 0.;
  for (int i = 0; i < SHAPEOP_INNER_ITERATIONS; ++i)
  {
	  cuda_scalar v = S0 * S1;
	  cuda_scalar f = v - CLAMP(v, rangeMin_, rangeMax_);
	  cuda_scalar g0 = S1, g1 = S0;
	  cuda_scalar dot_gd = DOT_2(g0, g1, d0, d1);
	  cuda_scalar dot_gg = DOT_2(g0, g1, g0, g1);
	  cuda_scalar frac = -((f - dot_gd) / dot_gg);
	  d0 = frac * g0; d1 = frac * g1;
	  S0 = SIG00 + d0; S1 = SIG11 + d1;
  }
  SIG00 = S0;
  SIG11 = S1;

  Matrix_Product_22_22(U00, U01,
	                   U10, U11,
	//
	      SIG00, (cuda_scalar)0,
	      (cuda_scalar)0, SIG11,
	//
	     tmp00, tmp01,
	     tmp10, tmp11);

Matrix_Product_22_22(tmp00, tmp01,
	                 tmp10, tmp11,

	                 V00, V10,
	                 V01, V11,

	                 F00, F01,
	                 F10, F11);

  cuda_scalar PF00, PF01,
	          PF10, PF11,
	          PF20, PF21;

  Matrix_Product_32_22(P00, P01,
	                   P10, P11,
	                   P20, P21,
	//
	                   F00, F01,
	                   F10, F11,
	//
	                   PF00, PF01,
	                   PF10, PF11,
	                   PF20, PF21);

	//printf("ido %d \n", idO_);

   projections_d[IDX(idO_ , 0, n_projected_points)] = weight_ * PF00; projections_d[IDX(idO_ , 0, n_projected_points)] = weight_ * PF01;
   projections_d[IDX(idO_ , 1, n_projected_points)] = weight_ * PF10; projections_d[IDX(idO_ , 1, n_projected_points)] = weight_ * PF11;
   projections_d[IDX(idO_ , 2, n_projected_points)] = weight_ * PF20; projections_d[IDX(idO_ , 2, n_projected_points)] = weight_ * PF21;
  /*
  if (idI0_ == 0 && idI1_ == 1 && idI2_ == 2) {
	  printf("force \n%f %f %f\n%f %f %f\n%f %f %f\n\________________  end frame _________________\n\n\n",
	  f_int_nonePD_d[IDX(idI0_, 0, n_points)], f_int_nonePD_d[IDX(idI0_, 1, n_points)], f_int_nonePD_d[IDX(idI0_, 2, n_points)],
	  f_int_nonePD_d[IDX(idI1_, 0, n_points)], f_int_nonePD_d[IDX(idI1_, 1, n_points)], f_int_nonePD_d[IDX(idI1_, 2, n_points)],
	  f_int_nonePD_d[IDX(idI2_, 0, n_points)], f_int_nonePD_d[IDX(idI2_, 1, n_points)], f_int_nonePD_d[IDX(idI2_, 2, n_points)]);
  }
 */

}
////////////////////////////////////////
__device__
void project_surface_material(int tid, int n_points, int n_constraints, int n_projected_points, int nb_cells,
	ShapeOpScalar *points_d, ShapeOpScalar *projections_d, cuda_scalar *f_int_nonePD_d,
	int *ConstraintType_d, int *idO_d, ShapeOpScalar *rangeMin_d, ShapeOpScalar *rangeMax_d, ShapeOpScalar *Scalar1_d, ShapeOpScalar *weight_d, ShapeOpScalar *E_nonePD_d,
	int *idI_d, ShapeOpScalar *vectorx_d, ShapeOpScalar *matrix22_d, ShapeOpScalar *matrix33_d, cuda_scalar *A_d, const float miu, const float lambda, const float kappa)
{
	//printf("called material \n");

	cuda_scalar rangeMin_ = rangeMin_d[tid];
	cuda_scalar rangeMax_ = rangeMax_d[tid];
	cuda_scalar weight_ = weight_d[tid];
	//int idO_ = idO_d[tid] + blockIdx.x*3*n_projected_points;;

	cuda_scalar rest00_, rest01_,
		rest10_, rest11_;
	int idI0_, idI1_, idI2_, cell_shift = blockIdx.x*3*n_points;

	// Involved points (indices)
	idI0_ = idI_d[IDX(tid, 0, n_constraints)] + cell_shift;
	idI1_ = idI_d[IDX(tid, 1, n_constraints)] + cell_shift;
	idI2_ = idI_d[IDX(tid, 2, n_constraints)] + cell_shift;

	// ColMajor Order
	rest00_ = matrix22_d[IDX(tid, 0, n_constraints)];
	rest10_ = matrix22_d[IDX(tid, 1, n_constraints)];
	rest01_ = matrix22_d[IDX(tid, 2, n_constraints)];
	rest11_ = matrix22_d[IDX(tid, 3, n_constraints)];

	cuda_scalar edges00, edges01,
		edges10, edges11,
		edges20, edges21;

	cuda_scalar P00, P01,
		P10, P11,
		P20, P21;

	edges00 = points_d[IDX(idI1_, 0, n_points)] - points_d[IDX(idI0_, 0, n_points)];
	edges10 = points_d[IDX(idI1_, 1, n_points)] - points_d[IDX(idI0_, 1, n_points)];
	edges20 = points_d[IDX(idI1_, 2, n_points)] - points_d[IDX(idI0_, 2, n_points)];
	// edges.col(1)
	edges01 = points_d[IDX(idI2_, 0, n_points)] - points_d[IDX(idI0_, 0, n_points)];
	edges11 = points_d[IDX(idI2_, 1, n_points)] - points_d[IDX(idI0_, 1, n_points)];
	edges21 = points_d[IDX(idI2_, 2, n_points)] - points_d[IDX(idI0_, 2, n_points)];

	//debug to be deleled!!
	//edges00 += 0.01;

	P00 = edges00;
	P10 = edges10;
	P20 = edges20;
	Normalize_3(P00, P10, P20);

	cuda_scalar dot = DOT_3(edges01, edges11, edges21, P00, P10, P20);
	P01 = edges01 - dot * P00;
	P11 = edges11 - dot * P10;
	P21 = edges21 - dot * P20;
	Normalize_3(P01, P11, P21);
	// P.col(1)

	cuda_scalar tmp00, tmp01,
		tmp10, tmp11;
	cuda_scalar F00, F01,
		F10, F11;

	Matrix_Product_23_32(P00, P10, P20,
		P01, P11, P21,
		//
		edges00, edges01,
		edges10, edges11,
		edges20, edges21,
		//
		tmp00, tmp01,
		tmp10, tmp11);

	Matrix_Product_22_22(tmp00, tmp01,
		tmp10, tmp11,
		//
		rest00_, rest01_,
		rest10_, rest11_,
		//
		F00, F01,
		F10, F11);


	// [U SIG V] = SVD(F)
	// SVD
	cuda_scalar U00, U01,
		U10, U11;

	cuda_scalar SIG00, SIG11;

	cuda_scalar V00, V01,
		V10, V11;

	svd_joel22(F00, F01,
		F10, F11,
		//
		U00, U01,
		U10, U11,
		//
		SIG00, SIG11,
		//
		V00, V01,
		V10, V11);

	/*
	if( idI0_ == 0 && idI1_ == 1 && idI2_ == 2){
	printf("rest %f %f \nrest %f %f \n", rest00_, rest01_, rest10_, rest11_);
	printf("edge %f %f \nedge %f %f  \nedge %f %f\n", edges00, edges01, edges10, edges11, edges20, edges21);
	printf("P %f %f \nP %f %f \nP %f %f  \n", P00, P01, P10, P11, P20, P21);
	printf("F %f %f \nF %f %f \n", F00, F01, F10, F11);
	printf("U %f %f \nU %f %f \n\nV %f %f \nV %f %f \n---------------------------------\n", U00, U01, U10, U11, V00, V01, V10, V11);
	printf("singular %f %f \n", SIG00, SIG11);
	}
	*/

	cuda_scalar A = -0.5*A_d[tid];
	//none linear business
	E_nonePD_d[blockIdx.x*n_constraints + tid] = A_d[tid] * (f_tr(SIG00, miu, lambda, kappa) + f_tr(SIG11, miu, lambda, kappa)) / 2;
	//printf("e_non_pd %f %f\n", E_nonePD_d[tid], A_d[tid]);
	/*
	if (idI0_ == 0 && idI1_ == 1 && idI2_ == 2) {
	printf("singular %f %f \n", SIG00, SIG11);
	}
	*/
	//F = U*S*VT = PIOLA
	Matrix_Product_22_22(U00, U01,
		U10, U11,
		//
		f_prime_tr(SIG00, miu, lambda, kappa), (cuda_scalar)0.,
		(cuda_scalar)0., f_prime_tr(SIG11, miu, lambda, kappa),
		//
		tmp00, tmp01,
		tmp10, tmp11);

	Matrix_Product_22_22(tmp00, tmp01,
		tmp10, tmp11,
		//
		V00, V10,
		V01, V11,
		//
		F00, F01,
		F10, F11);
	/*
	if (idI0_ == 0 && idI1_ == 1 && idI2_ == 2) {
	printf("Piola %f %f \nPiola %f %f A %f \n", F00, F01, F10, F11, A_d[tid]);
	}
	*/
	// tmp =  Piola*rest_T
	Matrix_Product_22_22(F00, F01,
		F10, F11,
		rest00_, rest10_,
		rest01_, rest11_,
		tmp00, tmp01,
		tmp10, tmp11);
	/*
	if (idI0_ == 0 && idI1_ == 1 && idI2_ == 2) {
	printf("Piola*rest %f %f \nPiola*rest %f %f\n", tmp00, tmp01, tmp10, tmp11);
	}
	*/
	cuda_scalar F20, F21;
	//F = P*Piola*rest_T;
	Matrix_Product_32_22(P00, P01,
		P10, P11,
		P20, P21,
		tmp00, tmp01,
		tmp10, tmp11,
		F00, F01,
		F10, F11,
		F20, F21);
	/*
	if (idI0_ == 0 && idI1_ == 1 && idI2_ == 2) {
	printf("P*Piola*rest %f %f \nP*Piola*rest %f %f\n P*Piola*rest %f %f\n", A*F00, A*F01, A*F10, A*F11, A*F20, A*F21);
	}
	*/
	atomicAdd(f_int_nonePD_d + IDX(idI0_, 0, n_points), A*(-F00 - F01));
	atomicAdd(f_int_nonePD_d + IDX(idI0_, 1, n_points), A*(-F10 - F11));
	atomicAdd(f_int_nonePD_d + IDX(idI0_, 2, n_points), A*(-F20 - F21));

	atomicAdd(f_int_nonePD_d + IDX(idI1_, 0, n_points), A*F00);
	atomicAdd(f_int_nonePD_d + IDX(idI1_, 1, n_points), A*F10);
	atomicAdd(f_int_nonePD_d + IDX(idI1_, 2, n_points), A*F20);

	atomicAdd(f_int_nonePD_d + IDX(idI2_, 0, n_points), A*F01);
	atomicAdd(f_int_nonePD_d + IDX(idI2_, 1, n_points), A*F11);
	atomicAdd(f_int_nonePD_d + IDX(idI2_, 2, n_points), A*F21);


	/*
	if (idI0_ == 0 && idI1_ == 1 && idI2_ == 2) {
	printf("force \n%f %f %f\n%f %f %f\n%f %f %f\n\________________  end frame _________________\n\n\n",
	f_int_nonePD_d[IDX(idI0_, 0, n_points)], f_int_nonePD_d[IDX(idI0_, 1, n_points)], f_int_nonePD_d[IDX(idI0_, 2, n_points)],
	f_int_nonePD_d[IDX(idI1_, 0, n_points)], f_int_nonePD_d[IDX(idI1_, 1, n_points)], f_int_nonePD_d[IDX(idI1_, 2, n_points)],
	f_int_nonePD_d[IDX(idI2_, 0, n_points)], f_int_nonePD_d[IDX(idI2_, 1, n_points)], f_int_nonePD_d[IDX(idI2_, 2, n_points)]);
	}
	*/

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__device__ void nearest_neighbor_linear_d(double *points, double *colid_points, double *colid_normals, double *nearest_points, 
										  double *nearest_normals, double *center, int debut, int fin, int n_points, int id,
                                          const float threshold_rep, const float threshold_nonRep)
{
	float x = points[id             ] + center[0];
	float y = points[id + n_points  ] + center[1];
	float z = points[id + 2*n_points] + center[2];

	float d = 0.0f;
	float min = 9999999.0f;
	int argmin = -1;

    double threshold;
    if (threshold_rep >= threshold_nonRep)
        threshold = threshold_rep;
    else
        threshold = threshold_nonRep;

	for (int i = 3*debut; i < 3*fin; i += 3)
    {
		d = Magnitude_3(x - colid_points[i], y - colid_points[i + 1], z - colid_points[i + 2]);

		if (d < min && d <= threshold)
        {
			argmin = i;
			min = d;
		}
	}

    if (argmin > -1)
    {
        if (Norm_3(colid_normals[argmin], colid_normals[argmin + 1], colid_normals[argmin + 2]) >= 2.0)
        {
            threshold = threshold_rep;

            // 0.5 is to decode the info of repulsion
            nearest_points[IDX(id, 0, n_points)] = colid_points[argmin    ] + threshold*0.5*colid_normals[argmin    ] - center[0];
            nearest_points[IDX(id, 1, n_points)] = colid_points[argmin + 1] + threshold*0.5*colid_normals[argmin + 1] - center[1];
            nearest_points[IDX(id, 2, n_points)] = colid_points[argmin + 2] + threshold*0.5*colid_normals[argmin + 2] - center[2];

            nearest_normals[IDX(id, 0, n_points)] = colid_normals[argmin    ];
            nearest_normals[IDX(id, 1, n_points)] = colid_normals[argmin + 1];
            nearest_normals[IDX(id, 2, n_points)] = colid_normals[argmin + 2];
        }
        else
        {
            threshold = threshold_nonRep;

            if (min <= threshold)
            {
                nearest_points[IDX(id, 0, n_points)] = colid_points[argmin    ] + threshold*colid_normals[argmin    ] - center[0];
                nearest_points[IDX(id, 1, n_points)] = colid_points[argmin + 1] + threshold*colid_normals[argmin + 1] - center[1];
                nearest_points[IDX(id, 2, n_points)] = colid_points[argmin + 2] + threshold*colid_normals[argmin + 2] - center[2];

                nearest_normals[IDX(id, 0, n_points)] = colid_normals[argmin    ];
                nearest_normals[IDX(id, 1, n_points)] = colid_normals[argmin + 1];
                nearest_normals[IDX(id, 2, n_points)] = colid_normals[argmin + 2];
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////////
__device__
void project_collision_d(int n_points, int nb_cells, double *points, double *nearest_points, 
						double *nearest_normals, cuda_scalar *f_int_nonePD_d, double *E_nonePD, const int id, const int n_constraints,
                        const float weight_col_rep, const float threshold_rep, const float weight_col_nonRep, const float threshold_nonRep, const float beta_morse)
{
	float dif[3];

	dif[0] = points[IDX(id, 0, n_points)] - nearest_points[IDX(id, 0, n_points)];
	dif[1] = points[IDX(id, 1, n_points)] - nearest_points[IDX(id, 1, n_points)];
	dif[2] = points[IDX(id, 2, n_points)] - nearest_points[IDX(id, 2, n_points)];

	if (DOT_3(dif[0], dif[1], dif[2], nearest_normals[IDX(id, 0, n_points)], nearest_normals[IDX(id, 1, n_points)], nearest_normals[IDX(id, 2, n_points)]) < 0.)
    {	
        if (Norm_3(nearest_normals[IDX(id, 0, n_points)], nearest_normals[IDX(id, 1, n_points)], nearest_normals[IDX(id, 2, n_points)]) >= 2.0)
        {
            // 0.5 is to decode the info of repulsion
            
            /* Morse potential
            f_int_nonePD_d[IDX(id, 0, n_points)] -= 2.0f*weight_col_rep*beta_morse*(expf(2.0f*beta_morse*dif[0]) - expf(beta_morse*dif[0]));
            f_int_nonePD_d[IDX(id, 1, n_points)] -= 2.0f*weight_col_rep*beta_morse*(expf(2.0f*beta_morse*dif[1]) - expf(beta_morse*dif[1]));
            f_int_nonePD_d[IDX(id, 2, n_points)] -= 2.0f*weight_col_rep*beta_morse*(expf(2.0f*beta_morse*dif[2]) - expf(beta_morse*dif[2]));

            E_nonePD[blockIdx.x*n_constraints + threadIdx.x] += weight_col_rep*(expf(2.0f*beta_morse*dif[0]) - 2.0f*expf(beta_morse*dif[0])) +
                                                                weight_col_rep*(expf(2.0f*beta_morse*dif[1]) - 2.0f*expf(beta_morse*dif[1])) +
                                                                weight_col_rep*(expf(2.0f*beta_morse*dif[2]) - 2.0f*expf(beta_morse*dif[2]));
            //*/

            /* Custom Potential
            // "beta_morse" 40.0 is a good value
            dif[0] = points[IDX(id, 0, n_points)] - (nearest_points[IDX(id, 0, n_points)] + (beta_morse - threshold_rep)*0.5*nearest_normals[IDX(id, 0, n_points)]);
            dif[1] = points[IDX(id, 1, n_points)] - (nearest_points[IDX(id, 1, n_points)] + (beta_morse - threshold_rep)*0.5*nearest_normals[IDX(id, 1, n_points)]);
            dif[2] = points[IDX(id, 2, n_points)] - (nearest_points[IDX(id, 2, n_points)] + (beta_morse - threshold_rep)*0.5*nearest_normals[IDX(id, 2, n_points)]);

            f_int_nonePD_d[IDX(id, 0, n_points)] -= weight_col_rep*dif[0];
            f_int_nonePD_d[IDX(id, 1, n_points)] -= weight_col_rep*dif[1];
            f_int_nonePD_d[IDX(id, 2, n_points)] -= weight_col_rep*dif[2];

            E_nonePD[blockIdx.x*n_constraints + threadIdx.x] += weight_col_rep*(dif[0] * dif[0] + dif[1] * dif[1] + dif[2] * dif[2]) / 2.0f;
            //*/

            ///* Classic Potential
            f_int_nonePD_d[IDX(id, 0, n_points)] -= weight_col_rep*dif[0];
            f_int_nonePD_d[IDX(id, 1, n_points)] -= weight_col_rep*dif[1];
            f_int_nonePD_d[IDX(id, 2, n_points)] -= weight_col_rep*dif[2];

            E_nonePD[blockIdx.x*n_constraints + threadIdx.x] += weight_col_rep*(dif[0] * dif[0] + dif[1] * dif[1] + dif[2] * dif[2]) / 2.0f;
            //*/
        }
        else
        {
            f_int_nonePD_d[IDX(id, 0, n_points)] -= weight_col_nonRep*dif[0];
            f_int_nonePD_d[IDX(id, 1, n_points)] -= weight_col_nonRep*dif[1];
            f_int_nonePD_d[IDX(id, 2, n_points)] -= weight_col_nonRep*dif[2];

            E_nonePD[blockIdx.x*n_constraints + threadIdx.x] += weight_col_nonRep*(dif[0] * dif[0] + dif[1] * dif[1] + dif[2] * dif[2]) / 2.0f;
        }
	}
}

}
}
///////////////////////////////////////////////////////////////////////////////
#endif // PROJECTIONS_GPU_H
///////////////////////////////////////////////////////////////////////////////
