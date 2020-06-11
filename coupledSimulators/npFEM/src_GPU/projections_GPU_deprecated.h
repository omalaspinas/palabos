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
#include "projections_GPU_MATH.h"
///////////////////////////////////////////////////////////////////////////////
namespace plb {
namespace npfem {

__device__ __forceinline__
void project_Area( int tid, int n_points, int n_constraints, int n_projected_points, int nb_cells,
                   ShapeOpScalar *points_d, ShapeOpScalar *projections_d, ShapeOpScalar *f_int_nonePD_d,
                   int *ConstraintType_d, int *idO_d, ShapeOpScalar *rangeMin_d, ShapeOpScalar *rangeMax_d, ShapeOpScalar *Scalar1_d, ShapeOpScalar *weight_d, ShapeOpScalar *E_nonePD_d,
                   int *idI_d, ShapeOpScalar *vectorx_d, ShapeOpScalar *matrix22_d, ShapeOpScalar *matrix33_d )
{
  ShapeOpScalar rangeMin_ = rangeMin_d[tid];
  ShapeOpScalar rangeMax_ = rangeMax_d[tid];
  ShapeOpScalar weight_   = weight_d[tid];
  int idO_                = idO_d[tid];

  ShapeOpScalar rest00_, rest01_,
                rest10_, rest11_;
  int idI0_, idI1_, idI2_;

  // Involved points (indices)
  idI0_ = idI_d[IDX(0, tid, 4)];
  idI1_ = idI_d[IDX(1, tid, 4)];
  idI2_ = idI_d[IDX(2, tid, 4)];

  // ColMajor Order
  rest00_ = matrix22_d[IDX(0, tid, 4)];
  rest10_ = matrix22_d[IDX(1, tid, 4)];
  rest01_ = matrix22_d[IDX(2, tid, 4)];
  rest11_ = matrix22_d[IDX(3, tid, 4)];

  ShapeOpScalar edges00, edges01,
                edges10, edges11,
                edges20, edges21;
  ShapeOpScalar P00, P01,
                P10, P11,
                P20, P21;
  ShapeOpScalar PT00, PT01, PT02,
                PT10, PT11, PT12;

  // edges.col(0)
  
  edges00 = points_d[IDX(0, idI1_, 3)] - points_d[IDX(0, idI0_, 3)];
  edges10 = points_d[IDX(1, idI1_, 3)] - points_d[IDX(1, idI0_, 3)];
  edges20 = points_d[IDX(2, idI1_, 3)] - points_d[IDX(2, idI0_, 3)];
  // edges.col(1)
  edges01 = points_d[IDX(0, idI2_, 3)] - points_d[IDX(0, idI0_, 3)];
  edges11 = points_d[IDX(1, idI2_, 3)] - points_d[IDX(1, idI0_, 3)];
  edges21 = points_d[IDX(2, idI2_, 3)] - points_d[IDX(2, idI0_, 3)];
  

  ShapeOpScalar Vec3_0_0, Vec3_0_1, Vec3_0_2;
  Vec3_0_0 = edges00;
  Vec3_0_1 = edges10;
  Vec3_0_2 = edges20;
  Normalize_3(Vec3_0_0, Vec3_0_1, Vec3_0_2);
  // P.col(0)
  P00 = Vec3_0_0;
  P10 = Vec3_0_1;
  P20 = Vec3_0_2;

  ShapeOpScalar Vec3_1_0, Vec3_1_1, Vec3_1_2;
  Vec3_1_0 = edges01;
  Vec3_1_1 = edges11;
  Vec3_1_2 = edges21;
  ShapeOpScalar Vec3_2_0, Vec3_2_1, Vec3_2_2;
  Vec3_2_0 = P00;
  Vec3_2_1 = P10;
  Vec3_2_2 = P20;
  
  ShapeOpScalar dot = DOT_3(Vec3_1_0, Vec3_1_1, Vec3_1_2, Vec3_2_0, Vec3_2_1, Vec3_2_2);
  Vec3_0_0 = edges01 - dot * P00;
  Vec3_0_1 = edges11 - dot * P10;
  Vec3_0_2 = edges21 - dot * P20;
  Normalize_3(Vec3_0_0, Vec3_0_1, Vec3_0_2);
  // P.col(1)
  P01 = Vec3_0_0;
  P11 = Vec3_0_1;
  P21 = Vec3_0_2;

  ShapeOpScalar tmp00, tmp01,
                tmp10, tmp11;
  ShapeOpScalar F00, F01,
                F10, F11;
  
  Matrix_Transpose_32_23(P00, P01,
                         P10, P11,
                         P20, P21,
                         //
                         PT00, PT01, PT02,
                         PT10, PT11, PT12);

  Matrix_Product_23_32(PT00, PT01, PT02,
                       PT10, PT11, PT12,
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
  ShapeOpScalar U00, U01,
                U10, U11;
  ShapeOpScalar SIG00, SIG01,
                SIG10, SIG11;  
  ShapeOpScalar V00, V01,
                V10, V11;

  SVD_2x2(F00, F01,
          F10, F11,
          //
          U00, U01,
          U10, U11,
          //
          SIG00, SIG01,
          SIG10, SIG11,
          //
          V00, V01,
          V10, V11);

  ShapeOpScalar S0 = SIG00, S1 = SIG11; 
  ShapeOpScalar d0 = 0., d1 = 0.;
  for (int i = 0; i < SHAPEOP_INNER_ITERATIONS; ++i)
  {
    ShapeOpScalar v = S0 * S1;
    ShapeOpScalar f = v - CLAMP(v, rangeMin_, rangeMax_);
    ShapeOpScalar g0 = S1, g1 = S0;
    ShapeOpScalar dot_gd = DOT_2(g0, g1, d0, d1);
    ShapeOpScalar dot_gg = DOT_2(g0, g1, g0, g1);
    ShapeOpScalar frac = -((f - dot_gd) / dot_gg);
    d0 = frac * g0; d1 = frac * g1;
    S0 = SIG00 + d0; S1 = SIG11 + d1;
  }
  SIG00 = S0; SIG01 = 0.;
  SIG10 = 0.; SIG11 = S1;

  ShapeOpScalar VT00, VT01,
                VT10, VT11;
  Matrix_Transpose_22(V00, V01,
                      V10, V11,
                      //
                      VT00, VT01,
                      VT10, VT11);

  Matrix_Product_22_22(U00, U01,
                       U10, U11,
                       //
                       SIG00, SIG01,
                       SIG10, SIG11,
                       //
                       tmp00, tmp01,
                       tmp10, tmp11);

  Matrix_Product_22_22(tmp00, tmp01,
                       tmp10, tmp11,
                       //
                       VT00, VT01,
                       VT10, VT11,
                       //
                       F00, F01,
                       F10, F11);

  ShapeOpScalar PF00, PF01,
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

  projections_d[IDX(0, idO_, 3)] = weight_ * PF00; projections_d[IDX(0, idO_ + 1, 3)] = weight_ * PF01;
  projections_d[IDX(1, idO_, 3)] = weight_ * PF10; projections_d[IDX(1, idO_ + 1, 3)] = weight_ * PF11;
  projections_d[IDX(2, idO_, 3)] = weight_ * PF20; projections_d[IDX(2, idO_ + 1, 3)] = weight_ * PF21;
  
  if (tid == 786) {
	  printf("projection load %f %f %f | %d \n",  points_d[IDX(0, idI2_, 3)], points_d[IDX(1, idI2_, 3)], points_d[IDX(2, idI2_, 3)], idI2_);
	  //printf("projection %f %f %f\n", projections_d[IDX(0, idO_, 3)], projections_d[IDX(1, idO_, 3)], projections_d[IDX(2, idO_, 3)]);
	  //printf("ido_ %d  idI0_ %d \n", idO_, idI0_);
  }
  
}
///////////////////////////////////////////////////////////////////////////////
__device__ __forceinline__
void project_Volume( int tid, int n_points, int n_constraints, int n_projected_points, int nb_cells,
                     ShapeOpScalar *points_d, ShapeOpScalar *projections_d, ShapeOpScalar *f_int_nonePD_d,
                     int *ConstraintType_d, int *idO_d, ShapeOpScalar *rangeMin_d, ShapeOpScalar *rangeMax_d, ShapeOpScalar *Scalar1_d, ShapeOpScalar *weight_d, ShapeOpScalar *E_nonePD_d,
                     int *idI_d, ShapeOpScalar *vectorx_d, ShapeOpScalar *matrix22_d, ShapeOpScalar *matrix33_d )
{
  ShapeOpScalar rangeMin_ = rangeMin_d[tid];
  ShapeOpScalar rangeMax_ = rangeMax_d[tid];
  ShapeOpScalar weight_   = weight_d[tid];
  int idO_                = idO_d[tid];

  int idI0_, idI1_, idI2_, idI3_;
  ShapeOpScalar rest00_, rest01_, rest02_,
                rest10_, rest11_, rest12_,
                rest20_, rest21_, rest22_;

  // ColMajor Order
  rest00_ = matrix33_d[IDX(0, tid, 9)];
  rest10_ = matrix33_d[IDX(1, tid, 9)];
  rest20_ = matrix33_d[IDX(2, tid, 9)];
  rest01_ = matrix33_d[IDX(3, tid, 9)];
  rest11_ = matrix33_d[IDX(4, tid, 9)];
  rest21_ = matrix33_d[IDX(5, tid, 9)];
  rest02_ = matrix33_d[IDX(6, tid, 9)];
  rest12_ = matrix33_d[IDX(7, tid, 9)];
  rest22_ = matrix33_d[IDX(8, tid, 9)];

  // Involved points (indices)
  idI0_ = idI_d[IDX(0, tid, 4)];
  idI1_ = idI_d[IDX(1, tid, 4)];
  idI2_ = idI_d[IDX(2, tid, 4)];
  idI3_ = idI_d[IDX(3, tid, 4)];

  // Matrix 3x3
  ShapeOpScalar edges00, edges01, edges02,
                edges10, edges11, edges12,
                edges20, edges21, edges22;

  // edges.col(0)
  edges00 = points_d[IDX(0, idI1_, 3)] - points_d[IDX(0, idI0_, 3)];
  edges10 = points_d[IDX(1, idI1_, 3)] - points_d[IDX(1, idI0_, 3)];
  edges20 = points_d[IDX(2, idI1_, 3)] - points_d[IDX(2, idI0_, 3)];
  // edges.col(1)
  edges01 = points_d[IDX(0, idI2_, 3)] - points_d[IDX(0, idI0_, 3)];
  edges11 = points_d[IDX(1, idI2_, 3)] - points_d[IDX(1, idI0_, 3)];
  edges21 = points_d[IDX(2, idI2_, 3)] - points_d[IDX(2, idI0_, 3)];
  // edges.col(2)
  edges02 = points_d[IDX(0, idI3_, 3)] - points_d[IDX(0, idI0_, 3)];
  edges12 = points_d[IDX(1, idI3_, 3)] - points_d[IDX(1, idI0_, 3)];
  edges22 = points_d[IDX(2, idI3_, 3)] - points_d[IDX(2, idI0_, 3)];
  
  ShapeOpScalar F00, F01, F02,
                F10, F11, F12,
                F20, F21, F22;
  Matrix_Product_33_33(edges00, edges01, edges02,
                       edges10, edges11, edges12,
                       edges20, edges21, edges22,
                       //
                       rest00_, rest01_, rest02_,
                       rest10_, rest11_, rest12_,
                       rest20_, rest21_, rest22_,
                       //
                       F00, F01, F02,
                       F10, F11, F12,
                       F20, F21, F22);

  // SVD
  ShapeOpScalar U00, U01, U02,
                U10, U11, U12,
                U20, U21, U22;
  ShapeOpScalar SIG00, SIG01, SIG02,
                SIG10, SIG11, SIG12,
                SIG20, SIG21, SIG22;
  ShapeOpScalar V00, V01, V02,
                V10, V11, V12,
                V20, V21, V22;

  SVD_3x3(F00, F01, F02,
          F10, F11, F12,
          F20, F21, F22,
          //
          U00, U01, U02,
          U10, U11, U12,
          U20, U21, U22,
          //
          SIG00, SIG01, SIG02,
          SIG10, SIG11, SIG12,
          SIG20, SIG21, SIG22,
          //
          V00, V01, V02,
          V10, V11, V12,
          V20, V21, V22);

  ShapeOpScalar S0 = SIG00, S1 = SIG11, S2 = SIG22;
  ShapeOpScalar d0 = 0., d1 = 0., d2 = 0.;
  for (int i = 0; i < SHAPEOP_INNER_ITERATIONS; ++i)
  {
    ShapeOpScalar v = S0 * S1 * S2;
    ShapeOpScalar f = v - CLAMP(v, rangeMin_, rangeMax_);
    ShapeOpScalar g0 = S1 * S2, g1 = S0 * S2, g2 = S0 * S1;
    ShapeOpScalar dot_gd = DOT_3(g0, g1, g2, d0, d1, d2);
    ShapeOpScalar dot_gg = DOT_3(g0, g1, g2, g0, g1, g2);
    ShapeOpScalar frac = -((f - dot_gd) / dot_gg);
    d0 = frac * g0; d1 = frac * g1; d2 = frac * g2;
    S0 = SIG00 + d0; S1 = SIG11 + d1; S2 = SIG22 + d2;
  }
  ShapeOpScalar det_U = Determinant_3(U00, U01, U02,
                                      U10, U11, U12,
                                      U20, U21, U22);
  ShapeOpScalar det_V = Determinant_3(V00, V01, V02,
                                      V10, V11, V12,
                                      V20, V21, V22);

  if ( det_U * det_V < 0. ) S2 = -S2;
  SIG00 = S0; SIG01 = 0.; SIG02 = 0.;
  SIG10 = 0.; SIG11 = S1; SIG12 = 0.;
  SIG20 = 0.; SIG21 = 0.; SIG22 = S2;

  ShapeOpScalar VT00, VT01, VT02,
                VT10, VT11, VT12,
                VT20, VT21, VT22;
  Matrix_Transpose_33(V00, V01, V02,
                      V10, V11, V12,
                      V20, V21, V22,
                      //
                      VT00, VT01, VT02,
                      VT10, VT11, VT12,
                      VT20, VT21, VT22);

  ShapeOpScalar tmp00, tmp01, tmp02,
                tmp10, tmp11, tmp12,
                tmp20, tmp21, tmp22;
  Matrix_Product_33_33(U00, U01, U02,
                       U10, U11, U12,
                       U20, U21, U22,
                       //
                       SIG00, SIG01, SIG02,
                       SIG10, SIG11, SIG12,
                       SIG20, SIG21, SIG22,
                       //
                       tmp00, tmp01, tmp02,
                       tmp10, tmp11, tmp12,
                       tmp20, tmp21, tmp22);

  Matrix_Product_33_33(tmp00, tmp01, tmp02,
                       tmp10, tmp11, tmp12,
                       tmp20, tmp21, tmp22,
                       //
                       VT00, VT01, VT02,
                       VT10, VT11, VT12,
                       VT20, VT21, VT22,
                       //
                       F00, F01, F02,
                       F10, F11, F12,
                       F20, F21, F22);

  projections_d[IDX(0, idO_, 3)] = weight_ * F00; projections_d[IDX(0, idO_ + 1, 3)] = weight_ * F01; projections_d[IDX(0, idO_ + 2, 3)] = weight_ * F02;
  projections_d[IDX(1, idO_, 3)] = weight_ * F10; projections_d[IDX(1, idO_ + 1, 3)] = weight_ * F11; projections_d[IDX(1, idO_ + 2, 3)] = weight_ * F12;
  projections_d[IDX(2, idO_, 3)] = weight_ * F20; projections_d[IDX(2, idO_ + 1, 3)] = weight_ * F21; projections_d[IDX(2, idO_ + 2, 3)] = weight_ * F22;
}
///////////////////////////////////////////////////////////////////////////////
__device__ __forceinline__
void project_Bending( int tid, int n_points, int n_constraints, int n_projected_points, int nb_cells,
                      ShapeOpScalar *points_d, ShapeOpScalar *projections_d, ShapeOpScalar *f_int_nonePD_d,
                      int *ConstraintType_d, int *idO_d, ShapeOpScalar *rangeMin_d, ShapeOpScalar *rangeMax_d, ShapeOpScalar *Scalar1_d, ShapeOpScalar *weight_d, ShapeOpScalar *E_nonePD_d,
                      int *idI_d, ShapeOpScalar *vectorx_d, ShapeOpScalar *matrix22_d, ShapeOpScalar *matrix33_d )
{
  ShapeOpScalar rangeMin_ = rangeMin_d[tid];
  ShapeOpScalar rangeMax_ = rangeMax_d[tid];
  ShapeOpScalar weight_   = weight_d[tid];
  int idO_                = idO_d[tid];
  ShapeOpScalar n_        = Scalar1_d[tid];
  
  ShapeOpScalar w0_, w1_, w2_, w3_;
  int idI0_, idI1_, idI2_, idI3_;

  w0_ = vectorx_d[IDX(0, tid, 4)];
  w1_ = vectorx_d[IDX(1, tid, 4)];
  w2_ = vectorx_d[IDX(2, tid, 4)];
  w3_ = vectorx_d[IDX(3, tid, 4)];

  // Involved points (indices)
  idI0_ = idI_d[IDX(0, tid, 4)];
  idI1_ = idI_d[IDX(1, tid, 4)];
  idI2_ = idI_d[IDX(2, tid, 4)];
  idI3_ = idI_d[IDX(3, tid, 4)];

  ShapeOpScalar e0, e1, e2;
  if (n_ > 1e-6)
  {
	
	  e0 = w0_ * points_d[IDX(0, idI0_, 3)] + w1_ * points_d[IDX(0, idI1_, 3)] + w2_ * points_d[IDX(0, idI2_, 3)] + w3_ * points_d[IDX(0, idI3_, 3)];
	  e1 = w0_ * points_d[IDX(1, idI0_, 3)] + w1_ * points_d[IDX(1, idI1_, 3)] + w2_ * points_d[IDX(1, idI2_, 3)] + w3_ * points_d[IDX(1, idI3_, 3)];
	  e2 = w0_ * points_d[IDX(2, idI0_, 3)] + w1_ * points_d[IDX(2, idI1_, 3)] + w2_ * points_d[IDX(2, idI2_, 3)] + w3_ * points_d[IDX(2, idI3_, 3)];

    ShapeOpScalar l = Norm_3(e0, e1, e2);

    if (l > 1e-6)
    {
      e0 /= l;
      e1 /= l;
      e2 /= l;
      l = n_ * CLAMP(l / n_, rangeMin_, rangeMax_);
      e0 *= l;
      e1 *= l;
      e2 *= l;
    }
  }

  projections_d[IDX(0, idO_, 3)] = weight_ * e0;
  projections_d[IDX(1, idO_, 3)] = weight_ * e1;
  projections_d[IDX(2, idO_, 3)] = weight_ * e2;
/*
  if (tid == 615 || tid == 622 || tid == 624) {
	  printf("projection bending %f %f %f\n", projections_d[IDX(0, idO_, 3)], projections_d[IDX(1, idO_, 3)], projections_d[IDX(2, idO_, 3)]);
	  printf("ido_ %d  idI0_ %d \n", idO_, idI0_);
  }
  */
}
///////////////////////////////////////////////////////////////////////////////
__device__ __forceinline__
void project_TriangleStrainLimiting( int tid, int n_points, int n_constraints, int n_projected_points, int nb_cells,
                                     ShapeOpScalar *points_d, ShapeOpScalar *projections_d, ShapeOpScalar *f_int_nonePD_d,
                                     int *ConstraintType_d, int *idO_d, ShapeOpScalar *rangeMin_d, ShapeOpScalar *rangeMax_d, ShapeOpScalar *Scalar1_d, ShapeOpScalar *weight_d, ShapeOpScalar *E_nonePD_d,
                                     int *idI_d, ShapeOpScalar *vectorx_d, ShapeOpScalar *matrix22_d, ShapeOpScalar *matrix33_d )
{
  ShapeOpScalar rangeMin_ = rangeMin_d[tid];
  ShapeOpScalar rangeMax_ = rangeMax_d[tid];
  ShapeOpScalar weight_   = weight_d[tid];
  int idO_                = idO_d[tid];

  ShapeOpScalar rest00_, rest01_,
                rest10_, rest11_;
  int idI0_, idI1_, idI2_;

  // Involved points (indices)
  idI0_ = idI_d[IDX(0, tid, 4)];
  idI1_ = idI_d[IDX(1, tid, 4)];
  idI2_ = idI_d[IDX(2, tid, 4)];

  // ColMajor Order
  rest00_ = matrix22_d[IDX(0, tid, 4)];
  rest10_ = matrix22_d[IDX(1, tid, 4)];
  rest01_ = matrix22_d[IDX(2, tid, 4)];
  rest11_ = matrix22_d[IDX(3, tid, 4)];

  ShapeOpScalar edges00, edges01,
                edges10, edges11,
                edges20, edges21;
  ShapeOpScalar P00, P01,
                P10, P11,
                P20, P21;
  ShapeOpScalar PT00, PT01, PT02,
                PT10, PT11, PT12;

  // edges.col(0)
  edges00 = points_d[IDX(0, idI1_, 3)] - points_d[IDX(0, idI0_, 3)];
  edges10 = points_d[IDX(1, idI1_, 3)] - points_d[IDX(1, idI0_, 3)];
  edges20 = points_d[IDX(2, idI1_, 3)] - points_d[IDX(2, idI0_, 3)];
  // edges.col(1)
  edges01 = points_d[IDX(0, idI2_, 3)] - points_d[IDX(0, idI0_, 3)];
  edges11 = points_d[IDX(1, idI2_, 3)] - points_d[IDX(1, idI0_, 3)];
  edges21 = points_d[IDX(2, idI2_, 3)] - points_d[IDX(2, idI0_, 3)];

  ShapeOpScalar Vec3_0_0, Vec3_0_1, Vec3_0_2;
  Vec3_0_0 = edges00;
  Vec3_0_1 = edges10;
  Vec3_0_2 = edges20;
  Normalize_3(Vec3_0_0, Vec3_0_1, Vec3_0_2);
  // P.col(0)
  P00 = Vec3_0_0;
  P10 = Vec3_0_1;
  P20 = Vec3_0_2;

  ShapeOpScalar Vec3_1_0, Vec3_1_1, Vec3_1_2;
  Vec3_1_0 = edges01;
  Vec3_1_1 = edges11;
  Vec3_1_2 = edges21;
  ShapeOpScalar Vec3_2_0, Vec3_2_1, Vec3_2_2;
  Vec3_2_0 = P00;
  Vec3_2_1 = P10;
  Vec3_2_2 = P20;
  
  ShapeOpScalar dot = DOT_3(Vec3_1_0, Vec3_1_1, Vec3_1_2, Vec3_2_0, Vec3_2_1, Vec3_2_2);
  Vec3_0_0 = edges01 - dot * P00;
  Vec3_0_1 = edges11 - dot * P10;
  Vec3_0_2 = edges21 - dot * P20;
  Normalize_3(Vec3_0_0, Vec3_0_1, Vec3_0_2);
  // P.col(1)
  P01 = Vec3_0_0;
  P11 = Vec3_0_1;
  P21 = Vec3_0_2;

  ShapeOpScalar tmp00, tmp01,
                tmp10, tmp11;
  ShapeOpScalar F00, F01,
                F10, F11;
  
  Matrix_Transpose_32_23(P00, P01,
                         P10, P11,
                         P20, P21,
                         //
                         PT00, PT01, PT02,
                         PT10, PT11, PT12);

  Matrix_Product_23_32(PT00, PT01, PT02,
                       PT10, PT11, PT12,
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
  ShapeOpScalar U00, U01,
                U10, U11;
  ShapeOpScalar SIG00, SIG01,
                SIG10, SIG11;  
  ShapeOpScalar V00, V01,
                V10, V11;

  SVD_2x2(F00, F01,
          F10, F11,
          //
          U00, U01,
          U10, U11,
          //
          SIG00, SIG01,
          SIG10, SIG11,
          //
          V00, V01,
          V10, V11);

  // Strain Limiting
  SIG00 = CLAMP(SIG00, rangeMin_, rangeMax_);
  SIG11 = CLAMP(SIG11, rangeMin_, rangeMax_);

  ShapeOpScalar VT00, VT01,
                VT10, VT11;
  Matrix_Transpose_22(V00, V01,
                      V10, V11,
                      //
                      VT00, VT01,
                      VT10, VT11);

  Matrix_Product_22_22(U00, U01,
                       U10, U11,
                       //
                       SIG00, SIG01,
                       SIG10, SIG11,
                       //
                       tmp00, tmp01,
                       tmp10, tmp11);

  Matrix_Product_22_22(tmp00, tmp01,
                       tmp10, tmp11,
                       //
                       VT00, VT01,
                       VT10, VT11,
                       //
                       F00, F01,
                       F10, F11);

  ShapeOpScalar PF00, PF01,
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

  projections_d[IDX(0, idO_, 3)] = weight_ * PF00; projections_d[IDX(0, idO_ + 1, 3)] = weight_ * PF01;
  projections_d[IDX(1, idO_, 3)] = weight_ * PF10; projections_d[IDX(1, idO_ + 1, 3)] = weight_ * PF11;
  projections_d[IDX(2, idO_, 3)] = weight_ * PF20; projections_d[IDX(2, idO_ + 1, 3)] = weight_ * PF21;
}

}
}
///////////////////////////////////////////////////////////////////////////////
#endif // PROJECTIONS_GPU_H
///////////////////////////////////////////////////////////////////////////////
