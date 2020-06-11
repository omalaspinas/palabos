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
#ifndef PROJECTIONS_GPU_MATH_H
#define PROJECTIONS_GPU_MATH_H
///////////////////////////////////////////////////////////////////////////////
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

#include "common.h"
#include "svd3_cuda.h"
///////////////////////////////////////////////////////////////////////////////
// Zero-based indexing for raw pointer arrays, where ld is the number of rows
// for ColMajor matrices
#define IDX(i,j,ld)                     ( ( (j) * (ld) ) + (i) )
#define SHAPEOP_INNER_ITERATIONS        4
#define DOT_2(a0, a1,     b0, b1)		    ((a0)*(b0) + (a1)*(b1))
#define DOT_3(a0, a1, a2, b0, b1, b2)   ((a0)*(b0) + (a1)*(b1) + (a2)*(b2))
#define SIGN(a) 			                  ((a)<0?-1:1)
#define CLAMP(a, l, h)                  (((a)>(h))?(h):(((a)<(l))?(l):(a)))

namespace plb {
namespace npfem {
///////////////////////////////////////////////////////////////////////////////
template <typename T>
__device__ __forceinline__
ShapeOpScalar Magnitude_2(T &x0,T &x1)
{
  return sqrt(DOT_2(x0, x1, x0, x1));
}
////////////////////////////////////////////////////////////////////////////////
template <typename T>
__device__ __forceinline__
void cross(T const x0, T const x1, T const x2, T const y0, T const y1, T const y2, T &out0, T &out1, T &out2)
{
	out0 = x1*y2 - x2*y1;
	out1 = x2*y0 - x0*y2;
	out2 = x0*y1 - x1*y0;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
__device__ __forceinline__
ShapeOpScalar cross_norm(float &x0, float &x1, float &x2, float &y0, float &y1, float &y2)
{
	float c0 = x1*y2 - x2*y1;
	float sum = c0*c0;
	c0 = x2*y0 - x0*y2;
	sum += c0*c0;
	c0 = x0*y1 - x1*y0;
	sum += c0*c0;

	return sqrt(sum);
}
////////////////////////////////////////////////////////////////////////////////
template <typename T>
__device__ __forceinline__
T Magnitude_3(T &&x0, T &&x1, T &&x2)
{
	return sqrt(DOT_3(x0, x1, x2, x0, x1, x2));
}
///////////////////////////////////////////////////////////////////////////////
template <typename T>
__device__ __forceinline__
T Magnitude_3(T &x0, T &x1, T &x2)
{
  return sqrt(DOT_3(x0, x1, x2, x0, x1, x2));
}
///////////////////////////////////////////////////////////////////////////////
template <typename T>
__device__ __forceinline__
ShapeOpScalar Norm_2(T &x0, T &x1)
{
  return sqrt(DOT_2(x0, x1, x0, x1));
}

template <typename T>
__device__ __forceinline__
ShapeOpScalar Norm_3(T &x0, T &x1,T &x2)
{
  return sqrt(DOT_3(x0, x1, x2, x0, x1, x2));
}
///////////////////////////////////////////////////////////////////////////////
template <typename T>
__device__ __forceinline__
void Normalize_2(T &x0, T &x1)
{
 T m = Magnitude_2(x0, x1);
 T inv_m = 1 / m;
  x0 *= inv_m;
  x1 *= inv_m;
}
///////////////////////////////////////////////////////////////////////////////
template <class T> 
__device__ __forceinline__ T Normalize_3(T &x0, T &x1, T &x2)
{
  T m = Magnitude_3(x0, x1, x2);
  T inv_m = 1 / m;
  x0 *= inv_m;
  x1 *= inv_m;
  x2 *= inv_m;
  return m;
}
///////////////////////////////////////////////////////////////////////////////
template <typename T>
__device__ __forceinline__
ShapeOpScalar Determinant_3(T &a,T &b,T &c,
                           T &d,T &e,T &f,
                           T &g,T &h,T &i)
{
  return a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h;
}
template <typename T>
__device__ 
void Matrix_invers3_3X3(const T m[9], T m_inv[9]) {
	T det = m[0] * m[4] * m[8] + m[1] * m[5] * m[6] + m[2] * m[3] * m[7] - m[6] * m[4] * m[2] - m[7] * m[5] * m[0] - m[8] * m[3] * m[1];
	__syncthreads();
	if (threadIdx.x == 0)m_inv[0] =  (m[4] * m[8] - m[5] * m[7]) / det;
	if (threadIdx.x == 1)m_inv[1] = -(m[1] * m[8] - m[2] * m[7]) / det;
	if (threadIdx.x == 2)m_inv[2] =  (m[1] * m[5] - m[2] * m[4]) / det;

	if (threadIdx.x == 3)m_inv[3] = -(m[3] * m[8] - m[5] * m[6]) / det;
	if (threadIdx.x == 4)m_inv[4] =  (m[0] * m[8] - m[2] * m[6]) / det;
	if (threadIdx.x == 5)m_inv[5] = -(m[0] * m[5] - m[2] * m[3]) / det;

	if (threadIdx.x == 6)m_inv[6] =  (m[3] * m[7] - m[4] * m[6]) / det;
	if (threadIdx.x == 7)m_inv[7] = -(m[0] * m[7] - m[1] * m[6]) / det;
	if (threadIdx.x == 8)m_inv[8] =  (m[0] * m[4] - m[1] * m[3]) / det;
}
template <typename T>
__device__ void vect_matrix_mult_3X3_par(const T A[9], const T x[3],T b[3]) {
	__shared__ T prod[9];

	if (threadIdx.x < 9) {
		prod[threadIdx.x] = A[threadIdx.x]*x[threadIdx.x%3];
	}
	__syncthreads();
	if (threadIdx.x < 3){
		b[threadIdx.x] = prod[3*threadIdx.x] +  prod[3*threadIdx.x + 1] + prod[3*threadIdx.x + 2];
	}
}

///////////////////////////////////////////////////////////////////////////////
// A:nx*ny  || B:ny*nz  || R:nx*nz
template <typename T>
__device__ __forceinline__
void Matrix_Product_22_22(const T A00, const T A01,
						  const T A10, const T A11,
                          //
						 const T B00, const T B01,
						 const T B10, const T B11,
                          //
                         T &R00,T &R01,
                         T &R10,T &R11)
{
  R00 = A00*B00 + A01*B10;
  R01 = A00*B01 + A01*B11;
  
  R10 = A10*B00 + A11*B10;
  R11 = A10*B01 + A11*B11;
}

//R = BA
template <typename T>
__device__ __forceinline__
void Matrix_Product_33_33(const T A00, const T A01, const T A02,
						  const T A10, const T A11, const T A12,
						  const T A20, const T A21, const T A22,
                          //
						  const T B00, const T B01, const T &B02,
						  const T B10, const T B11, const T &B12,
						  const T B20, const T B21, const T &B22,
                          //
                          T &R00,T &R01,T &R02,
                          T &R10,T &R11,T &R12,
                          T &R20,T &R21,T &R22)
{
  R00 = A00*B00 + A01*B10 + A02*B20;
  R01 = A00*B01 + A01*B11 + A02*B21;
  R02 = A00*B02 + A01*B12 + A02*B22;

  R10 = A10*B00 + A11*B10 + A12*B20;
  R11 = A10*B01 + A11*B11 + A12*B21;
  R12 = A10*B02 + A11*B12 + A12*B22;

  R20 = A20*B00 + A21*B10 + A22*B20;
  R21 = A20*B01 + A21*B11 + A22*B21;
  R22 = A20*B02 + A21*B12 + A22*B22;
}
template <typename T>
__device__ __forceinline__
void Matrix_Product_23_32(T &A00,T &A01,T &A02,
                         T &A10,T &A11,T &A12,
                          //
                         T &B00,T &B01,
                         T &B10,T &B11,
                         T &B20,T &B21,
                          //
                         T &R00,T &R01,
                         T &R10,T &R11)
{
  R00 = A00*B00 + A01*B10 + A02*B20;
  R01 = A00*B01 + A01*B11 + A02*B21;
  
  R10 = A10*B00 + A11*B10 + A12*B20;
  R11 = A10*B01 + A11*B11 + A12*B21;
}

template <typename T>
__device__ __forceinline__
void Matrix_Product_32_23(T &A00,T &A01,
                         T &A10,T &A11,
                         T &A20,T &A21,
                          //
                         T &B00,T &B01,T &B02,
                         T &B10,T &B11,T &B12,
                          //
                         T &R00,T &R01,T &R02,
                         T &R10,T &R11,T &R12,
                         T &R20,T &R21,T &R22)
{
  R00 = A00*B00 + A01*B10;
  R01 = A00*B01 + A01*B11;
  R02 = A00*B02 + A01*B12;

  R10 = A10*B00 + A11*B10;
  R11 = A10*B01 + A11*B11;
  R12 = A10*B02 + A11*B12;

  R20 = A20*B00 + A21*B10;
  R21 = A20*B01 + A21*B11;
  R22 = A20*B02 + A21*B12;
}

template <typename T>
__device__ __forceinline__
void Matrix_Product_32_22(T &A00,T &A01,
                         T &A10,T &A11,
                         T &A20,T &A21,
                          //
                         T &B00,T &B01,
                         T &B10,T &B11,
                          //
                         T &R00,T &R01,
                         T &R10,T &R11,
                         T &R20,T &R21)
{
  R00 = A00*B00 + A01*B10; R01 = A00*B01 + A01*B11;

  R10 = A10*B00 + A11*B10; R11 = A10*B01 + A11*B11;

  R20 = A20*B00 + A21*B10; R21 = A20*B01 + A21*B11;
}
///////////////////////////////////////////////////////////////////////////////
template <typename T>
__device__ __forceinline__
void Matrix_Transpose_22(T &A00,T &A01,
                        T &A10,T &A11,
                         //
                        T &R00,T &R01,
                        T &R10,T &R11)
{
  R00 = A00;
  R01 = A10;

  R10 = A01;
  R11 = A11;
}

template <typename T>
__device__ __forceinline__
void Matrix_Transpose_33(T &A00,T &A01,T &A02,
                        T &A10,T &A11,T &A12,
                        T &A20,T &A21,T &A22,
                         //
                        T &R00,T &R01,T &R02,
                        T &R10,T &R11,T &R12,
                        T &R20,T &R21,T &R22)
{
  R00 = A00;
  R01 = A10;
  R02 = A20;

  R10 = A01;
  R11 = A11;
  R12 = A21;

  R20 = A02;
  R21 = A12;
  R22 = A22;  
}

template <typename T>
__device__ __forceinline__
void Matrix_Transpose_32_23(T &A00,T &A01,
                           T &A10,T &A11,
                           T &A20,T &A21,
                            //
                           T &R00,T &R01,T &R02,
                           T &R10,T &R11,T &R12)
{
  R00 = A00;
  R01 = A10;
  R02 = A20;

  R10 = A01;
  R11 = A11;
  R12 = A21; 
}
//////////////////////////////////////////////////////////////////////////////

template <typename T>
__device__ __forceinline__
void svd_joel22(const T A00, const T A01,
				const T A10, const T A11,
	//
				T &U00, T &U01,
				T &U10, T &U11,
	//
				T &SIG00, T &SIG11,
	//
				T &V00, T &V01,
				T &V10, T &V11) {

	cuda_scalar E = (A00 + A11)/2;
	cuda_scalar F = (A00 - A11)/2;
	cuda_scalar G = (A10 + A01)/2;
	cuda_scalar H = (A10 - A01)/2;

	cuda_scalar Q = sqrt(E*E + H*H);
	cuda_scalar R = sqrt(F*F + G*G);
	SIG00 = Q + R;
	SIG11 = Q - R;
	cuda_scalar a1 = atan2(G, F);
	cuda_scalar a2 = atan2(H, E);

	cuda_scalar theta = (a2 - a1)/2;
	cuda_scalar beta  = (a2 + a1)/2;
	
	U00 =cos(beta); U01 = -sin(beta);
	U10 = -U01; U11 = U00;

	V00 = cos(theta); V01 = sin(theta);
	V10 = -V01; V11 = V00;
}



///////////////////////////////////////////////////////////////////////////////
// A = U * SIG * VT
template <typename T>
__device__ __forceinline__
void SVD_2x2(const T A00, const T A01,
			 const T A10, const T A11,
             //
			 T &U00, T &U01,
			 T &U10, T &U11,
             //
             T &SIG00,T &SIG11,
             //
             T &V00,T &V01,
             T &V10,T &V11)
{
 T Su00, Su01,
   Su10, Su11;

  Matrix_Product_22_22(A00, A01,
                       A10, A11,
                       //
					   A00, A10,
					   A01, A11,
                       //
                       Su00, Su01,
                       Su10, Su11);

 T phi = 0.5 * atan2(Su01 + Su10, Su00 - Su11);
 T Cphi = cos(phi);
 T Sphi = sin(phi);

  U00 = Cphi; U01 = -Sphi;
  U10 = Sphi; U11 = Cphi;

 T Sw00, Sw01,
   Sw10, Sw11;

  Matrix_Product_22_22(A00, A10,
                       A01, A11,
                       //
                       A00, A01,
                       A10, A11,
                       //
                       Sw00, Sw01,
                       Sw10, Sw11);
  
 T theta = 0.5 * atan2(Sw01 + Sw10, Sw00 - Sw11);
 T Ctheta = cos(theta);
 T Stheta = sin(theta);

 T W00, W01,
   W10, W11;
  W00 = Ctheta; W01 = -Stheta;
  W10 = Stheta; W11 = Ctheta;

 T SUsum  = Su00 + Su11;
 T SUdiff = sqrt((Su00 - Su11)*(Su00 - Su11) + 4 * Su01 * Su10);
  
 T svals0, svals1;
  svals0 = sqrt((SUsum + SUdiff) / 2);
  svals1 = sqrt((SUsum - SUdiff) / 2);

  SIG00 = svals0; 
  SIG11 = svals1;
 
 T tmp00, tmp01,
   tmp10, tmp11;
  Matrix_Product_22_22(U00, U10,
                       U01, U11,
                       //
                       A00, A01,
                       A10, A11,
                       //
                       tmp00, tmp01,
                       tmp10, tmp11);
  
 T S00, S01,
   S10, S11;
  Matrix_Product_22_22(tmp00, tmp01,
                       tmp10, tmp11,
                       //
                       W00, W01,
                       W10, W11,
                       //
                       S00, S01,
                       S10, S11);

  tmp00 = SIGN(S00); tmp01 = 0.;
  tmp10 = 0.       ; tmp11 = SIGN(S11);

  Matrix_Product_22_22(W00, W01,
                       W10, W11,
                       //
                       tmp00, tmp01,
                       tmp10, tmp11,
                       //
                       V00, V01,
                       V10, V11);
}
 ///////////////////////////////////////////////////////////////////////////////
__device__ cuda_scalar f_tr(cuda_scalar x, cuda_scalar miu, cuda_scalar lambda, cuda_scalar kappa) {

	cuda_scalar x_un_p4 = x - 1;
	x_un_p4*= x_un_p4;
	x_un_p4*= x_un_p4;

	cuda_scalar x_p2 = x*x;
	cuda_scalar x_p4 = x_p2*x_p2;
	//TODO ask christos why energie can be negative
	return miu/4* x_un_p4 + lambda/8*x_p4 - lambda/4 * x_p2;
}
///////////////////////////////////////////////////////////////////////////////
__device__ cuda_scalar f_prime_tr(cuda_scalar x, cuda_scalar miu, cuda_scalar lambda, cuda_scalar kappa) {
	
	cuda_scalar x_un = x - 1;

	return miu*x_un*x_un*x_un + lambda/2*x*x*x - lambda/2*x;
}
//////////////////////////////////////////////////////////////////////////////
__device__ cuda_scalar g_tr(cuda_scalar x, cuda_scalar miu, cuda_scalar lambda, cuda_scalar kappa) {

	return 0;
}
//////////////////////////////////////////////////////////////////////////////
__device__ cuda_scalar g_prime_tr(cuda_scalar x, cuda_scalar miu, cuda_scalar lambda, cuda_scalar kappa) {

	return 0;
}

}
}
//////////////////////////////////////////////////////////////////////////////
#endif // PROJECTIONS_GPU_MATH_H
///////////////////////////////////////////////////////////////////////////////