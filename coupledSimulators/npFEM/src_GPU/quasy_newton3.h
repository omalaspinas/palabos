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
#ifndef QUASY3
#define QUASY3

#include "quasy_newton.h"

namespace plb {
namespace npfem {

__global__ void  test_convergence_d(double *energies, double *energies_prev, int *has_converged, double *descent_direction_dot_gradient, double scalar, int ncell) {

	//energies est energies_prev n'ont pas la meme structure, energies ne contient que l'energie par composant, cell-ci doit etre somme.
	double E = energies[threadIdx.x] + energies[threadIdx.x + ncell] + energies[threadIdx.x + 2 * ncell];
	#if DEBUG
		if (threadIdx.x == 0)printf("%f \n", E);
	#endif
	has_converged[blockIdx.x*blockDim.x + threadIdx.x] = E <= energies_prev[blockIdx.x*blockDim.x + threadIdx.x] + scalar*descent_direction_dot_gradient[blockIdx.x*blockDim.x + threadIdx.x];
	energies_prev[blockIdx.x*blockDim.x + threadIdx.x] = E;
}

__device__ void  center_d(double *points, int n_points, cuda_scalar *center, double *center_output, int point_id) {

	center[threadIdx.x             ] = points[point_id             ];
	center[threadIdx.x +   n_points] = points[point_id +   n_points];
	center[threadIdx.x + 2*n_points] = points[point_id + 2*n_points];
	//TODO aaargg big wast of thread
	//Now I use lest thread, is it faster?
	__syncthreads();
	triple_sum(center, center + n_points, center + 2*n_points, 256, n_points, point_id);
	__syncthreads();

	center_output[0] = center[0         ]/n_points;
	center_output[1] = center[  n_points]/n_points;
	center_output[2] = center[2*n_points]/n_points;

	points[point_id             ] -= center[0         ]/n_points;
	points[point_id +   n_points] -= center[n_points  ]/n_points;
	points[point_id + 2*n_points] -= center[2*n_points]/n_points;

}

__device__ void printL(sparse_matrix_cuda L) {
	printf("L %d \n", L.degree);
}

//trace((x-y)'*M*(x-y)) + trace(x' L x) - trace (x' J p ) + p*p/2
__device__ double eval_objectif_and_grandient3(const sparse_matrix_cuda L, const sparse_matrix_cuda J, const sparse_matrix_cuda M_star,
								 const double *m, const double *x, const double *y, const double *p, 
								 const double *E_nonePD, cuda_scalar *f_nonePD, const double h2, const int x_n, const int p_n, const int n_const, const int nb_sum_thread, double *gradient, double *buffer, cuda_scalar *buffer2, int id) {

	const int x_adress_shift = blockIdx.x*x_n*3;
	const int p_adress_shift = blockIdx.x*p_n*3;

	double gradientx = 0;
	double gradienty = 0;
	double gradientz = 0;

	double tpx = 0;
	double tpy = 0;
	double tpz = 0;

	double matrix_value = 0;
	int j = 0, column_id;

	double xx = x[id        ];
	double xy = x[id +   x_n];
	double xz = x[id + 2*x_n];

	double x_yx = xx - y[id        ];
	double x_yy = xy - y[id +   x_n];
	double x_yz = xz - y[id + 2*x_n];

	//G = M_star(x-y)
	for (j = 0; j< M_star.degree; j++) {
		column_id    = M_star.index[threadIdx.x + x_n*j] + x_adress_shift;
		matrix_value = M_star.value[threadIdx.x + x_n*j];
		tpx += (x[column_id        ] - y[column_id        ])*matrix_value;
		tpy += (x[column_id +   x_n] - y[column_id +   x_n])*matrix_value;
		tpz += (x[column_id + 2*x_n] - y[column_id + 2*x_n])*matrix_value; // /h2/2;
	}	
	//E = (x-y)M_star(x-y)/2
	buffer[threadIdx.x] = (tpx*x_yx + tpy*x_yy + tpz*x_yz)/ 2;

	gradientx = tpx;
	gradienty = tpy;
	gradientz = tpz;
	
	tpx = 0;
	tpy = 0;
	tpz = 0;

	//G += L*x
	for (j = 0; j< L.degree; j++) {
		column_id    = L.index[threadIdx.x + x_n*j] + x_adress_shift;
		matrix_value = L.value[threadIdx.x + x_n*j];
		tpx += x[column_id          ]*matrix_value;
		tpy += x[column_id +     x_n]*matrix_value;
		tpz += x[column_id + 2 * x_n]*matrix_value;	
	}

	//E += x'*L*x/2
	buffer[threadIdx.x] += (tpx*xx + tpy*xy + tpz*xz)/2;

	gradientx += tpx;
	gradienty += tpy;
	gradientz += tpz;
	tpx = 0;
	tpy = 0;
	tpz = 0;

    //G -= J*p
	for (j = 0; j< J.degree; j++) {
		column_id    = J.index[threadIdx.x + x_n*j] + p_adress_shift;
		matrix_value = J.value[threadIdx.x + x_n*j];
		tpx += p[column_id          ]*matrix_value;
		tpy += p[column_id +     p_n]*matrix_value;
		tpz += p[column_id + 2 * p_n]*matrix_value;
	}

	__syncthreads();
	//E-= x'*J*p
	buffer[threadIdx.x] -= (tpx*xx + tpy*xy + tpz*xz);
	
	gradient[id        ] = gradientx - tpx - f_nonePD[id        ];
	gradient[id +   x_n] = gradienty - tpy - f_nonePD[id +   x_n];
	gradient[id + 2*x_n] = gradientz - tpz - f_nonePD[id + 2*x_n];
	
	//||p||^2/2 mmm pas sure a voir
	tpx = 0;
	for (int i = threadIdx.x; i < p_n; i += blockDim.x) {
		column_id = p_adress_shift + i;
		tpx += p[column_id        ] * p[column_id        ];
		tpx += p[column_id +   p_n] * p[column_id +   p_n];
		tpx += p[column_id + 2*p_n] * p[column_id + 2*p_n];
	}
	
	tpx /= 2;
	for (int i = blockIdx.x*n_const + threadIdx.x; i < (blockIdx.x + 1)*n_const; i += x_n) {
		tpx += E_nonePD[i];
	}
	buffer[threadIdx.x] += tpx;
	buffer2[threadIdx.x] = gradient[id]*gradient[id] + gradient[id + x_n]*gradient[id + x_n] + gradient[id + 2*x_n]*gradient[id + 2*x_n];

	__syncthreads();
	//TODO BUG buffer and buffer
	double_sum(buffer, buffer2, nb_sum_thread, x_n, threadIdx.x);
	//sum(buffer, nb_sum_thread, x_n, threadIdx.x);
	__syncthreads();
	//if(threadIdx.x == 0)printf("BUFFER 0 %f \n", buffer[0]);
	return buffer[0];
}

//DEAD CODE, I dont use it anymore but for some reason I don't feel like deleting it
__device__ void eval_gradient3(const sparse_matrix_cuda L, const sparse_matrix_cuda J, const sparse_matrix_cuda M_star, const double *x, const double *y, const double *p, const double*m,
                               const cuda_scalar *force, double *result, const double h2, const int x_n, const int p_n) {

	int x_adress_shift = blockIdx.x*x_n * 3;
	int p_adress_shift = blockIdx.x*p_n * 3;

	double gx, gy, gz, matrix_value;

	//M(x - y)/h2
	gx = 0;// m[threadIdx.x] * (x[x_adress_shift + threadIdx.x          ] - y[x_adress_shift + threadIdx.x          ])/h2;
	gy = 0;// m[threadIdx.x] * (x[x_adress_shift + threadIdx.x +     x_n] - y[x_adress_shift + threadIdx.x +     x_n])/h2;
	gz = 0;// m[threadIdx.x] * (x[x_adress_shift + threadIdx.x + 2 * x_n] - y[x_adress_shift + threadIdx.x + 2 * x_n])/h2;

	int j = 0, column_id;

	//M_start(x-y)
	
	for (j = 0; j < M_star.degree; j++) {
		column_id    = M_star.index[threadIdx.x + x_n*j];
		matrix_value = M_star.value[threadIdx.x + x_n*j];
		//printf("%d \n", column_id);
		gx += matrix_value * (x[x_adress_shift + column_id          ] - y[x_adress_shift + column_id          ]);
		gy += matrix_value * (x[x_adress_shift + column_id +     x_n] - y[x_adress_shift + column_id +     x_n]);
		gz += matrix_value * (x[x_adress_shift + column_id + 2 * x_n] - y[x_adress_shift + column_id + 2 * x_n]);
		
	}
	
	//L*x 
	for (j = 0; j < L.degree; j++) {
		column_id    = L.index[threadIdx.x + x_n*j];
		matrix_value = L.value[threadIdx.x + x_n*j];
		gx += matrix_value * x[x_adress_shift + column_id          ];
		gy += matrix_value * x[x_adress_shift + column_id +     x_n];
		gz += matrix_value * x[x_adress_shift + column_id + 2 * x_n];
	}
	
	//-J*p 
	for (j = 0; j < J.degree; j++) {
		column_id    = J.index[threadIdx.x + x_n*j];
		matrix_value = J.value[threadIdx.x + x_n*j];
		gx -= matrix_value * p[p_adress_shift + column_id          ];
		gy -= matrix_value * p[p_adress_shift + column_id +     p_n];
		gz -= matrix_value * p[p_adress_shift + column_id + 2 * p_n];

		//if (threadIdx.x == 109)printf("proj %f | %f %f %f \n", matrix_value, p[p_adress_shift + column_id], p[p_adress_shift + column_id + p_n], p[p_adress_shift + column_id + 2 * p_n]);
	}
	
	result[x_adress_shift + threadIdx.x        ] = gx - force[x_adress_shift + threadIdx.x        ];
	result[x_adress_shift + threadIdx.x +   x_n] = gy - force[x_adress_shift + threadIdx.x +   x_n];
	result[x_adress_shift + threadIdx.x + 2*x_n] = gz - force[x_adress_shift + threadIdx.x + 2*x_n];
	
	//if (threadIdx.x == 0)printf("x %f %f %f grad %f %f %f force %f %f %f\n", x[x_adress_shift + threadIdx.x], x[x_adress_shift + threadIdx.x + x_n], x[x_adress_shift + threadIdx.x + 2*x_n],
	//	gx, gy, gz, force[x_adress_shift + threadIdx.x], force[x_adress_shift + threadIdx.x + x_n], force[x_adress_shift + threadIdx.x + 2 * x_n]);
}

//L_BFG
__device__ void comput_s_t_rho_d(const double *point, const double *prev_point, const double *gradient, const double *prev_gradient, 
								 double *s, double *t, double *rho, const int head, const int n, const int nb_cell, const int nb_sum_thread, double *buffer, int id) {


	int s_id = CIRCULAR_ID(head, MEM_SIZE) * 3 * n*nb_cell + id;

	double sx = point[id        ] - prev_point[id        ];
	double sy = point[id +     n] - prev_point[id +     n];
	double sz = point[id + 2 * n] - prev_point[id + 2 * n];

	double tx = gradient[id        ] - prev_gradient[id        ];
	double ty = gradient[id +     n] - prev_gradient[id +     n];
	double tz = gradient[id + 2 * n] - prev_gradient[id + 2 * n];
  
	//est ce que c'est dans la cache
	buffer[threadIdx.x] =  sx*tx;
	buffer[threadIdx.x] += sy*ty;
	buffer[threadIdx.x] += sz*tz;

	s[s_id        ] = sx;//point[point_id      ] - prev_point[point_id      ];
	s[s_id +     n] = sy;//point[point_id +   n] - prev_point[point_id +   n];
	s[s_id + 2 * n] = sz;//point[point_id + 2*n] - prev_point[point_id + 2*n];

	t[s_id        ] = tx;//gradient[point_id      ] - prev_gradient[point_id      ];
	t[s_id +     n] = ty;//gradient[point_id +   n] - prev_gradient[point_id +   n];
	t[s_id + 2 * n] = tz;//gradient[point_id + 2*n] - prev_gradient[point_id + 2*n];

	//sommation
	__syncthreads();

	sum(buffer, nb_sum_thread, n, threadIdx.x);

	__syncthreads();

	if (threadIdx.x == 0) {
		rho[CIRCULAR_ID(head, MEM_SIZE)*nb_cell + blockIdx.x] = buffer[0];
		//printf("rho %f \n", rho[CIRCULAR_ID(head, MEM_SIZE)*nb_cell + blockIdx.x]);
	}
}

__device__ void x__plus_equal_alpha_time_y3(const double *y, double *x, const double alpha, const int n, const int id) {
	//if (threadIdx.x == 0)printf("x %f %f %f y %f %f %f alpha %f\n", x[id], x[id + n], x[id + 2 * n], y[id], y[id + n], y[id + 2 * n], alpha);

	x[id      ] += alpha*y[id      ];
	x[id +   n] += alpha*y[id +   n];
	x[id + 2*n] += alpha*y[id + 2*n];
}

__device__ void l_bfg_d(const double *s, const double *t, const double *rho, double *alpha_global, 
                        double *q, int tail, int head, const int n, const int nb_cell, const int nb_sum_thread, double *alpha, int id) {

	double qx = q[id        ];
	double qy = q[id +     n];
	double qz = q[id + 2 * n];

	//if(threadIdx.x == 0)printf("q %f %f %f \n", qx, qy, qz);

	for (int i = head - 1; i >= tail; i--) {
		int s_id = CIRCULAR_ID(i, MEM_SIZE)*3*n*nb_cell + id;

		//double local_rho = rho[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x];

		alpha[threadIdx.x] =  (s[s_id        ])*qx;
		alpha[threadIdx.x] += (s[s_id +     n])*qy;
		alpha[threadIdx.x] += (s[s_id + 2 * n])*qz;

		__syncthreads();

		sum(alpha, nb_sum_thread, n, threadIdx.x);

		__syncthreads();
		double local_alpha = 0;
		if( rho[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x] != 0){
			local_alpha = alpha[0] / rho[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x];
		}
		alpha_global[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x] = local_alpha;
		//if(id == 0)printf("alpha %.16g %.16g %.16g\n", alpha[0]/rho[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x], alpha[0], rho[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x]);
		//if(id == 0)printf("alpha %.16g  rho %f\n", local_alpha, rho[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x] );

		qx -= local_alpha*t[s_id        ];
		qy -= local_alpha*t[s_id +     n];
		qz -= local_alpha*t[s_id + 2 * n];

		__syncthreads();

	}
	q[id        ] = qx;
	q[id +     n] = qy;
	q[id + 2 * n] = qz;
}

__device__ void l_bfg2_d(const double *s, const double *t, const double *rho, const double *alpha, double *r, int tail, int head, const int n, const int nb_cell, const int nb_sum_thread, double *nu, int id) {

	double rx = r[id        ];
	double ry = r[id +     n];
	double rz = r[id + 2 * n];

	for (int i = tail; i < head; i++) {
		int s_id = CIRCULAR_ID(i, MEM_SIZE) * 3 * n*nb_cell + id;

		//double local_rho = rho[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x];
		nu[threadIdx.x] =  (t[s_id        ])*rx;
		nu[threadIdx.x] += (t[s_id +     n])*ry;
		nu[threadIdx.x] += (t[s_id + 2 * n])*rz;

		__syncthreads();

		sum(nu, nb_sum_thread, n, threadIdx.x);

		__syncthreads();
		double local_nu = 0;
	
		if (rho[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x] != 0) {
			local_nu = alpha[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x] - nu[0] / rho[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x];
		}
		rx += local_nu*s[s_id        ];
		ry += local_nu*s[s_id +     n];
		rz += local_nu*s[s_id + 2 * n];

		__syncthreads();
	}

	r[id        ] = rx;
	r[id +     n] = ry;
	r[id + 2 * n] = rz;
}


__device__ void  mat_vect_mult3(const double *A, const double *x, double *result, const int n, int id) {

	double Ax_x = 0, Ax_y = 0, Ax_z = 0, A_coef;
	int tp_index;
	for (int j = 0; j< n; j++) {
		A_coef = A[ID_COL(threadIdx.x, j, n)];
		tp_index = blockIdx.x*3*n + j;
		Ax_x += x[tp_index      ] * A_coef;
		Ax_y += x[tp_index +   n] * A_coef;
		Ax_z += x[tp_index + 2*n] * A_coef;
	}
	result[id      ] = Ax_x;
	result[id +   n] = Ax_y;
	result[id + 2*n] = Ax_z;
}

__device__ void curvature_test_d(const double *gradient, const double *descent_direction, const double *descent_direction_dot_gradient, 
								int *has_converged, const double gamma2, const int n, const int nb_sum_thread, const int id, double *buffer) {

	if (threadIdx.x < n) {
		buffer[threadIdx.x] = gradient[id      ] * descent_direction[id      ] +
							  gradient[id +   n] * descent_direction[id +   n] +
							  gradient[id + 2*n] * descent_direction[id + 2*n];
	}
	__syncthreads();
	//compute gradient.dot(desent direction)
	sum(buffer, nb_sum_thread, n, threadIdx.x);

	__syncthreads();

	if (threadIdx.x == 0) {
		has_converged[blockIdx.x] = (buffer[0] >= gamma2*descent_direction_dot_gradient[blockIdx.x]);
	}
}


__device__ void momentum_points_first_guess3(double *points, const double *velocities, const double *M_tild_inv, const double *force_extern, double *momentum, const  double *M, const double h, const int n, const int id) {
	
	//TODO make it multi_cell
	double Ax_x = 0, Ax_y = 0, Ax_z = 0, A_coef;

	for (int j = 0; j<n; j++) {
		A_coef = M_tild_inv[ID_COL(threadIdx.x, j, n)];
		int local_id = j + 3*n*blockIdx.x;
		//if (threadIdx.x == 0)printf("m_tild %f h %f force %f \n", A_coef, h, force_extern[local_id]);
		Ax_x += A_coef *(M[j] * h*velocities[local_id      ] + h*h*force_extern[local_id      ]);
		Ax_y += A_coef *(M[j] * h*velocities[local_id +   n] + h*h*force_extern[local_id +   n]);
		Ax_z += A_coef *(M[j] * h*velocities[local_id + 2*n] + h*h*force_extern[local_id + 2*n]);
	}
	
	momentum[id      ] = points[id      ] + Ax_x;
	momentum[id +   n] = points[id +   n] + Ax_y;
	momentum[id + 2*n] = points[id + 2*n] + Ax_z;

	points[id      ] = momentum[id      ];
	points[id +   n] = momentum[id +   n];
	points[id + 2*n] = momentum[id + 2*n];
}



/*
__device__ void momentum_points_first_guess3(double *points, 
											 const double *velocities, const double *M_tild_inv, double *force_extern,  double *momentum, const  double *M, const double h, const int n, const int id) {

	momentum[id      ] = points[id      ] + h*(velocities[id      ] + h*force_extern[id      ] / M[id]);
	momentum[id +   n] = points[id +   n] + h*(velocities[id +   n] + h*force_extern[id +   n] / M[id]);
	momentum[id + 2*n] = points[id + 2*n] + h*(velocities[id + 2*n] + h*force_extern[id + 2*n] / M[id]);

	//if(id == 109)printf("momentum %f %f %f M %f \n", momentum[id], momentum[id + n], momentum[id + 2 * n], M[id]);
 
	points[id      ] = momentum[id      ];
	points[id +   n] = momentum[id +   n];
	points[id + 2*n] = momentum[id + 2*n];
}
*/

__device__ void copy_gradient_and_points3(const double *points, double *points_prev, const double *gradient, double *gradient_prev, double *q, const int n, const int id) {
	
	points_prev[id      ] = points[id      ];
	points_prev[id +   n] = points[id +   n];
	points_prev[id + 2*n] = points[id + 2*n];

	gradient_prev[id      ] = gradient[id      ];
	gradient_prev[id   + n] = gradient[id +   n];
	gradient_prev[id + 2*n] = gradient[id + 2*n];

	q[id      ]= -gradient[id      ];
	q[id   + n]= -gradient[id   + n];
	q[id + 2*n]= -gradient[id + 2*n];
}

__device__ void velocities_equals_points_minus_points_prev3(const double *points, const double *points_prev, double *velocities, const double one_over_h, const int n, const int id) {

	velocities[id      ] = (points[id      ] - points_prev[id      ])*one_over_h;
	velocities[id   + n] = (points[id +   n] - points_prev[id +   n])*one_over_h;
	velocities[id + 2*n] = (points[id + 2*n] - points_prev[id + 2*n])*one_over_h;
	//if (threadIdx.x == 0)printf("velocitiy %f %f %f \n", velocities[id], velocities[id + n], velocities[id + 2 * n]);
}

__device__ void damping(const double *points, const double *points_prev, double *velocities, 
						const double one_over_h, const double h2, const double *m, const double mass_total, const int n, 
						const int nb_sum_thread, const int id, cuda_scalar *buffer, double center[3], const float calpha, const float vel_cap, const float vel_cap_fin, const float dx_p) {

	velocities_equals_points_minus_points_prev3(points, points_prev, velocities, one_over_h, n, id);
	
	cuda_scalar mi = m[threadIdx.x];
	
	buffer[threadIdx.x      ] = points[id      ]*mi;
	buffer[threadIdx.x +   n] = points[id +   n]*mi;
	buffer[threadIdx.x + 2*n] = points[id + 2*n]*mi;

	__syncthreads();
	triple_sum(buffer, buffer + n, buffer + 2*n, nb_sum_thread, n, threadIdx.x);
	__syncthreads();

	//if(threadIdx.x == 0)printf("mi %f mtotal %f \n", mi, mass_total);

	cuda_scalar Xcmx = buffer[  0]/mass_total;
	cuda_scalar Xcmy = buffer[  n]/mass_total;
	cuda_scalar Xcmz = buffer[2*n]/mass_total;

	if (threadIdx.x == 0) {
		center[0] = Xcmx;
		center[1] = Xcmy;
		center[2] = Xcmz;
	}
	
	__syncthreads();
	//if (threadIdx.x == 0)printf("Xcmx %f Xcmy %f Xcmz %f | xi %f %f %f\n", Xcmx, Xcmy, Xcmz, points[id], points[id+n], points[id+2*n]);

	buffer[threadIdx.x      ] = velocities[id      ]*mi;
	buffer[threadIdx.x +   n] = velocities[id +   n]*mi;
	buffer[threadIdx.x + 2*n] = velocities[id + 2*n]*mi;

	__syncthreads();
	triple_sum(buffer, buffer + n, buffer + 2*n, nb_sum_thread, n, threadIdx.x);
	__syncthreads();

	cuda_scalar Vcmx = buffer[  0]/mass_total;
	cuda_scalar Vcmy = buffer[  n]/mass_total;
	cuda_scalar Vcmz = buffer[2*n]/mass_total;

	//if (threadIdx.x == 0)printf("Vcmx %f Vcmy %f Vcmz %f | vi %f %f %f\n", Vcmx, Vcmy, Vcmz, velocities[id], velocities[id + n], velocities[id + 2 * n]);

	cuda_scalar rx = points[id      ] - Xcmx;
	cuda_scalar ry = points[id +   n] - Xcmy;
	cuda_scalar rz = points[id + 2*n] - Xcmz;
	
	//if (threadIdx.x == 0)printf("rx %f ry %f rz %f | %f %f %f\n", rx, ry, rz, velocities[id] * mi, velocities[id + n] * mi, velocities[id + 2 * n] * mi);

	cross(rx, ry, rz,
		  (cuda_scalar)velocities[id]*mi, (cuda_scalar)velocities[id + n]*mi, (cuda_scalar)velocities[id + 2*n]*mi,
		  buffer[threadIdx.x], buffer[threadIdx.x + n], buffer[threadIdx.x + 2*n]);

	//if (threadIdx.x == 0)printf("rx X Vimi %f %f %f \n", buffer[threadIdx.x], buffer[threadIdx.x + n], buffer[threadIdx.x + 2*n]);

	__syncthreads();
	triple_sum(buffer, buffer + n, buffer + 2*n, nb_sum_thread, n, threadIdx.x);
	__syncthreads();

	__shared__ cuda_scalar L[3];

	if(threadIdx.x == 0){
		L[0] = buffer[  0];
		L[1] = buffer[  n];
		L[2] = buffer[2*n];
	}
	cuda_scalar r[9];
	__shared__ cuda_scalar I[9];
	__shared__ cuda_scalar I_inv[9];
	__shared__ cuda_scalar omega[3];

	Matrix_Product_33_33((cuda_scalar)0., -rz,  ry,
						 rz, (cuda_scalar)0. , -rx,
						-ry, rx, (cuda_scalar)0.,
		//
						(cuda_scalar)0., rz, -ry,
						-rz, (cuda_scalar)0., rx,
						 ry, -rx, (cuda_scalar)0.,
		//
						r[0], r[1], r[2],
						r[3], r[4], r[5],
						r[6], r[7], r[8]);

	//if (threadIdx.x == 0)printf("\n I \n %f %f %f\n %f %f %f\n %f %f %f\n", r[0], r[1], r[2],r[3], r[4], r[5],r[6], r[7], r[8]);

	buffer[threadIdx.x      ] = r[0]*mi;
	buffer[threadIdx.x +   n] = r[1]*mi;
	buffer[threadIdx.x + 2*n] = r[2]*mi;
	__syncthreads();
	triple_sum(buffer, buffer + n, buffer + 2 * n, nb_sum_thread, n, threadIdx.x);
	__syncthreads();
	if (threadIdx.x == 0) {
		I[0] = buffer[  0];
		I[1] = buffer[  n];
		I[2] = buffer[2*n];
	}
	__syncthreads();
	buffer[threadIdx.x      ] = r[3]*mi;
	buffer[threadIdx.x +   n] = r[4]*mi;
	buffer[threadIdx.x + 2*n] = r[5]*mi;
	__syncthreads();
	triple_sum(buffer, buffer + n, buffer + 2 * n, nb_sum_thread, n, threadIdx.x);
	__syncthreads();
	if (threadIdx.x == 0) {
		I[3] = buffer[  0];
		I[4] = buffer[  n];
		I[5] = buffer[2*n];
	}
	__syncthreads();
	buffer[threadIdx.x      ] = r[6]*mi;
	buffer[threadIdx.x +   n] = r[7]*mi;
	buffer[threadIdx.x + 2*n] = r[8]*mi;
	__syncthreads();
	triple_sum(buffer, buffer + n, buffer + 2 * n, nb_sum_thread, n, threadIdx.x);
	__syncthreads();
	if (threadIdx.x == 0) {
		I[6] = buffer[  0];
		I[7] = buffer[  n];
		I[8] = buffer[2*n];
		
	}
	//if (threadIdx.x == 0)printf("\n I \n %f %f %f\n %f %f %f\n %f %f %f\n", I[0], I[1], I[2], I[3], I[4], I[5], I[6], I[7], I[8]);
	__syncthreads();
	Matrix_invers3_3X3(I, I_inv);
	__syncthreads();
	//if (threadIdx.x == 0)printf("\n L \n %f %f %f\n", L[0], L[1], L[2]);
	//if (threadIdx.x == 0)printf("\n Iinv \n %f %f %f\n %f %f %f\n %f %f %f\n", I_inv[0], I_inv[1], I_inv[2], I_inv[3], I_inv[4], I_inv[5], I_inv[6], I_inv[7], I_inv[8]);
	vect_matrix_mult_3X3_par(I_inv, L, omega);
	__syncthreads();
	//if (threadIdx.x == 0)printf("\n omeg = I_inv*L\n %f %f %f\n ", omega[0], omega[1], omega[2]);
	cuda_scalar deltaV[3];
	cross(omega[0], omega[1], omega[2],
		  rx, ry, rz,
		  deltaV[0], deltaV[1], deltaV[2]);

	deltaV[0] = Vcmx + deltaV[0] - velocities[id      ];
	deltaV[1] = Vcmy + deltaV[1] - velocities[id +   n];
	deltaV[2] = Vcmz + deltaV[2] - velocities[id + 2*n];


	velocities[id      ] += calpha*deltaV[0];
	velocities[id +   n] += calpha*deltaV[1];
	velocities[id + 2*n] += calpha*deltaV[2];

#ifndef NPFEM_SA
	//if(threadIdx.x == 0)printf("1/h %f  dx_P %f", one_over_h*h2, DX_P);
	ShapeOpScalar vel_lattx = velocities[id      ]/(dx_p*one_over_h*h2);
	ShapeOpScalar vel_latty = velocities[id +   n]/(dx_p*one_over_h*h2);
	ShapeOpScalar vel_lattz = velocities[id + 2*n]/(dx_p*one_over_h*h2);
	cuda_scalar norm_vel =  Magnitude_3(vel_lattx, vel_latty, vel_lattz);
	if (norm_vel >= vel_cap) {
		velocities[id      ] = vel_lattx/norm_vel*vel_cap_fin*(dx_p*one_over_h*h2);
		velocities[id +   n] = vel_latty/norm_vel*vel_cap_fin*(dx_p*one_over_h*h2);
		velocities[id + 2*n] = vel_lattz/norm_vel*vel_cap_fin*(dx_p*one_over_h*h2);
	}
#endif
}


__device__ void trace_x_time_y3(const double *x, const double *y, double *result, const int n, const int id, const int nb_sum_thread, double *buffer) {

	if (threadIdx.x < n) {
		buffer[threadIdx.x] = x[id] * y[id] + x[id + n] * y[id + n] + x[id + 2 * n] * y[id + 2 * n];
	}
	__syncthreads();

	sum(buffer, nb_sum_thread, n, threadIdx.x);

	__syncthreads();

	if (threadIdx.x == 0) {
		//printf("blockIdx.x %d sum %f \n", blockIdx.x, buffer[0]);
		*result = buffer[0];
	}
}

}
}

#endif //QUASY3
