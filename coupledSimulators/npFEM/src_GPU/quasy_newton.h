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
#ifndef QUASY
#define QUASY

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "common.h"
#include "sum_cuda.h"
#include "sparse_matrix.h"


#define CIRCULAR_ID(id, n) (id)%(n)
#define ID_COL(i,j,N) ((N*(j)) + (i))
#define ID_ROW(i,j,N) ((N*(i)) + (j))

namespace plb {
namespace npfem {

#define OTHER_ERROR if (cudaPeekAtLastError())printf(" error %s \n", cudaGetErrorString(cudaPeekAtLastError()))
#define HANDLE_KERNEL_ERROR(...) \
do {                                         \
    __VA_ARGS__;                             \
/*    HANDLE_ERROR( cudaPeekAtLastError() );  */ \
    HANDLE_ERROR( cudaDeviceSynchronize() ); \
} while(0)

#define HANDLE_ERROR(ans) (handleError((ans), __FILE__, __LINE__))
inline void handleError(cudaError_t code, const char *file, int line)
{
	if (code != cudaSuccess) {
		fprintf(stderr, "CUDA assert: %s %s %d\n", cudaGetErrorString(code), file, line);
		exit(EXIT_FAILURE);
	}
}

//trace((x-y)'*M*(x-y)) + trace(x' L x) - trace (x' J p )
__global__ void  eval_objectif(sparse_matrix_cuda L, sparse_matrix_cuda J, double *m, double *x, double *y, double *p, double *energies, double h2, int x_n, int p_n, int ncell, int nb_sum_thread){

	int id = threadIdx.x;
	int vector_id = blockIdx.x;
	int x_adress_shift = vector_id*x_n;
	int p_adress_shift = vector_id*p_n;

	extern __shared__ double buffer[];

	buffer[id] = 0;

	//(x - y)'*(x - y)*M/(2*h^2)
	double local_x = x[x_adress_shift + id];
	double x_minus_y = local_x - y[x_adress_shift + id];

	buffer[id] += m[id]/h2 * x_minus_y*x_minus_y;

	int j = 0;
	//x'*L*x/2
	double sum_line = 0;
	for (j = 0; j< L.degree; j++) {
		sum_line += x[x_adress_shift + L.index[id + x_n*j]] * L.value[id + x_n*j];
	}
	buffer[id] += sum_line*local_x / 2;

	//-x'*J*p
	sum_line = 0;
	for (j = 0; j< J.degree; j++) {
		sum_line -= p[p_adress_shift + J.index[id + x_n*j]] * J.value[id + x_n*j];
	}

	buffer[id] += sum_line*local_x;
    //||p||^2/2
	sum_line = 0;
	for (int i = id; i < p_n; i += x_n) {
		sum_line += p[p_adress_shift + i] * p[p_adress_shift + i];
	}
	buffer[id] += sum_line*0.5;

	__syncthreads();
	sum(buffer, nb_sum_thread, x_n, id);
	__syncthreads();

	if (id == 0) {
		energies[ncell*(vector_id % 3) + vector_id / 3] = buffer[0];
		//if(blockIdx.x == 0)
		//printf("E = %f \n", buffer[0]);
	}
}

__global__ void  center(double *points, int n_points) {

	extern __shared__ double buffer[];
	
	buffer[threadIdx.x] = points[threadIdx.x + blockIdx.x*n_points];

	__syncthreads();
		sum(buffer, 256, n_points, threadIdx.x);
	__syncthreads();
	
	points[threadIdx.x + blockIdx.x*n_points] -= buffer[0]/ n_points;

}

__global__ void  test_convergence(double *energies, double *energies_prev, int *has_converged, double *descent_direction_dot_gradient, double scalar, int ncell) {

	//energies est energies_prev n'ont pas la meme structure, energies ne contient que l'energie par composant, cell-ci doit etre somme.
	double E = energies[threadIdx.x] + energies[threadIdx.x + ncell] + energies[threadIdx.x + 2 * ncell];
	#if DEBUG
		if (threadIdx.x == 0)printf("%f \n", E);
	#endif
	has_converged[blockIdx.x*blockDim.x + threadIdx.x] = E <= energies_prev[blockIdx.x*blockDim.x + threadIdx.x] + scalar*descent_direction_dot_gradient[blockIdx.x*blockDim.x + threadIdx.x];
	energies_prev[blockIdx.x*blockDim.x + threadIdx.x] = E;
}

//eval gradient 1/h^2M(x - y) + Lx -Jp 
__global__ void eval_gradient(sparse_matrix_cuda L, sparse_matrix_cuda J, double *x, double *y, double *p, double*m, cuda_scalar *force, double *result, double h2, int x_n, int p_n) {

	int id = threadIdx.x;
	int x_adress_shift = blockIdx.x*x_n;
	int p_adress_shift = blockIdx.x*p_n;

	double sum_line = 0;

	//M(x - y)/h2
	sum_line += m[id] * (x[x_adress_shift + id] - y[x_adress_shift + id]);
	sum_line /= h2;

	int j = 0;
	//L*x 
	for (j = 0; j < L.degree; j++) {
		sum_line += L.value[id + x_n*j] * x[x_adress_shift + L.index[id + x_n*j]];
	}

	//-J*p 
	for (j = 0; j < J.degree; j++) {
		sum_line -= J.value[id + x_n*j] * p[p_adress_shift + J.index[id + x_n*j]];
	}
	result[x_adress_shift + id] = sum_line - force[x_adress_shift + id];
	force[x_adress_shift + id] = 0;
}


__global__ void eval_gradient_and_copy(sparse_matrix_cuda L, sparse_matrix_cuda J, double *x, double *y, double *p, double*m, cuda_scalar *force, double *gradient, double *gradient_prev, double *descent_direction, double h2, int x_n, int p_n) {

	int id = threadIdx.x;
	int vector_id = blockIdx.x;
	int x_adress_shift = vector_id*x_n;
	int p_adress_shift = vector_id*p_n;

	double sum_line = 0;

	//M(x - y)/h2
	sum_line += m[id] * (x[x_adress_shift + id] - y[x_adress_shift + id]);
	sum_line /= h2;
	
	int j = 0;
	//L*x 
	for (j = 0; j < L.degree; j++) {
		//printf("l %f \n", L.value[id + x_n*j]);
		sum_line += L.value[id + x_n*j]* x[x_adress_shift + L.index[id + x_n*j]];
	}

	//-J*p 
	for (j = 0; j < J.degree; j++) {
		sum_line -= J.value[id + x_n*j]*p[p_adress_shift + J.index[id + x_n*j]];	
	}
	
	//printf("force in gradient %f %d\n", force[x_adress_shift + id], threadIdx.x);

	gradient_prev[x_adress_shift + id] = gradient[x_adress_shift + id];
	gradient[x_adress_shift + id] = sum_line - force[x_adress_shift + id];
	descent_direction[x_adress_shift + id] = - gradient[x_adress_shift + id];
	force[x_adress_shift + id] = 0;

	//if (threadIdx.x == 0)printf("grad %f \n", gradient[x_adress_shift + id]);
}


//L_BFG
__global__ void comput_s_t_rho(double *point, double *prev_point, double *gradient, double *prev_gradient, double *s, double *t, double *rho, int head, int n, int nb_cell,  int nb_sum_thread) {

	int id = threadIdx.x;
	//mouai a factorise

	int point_id = blockIdx.x * 3 * n + id;
	int s_id = CIRCULAR_ID(head, MEM_SIZE) * 3 * n*nb_cell + point_id;

	extern __shared__ double buffer[];
	double sx = point[point_id        ] - prev_point[point_id        ];
	double sy = point[point_id +     n] - prev_point[point_id +     n];
	double sz = point[point_id + 2 * n] - prev_point[point_id + 2 * n];

	double tx = gradient[point_id        ] - prev_gradient[point_id        ];
	double ty = gradient[point_id +     n] - prev_gradient[point_id +     n];
	double tz = gradient[point_id + 2 * n] - prev_gradient[point_id + 2 * n];
    /*
	if (threadIdx.x == 0) {
		printf("%f %f %f \n", gradient[point_id], gradient[point_id + n], gradient[point_id + 2 * n]);
		printf("%T %f %f \n",tx, ty, tz);
	}
	*/
	//est ce que c'est dans la cache
	buffer[id] =  sx*tx;
	buffer[id] += sy*ty;
	buffer[id] += sz*tz;

	s[s_id        ] = sx;//point[point_id      ] - prev_point[point_id      ];
	s[s_id +     n] = sy;//point[point_id +   n] - prev_point[point_id +   n];
	s[s_id + 2 * n] = sz;//point[point_id + 2*n] - prev_point[point_id + 2*n];

	t[s_id        ] = tx;//gradient[point_id      ] - prev_gradient[point_id      ];
	t[s_id +     n] = ty;//gradient[point_id +   n] - prev_gradient[point_id +   n];
	t[s_id + 2 * n] = tz;//gradient[point_id + 2*n] - prev_gradient[point_id + 2*n];

	//sommation
	__syncthreads();

	sum(buffer, nb_sum_thread, n, id);

	__syncthreads();


	if (id == 0) {
		rho[CIRCULAR_ID(head, MEM_SIZE)*nb_cell + blockIdx.x] = buffer[0];
		//printf("rho %f \n", rho[CIRCULAR_ID(head, MEM_SIZE)*nb_cell + blockIdx.x]);
	}

}

__global__ void x_eq_y(double *y, double *x, int  n) {

	int id = blockIdx.x*n + threadIdx.x;
	x[id]  = y[id];
	//if (threadIdx.x == 0)printf("descente %f \n", x[id]);
}


__global__ void x__plus_equal_alpha_time_y(double *y, double *x, double alpha, int  n) {

	int id = blockIdx.x*n + threadIdx.x;
	x[id] += alpha*y[id];
	//if (threadIdx.x == 0)printf("descente %f point %f\n", alpha*y[id], x[id]);
}


__global__ void l_bfg(double *s, double *t, double *rho, double *alpha_global, double *q, int tail, int head, int n, int nb_cell, int nb_sum_thread) {

	int id = threadIdx.x;

	int point_id = (blockIdx.x) * 3 * n + id;

	extern __shared__ double alpha[];

	double qx = q[point_id        ];
	double qy = q[point_id +     n];
	double qz = q[point_id + 2 * n];

	//if(threadIdx.x == 0)printf("q %f %f %f \n", qx, qy, qz);

	for (int i = head - 1; i >= tail; i--) {
		int s_id = CIRCULAR_ID(i, MEM_SIZE) * 3 * n*nb_cell + point_id;

		//double local_rho = rho[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x];

		alpha[id] =  (s[s_id        ])*qx;
		alpha[id] += (s[s_id +     n])*qy;
		alpha[id] += (s[s_id + 2 * n])*qz;

		__syncthreads();

		sum(alpha, nb_sum_thread, n, id);

		__syncthreads();

		double local_alpha = alpha[0] / rho[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x];
		alpha_global[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x] = local_alpha;
		//if(id == 0)printf("alpha %.16g %.16g %.16g\n", alpha[0]/rho[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x], alpha[0], rho[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x]);
		//if(id == 0)printf("alpha %.16g  rho %f\n", local_alpha, rho[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x] );

		qx -= local_alpha*t[s_id        ];
		qy -= local_alpha*t[s_id +     n];
		qz -= local_alpha*t[s_id + 2 * n];

		__syncthreads();

	}

	q[point_id        ] = qx;
	q[point_id +     n] = qy;
	q[point_id + 2 * n] = qz;
}

__global__ void l_bfg2(double *s, double *t, double *rho, double *alpha, double *r, int tail, int head, int n, int nb_cell,  int nb_sum_thread) {

	int id = threadIdx.x;

	int point_id = (blockIdx.x) * 3 * n + id;

	extern __shared__ double nu[];

	double rx = r[point_id        ];
	double ry = r[point_id +     n];
	double rz = r[point_id + 2 * n];

	for (int i = tail; i < head; i++) {
		int s_id = CIRCULAR_ID(i, MEM_SIZE) * 3 * n*nb_cell + point_id;

		//double local_rho = rho[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x];

		nu[id] =  (t[s_id        ])*rx;
		nu[id] += (t[s_id +     n])*ry;
		nu[id] += (t[s_id + 2 * n])*rz;

		__syncthreads();

		sum(nu, nb_sum_thread, n, id);

		__syncthreads();

		double local_nu = alpha[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x] - nu[0] / rho[CIRCULAR_ID(i, MEM_SIZE)*nb_cell + blockIdx.x];

		rx += local_nu*s[s_id        ];
		ry += local_nu*s[s_id +     n];
		rz += local_nu*s[s_id + 2 * n];

		__syncthreads();

	}

	r[point_id        ] = rx;
	r[point_id +     n] = ry;
	r[point_id + 2 * n] = rz;
}

//leger overhead pour 10 matrix en parallel
__global__ void  mat_vect_mult_(double *A, double *x, double *result, int n) {
	int id = threadIdx.x;
	//int matid = blockIdx.x;

	double sum = 0;
	for (int j = 0; j<n; j++) {
        //if (j == 0 && blockIdx.x == 0)printf("%f |%d\n", A[ID_COL(id, j, n)], ID_COL(id, j, n));
		sum += x[blockIdx.x*n + j]*A[ID_COL(id, j, n)];
	}
	
	//if(id == 1 && blockIdx.x == 0)printf("\n id %d sum %f\n", id, sum);
	result[blockIdx.x*n + id] = sum;
}

//compute sum(x*y) = trace(x'*y) = trace(y'*x)
__device__ void trace_x_time_y(double *x, double *y, double *result, int n, int id, int  nb_sum_thread) {

	extern __shared__ double buffer[];

	if (threadIdx.x < n) {
		buffer[threadIdx.x] = x[id] * y[id] + x[id + n] * y[id + n] + x[id + 2 * n] * y[id + 2 * n];
	}
	__syncthreads();

	sum(buffer, nb_sum_thread, n, threadIdx.x);

	__syncthreads();

	if (threadIdx.x == 0) {
		//printf("blockIdx.x %d sum %f \n", blockIdx.x, buffer[0]);
		result[blockIdx.x] = buffer[0];
	}
}

__global__ void reduce(double *x, double *y, double *result, int n, int  nb_sum_thread) {

	trace_x_time_y(x, y, result, n, blockIdx.x * 3 * n + threadIdx.x, nb_sum_thread);
}

__global__ void curvature_test(double *gradient, double *descent_direction, double *descent_direction_dot_gradient, int *has_converged, double gamma2, int n, int nb_sum_thread) {

	extern __shared__ double descent_direction_dot_new_gradient[];
	int id = blockIdx.x * 3 * n + threadIdx.x;

	if (threadIdx.x < n) {
		descent_direction_dot_new_gradient[threadIdx.x] = gradient[id] * descent_direction[id] +
			gradient[id + n] * descent_direction[id + n] +
			gradient[id + 2 * n] * descent_direction[id + 2 * n];
	}
	__syncthreads();

	sum(descent_direction_dot_new_gradient, nb_sum_thread, n, threadIdx.x);

	__syncthreads();

	if (threadIdx.x == 0) {
		has_converged[blockIdx.x] = (descent_direction_dot_new_gradient[0] >= gamma2*descent_direction_dot_gradient[blockIdx.x]);
	}

}

__global__ void momentum_points_first_guess(double *points, double *points_prev, double *velocities, double *force_extern,  double *momentum, double *M, double h, int n) {

	int id = threadIdx.x;
	int vector_id = blockIdx.x*n;
	
	momentum[vector_id + id] = points[vector_id + id] + h*(velocities[vector_id + id] + h*force_extern[vector_id + id] / M[id]);
	points_prev[vector_id + id] = points[vector_id + id];
	points[vector_id + id] = momentum[vector_id + id];
	//if (threadIdx.x == 0)printf("velocities %f force %f momentum %f\n", velocities[vector_id + id], force_extern[vector_id + id], momentum[vector_id + id]);
}

__global__ void copy_gradient_and_points(double *points, double *points_prev, double *gradient, double *gradient_prev, double *q, int n) {
	int id = blockIdx.x*n + threadIdx.x;

	points_prev[id] = points[id];
	gradient_prev[id] = gradient[id];
	q[id] = -gradient[id];
}

__global__ void velocities_equals_points_minus_points_prev(double *points, double *points_prev, double *velocities, double one_over_h, int n) {
	int id = blockIdx.x*n + threadIdx.x;
	velocities[id] = (points[id] - points_prev[id])*one_over_h;
}

}
}

#endif //QUASY
