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
#ifndef DEVICE_UTILITIES_CU
#define DEVICE_UTILITIES_CU
///////////////////////////////////////////////////////////////////////////////
#include <chrono>

#include "device_utilities.h"
#include "projections_GPU_soa.h"
#include "GPU_data.h"
#include "quasy_newton3.h"
///////////////////////////////////////////////////////////////////////////////
namespace plb {
namespace npfem {
///////////////////////////////////////////////////////////////////////////////
std::chrono::duration<ShapeOpScalar> time_elapsed_area;
std::chrono::duration<ShapeOpScalar> time_elapsed_volume;
std::chrono::duration<ShapeOpScalar> time_elapsed_bending;
std::chrono::duration<ShapeOpScalar> time_elapsed_triangle;
double time_elapsed_area_d = 0;
double time_elapsed_volume_d = 0;
double time_elapsed_bending_d = 0;
double time_elapsed_triangle_d = 0;
///////////////////////////////////////////////////////////////////////////////
// Handle CUDA error messages
inline void CUDA_HandleError(cudaError_t err, const char *file, int line)
{
  if (err != cudaSuccess) {
    std::cout << cudaGetErrorString(err) << " in " << file << " at line " << line << std::endl;
    exit(EXIT_FAILURE);
  }
}
///////////////////////////////////////////////////////////////////////////////
void GPU_Mem_check()
{
  size_t free_byte;
  size_t total_byte;

  cudaMemGetInfo(&free_byte, &total_byte);

  ShapeOpScalar free_db = (ShapeOpScalar)free_byte;
  ShapeOpScalar total_db = (ShapeOpScalar)total_byte;

  ShapeOpScalar used_db = total_db - free_db;

  //std::cout << "GPU memory usage: used = " << used_db / 1024.0 / 1024.0 << ", free = " << free_db / 1024.0 / 1024.0 << " MB, total = " << total_db / 1024.0 / 1024.0 << " MB" << std::endl;
  //printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n", used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);
}
///////////////////////////////////////////////////////////////////////////////
template <typename T>
__global__ void set_to_zero(T *data, int nb) {
	int id = blockIdx.x*blockDim.x + threadIdx.x;
	if (id < nb)data[id] = 0;
}
//////////////////////////////////////////////////////////////////////////////////////////

__device__ void swap(double *a, double *b){
	double tp = *a;
	*a = *b;
	*b = tp;
}

template <typename T>
__device__ T kselect(T *seq, int n, int k) {
	int i = 0;
	int start_n = n;
	//randome pivot
	//swap(seq[0], seq[(int)(((double)rand())/MAXINT*(n-1))]);

	double p = seq[0];
	int t = 0;
	int i_prev = 0;

	while (t < start_n*start_n) {
		t++;
		i_prev = i;
		//swap(seq + i, seq + i + (int)(((double)rand())/MAXINT*(n - 1 - i)));
		p = seq[i];
		//is seq[i]=p
		while (i + 1 < n && p >seq[i + 1]) {
			t++;
			swap(&seq[i], &seq[i + 1]);
			i++;
		}
		int j = i + 1;

		while (j<n) {

			if (p >= seq[j]) {
				t++;
				//printf("%d > %d\n", p, seq[j]);
				swap(&seq[i + 1], &seq[j]);
				swap(&seq[i + 1], &seq[i]);
				i += 1;
			}
			j++;
		}

		if (i == k) {
			//printf("temps select Vector %d total %d  i%d | %f\n", t, start_n, i, seq[i]);
			return seq[i];
		}
		else if (i < k) {
			i++;
		}
		else {
			n = i;
			i = i_prev;
		}
	}
	return 0;
	//printf("temps select Vector END %d \n", t);	
}
////////////////////////////////////////////////////////////////////////////////////////
__launch_bounds__(258, 1)
__global__ void median_filter_g(graph_data graph, ShapeOpScalar *force, int nb){
	
	double x[10];
	double y[10];
	double z[10];
	int k = 0;
	
	for(int i = graph.indexs[threadIdx.x]; i < graph.indexs[threadIdx.x + 1]; i++){
		int id = graph.edges[i];
		//printf("id %d %d f \n", id, i );

		x[k] = force[blockIdx.x*nb*3 + id       ];
		y[k] = force[blockIdx.x*nb*3 + id +   nb];
		z[k] = force[blockIdx.x*nb*3 + id + 2*nb];
		k++;
	}
	
	int size = graph.indexs[threadIdx.x + 1] - graph.indexs[threadIdx.x];

	if(size > 10)printf("WTF !\n");
	
	ShapeOpScalar mx=0, my=0, mz=0;
	if (size % 2) {
		mx = kselect(x, size, size/2);
		my = kselect(y, size, size/2);
		mz = kselect(z, size, size/2);
	}else{
		mx = (kselect(x, size, size/2) + kselect(x, size, size/2 - 1))/2;
		my = (kselect(y, size, size/2) + kselect(y, size, size/2 - 1))/2;
		mz = (kselect(z, size, size/2) + kselect(z, size, size/2 - 1))/2;
	}

	__syncthreads();
	force[blockIdx.x*nb*3 + threadIdx.x       ] = mx;
	force[blockIdx.x*nb*3 + threadIdx.x +   nb] = my;
	force[blockIdx.x*nb*3 + threadIdx.x + 2*nb] = mz;
	
}
///////////////////////////////////////////////////////////////////////////////////////////
void median_filter(graph_data graph, ShapeOpScalar *force, int nb_cells, int nb){
	HANDLE_KERNEL_ERROR(median_filter_g<<< nb_cells, nb >>>(graph, force,  nb));
}

//////////////////////////////////////////////////////////////////////////////////////////
void send_graph(graph_data *graph_h, graph_data *graph_d, int n)
{
	CUDA_HANDLE_ERROR(cudaMalloc(&graph_d->indexs, n * sizeof(int)));
	CUDA_HANDLE_ERROR(cudaMalloc(&graph_d->edges, graph_h->nb_edges * sizeof(int)));

	CUDA_HANDLE_ERROR(cudaMemcpy(graph_d->indexs, graph_h->indexs, n * sizeof(int), cudaMemcpyHostToDevice));
	CUDA_HANDLE_ERROR(cudaMemcpy(graph_d->edges, graph_h->edges, graph_h->nb_edges * sizeof(int), cudaMemcpyHostToDevice));
}
//////////////////////////////////////////////////////////////////////////////////////////
__global__ void clone_cells_points(Mesh_info info, Simulation_input input, Simulation_data sim) {

	double px = input.points[threadIdx.x                  ];
	double py = input.points[threadIdx.x +   info.n_points];
	double pz = input.points[threadIdx.x + 2*info.n_points];

	//if (threadIdx.x == 0) printf("points  %f %f %f | %d\n", input.points[threadIdx.x], input.points[threadIdx.x + info.n_points], input.points[threadIdx.x + 2 * info.n_points], blockIdx.x);
	int id = (blockIdx.x+1)*info.n_points*3 + threadIdx.x;

	input.points[id                  ] = px + sim.center[3*(blockIdx.x+1)    ];
	input.points[id +   info.n_points] = py + sim.center[3*(blockIdx.x+1) + 1];
	input.points[id + 2*info.n_points] = pz + sim.center[3*(blockIdx.x+1) + 2];

	#ifdef DEBUG
		if(threadIdx.x == 0) printf("center  %f %f %f | %d\n", sim.center[3 *blockIdx.x], sim.center[3 * blockIdx.x  + 1], sim.center[3*blockIdx.x  + 2], blockIdx.x);
	#endif
}
///////////////////////////////////////////////////////////////////////////////////////////////
__global__ void move_first(Mesh_info info, Simulation_input input, Simulation_data sim) {
	input.points[threadIdx.x                  ] += sim.center[0];
	input.points[threadIdx.x +   info.n_points] += sim.center[1];
	input.points[threadIdx.x + 2*info.n_points] += sim.center[2];
}
///////////////////////////////////////////////////////////////////////////////////////////////
__global__ void shift_cells_points_g(ShapeOpScalar *points, ShapeOpScalar *center, int n) {

	int id = blockIdx.x*n + threadIdx.x;
	points[id      ] += center[3*(blockIdx.x + 1)    ];
	points[id +   n] += center[3*(blockIdx.x + 1) + 1];
	points[id + 2*n] += center[3*(blockIdx.x + 1) + 2];
}
///////////////////////////////////////////////////////////////////////////////////////////
__global__ void first_centerG(Mesh_info info,  Simulation_input input, Simulation_data sim, double *center){

	__shared__ extern double buffer[];

	int n = info.n_points;

	int id = blockIdx.x*3*n + threadIdx.x;

	buffer[threadIdx.x      ] = input.points[id      ];
	buffer[threadIdx.x +   n] = input.points[id +   n];
	buffer[threadIdx.x + 2*n] = input.points[id + 2*n];

	__syncthreads();
	triple_sum(buffer, buffer + n, buffer + 2 * n, info.sum_points, n, threadIdx.x);
	__syncthreads();
	if (threadIdx.x == 0) {
		center[ blockIdx.x*3    ] = buffer[0  ]/n;
		center[ blockIdx.x*3 + 1] = buffer[  n]/n;
		center[ blockIdx.x*3 + 2] = buffer[2*n]/n;
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////
void first_center(Mesh_info info, Simulation_input input, Simulation_data sim, double *center_d) {
	HANDLE_KERNEL_ERROR(first_centerG<<<info.nb_cells, info.n_points, 3 * info.n_points * sizeof(double)>>>(info, input, sim, center_d));
}
///////////////////////////////////////////////////////////////////////////////////////////////
void set_cells_initial_position(const double *center_h, double *center_d,const int nb_cells, Mesh_info info, Mesh_data mesh, Simulation_input input, Simulation_data sim) {
	CUDA_HANDLE_ERROR(cudaMemcpy(center_d, center_h, nb_cells*3*sizeof(double), cudaMemcpyHostToDevice));
	//std::cout << "cloning point "<< info.nb_cells << " " << info.n_points  <<std::endl;
	int nb_clone = info.nb_cells - 1;
	if(nb_clone){
		HANDLE_KERNEL_ERROR(clone_cells_points <<< nb_clone, info.n_points>>>(info, input, sim));
	}
	HANDLE_KERNEL_ERROR(move_first<<<1, info.n_points >>>(info, input, sim));
	HANDLE_KERNEL_ERROR(first_centerG <<<info.nb_cells, info.n_points, 3*info.n_points*sizeof(double) >>>(info, input, sim, center_d));
}
////////////////////////////////////////////////////////////////////////////////////////////////
void reset_position_d(const double *center_h, double *center_d, ShapeOpScalar *points, Mesh_info info, Mesh_data mesh, Simulation_input input_d, Simulation_data sim) {
	std::vector <ShapeOpScalar> tp_points(3*info.n_points);

	for (int i = 0; i < 3*info.n_points; i+=3) {
		//printf( "%f %f %f %d \n", points[i] , points[i+1], points[i + 2], i);
		tp_points[i/3                  ] = points[i    ];
		tp_points[i/3 +   info.n_points] = points[i + 1];
		tp_points[i/3 + 2*info.n_points] = points[i + 2];
	}
	CUDA_HANDLE_ERROR(cudaMemcpy(input_d.points, &tp_points[0],  3*info.n_points * sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
	set_cells_initial_position(center_h, center_d, info.nb_cells, info, mesh,  input_d,  sim);
}
///////////////////////////////////////////////////////////////////////////////////////////////
void shift_cells_points_g(ShapeOpScalar *points, ShapeOpScalar *center_h, ShapeOpScalar *center_d, int n, int nb_cells) {
	CUDA_HANDLE_ERROR(cudaMemcpy(center_d, center_h, nb_cells*3*sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
	HANDLE_KERNEL_ERROR(shift_cells_points_g<<<nb_cells,  n>>>(points, center_d, n));
}
///////////////////////////////////////////////////////////////////////////////////////////////
void GPU_Init(Mesh_info        *mesh_info,
			  Mesh_data        *mesh_data_d,        Mesh_data        *mesh_data_h,
              Simulation_input *simulation_input_d, Simulation_input *simulation_input_h,
              Simulation_data  *simulation_data_d ,
              Collision_data   *collision_data_d  , Collision_data   *collision_data_h, cuda_scalar **matrices_d)
{
  GPU_Mem_check();

  int n_points = mesh_info->n_points;
  int nb_cells = mesh_info->nb_cells;
  int n_constraints = mesh_info->n_constraints;
  int n_triangles = mesh_info->n_triangles;
  int idO = mesh_info->idO;
  //not for joel check palabos force
  // ALLOCATE SPACE ON GPU
  CUDA_HANDLE_ERROR(cudaMalloc(matrices_d, 9*nb_cells*sizeof(cuda_scalar)));

  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_input_d->forces_ex,        3 * n_points * nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_input_d->points,           3 * n_points * nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_input_d->velocities,       3 * n_points * nb_cells * sizeof(ShapeOpScalar)));

  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->force_npd,         3 * n_points * nb_cells * sizeof(cuda_scalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->oldPoints,         3 * n_points * nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->momentum,          3 * n_points * nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->projections,            3 * idO * nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->points_last_iter,  3 * n_points * nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->nearest_normals,   3 * n_points * nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->nearest_points,    3 * n_points * nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->gradient,          3 * n_points * nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->gradient_prev,     3 * n_points * nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->q,                 3 * n_points * nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->descent_direction, 3 * n_points * nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->energies,                         nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->energies_prev,                    nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->has_converged,                    nb_cells * sizeof(int)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->descent_direction_dot_gradient,   nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->E_nonePD,         n_constraints * nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->s,       3 * n_points* MEM_SIZE * nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->t,       3 * n_points* MEM_SIZE * nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->rho,                   MEM_SIZE * nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->alpha,                 MEM_SIZE * nb_cells * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&simulation_data_d->center,						  3* nb_cells * sizeof(ShapeOpScalar)));
  //collision data
  CUDA_HANDLE_ERROR(cudaMalloc(&collision_data_d->nb_neighbours,  (nb_cells + 1)*sizeof(int)));
  // Matrices
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->N_inv_dense, n_points * n_points * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->M_tild_inv,  n_points * n_points * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->mass,        n_points * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->M_star.index, mesh_data_d->M_star.degree*n_points*sizeof(int)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->M_star.value, mesh_data_d->M_star.degree*n_points*sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->J.index,      mesh_data_d->J.degree*n_points*sizeof(int)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->J.value,      mesh_data_d->J.degree*n_points*sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->L.index,      mesh_data_d->L.degree*n_points*sizeof(int)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->L.value,      mesh_data_d->L.degree*n_points*sizeof(ShapeOpScalar)));
  
  // Allocate constraint data structures, this data is relative to the mesh not its instances
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->triangles,   3 *n_triangles   * sizeof(short)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->A,              n_constraints * sizeof(cuda_scalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->ConstraintType, n_constraints * sizeof(int)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->idO,            n_constraints * sizeof(int)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->rangeMin,       n_constraints * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->rangeMax,       n_constraints * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->Scalar1,        n_constraints * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->weight,         n_constraints * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->idI,            n_constraints * 4 * sizeof(int)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->vectorx,        n_constraints * 4 * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->matrix22,       n_constraints * 4 * sizeof(ShapeOpScalar)));
  CUDA_HANDLE_ERROR(cudaMalloc(&mesh_data_d->matrix33,       n_constraints * 9 * sizeof(ShapeOpScalar)));

  //printf("gpu alloc done \n");

  //printf("pointer %p %d\n", simulation_input_d->points, 3 * n_points * nb_cells * sizeof(ShapeOpScalar));
  // Fill memory ON GPU
  CUDA_HANDLE_ERROR(cudaMemcpy(simulation_input_d->forces_ex , simulation_input_h->forces_ex ,          3*n_points*sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
  CUDA_HANDLE_ERROR(cudaMemcpy(simulation_input_d->points    , simulation_input_h->points    , nb_cells*3*n_points*sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
  CUDA_HANDLE_ERROR(cudaMemcpy(simulation_input_d->velocities, simulation_input_h->velocities,          3*n_points*sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));

  //printf("pointer l %p %p \n", &mesh_data_d->L.value,  mesh_data_h->L.value);
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->triangles,    mesh_data_h->triangles, 3 * n_triangles * sizeof(short), cudaMemcpyHostToDevice));
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->L.index,      mesh_data_h->L.index,      mesh_data_h->L.degree*n_points*sizeof(int), cudaMemcpyHostToDevice));
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->L.value,      mesh_data_h->L.value,      mesh_data_h->L.degree*n_points*sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->J.index,      mesh_data_h->J.index,      mesh_data_h->J.degree*n_points*sizeof(int), cudaMemcpyHostToDevice));
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->J.value,      mesh_data_h->J.value,      mesh_data_h->J.degree*n_points*sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->M_star.index, mesh_data_h->M_star.index, mesh_data_h->M_star.degree*n_points*sizeof(int), cudaMemcpyHostToDevice));
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->M_star.value, mesh_data_h->M_star.value, mesh_data_h->M_star.degree*n_points*sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
  
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->M_tild_inv , mesh_data_h->M_tild_inv,  n_points * n_points * sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->N_inv_dense, mesh_data_h->N_inv_dense, n_points * n_points * sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->mass       , mesh_data_h->mass,        n_points * sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
  // Constraint data
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->A,              mesh_data_h->A,              n_constraints * sizeof(cuda_scalar), cudaMemcpyHostToDevice));
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->ConstraintType, mesh_data_h->ConstraintType, n_constraints * sizeof(int), cudaMemcpyHostToDevice));
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->idO,            mesh_data_h->idO,            n_constraints * sizeof(int), cudaMemcpyHostToDevice));
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->rangeMin,       mesh_data_h->rangeMin,       n_constraints * sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
  
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->rangeMax,       mesh_data_h->rangeMax,       n_constraints * sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->Scalar1,        mesh_data_h->Scalar1,        n_constraints * sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->weight,         mesh_data_h->weight,         n_constraints * sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
  //CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->E_nonePD,       mesh_data_h->E_nonePD,       n_constraints * sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
 
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->idI,            mesh_data_h->idI,            n_constraints * 4 * sizeof(int), cudaMemcpyHostToDevice));
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->vectorx,        mesh_data_h->vectorx,        n_constraints * 4 * sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->matrix22,       mesh_data_h->matrix22,       n_constraints * 4 * sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
  CUDA_HANDLE_ERROR(cudaMemcpy(mesh_data_d->matrix33,       mesh_data_h->matrix33,       n_constraints * 9 * sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
 
  //CUDA_HANDLE_ERROR(cudaDeviceSynchronize());
  //printf("copy stuff 2 done \n");
  //I m not even sure this is needed, but I prefer to  be safe.
  HANDLE_KERNEL_ERROR(set_to_zero<<<nb_cells*n_constraints/1025 + 1,1024>>>(simulation_data_d->E_nonePD, nb_cells * n_constraints));
  HANDLE_KERNEL_ERROR(set_to_zero<<<nb_cells, 3>>>(simulation_data_d->center, 3*nb_cells));
  HANDLE_KERNEL_ERROR(set_to_zero<<<3*nb_cells, n_points >>>(simulation_input_d->forces_ex,3 * n_points * nb_cells));
  HANDLE_KERNEL_ERROR(set_to_zero<<<3*nb_cells, n_points >>>(simulation_data_d->force_npd, 3 * n_points * nb_cells));
  HANDLE_KERNEL_ERROR(set_to_zero<<<3*nb_cells, n_points >>>(simulation_data_d->oldPoints, 3 * n_points * nb_cells));
  //HANDLE_KERNEL_ERROR(set_to_zero << <3 * nb_cells, n_points >> >(simulation_input_d->velocities, 3 * n_points * nb_cells));
  HANDLE_KERNEL_ERROR(set_to_zero<<<3*nb_cells, n_points >>>(simulation_data_d->points_last_iter, 3 * n_points * nb_cells));
  GPU_Mem_check();
}
///////////////////////////////////////////////////////////////////////////////

__device__ void project_d(int n_points, int n_constraints, int n_projected_points, int nb_cells,
		ShapeOpScalar *points_d, ShapeOpScalar *projections_d, cuda_scalar *f_int_nonePD_d,
		int *ConstraintType_d, int *idO_d, ShapeOpScalar *rangeMin_d, ShapeOpScalar *rangeMax_d, ShapeOpScalar *Scalar1_d, ShapeOpScalar *weight_d, ShapeOpScalar *E_nonePD_d,
		int *idI_d, ShapeOpScalar *vectorx_d, ShapeOpScalar *matrix22_d, ShapeOpScalar *matrix33_d, cuda_scalar *A_d, const float miu, const float lambda, const float kappa){
#ifdef SEQ

	int tid = 0;
	if (threadIdx.x == 0) {
#else
	int tid = threadIdx.x;
#endif
	while (tid < n_constraints) {
		//if(threadIdx.x == n_points - 1)printf("tid %d \n", ConstraintType_d[tid]);

		if (ConstraintType_d[tid] == 1) {
			project_Area(tid, n_points, n_constraints, n_projected_points, nb_cells,
						 points_d, projections_d, f_int_nonePD_d,
						 ConstraintType_d, idO_d, rangeMin_d, rangeMax_d, Scalar1_d, weight_d, E_nonePD_d,
						 idI_d, vectorx_d, matrix22_d, matrix33_d);

		} else if (ConstraintType_d[tid] == 3) {
			project_Bending(tid, n_points, n_constraints, n_projected_points, nb_cells,
							points_d, projections_d, f_int_nonePD_d,
							ConstraintType_d, idO_d, rangeMin_d, rangeMax_d, Scalar1_d, weight_d, E_nonePD_d,
							idI_d, vectorx_d, matrix22_d, matrix33_d);

		} else if (ConstraintType_d[tid] == 5) {
			project_surface_material(tid, n_points, n_constraints, n_projected_points, nb_cells,
									 points_d, projections_d, f_int_nonePD_d,
									 ConstraintType_d, idO_d, rangeMin_d, rangeMax_d, Scalar1_d, weight_d, E_nonePD_d,
									 idI_d, vectorx_d, matrix22_d, matrix33_d, A_d, miu, lambda, kappa);
		}																						
																								 
#ifdef SEQ
		tid++;
#else
		tid += blockDim.x;
#endif
	}
#ifdef SEQ
	}
#endif
}

///////////////////////////////////////////////////////////////////////////////
//TODO make it multicell
void points_from_Host_to_Device(int n_points, ShapeOpScalar *points_d, ShapeOpScalar *points_h, int cell){
	CUDA_HANDLE_ERROR(cudaMemcpy(points_d + 3*n_points*cell, points_h, 3 * n_points * sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
	//CUDA_HANDLE_ERROR(cudaMemcpyAsync(Palabos_Forces_d, Palabos_Forces_h, 3 * n_points * sizeof(ShapeOpScalar), cudaMemcpyHostToDevice, stream));
}
void external_forces_from_Host_to_Device(int n_points, ShapeOpScalar *Palabos_Forces_h, ShapeOpScalar *Palabos_Forces_d, cudaStream_t stream){
	 CUDA_HANDLE_ERROR(cudaMemcpy(Palabos_Forces_d, Palabos_Forces_h, 3 * n_points * sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
	 //CUDA_HANDLE_ERROR(cudaMemcpyAsync(Palabos_Forces_d, Palabos_Forces_h, 3 * n_points * sizeof(ShapeOpScalar), cudaMemcpyHostToDevice, stream));
}
///////////////////////////////////////////////////////////////////////////////
void points_from_Device_to_Host(int n_points, ShapeOpScalar *points_d, ShapeOpScalar *points_h, cudaStream_t stream){
	CUDA_HANDLE_ERROR(cudaMemcpy(points_h, points_d, 3 * n_points * sizeof(ShapeOpScalar), cudaMemcpyDeviceToHost));
	//CUDA_HANDLE_ERROR(cudaMemcpyAsync(points_h, points_d, 3 * n_points * sizeof(ShapeOpScalar), cudaMemcpyDeviceToHost, stream));
}
//////////////////////////////////////////////////////////////////////////////

void send_GPU_collinding_points(ShapeOpScalar *points, ShapeOpScalar **colid_points_d, ShapeOpScalar *normals, ShapeOpScalar **colid_normals_d, int *n_per_cell, int *n_per_cell_d, int nb, int nb_cell)
{
	//printf("nb_cell %d\n", nb_cell);
	CUDA_HANDLE_ERROR(cudaMemcpy(n_per_cell_d, n_per_cell, (nb_cell + 1)*sizeof(int), cudaMemcpyHostToDevice));

	CUDA_HANDLE_ERROR(cudaMalloc(colid_points_d, 3*nb*sizeof(ShapeOpScalar)));
	CUDA_HANDLE_ERROR(cudaMemcpy(*colid_points_d, points, 3*nb*sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));

	CUDA_HANDLE_ERROR(cudaMalloc(colid_normals_d, 3 * nb * sizeof(ShapeOpScalar)));
	CUDA_HANDLE_ERROR(cudaMemcpy(*colid_normals_d, normals, 3 * nb * sizeof(ShapeOpScalar), cudaMemcpyHostToDevice));
}
///////////////////////////////////////////////////////////////////////////////////////////
void free_GPU_pointer(void *pointer)
{
	cudaFree(pointer);
}
///////////////////////////////////////////////////////////////////////////////////////////
__launch_bounds__(258, 1)
__global__ void compute_next_frame_rbc_g(Mesh_info info, Mesh_data mesh, Simulation_input input, Simulation_data sim, Collision_data coll, ShapeOpScalar h, ShapeOpScalar h2, int it_max) {
	
	int line_search_it = 0, it = 0, tail = 0, head = 0, total_conv = 1, head_o = 0, x_n = info.n_points;

	int point_id = blockIdx.x*x_n*3 + threadIdx.x;

	__shared__ extern ShapeOpScalar buffer[];
	__shared__ ShapeOpScalar objectif[CONV_WINDOWS];
	__shared__ ShapeOpScalar objectif_sqr[CONV_WINDOWS];
	__shared__ ShapeOpScalar directional_derivative_prev, directional_derivative;

	//__shared__ ShapeOpScalar center[3];
	ShapeOpScalar newE = 0;
	cuda_scalar gradient_norm = 100;
	cuda_scalar *volume_der = (cuda_scalar*)(buffer + x_n);

	sim.points_last_iter[point_id        ] = input.points[point_id        ];
	sim.points_last_iter[point_id +   x_n] = input.points[point_id +   x_n];
	sim.points_last_iter[point_id + 2*x_n] = input.points[point_id + 2*x_n];

	#ifdef DEBUG
	if (threadIdx.x == 102 && blockIdx.x == BL) {
		printf("points %f %f %f | %d h %f\n", input.points[point_id], input.points[point_id + x_n], input.points[point_id + 2 * x_n], blockIdx.x, h);
	}
	#endif
	//center_d(input.points, x_n, volume_der, center, point_id);
	input.points[point_id        ] -= sim.center[3*blockIdx.x    ];
	input.points[point_id +   x_n] -= sim.center[3*blockIdx.x + 1];
	input.points[point_id + 2*x_n] -= sim.center[3*blockIdx.x + 2];

	#ifdef DEBUG
		if (threadIdx.x == 102 && blockIdx.x == BL) {
			printf("forces %f %f %f | %d \n", input.forces_ex[point_id], input.forces_ex[point_id + x_n], input.forces_ex[point_id + 2 * x_n], blockIdx.x);
			printf("velocity %f %f %f | %d \n", input.velocities[point_id], input.velocities[point_id + x_n], input.velocities[point_id + 2 * x_n], blockIdx.x);
			printf("points %f %f %f | %d \n", input.points[point_id], input.points[point_id + x_n], input.points[point_id + 2 * x_n], blockIdx.x);
		}
	#endif

	//printf("coll.total_neighbours %f \n", coll.total_neighbours);
	//if(threadIdx.x == 0)printf("nb voisin %d \n", coll.nb_neighbours[blockIdx.x + 1] - coll.nb_neighbours[blockIdx.x]);
	if(coll.nb_neighbours[blockIdx.x + 1] - coll.nb_neighbours[blockIdx.x]  > 0){
		nearest_neighbor_linear_d(input.points, coll.colid_points, coll.colid_normals, sim.nearest_points, sim.nearest_normals, sim.center + 3*blockIdx.x, 
								  coll.nb_neighbours[blockIdx.x], coll.nb_neighbours[blockIdx.x + 1], x_n, point_id, info.threshold_rep, info.threshold_nonRep);
	}
	__syncthreads();
	//compute momentum, set point = momentum
	//printf("%f %f %f |%d\n", input.points[point_id], input.points[point_id + x_n], input.points[point_id + 2 * x_n], point_id);
	momentum_points_first_guess3(input.points, input.velocities, mesh.M_tild_inv, input.forces_ex, sim.momentum, mesh.mass, h, x_n, point_id);

	#ifdef DEBUG
		if (threadIdx.x == 102 && blockIdx.x == BL)printf("momentum %f %f %f \n", sim.momentum[threadIdx.x], sim.momentum[threadIdx.x + x_n], sim.momentum[threadIdx.x + 2*x_n]);
	#endif
	__syncthreads();
	

	sim.force_npd[point_id        ] = 0;
	sim.force_npd[point_id +   x_n] = 0;
	sim.force_npd[point_id + 2*x_n] = 0;
	sim.E_nonePD[threadIdx.x + blockIdx.x*info.n_constraints] = 0;
	
	project_d(x_n, info.n_constraints, info.idO, 1, input.points, sim.projections, sim.force_npd,
			  mesh.ConstraintType, mesh.idO, mesh.rangeMin, mesh.rangeMax, mesh.Scalar1, mesh.weight, sim.E_nonePD,
			  mesh.idI, mesh.vectorx, mesh.matrix22, mesh.matrix33, mesh.A, info.miu, info.lambda, info.kappa);

	project_volume_d(x_n, info.n_triangles, info.n_constraints, input.points, mesh.triangles, sim.force_npd, sim.E_nonePD, info.volume, info.sum_points, buffer, volume_der, point_id, info.volume_weight);
	project_collision_d(x_n, info.nb_cells, input.points, sim.nearest_points, sim.nearest_normals, sim.force_npd, sim.E_nonePD, point_id, info.n_constraints, info.weight_col_rep, info.threshold_rep, info.weight_col_nonRep, info.threshold_nonRep, info.beta_morse);

	__syncthreads();

	objectif[head_o] = eval_objectif_and_grandient3(mesh.L, mesh.J, mesh.M_star, mesh.mass, input.points, sim.momentum, sim.projections, sim.E_nonePD, sim.force_npd,
													h*h, x_n, info.idO, info.n_constraints, info.sum_points, sim.gradient, buffer, volume_der, point_id);

	gradient_norm = volume_der[0];
	//todo test gradient norm
	objectif_sqr[head_o] = objectif[head_o]*objectif[head_o];

	sim.q[point_id        ] = -sim.gradient[point_id        ];
	sim.q[point_id +   x_n] = -sim.gradient[point_id +   x_n];
	sim.q[point_id + 2*x_n] = -sim.gradient[point_id + 2*x_n];
	
	sim.gradient_prev[point_id        ] = -sim.gradient[point_id        ];
	sim.gradient_prev[point_id +   x_n] = -sim.gradient[point_id +   x_n];
	sim.gradient_prev[point_id + 2*x_n] = -sim.gradient[point_id + 2*x_n];

	int attempt_stagn = 0;
	bool stagn = false;
//main loop
	while (it < it_max) {
		__syncthreads();

		sim.E_nonePD[threadIdx.x + blockIdx.x*info.n_constraints] = 0;

		if (gradient_norm <= TOL*TOL){
			//printf("grad == 0");
			break;
		}

#ifdef DEBUG
		if (threadIdx.x == 102 && blockIdx.x == BL) {
			printf("\nforce  npd %f %f %f | %d \n", sim.force_npd[point_id], sim.force_npd[point_id + x_n], sim.force_npd[point_id + 2 * x_n], blockIdx.x);
			printf("gradient %.17f %.17f %.17f \n", sim.gradient[point_id], sim.gradient[point_id + x_n], sim.gradient[point_id + 2 * x_n]);
			printf("Energy %.17f  it %d | %d \n", objectif[head_o], it, blockIdx.x);
		}
#endif
		if (total_conv == CONV_WINDOWS){
			ShapeOpScalar mean = 0, std = 0;
			for (int i = 0; i < CONV_WINDOWS; i++) {
				mean += objectif[i];
				std += objectif_sqr[i];
			}
			mean = mean/CONV_WINDOWS;
			std = (std + (5 - 2 * CONV_WINDOWS)*mean*mean)/(CONV_WINDOWS - 1);
			std = (std <= 0 ? 0 : sqrt((double)std));
#ifdef DEBUG
			if (threadIdx.x == 102 && blockIdx.x == BL)printf("Energy %f mean %f std %f  |%d\n", objectif[head_o], mean, std, it);
#endif
			if (abs((double)(std / mean)) <= TOL)break;
		}

		if (it > 0){
			comput_s_t_rho_d(input.points, sim.oldPoints, sim.gradient, sim.gradient_prev, sim.s, sim.t, sim.rho, head, x_n, info.nb_cells, info.sum_points, buffer, point_id);
			head++;
			if (head - tail > MEM_SIZE) {
				tail++;
			}
		}

		__syncthreads();
		copy_gradient_and_points3(input.points, sim.oldPoints, sim.gradient, sim.gradient_prev, sim.q, x_n, point_id);
		//LBFG optimisation
		//comput q according to s, t and rho
		__syncthreads();
		l_bfg_d(sim.s, sim.t, sim.rho, sim.alpha, sim.q, tail, head, x_n, info.nb_cells, info.sum_points, buffer, point_id);
		//descent_direction = inv(H)*q
		__syncthreads();
		mat_vect_mult3(mesh.N_inv_dense, sim.q, sim.descent_direction, x_n, point_id);
		//correct descent_direction accordind to s, t rho and alpha
		__syncthreads();
		l_bfg2_d(sim.s, sim.t, sim.rho, sim.alpha, sim.descent_direction, tail, head, x_n, info.nb_cells, info.sum_points, buffer, point_id);
		//armijo condition gradient = gradient(xk)
		__syncthreads();
		trace_x_time_y3(sim.gradient, sim.descent_direction, &directional_derivative_prev, x_n, point_id, info.sum_points, buffer);

		ShapeOpScalar a = 2;
		line_search_it = 0;
		while (line_search_it < MAX_LINE_SEARCH) {

			//if(threadIdx.x==0)printf("line %d \n", line_search_it);
			__syncthreads();
			a /= 2.;
	#ifdef DEBUG
			if (threadIdx.x == 102 && blockIdx.x == BL)printf("descent %.17f %.17f %.17f \n", sim.descent_direction[point_id], sim.descent_direction[point_id + x_n], sim.descent_direction[point_id + 2 * x_n]);
	#endif
			x__plus_equal_alpha_time_y3(sim.descent_direction, input.points, a, x_n, point_id);

			sim.force_npd[point_id          ] = 0;
			sim.force_npd[point_id +     x_n] = 0;
			sim.force_npd[point_id + 2 * x_n] = 0;
			sim.E_nonePD[threadIdx.x + blockIdx.x*info.n_constraints] = 0;

			__syncthreads();

			project_d(x_n, info.n_constraints, info.idO, 1,
				input.points, sim.projections, sim.force_npd,
				mesh.ConstraintType, mesh.idO, mesh.rangeMin, mesh.rangeMax, mesh.Scalar1, mesh.weight, sim.E_nonePD,
				mesh.idI, mesh.vectorx, mesh.matrix22, mesh.matrix33, mesh.A, info.miu, info.lambda, info.kappa);

			project_volume_d(x_n, info.n_triangles, info.n_constraints, input.points, mesh.triangles, sim.force_npd, sim.E_nonePD, info.volume, info.sum_points, buffer, volume_der, point_id, info.volume_weight);
			project_collision_d(x_n, info.nb_cells, input.points, sim.nearest_points, sim.nearest_normals, sim.force_npd, sim.E_nonePD, point_id, info.n_constraints, info.weight_col_rep, info.threshold_rep, info.weight_col_nonRep, info.threshold_nonRep, info.beta_morse);
			__syncthreads();
			//gradient = gradient(xk + 1);

			//armijo condition
			newE = eval_objectif_and_grandient3(mesh.L, mesh.J, mesh.M_star, mesh.mass, input.points, sim.momentum, sim.projections, sim.E_nonePD, sim.force_npd,
												h*h, x_n, info.idO, info.n_constraints, info.sum_points, sim.gradient, buffer, volume_der, point_id);
			gradient_norm = volume_der[0];

			//curvature_condition
			//curvature_test_d(sim.gradient, sim.descent_direction, sim.descent_direction_dot_gradient, sim.has_converged, 0.3, x_n, info.sum_points, point_id, buffer);
			trace_x_time_y3(sim.gradient, sim.descent_direction, &directional_derivative, x_n, point_id, info.sum_points, buffer);
			//if (threadIdx.x == 0)printf("line search criteria gamma %f  dot %f newE %f  objectif[head_o] %f %d\n", GAMMA, directional_derivative, newE, objectif[head_o], line_search_it);
			__syncthreads();

			line_search_it++;
		
			if (newE <= objectif[head_o] + a*GAMMA*directional_derivative_prev  && directional_derivative >= GAMMA2*directional_derivative_prev) { // 
				head_o = (head_o + 1) % CONV_WINDOWS;
				total_conv = (total_conv < CONV_WINDOWS ? total_conv + 1 : CONV_WINDOWS);
				objectif[head_o] = newE;
				objectif_sqr[head_o] = newE*newE;
				break;
			}
			#ifdef DEBUG
			if (threadIdx.x == 102 && blockIdx.x == BL) {
				printf("objectif tmp %.17f | directional_derivative_prev %.17f \n", newE, sim.descent_direction_dot_gradient[blockIdx.x]);
			}
			#endif
			
			if (line_search_it == MAX_LINE_SEARCH) {
				
				head_o = (head_o + 1) % CONV_WINDOWS;
				total_conv = (total_conv < CONV_WINDOWS ? total_conv + 1 : CONV_WINDOWS);
				objectif[head_o] = newE;
				objectif_sqr[head_o] = newE*newE;
				attempt_stagn++;
				break;
				//#ifdef DEBUG
				//if (threadIdx.x == 0 && line_search_it == MAX_LINE_SEARCH)printf("max Line srch\n");
				//#endif
				//goto end_loop;
			}
			
			input.points[point_id        ] = sim.oldPoints[point_id        ]; 
			input.points[point_id +   x_n] = sim.oldPoints[point_id +   x_n];
			input.points[point_id + 2*x_n] = sim.oldPoints[point_id + 2*x_n];

			sim.gradient[point_id        ] = sim.gradient_prev[point_id        ];
			sim.gradient[point_id +   x_n] = sim.gradient_prev[point_id +   x_n];
			sim.gradient[point_id + 2*x_n] = sim.gradient_prev[point_id + 2*x_n];
			/*
			if (line_search_it == MAX_LINE_SEARCH) {
				goto end_loop;
			}
			*/
		}
		it++;
		if(attempt_stagn == MAX_ATTEMPT) {
			stagn = true;
			break;
		}
	}
	end_loop:

	#ifdef DEBUG
		if(threadIdx.x == 0 && blockIdx.x == BL)printf("nb iter %d %d | %d \n", it, line_search_it, blockIdx.x);
	#endif
	//input.points[point_id       ] += center[0];
	//input.points[point_id +  x_n] += center[1];
	//input.points[point_id +2*x_n] += center[2];
	
	#ifdef DEBUG
		if (threadIdx.x == 102)printf("point GPU %f %f %f \n", input.points[point_id], input.points[point_id + x_n], input.points[point_id + 2*x_n]);
	#endif

	input.points[point_id       ] += sim.center[3*blockIdx.x    ];
	input.points[point_id +  x_n] += sim.center[3*blockIdx.x + 1];
	input.points[point_id +2*x_n] += sim.center[3*blockIdx.x + 2];
	
	sim.nearest_normals[point_id       ] = 0;
	sim.nearest_normals[point_id +  x_n] = 0;
	sim.nearest_normals[point_id +2*x_n] = 0;
	
	//velocities_equals_points_minus_points_prev3(input.points, sim.points_last_iter, input.velocities, 1 / h, x_n, point_id);
	damping( input.points, sim.points_last_iter, input.velocities, 1/h, h2, mesh.mass, 
			info.mass_total, x_n, info.sum_points, point_id, volume_der, sim.center + 3*blockIdx.x, info.calpha, info.vel_cap, info.vel_cap_fin, info.dx_p);

	input.points[point_id        ] = sim.points_last_iter[point_id        ] + input.velocities[point_id        ]*h;
	input.points[point_id +   x_n] = sim.points_last_iter[point_id +   x_n] + input.velocities[point_id +   x_n]*h;
	input.points[point_id + 2*x_n] = sim.points_last_iter[point_id + 2*x_n] + input.velocities[point_id + 2*x_n]*h;

	//DEBUG
	//if(threadIdx.x == 0)printf("vel GPU in kernel:\n%f \n%f \n%f \n", input.velocities[point_id], input.velocities[point_id + x_n] , input.velocities[point_id + 2*x_n]);
	
}

//////////////////////////////////////////////////////////////////////////////////////////
__global__ void points_time_mat_3X3d(const cuda_scalar *mat, double *points, const int n) {
	const cuda_scalar *m = mat + 9*blockIdx.x;
	double px, py, pz;
	px = points[threadIdx.x      ];
	py = points[threadIdx.x +   n];
	pz = points[threadIdx.x + 2*n];

	//if(threadIdx.x == 0)printf("ROTATION\n%f %f %f \n%f %f %f \n%f %f %f \n", m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]);

	points[threadIdx.x     ] = px*m[0] + py*m[3] + pz*m[6];
	points[threadIdx.x +  n] = px*m[1] + py*m[4] + pz*m[7];
	points[threadIdx.x +2*n] = px*m[2] + py*m[5] + pz*m[8];
}
//////////////////////////////////////////////////////////////////////////////////////////
void points_time_mat_3X3(cuda_scalar *mat_d, const cuda_scalar *mat_h, double *points, const int nb_cells, const int n) {
	CUDA_HANDLE_ERROR(cudaMemcpy(mat_d, mat_h, nb_cells*9*sizeof(cuda_scalar), cudaMemcpyHostToDevice));
	points_time_mat_3X3d <<<1, n >>>(mat_d, points, n);
	OTHER_ERROR;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
void __global__ make_periodic_g( ShapeOpScalar *points, ShapeOpScalar *center, float nx, float ny, float nz, int npoints, int it) {

		int dim = blockIdx.x%3;
		float limit = nx;
		if(dim == 1)limit = ny;
		if(dim == 2)limit = nz;

		ShapeOpScalar c = center[blockIdx.x];

		//if (blockIdx.x/3 == 13 && threadIdx.x==0)printf(" c %f | blockIdx %d \n", c, blockIdx.x);

		double s = 0;

		if(c > limit){
			s = -limit;
			//if(threadIdx.x == 0)printf("dim %d| cells %d |blockId %d|center %f | nx %f | befor %f after %f | %d\n",dim, blockIdx.x/3, blockIdx.x,  c, limit, points[blockIdx.x*npoints + threadIdx.x], points[blockIdx.x*npoints + threadIdx.x] + s, it);
		}
		if(c < 0)s = limit;
		points[blockIdx.x*npoints + threadIdx.x] += s;
		
		__syncthreads();
		if(threadIdx.x == 0)center[blockIdx.x]  = c + s; 

}
////////////////////////////////////////////////////////////////////////////////////////////////////////
void device_make_periodic(ShapeOpScalar * points_d, ShapeOpScalar *center, float nx, float ny, float nz, int ncells, int npoints, int it){
	HANDLE_KERNEL_ERROR(make_periodic_g<<<3*ncells, npoints>>>(points_d, center, nx, ny, nz, npoints, it));
}

//////////////////////////////////////////////////////////////////////////////////////////
void compute_next_frame_rbc(Mesh_info  *info, Mesh_data *mesh, Simulation_input *input,
	                        Simulation_data *sim, Collision_data *coll, ShapeOpScalar h, ShapeOpScalar h_2, int it_max, cudaStream_t stream){

	//std::cout << "wtf " << info->nb_cells << " "<< info->n_points  << " " << h << " " << h_2  << " it max " << it_max << std::endl;
	HANDLE_KERNEL_ERROR(compute_next_frame_rbc_g<<< info->nb_cells, info->n_points, info->n_points*(sizeof(double) + 3*sizeof(cuda_scalar)), stream >>>(*info,*mesh,*input, *sim, *coll, h, h_2, it_max));

}

void debug_matrix_from_gpu(double *mat_d, int n) {
    double *matrix = new double[n*n];

	CUDA_HANDLE_ERROR(cudaMemcpy( matrix, mat_d, n * n * sizeof(double), cudaMemcpyDeviceToHost));
	for (int i= 0; i < n; i ++){
		for (int j = 0; j < 15; j++) {
				printf("%f ", matrix[i + j*n]);
		}
		printf(" \n");
	}
    delete[] matrix;
	
}

}
}
///////////////////////////////////////////////////////////////////////////////
#endif // DEVICE_UTILITIES_CU
///////////////////////////////////////////////////////////////////////////////