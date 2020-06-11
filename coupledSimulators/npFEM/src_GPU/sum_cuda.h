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
#ifndef CUDA_SUM
#define CUDA_SUM

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

namespace plb {
namespace npfem {

template <class T>
__forceinline__ __device__ void sum(T *buffer, int n_sum, int n, int id) {


	for (int i = n_sum; i < n - id; i += n_sum) {
			buffer[id] += buffer[id + i];
	}
	int r = n_sum / 2;
	__syncthreads();

	while (r) {
		__syncthreads();
		if (id < r) {
			buffer[id] += buffer[id + r];
		}
		r /= 2;
	}
}

template <class T,  class T2>
__forceinline__ __device__ void double_sum(T *buffer, T2 *buffer2, int n_sum, int n, int id) {

	for (int i = n_sum; i < n - id; i += n_sum) {
		buffer[id] += buffer[id + i];
		buffer2[id] += buffer2[id + i];
	}
	int r = n_sum / 2;
	__syncthreads();
	while (r) {
		__syncthreads();
		if (id < r ) {
			buffer[id  ] += buffer[id + r ];

		}else if (id < 2*r) {
			buffer2[id - r] += buffer2[id];
		}
		r /= 2;
	}
	
}

template <class T>
__forceinline__ __device__ void triple_sum(T *buffer, T *buffer2,T *buffer3, int n_sum, int n, int id) {

	for (int i = n_sum; i < n - id; i += n_sum) {
			buffer[id]  += buffer[id + i];
			buffer2[id] += buffer2[id + i];
			buffer3[id] += buffer3[id + i];
	}
	int r = n_sum / 2;
	__syncthreads();
	if (id < r) {
		buffer[id] += buffer[id + r];
	}
	else if (id < 2 * r) {
		buffer2[id - r] += buffer2[id];
	}
	r /= 2;
	while (r) {
		__syncthreads();
		if (id < r) {
			buffer[id] += buffer[id + r];
		}else if (id < 2*r) {
			buffer2[id - r] += buffer2[id];
		}else if (id < 4*r) {
			buffer3[id - 2*r] += buffer3[id];
		}
		r /= 2;
	}
	__syncthreads();
	if (id == 0){
		buffer3[0] += buffer3[1];
	}
}

__global__ void test_sum(){

	__shared__ double buffer[516];
	buffer[threadIdx.x] = 1;
	buffer[threadIdx.x + 258] = 2;
	__syncthreads();
	double_sum(buffer, buffer + 258, 256, 258, threadIdx.x);
	__syncthreads();
	if(threadIdx.x == 0) printf("test sum %f %f",buffer[0], buffer[258]);
}

__global__ void test_sum2() {
	__shared__ double buffer1[258];
	__shared__ double buffer2[258];
	__shared__ double buffer3[258];
	buffer1[threadIdx.x] = threadIdx.x;
	buffer2[threadIdx.x] = threadIdx.x;
	buffer3[threadIdx.x] = threadIdx.x;
	__syncthreads();
	triple_sum(buffer1, buffer2, buffer3, 256, 258, threadIdx.x);
	__syncthreads();
	if (threadIdx.x == 0) printf("test sum %f %f %f", buffer1[0], buffer2[0], buffer3[0]);
}

}
}

#endif

