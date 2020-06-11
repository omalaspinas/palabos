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
#include <stdio.h>

#include "sparse_matrix.h"

namespace plb {
namespace npfem {

void print_mat_sparse(sparse_matrix_cuda mat, int l, int n, int size) {

	for (int j = 0; j<l; j++) {
		for (int i = 0; i < size; i++) {
			printf("[%2d %2d: %.6f] ", i, mat.index[i + j*n], mat.value[i + j*n]);
		}
		printf("\n");
	}
}

sparse_matrix_cuda make_sparse_from_full(double *mat, int rows, int cols) {

	sparse_matrix_cuda out;

	double tol = 0.000;
	int degree = 0;

	for (int i = 0; i < rows; i++) {
		int k = 0;
		for (int j = 0; j < cols; j++) {
			int id = j*rows + i;
			if (mat[id] * mat[id] > tol) {
				//printf("%f ", mat[id]);
				k++;
			}
		}
		//printf("%d \n", k);

		if (k > degree) {
			degree = k;
		}
	}

	out.degree = degree;
	//printf("degree %d rows %d\n", out.degree, rows);
	out.value = new double[rows*degree]();
	out.index = new int[rows*degree]();

	for (int i = 0; i < rows; i++) {
		int k = 0;
		for (int j = 0; j < cols; j++) {
			int id = j*rows + i;
			//printf("j %d  \n", j);

			if (mat[id] * mat[id] > tol) {

				out.value[i + k*rows] = mat[id];
				out.index[i + k*rows] = j;
				//if(i==0)printf("mat %f | %d %d |  %d \n", out.value[i + k*n], i, out.index[i + k*n], i + k*n);
				k++;
				if (k >= degree) {
					break;
				}
			}
		}
	}
	return out;
}

}
}