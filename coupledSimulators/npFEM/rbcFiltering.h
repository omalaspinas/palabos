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
#pragma once

#include "palabos3D.h"
#include "palabos3D.hh"

namespace plb {
namespace npfem {

// Convert a matrix to a list (std::vector) of row vectors of the same size
//
// Template:
//   Mat  Matrix type, must implement:
//     .resize(m,n)
//     .row(i) = Row
//   T  type that can be safely cast to type in Mat via '='
// Inputs:
//   M  an m by n matrix
// Outputs:
//   V  a m-long list of vectors of size n
//
// See also: list_to_matrix
template <typename DerivedM>
inline void matrix_to_list(const Eigen::DenseBase<DerivedM>& M,
    std::vector<std::vector<typename DerivedM::Scalar>>& V)
{
    using namespace std;
    V.resize(M.rows(), vector<typename DerivedM::Scalar>(M.cols()));
    // loop over rows
    for (int i = 0; i < M.rows(); i++) {
        // loop over cols
        for (int j = 0; j < M.cols(); j++) {
            V[i][j] = M(i, j);
        }
    }
}

// Convert a matrix to a list (std::vector) of elements in column-major
// ordering.
//
// Inputs:
//    M  an m by n matrix
// Outputs:
//    V  an m*n list of elements
template <typename DerivedM>
inline void matrix_to_list(const Eigen::DenseBase<DerivedM>& M,
    std::vector<typename DerivedM::Scalar>& V)
{
    using namespace std;
    V.resize(M.size());
    // loop over cols then rows
    for (int j = 0; j < M.cols(); j++) {
        for (int i = 0; i < M.rows(); i++) {
            V[i + j * M.rows()] = M(i, j);
        }
    }
}

// Return wrapper
template <typename DerivedM>
inline std::vector<typename DerivedM::Scalar> matrix_to_list(
    const Eigen::DenseBase<DerivedM>& M)
{
    std::vector<typename DerivedM::Scalar> V;
    matrix_to_list(M, V);
    return V;
}
//////////////////////////////////////////////
void swap_local(double *a, double *b) {
	double tp = *a;
	*a = *b;
	*b = tp;
}

template <typename T>
T kselect(T *seq, int n, int k) {
	int i = 0;
	int start_n = n;
	//randome pivot
	//swap_local(seq[0], seq[(int)(((double)rand())/MAXINT*(n-1))]);

	double p = seq[0];
	int t = 0;
	int i_prev = 0;

	while (t < start_n*start_n) {
		t++;
		i_prev = i;
		//swap_local(seq + i, seq + i + (int)(((double)rand())/MAXINT*(n - 1 - i)));
		p = seq[i];
		//is seq[i]=p
		while (i + 1 < n && p >seq[i + 1]) {
			t++;
			swap_local(&seq[i], &seq[i + 1]);
			i++;
		}
		int j = i + 1;

		while (j<n) {

			if (p >= seq[j]) {
				t++;
				//printf("%d > %d\n", p, seq[j]);
				swap_local(&seq[i + 1], &seq[j]);
				swap_local(&seq[i + 1], &seq[i]);
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
double median_j(double *vector, int n) {
	if (n%2){
		return kselect(vector, n, n/2);
	}else{
		return (kselect(vector, n, n/2) + kselect(vector, n, n/2 - 1))/2;
	}
}


//////////////////////////////////////////////
// Compute the median of an eigen vector
//
// Inputs:
//   V  #V list of unsorted values
// Outputs:
//   m  median of those values
// Returns true on success, false on failure
template <typename DerivedV, typename mType>
inline bool median(const Eigen::MatrixBase<DerivedV>& V, mType& m)
{
    using namespace std;
    if (V.size() == 0) {
        return false;
    }
    vector<typename DerivedV::Scalar> vV;
    matrix_to_list(V, vV);
    // http://stackoverflow.com/a/1719155/148668
    size_t n = vV.size() / 2;
    nth_element(vV.begin(), vV.begin() + n, vV.end());
    if (vV.size() % 2 == 0) {
        nth_element(vV.begin(), vV.begin() + n - 1, vV.end());
        m = 0.5 * (vV[n] + vV[n - 1]);
    } else {
        m = vV[n];
    }
    return true;
}

}
}