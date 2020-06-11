///////////////////////////////////////////////////////////////////////////////
/* This file is distributed as part of the Palabos library.
 *
 * It has been adapted from a file of the ShapeOp library.
 * The ShapeOp library can be downloaded at the address https://www.shapeop.org.
 * It is governed by the terms of the Mozilla Public License v. 2.0.
 *
 * This file is subject to the terms of the Mozilla Public License v. 2.0.
 * If a copy of the MPL was not distributed with this file, 
 * you can obtain one at http://mozilla.org/MPL/2.0/.
 * 
 * Contact:
 * Christos Kotsalos
 * kotsaloscv@gmail.com
 * Computer Science Department
 * University of Geneva
 * 
 * The most recent release of Palabos can be downloaded at 
 * <https://palabos.unige.ch/>
*/
///////////////////////////////////////////////////////////////////////////////
// Remark on RowMajor:
// I have tried to use RowMajor as well, but there is an assertion failure
// for column vectors. Essentially, column vectors are not allowed by Eigen
// to be RowMajor (which makes sense).
// There is a work-around by using the ternary operator in the Options, but
// for now we go with the Eigen default which is ColMajor
///////////////////////////////////////////////////////////////////////////////
#ifndef TYPES_H
#define TYPES_H
///////////////////////////////////////////////////////////////////////////////
#include <Eigen3/Core>
#include <Eigen3/Dense>
#include <Eigen3/Sparse>
#include <Eigen3/Eigenvalues>
#include "nanoflann/nanoflann.hpp"

#include "common.h"
///////////////////////////////////////////////////////////////////////////////
/**
This file redefines EIGEN types using the scalar type ::ShapeOpScalar defined in
Common.h.*/
///////////////////////////////////////////////////////////////////////////////
/** Defines Eigen Alignment type.*/
#ifdef SHAPEOP_DONT_ALIGN
#define SHAPEOP_ALIGNMENT Eigen::DontAlign
#else
#define SHAPEOP_ALIGNMENT Eigen::AutoAlign
#endif
///////////////////////////////////////////////////////////////////////////////
/** \brief Namespace of the ShapeOp library.*/
namespace plb {
namespace npfem {
typedef ShapeOpScalar Scalar; // A scalar type, double or float, as defined in
// ::ShapeOpScalar in Common.h.
// Dense
// Biwise OR |: In our case is 0 | 0 => 0 => ColMajor
template <int Rows, int Cols,
    int Options = (Eigen::ColMajor | SHAPEOP_ALIGNMENT)>
using MatrixT = Eigen::Matrix<Scalar, Rows, Cols,
    Options>; // A typedef of the dense matrix of Eigen.
typedef MatrixT<2, 1> Vector2; // A 2d column vector.
typedef MatrixT<2, 2> Matrix22; // A 2 by 2 matrix.
typedef MatrixT<2, 3> Matrix23; // A 2 by 3 matrix.
typedef MatrixT<3, 1> Vector3; // A 3d column vector.
typedef MatrixT<3, 2> Matrix32; // A 3 by 2 matrix.
typedef MatrixT<3, 3> Matrix33; // A 3 by 3 matrix.
typedef MatrixT<3, 4> Matrix34; // A 3 by 4 matrix.
typedef MatrixT<4, 1> Vector4; // A 4d column vector.
typedef MatrixT<4, 4> Matrix44; // A 4 by 4 matrix.
typedef MatrixT<3, Eigen::Dynamic> Matrix3X; // A 3 by n matrix.
typedef MatrixT<Eigen::Dynamic, 3> MatrixX3; // A n by 3 matrix.
typedef MatrixT<Eigen::Dynamic, 1> VectorX; // A nd column vector.
typedef MatrixT<Eigen::Dynamic, Eigen::Dynamic> MatrixXX; // A n by m matrix.
typedef MatrixT<Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXXCuda; // A n by m matrix.
// Sparse
template <int Options = Eigen::ColMajor>
using SparseMatrixT = Eigen::SparseMatrix<Scalar,
    Options>; // A typedef of the sparse matrix of Eigen.
typedef SparseMatrixT<> SparseMatrix; // The default sparse matrix of Eigen.
typedef Eigen::Triplet<Scalar>
    Triplet; // A triplet, used in the sparse triplet representation for
// matrices.
///////////////////////////////////////////////////////////////////////////////
} // namespace npfem
} // namespace plb
///////////////////////////////////////////////////////////////////////////////
#endif // TYPES_H
///////////////////////////////////////////////////////////////////////////////
