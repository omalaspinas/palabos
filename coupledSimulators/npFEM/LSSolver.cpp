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
#ifndef LSSOLVER_CPP
#define LSSOLVER_CPP
///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "LSSolver.h"
///////////////////////////////////////////////////////////////////////////////
namespace plb {
namespace npfem {
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void SimplicialLDLTSolver::initialize(
    const SparseMatrix& A, unsigned int iteration)
{
    solver_.compute(A);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE VectorX SimplicialLDLTSolver::solve(
    const VectorX& b, const VectorX& x0) const
{
    return solver_.solve(b);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE Eigen::ComputationInfo SimplicialLDLTSolver::info() const
{
    return solver_.info();
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void CGSolver::initialize(
    const SparseMatrix& A, unsigned int iteration)
{
    solver_.compute(A);
    solver_.setMaxIterations(iteration);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE VectorX CGSolver::solve(
    const VectorX& b, const VectorX& x0) const
{
    return solver_.solveWithGuess(b, x0);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE Eigen::ComputationInfo CGSolver::info() const
{
    return solver_.info();
}
///////////////////////////////////////////////////////////////////////////////
/*
SHAPEOP_INLINE void MINRESSolver::initialize(
    const SparseMatrix& A, unsigned int iteration)
{
    solver_.compute(A);
    solver_.setMaxIterations(iteration);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE VectorX MINRESSolver::solve(
    const VectorX& b, const VectorX& x0) const
{
    return solver_.solveWithGuess(b, x0);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE Eigen::ComputationInfo MINRESSolver::info() const
{
    return solver_.info();
}
*/
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE SORSolver::SORSolver(plb::npfem::Scalar relaxation)
    : relaxation_(relaxation)
{
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void SORSolver::initialize(
    const SparseMatrix& A, unsigned int iteration)
{
    A_ = A;
    iteration_ = iteration;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE VectorX SORSolver::solve(
    const VectorX& b, const VectorX& x0) const
{
    VectorX x = x0;
    for (unsigned int k = 0; k < iteration_; ++k)
        for (int i = 0; i < A_.rows(); ++i) {
            x(i) += (relaxation_ / A_.coeff(i, i)) * (b(i) - A_.row(i).dot(x));
        }
    return x;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE Eigen::ComputationInfo SORSolver::info() const
{
    return Eigen::Success; // TODO this needs to be modified
}
///////////////////////////////////////////////////////////////////////////////
} // namespace npfem
} // namespace plb
///////////////////////////////////////////////////////////////////////////////
#endif // LSSOLVER_CPP
///////////////////////////////////////////////////////////////////////////////
