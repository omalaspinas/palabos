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
#ifndef LSSOLVER_H
#define LSSOLVER_H
///////////////////////////////////////////////////////////////////////////////
#include "Types.h"
///////////////////////////////////////////////////////////////////////////////
/* This file contains all the linear system solvers of the ShapeOp library. */
///////////////////////////////////////////////////////////////////////////////
namespace plb {
namespace npfem {
///////////////////////////////////////////////////////////////////////////////
/* 
  Base class of any sparse linear system solver. This class defines the
  main functionalities of the ShapeOp sparse linear system solvers (Ax = b).
*/
class SHAPEOP_API LSSolver {
public:
    virtual ~LSSolver(){};
    /* Initialize the linear system solver using the sparse matrix A. */
    virtual void initialize(const SparseMatrix& A, unsigned int iteration = 1)
        = 0;
    /* Solve the linear system Ax = b. */
    virtual VectorX solve(const VectorX& b, const VectorX& x0) const = 0;
    /* Reports whether previous computation was successful. */
    virtual Eigen::ComputationInfo info() const = 0;
};
///////////////////////////////////////////////////////////////////////////////
/* 
Sparse linear system solver based on Cholesky. This class implements
a sparse linear system solver based on the Cholesky LDL^T algorithm from
Eigen.
*/
class SHAPEOP_API SimplicialLDLTSolver : public LSSolver {
public:
    virtual ~SimplicialLDLTSolver(){};
    /* Prefactorize the sparse matrix (A = LDL^T). */
    virtual void initialize(
        const SparseMatrix& A, unsigned int iteration) override final;
    /* Solve the linear system by applying twice backsubstitution. */
    virtual VectorX solve(
        const VectorX& b, const VectorX& x0) const override final;
    /* Reports whether previous computation was successful. */
    virtual Eigen::ComputationInfo info() const override final;
private:
    Eigen::SimplicialLDLT<SparseMatrix> solver_;
};
///////////////////////////////////////////////////////////////////////////////
/*
Sparse linear system solver based on CG. This class implements a
sparse linear system solver based on the CG algorithm from Eigen.
*/
class SHAPEOP_API CGSolver : public LSSolver {
public:
    virtual ~CGSolver(){};
    /* Initialize PCG. */
    virtual void initialize(
        const SparseMatrix& A, unsigned int iteration) override final;
    /* Solve the linear system by applying CG. */
    virtual VectorX solve(
        const VectorX& b, const VectorX& x0) const override final;
    /* Reports whether previous computation was successful. */
    virtual Eigen::ComputationInfo info() const override final;
private:
    Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower,
        Eigen::IncompleteLUT<Scalar>>
        solver_;
};
///////////////////////////////////////////////////////////////////////////////
/*
Sparse linear system solver based on MINRES. This class implements a
sparse linear system solver based on the MINRES algorithm from Eigen.

/*
class SHAPEOP_API MINRESSolver : public LSSolver {
public:
    virtual ~MINRESSolver(){};
    // Initialize MINRES.
    virtual void initialize(
    const SparseMatrix& A, unsigned int iteration) override final;
    //Solve the linear system by applying MINRES.
    virtual VectorX solve(
    const VectorX& b, const VectorX& x0) const override final;
    // Reports whether previous computation was successful.
    virtual Eigen::ComputationInfo info() const override final;
private:
    Eigen::MINRES<SparseMatrix, Eigen::Lower, Eigen::IncompleteLUT<Scalar>>
    solver_;
};
*/
///////////////////////////////////////////////////////////////////////////////
/* Sparse linear system solver based on successive over-relaxation (SOR). */
class SHAPEOP_API SORSolver : public LSSolver {
public:
    SORSolver(plb::npfem::Scalar relaxation = 1.6);
    virtual ~SORSolver(){};
    /* Initialize SOR. */
    virtual void initialize(
        const SparseMatrix& A, unsigned int iteration) override final;
    /* Solve the linear system by applying SOR. */
    virtual VectorX solve(
        const VectorX& b, const VectorX& x0) const override final;
    /* Reports whether previous computation was successful. */
    virtual Eigen::ComputationInfo info() const override final;
private:
    SparseMatrixT<Eigen::RowMajor> A_;
    plb::npfem::Scalar relaxation_;
    unsigned int iteration_;
};
///////////////////////////////////////////////////////////////////////////////
} // namespace npfem
} // namespace plb
///////////////////////////////////////////////////////////////////////////////
#ifdef SHAPEOP_HEADER_ONLY
#include "LSSolver.cpp"
#endif
///////////////////////////////////////////////////////////////////////////////
#endif // LSSOLVER_H
///////////////////////////////////////////////////////////////////////////////
