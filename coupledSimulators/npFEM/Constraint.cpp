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
#ifndef CONSTRAINT_CPP
#define CONSTRAINT_CPP
///////////////////////////////////////////////////////////////////////////////
#include <cassert>
#include <algorithm>
#include <iostream>

#include "Constraint.h"
#include "Solver.h"
#include "hyperelasticity.h"
///////////////////////////////////////////////////////////////////////////////
#define SHAPEOP_INNER_ITERATIONS 4
///////////////////////////////////////////////////////////////////////////////
#ifdef SHAPEOP_OPENMP
#ifdef SHAPEOP_MSVC
#define SHAPEOP_OMP_CRITICAL __pragma(omp critical)
#else
#define SHAPEOP_OMP_CRITICAL _Pragma("omp critical")
#endif
#else
#define SHAPEOP_OMP_CRITICAL
#endif
///////////////////////////////////////////////////////////////////////////////
namespace plb {
namespace npfem {
///////////////////////////////////////////////////////////////////////////////
/* Clamps v to lie between vMin and vMax.*/
SHAPEOP_INLINE Scalar clamp(Scalar v, Scalar vMin, Scalar vMax)
{
    Scalar result = v > vMin ? v : vMin;
    return result > vMax ? vMax : result;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE std::shared_ptr<Constraint> Constraint::shapeConstraintFactory(
    const std::string& constraintType, const std::vector<int>& idI,
    Scalar weight, const Matrix3X& positions)
{
    std::size_t n = idI.size();
    std::shared_ptr<Constraint> c;

    if (constraintType.compare("SurfaceMaterial") == 0) {
        if (n != 3) {
            return c;
        }
        return std::make_shared<SurfaceMaterialConstraint>(
            idI, weight, positions);
    }
    if (constraintType.compare("VolumeMaterial") == 0) {
        if (n != 4) {
            return c;
        }
        return std::make_shared<VolumeMaterialConstraint>(
            idI, weight, positions);
    }
    if (constraintType.compare("VolumeDamping") == 0) {
        if (n != 4) {
            return c;
        }
        return std::make_shared<VolumeDampingConstraint>(
            idI, weight, positions);
    }
    if (constraintType.compare("Collision") == 0) {
        if (n != 1) {
            return c;
        }
        return std::make_shared<CollisionConstraint>(idI, weight, positions);
    }
    if (constraintType.compare("TriangleARAP") == 0) {
        if (n != 3) {
            return c;
        }
        return std::make_shared<TriangleARAPConstraint>(idI, weight, positions);
    }
    if (constraintType.compare("TetrahedronARAP") == 0) {
        if (n != 4) {
            return c;
        }
        return std::make_shared<TetrahedronARAPConstraint>(
            idI, weight, positions);
    }
    if (constraintType.compare("EdgeStrainLimiting") == 0) {
        if (n != 2) {
            return c;
        }
        return std::make_shared<EdgeStrainLimitingConstraint>(
            idI, weight, positions);
    }
    if (constraintType.compare("TriangleStrainLimiting") == 0) {
        if (n != 3) {
            return c;
        }
        return std::make_shared<TriangleStrainLimitingConstraint>(
            idI, weight, positions);
    }
    if (constraintType.compare("TetrahedronStrainLimiting") == 0) {
        if (n != 4) {
            return c;
        }
        return std::make_shared<TetrahedronStrainLimitingConstraint>(
            idI, weight, positions);
    }
    if (constraintType.compare("Area") == 0) {
        if (n != 3) {
            return c;
        }
        return std::make_shared<AreaConstraint>(idI, weight, positions);
    }
    if (constraintType.compare("Volume") == 0) {
        if (n != 4) {
            return c;
        }
        return std::make_shared<VolumeConstraint>(idI, weight, positions);
    }
    if (constraintType.compare("Bending") == 0) {
        if (n != 4) {
            return c;
        }
        return std::make_shared<BendingConstraint>(idI, weight, positions);
    }
    if (constraintType.compare("Closeness") == 0) {
        if (n != 1) {
            return c;
        }
        return std::make_shared<ClosenessConstraint>(idI, weight, positions);
    }

    return c;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE Constraint::Constraint(const std::vector<int>& idI, Scalar weight)
    : idI_(idI)
    , weight_(std::sqrt(weight)) // It is inside the Frobenius norm
    , idO_(0)
    , E_nonePD_(0.)
    , PDSys_Build(true)
{
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE SurfaceMaterialConstraint::SurfaceMaterialConstraint(
    const std::vector<int>& idI, Scalar weight, const Matrix3X& positions,
    Scalar rangeMin, Scalar rangeMax, Scalar miu, Scalar lambda, Scalar kappa)
    : Constraint(idI, weight)
    , rangeMin_(rangeMin)
    , rangeMax_(rangeMax)
    , miu_(miu)
    , lambda_(lambda)
    , kappa_(kappa)
{
    assert(idI.size() == 3);
    ConstraintType_ = "SurfaceMaterial";
    Matrix32 edges, P;
    edges.col(0) = positions.col(idI_[1]) - positions.col(idI_[0]);
    edges.col(1) = positions.col(idI_[2]) - positions.col(idI_[0]);
    P.col(0) = edges.col(0).normalized();
    P.col(1)
        = (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
    // Reduce the problem to a plane (given the axes of the new system)
    rest_ = (P.transpose() * edges).inverse();
    A_ = std::abs((P.transpose() * edges).determinant() / 2.);
    weight_ *= std::sqrt(A_);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void SurfaceMaterialConstraint::calculateArea(
    const Matrix3X& positions, Scalar& area)
{
    Matrix32 edges, P;
    edges.col(0) = positions.col(idI_[1]) - positions.col(idI_[0]);
    edges.col(1) = positions.col(idI_[2]) - positions.col(idI_[0]);
    P.col(0) = edges.col(0).normalized();
    P.col(1) = (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
    Scalar A = std::abs((P.transpose() * edges).determinant() / 2.);

    SHAPEOP_OMP_CRITICAL { area += A; }
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void SurfaceMaterialConstraint::mass_lumping(
    const Matrix3X& positions, std::vector<Triplet>& triplets)
{
    Matrix32 edges, P;
    edges.col(0) = positions.col(idI_[1]) - positions.col(idI_[0]);
    edges.col(1) = positions.col(idI_[2]) - positions.col(idI_[0]);
    P.col(0) = edges.col(0).normalized();
    P.col(1)
        = (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
    Scalar A = std::abs((P.transpose() * edges).determinant() / 2.);

    // RBC membrane thickness (micro-m)
    Scalar thickness = 0.005;

    for (int i = 0; i < 3; ++i)
        triplets.push_back(Triplet(idI_[i], idI_[i], A * thickness / 3.));
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void SurfaceMaterialConstraint::project(
    const Matrix3X& positions, Matrix3X& projections, Matrix3X& f_int_nonePD,
    const Matrix3X& oldPositions)
{
    Matrix32 edges, P;
    edges.col(0) = (positions.col(idI_[1]) - positions.col(idI_[0]));
    edges.col(1) = (positions.col(idI_[2]) - positions.col(idI_[0]));
    P.col(0) = edges.col(0).normalized();
    P.col(1)
        = (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
    Matrix22 F = P.transpose() * edges * rest_;

    Eigen::JacobiSVD<Matrix22> svd(
        F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Vector2 S = svd.singularValues();
    Matrix22 U = svd.matrixU();
    Matrix22 V = svd.matrixV();

    // Guarantees that U & V are SO(2) (see my notes for more)
    Scalar detU = U.determinant();
    Scalar detV = V.determinant();
    if (detU < 0) {
        U.block<2, 1>(0, 1) *= -1.;
        S(1) *= -1.;
    }
    if (detV < 0) {
        V.block<2, 1>(0, 1) *= -1.;
        S(1) *= -1.;
    }

    // Principal Stretches
    Scalar l1 = S(0), l2 = S(1);
    Scalar dPsi_dl1, dPsi_dl2;

    dPsi_dl1 = f_prime_tr(l1, miu_, lambda_, kappa_)
        + g_prime_tr(l1 * l2, miu_, lambda_, kappa_) * l2;
    dPsi_dl2 = f_prime_tr(l2, miu_, lambda_, kappa_)
        + g_prime_tr(l1 * l2, miu_, lambda_, kappa_) * l1;

    Matrix22 Piola_hat = Matrix22::Zero();
    Piola_hat(0, 0) = dPsi_dl1;
    Piola_hat(1, 1) = dPsi_dl2;
    Matrix22 Piola = U * Piola_hat * V.transpose();

    // Check Bender 2015 & Irving 2004
    // In PD the weight of an energy is A or V * k, where k is user defined
    // In none PD energies, the energy comes from the integation of the density
    // function and thus the A or V comes naturally. The k is defined through
    // the fitting as described by Liu.
    Matrix32 H = -0.5 * A_ * P * Piola * rest_.transpose();

    SHAPEOP_OMP_CRITICAL
    {
        f_int_nonePD.col(idI_[0]) += -H.col(0) - H.col(1);
        f_int_nonePD.col(idI_[1]) += H.col(0);
        f_int_nonePD.col(idI_[2]) += H.col(1);
    }

    // Energy calculation
    Scalar Psi = f_tr(l1, miu_, lambda_, kappa_)
        + f_tr(l2, miu_, lambda_, kappa_)
        + g_tr(l1 * l2, miu_, lambda_, kappa_);
    E_nonePD_ = 0.5 * A_ * Psi;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void SurfaceMaterialConstraint::addConstraint(
    std::vector<Triplet>& triplets, int& idO) const
{
    // None PD Energy
    // idO_ = idO;
    int n = 2;
    for (int i = 0; i < n; ++i) {
        triplets.push_back(
            Triplet(idO + i, idI_[0], -weight_ * (rest_(0, i) + rest_(1, i))));
        triplets.push_back(Triplet(idO + i, idI_[1], weight_ * rest_(0, i)));
        triplets.push_back(Triplet(idO + i, idI_[2], weight_ * rest_(1, i)));
    }
    idO += n;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE VolumeMaterialConstraint::VolumeMaterialConstraint(
    const std::vector<int>& idI, Scalar weight, const Matrix3X& positions,
    Scalar rangeMin, Scalar rangeMax, Scalar miu, Scalar lambda, Scalar kappa)
    : Constraint(idI, weight)
    , rangeMin_(rangeMin)
    , rangeMax_(rangeMax)
    , miu_(miu)
    , lambda_(lambda)
    , kappa_(kappa)
{
    assert(idI_.size() == 4);
    ConstraintType_ = "VolumeMaterial";
    Matrix33 edges;
    for (int i = 0; i < 3; ++i)
        edges.col(i) = positions.col(idI_[i + 1]) - positions.col(idI_[0]);
    rest_ = edges.inverse();
    V_ = std::abs((edges).determinant() / 6.);
    weight_ *= std::sqrt(V_);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void VolumeMaterialConstraint::calculateVolume(
    const Matrix3X& positions, Scalar& volume)
{
    Matrix33 edges;
    for (int i = 0; i < 3; ++i)
        edges.col(i) = positions.col(idI_[i + 1]) - positions.col(idI_[0]);
    Scalar V = std::abs((edges).determinant() / 6.);

    SHAPEOP_OMP_CRITICAL { volume += V; }
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void VolumeMaterialConstraint::mass_lumping(
    const Matrix3X& positions, std::vector<Triplet>& triplets)
{
    Matrix33 edges;
    for (int i = 0; i < 3; ++i)
        edges.col(i) = positions.col(idI_[i + 1]) - positions.col(idI_[0]);
    Scalar V = std::abs((edges).determinant() / 6.);

    for (int i = 0; i < 4; ++i)
        triplets.push_back(Triplet(idI_[i], idI_[i], V / 4.));
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void VolumeMaterialConstraint::project(const Matrix3X& positions,
    Matrix3X& projections, Matrix3X& f_int_nonePD, const Matrix3X& oldPositions)
{
    Matrix33 edges;
    for (int i = 0; i < 3; ++i)
        edges.col(i) = positions.col(idI_[i + 1]) - positions.col(idI_[0]);
    Matrix33 F = edges * rest_;

    Eigen::JacobiSVD<Matrix33> svd(
        F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Vector3 S = svd.singularValues();
    Matrix33 U = svd.matrixU();
    Matrix33 V = svd.matrixV();

    // For rangeMin_ = rangeMax_ = 1., d stays zero vector and the projection is
    // simply the deformation gradient!
    Vector3 d(0., 0., 0.);
    for (int i = 0; i < SHAPEOP_INNER_ITERATIONS; ++i) {
        // To keep the volume constant, you simply set l1*l2*l3 = 1 (li:
        // principal stretches)
        Scalar v = S(0) * S(1) * S(2);
        Scalar f = v - clamp(v, rangeMin_, rangeMax_);
        Vector3 g(S(1) * S(2), S(0) * S(2), S(0) * S(1));
        d = -((f - g.dot(d)) / g.dot(g)) * g;
        S = svd.singularValues() + d;
    }
    if (svd.matrixU().determinant() * svd.matrixV().determinant() < 0.)
        S(2) = -S(2);

    F = U * S.asDiagonal() * V.transpose();

    projections.block<3, 3>(0, idO_) = weight_ * F;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void VolumeMaterialConstraint::addConstraint(
    std::vector<Triplet>& triplets, int& idO) const
{
    // Firstly, we build the PD system
    if (PDSys_Build) {
        idO_ = idO;
        PDSys_Build = false;
    }
    
    int n = 3;
    for (int i = 0; i < n; ++i) {
        triplets.push_back(Triplet(idO + i, idI_[0],
            -weight_ * (rest_(0, i) + rest_(1, i) + rest_(2, i))));
        triplets.push_back(Triplet(idO + i, idI_[1], weight_ * rest_(0, i)));
        triplets.push_back(Triplet(idO + i, idI_[2], weight_ * rest_(1, i)));
        triplets.push_back(Triplet(idO + i, idI_[3], weight_ * rest_(2, i)));
    }
    idO += n;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE VolumeDampingConstraint::VolumeDampingConstraint(
    const std::vector<int>& idI, Scalar weight, const Matrix3X& positions,
    Scalar rangeMin, Scalar rangeMax, Scalar miu, Scalar lambda, Scalar kappa)
    : Constraint(idI, weight)
    , rangeMin_(rangeMin)
    , rangeMax_(rangeMax)
    , miu_(miu)
    , lambda_(lambda)
    , kappa_(kappa)
{
    assert(idI_.size() == 4);
    ConstraintType_ = "VolumeDamping";
    Matrix33 edges;
    for (int i = 0; i < 3; ++i)
        edges.col(i) = positions.col(idI_[i + 1]) - positions.col(idI_[0]);
    rest_ = edges.inverse();
    V_ = std::abs((edges).determinant() / 6.);
    weight_ *= std::sqrt(V_);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void VolumeDampingConstraint::project(const Matrix3X& positions,
    Matrix3X& projections, Matrix3X& f_int_nonePD, const Matrix3X& oldPositions)
{
    Matrix33 edges, rest;
    for (int i = 0; i < 3; ++i)
        edges.col(i)
            = oldPositions.col(idI_[i + 1]) - oldPositions.col(idI_[0]);
    rest = edges.inverse();
    Scalar Volume = std::abs((edges).determinant() / 6.);

    for (int i = 0; i < 3; ++i)
        edges.col(i) = positions.col(idI_[i + 1]) - positions.col(idI_[0]);

    Matrix33 F = edges * rest;

    Eigen::JacobiSVD<Matrix33> svd(
        F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Vector3 S = svd.singularValues();
    Matrix33 U = svd.matrixU();
    Matrix33 V = svd.matrixV();

    // Guarantees that U & V are SO(3) and inversions are handle correctly
    Scalar detU = U.determinant();
    Scalar detV = V.determinant();
    if (detU < 0) {
        U.block<3, 1>(0, 2) *= -1;
        S(2) *= -1;
    }
    if (detV < 0) {
        V.block<3, 1>(0, 2) *= -1;
        S(2) *= -1;
    }

    Scalar l1 = S(0), l2 = S(1), l3 = S(2);
    Scalar dPsi_dl1, dPsi_dl2, dPsi_dl3;

    dPsi_dl1 = f_prime_tet(l1, miu_, lambda_, kappa_)
        + g_prime_tet(l1 * l2, miu_, lambda_, kappa_) * l2
        + g_prime_tet(l1 * l3, miu_, lambda_, kappa_) * l3
        + h_prime_tet(l1 * l2 * l3, miu_, lambda_, kappa_) * l2 * l3;

    dPsi_dl2 = f_prime_tet(l2, miu_, lambda_, kappa_)
        + g_prime_tet(l1 * l2, miu_, lambda_, kappa_) * l1
        + g_prime_tet(l2 * l3, miu_, lambda_, kappa_) * l3
        + h_prime_tet(l1 * l2 * l3, miu_, lambda_, kappa_) * l1 * l3;

    dPsi_dl3 = f_prime_tet(l3, miu_, lambda_, kappa_)
        + g_prime_tet(l2 * l3, miu_, lambda_, kappa_) * l2
        + g_prime_tet(l1 * l3, miu_, lambda_, kappa_) * l1
        + h_prime_tet(l1 * l2 * l3, miu_, lambda_, kappa_) * l1 * l2;

    Matrix33 Piola_hat = Matrix33::Zero();
    Piola_hat(0, 0) = dPsi_dl1;
    Piola_hat(1, 1) = dPsi_dl2;
    Piola_hat(2, 2) = dPsi_dl3;
    Matrix33 Piola = U * Piola_hat * V.transpose();

    Matrix33 H = -0.5 * Volume * Piola * rest.transpose();

    SHAPEOP_OMP_CRITICAL
    {
        f_int_nonePD.col(idI_[0]) += -H.col(0) - H.col(1) - H.col(2);
        f_int_nonePD.col(idI_[1]) += H.col(0);
        f_int_nonePD.col(idI_[2]) += H.col(1);
        f_int_nonePD.col(idI_[3]) += H.col(2);
    }

    // Energy calculation
    Scalar Psi = f_tet(l1, miu_, lambda_, kappa_)
        + f_tet(l2, miu_, lambda_, kappa_) + f_tet(l3, miu_, lambda_, kappa_)
        + g_tet(l1 * l2, miu_, lambda_, kappa_)
        + g_tet(l2 * l3, miu_, lambda_, kappa_)
        + g_tet(l3 * l1, miu_, lambda_, kappa_)
        + h_tet(l1 * l2 * l3, miu_, lambda_, kappa_);

    E_nonePD_ = 0.5 * Volume * Psi;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void VolumeDampingConstraint::addConstraint(
    std::vector<Triplet>& triplets, int& idO) const
{
    int n = 3;
    for (int i = 0; i < n; ++i) {
        triplets.push_back(Triplet(idO + i, idI_[0],
            -weight_ * (rest_(0, i) + rest_(1, i) + rest_(2, i))));
        triplets.push_back(Triplet(idO + i, idI_[1], weight_ * rest_(0, i)));
        triplets.push_back(Triplet(idO + i, idI_[2], weight_ * rest_(1, i)));
        triplets.push_back(Triplet(idO + i, idI_[3], weight_ * rest_(2, i)));
    }
    idO += n;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE CollisionConstraint::CollisionConstraint(
    const std::vector<int>& idI, Scalar weight, const Matrix3X& positions)
    : Constraint(idI, weight)
{
    assert(idI.size() == 1);
    ConstraintType_ = "Collision";
    thisWeight_ = weight;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void CollisionConstraint::setCollindingPoint(
    const Vector3& position){ collidingPoint_ = position; }
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void CollisionConstraint::setCollindingPointNormal(
    const Vector3& normal) { collidingPointNormal_ = normal; }
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void CollisionConstraint::setParams(const Scalar& collisions_threshold_rep, const Scalar& collisionsWeight_rep,
    const Scalar& collisions_threshold_nonRep, const Scalar& collisionsWeight_nonRep,
    const Scalar& beta_morse)
{
    collisions_threshold_rep_ = collisions_threshold_rep;
    collisionsWeight_rep_ = collisionsWeight_rep;
    collisions_threshold_nonRep_ = collisions_threshold_nonRep;
    collisionsWeight_nonRep_ = collisionsWeight_nonRep;
    beta_morse_ = beta_morse;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void CollisionConstraint::project(const Matrix3X& positions,
    Matrix3X& projections, Matrix3X& f_int_nonePD, const Matrix3X& oldPositions)
{
    E_nonePD_ = 0.;
    Vector3 point = positions.col(idI_[0]);

    Vector3 CP = collidingPoint_ - point;
    Vector3 PC = point - collidingPoint_;

    if (PC.dot(collidingPointNormal_) < 0.)
    {
        // Repulsion
        if (collidingPointNormal_.norm() >= 2.)
        {
            SHAPEOP_OMP_CRITICAL
            {
                f_int_nonePD.col(idI_[0]) += collisionsWeight_rep_ * CP;
            }

            E_nonePD_ = 0.5 * collisionsWeight_rep_ * CP.squaredNorm();
        }
        // NO Repulsion
        else
        {
            SHAPEOP_OMP_CRITICAL
            {
                f_int_nonePD.col(idI_[0]) += collisionsWeight_nonRep_ * CP;
            }

            E_nonePD_ = 0.5 * collisionsWeight_nonRep_ * CP.squaredNorm();
        }
    }
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void CollisionConstraint::addConstraint(
    std::vector<Triplet>& triplets, int& idO) const
{
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* Below are the original ShapeOp Energies (Quadratic)                       */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE TriangleARAPConstraint::TriangleARAPConstraint(
    const std::vector<int>& idI, Scalar weight, const Matrix3X& positions)
    : Constraint(idI, weight)
{
    assert(idI.size() == 3);
    ConstraintType_ = "TriangleARAP";
    Matrix32 edges, P;
    edges.col(0) = positions.col(idI_[1]) - positions.col(idI_[0]);
    edges.col(1) = positions.col(idI_[2]) - positions.col(idI_[0]);
    P.col(0) = edges.col(0).normalized();
    P.col(1)
        = (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
    rest_ = (P.transpose() * edges).inverse();
    Scalar A = (P.transpose() * edges).determinant() / 2.;
    weight_ *= std::sqrt(std::abs(A));
}
///////////////////////////////////////////////////////////////////////////////
/* Aq = F & Bp = R( SO(3) ) */
SHAPEOP_INLINE void TriangleARAPConstraint::project(const Matrix3X& positions,
    Matrix3X& projections, Matrix3X& f_int_nonePD, const Matrix3X& oldPositions)
{
    Matrix32 edges, P;
    edges.col(0) = (positions.col(idI_[1]) - positions.col(idI_[0]));
    edges.col(1) = (positions.col(idI_[2]) - positions.col(idI_[0]));
    P.col(0) = edges.col(0).normalized();
    P.col(1)
        = (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
    Matrix22 F = P.transpose() * edges * rest_;
    Eigen::JacobiSVD<Matrix22> svd(
        F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    projections.block<3, 2>(0, idO_)
        = (weight_ * P * svd.matrixU() * svd.matrixV().transpose());
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void TriangleARAPConstraint::addConstraint(
    std::vector<Triplet>& triplets, int& idO) const
{
    // Firstly, we build the PD system
    if (PDSys_Build) {
        idO_ = idO;
        PDSys_Build = false;
    }

    int n = 2;
    for (int i = 0; i < n; ++i) {
        triplets.push_back(
            Triplet(idO + i, idI_[0], -weight_ * (rest_(0, i) + rest_(1, i))));
        triplets.push_back(Triplet(idO + i, idI_[1], weight_ * rest_(0, i)));
        triplets.push_back(Triplet(idO + i, idI_[2], weight_ * rest_(1, i)));
    }
    idO += n;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE TetrahedronARAPConstraint::TetrahedronARAPConstraint(
    const std::vector<int>& idI, Scalar weight, const Matrix3X& positions)
    : Constraint(idI, weight)
{
    assert(idI.size() == 4);
    ConstraintType_ = "TetrahedronARAP";
    Matrix33 edges;
    for (int i = 0; i < 3; ++i)
        edges.col(i) = positions.col(idI_[i + 1]) - positions.col(idI_[0]);
    rest_ = edges.inverse();
    Scalar V = (edges).determinant() / 6.;
    weight_ *= std::sqrt(std::abs(V));
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void TetrahedronARAPConstraint::project(
    const Matrix3X& positions, Matrix3X& projections, Matrix3X& f_int_nonePD,
    const Matrix3X& oldPositions)
{
    Matrix33 edges;
    for (int i = 0; i < 3; ++i)
        edges.col(i) = positions.col(idI_[i + 1]) - positions.col(idI_[0]);
    Matrix33 F = edges * rest_;
    Eigen::JacobiSVD<Matrix33> svd(
        F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    projections.block<3, 3>(0, idO_)
        = (weight_ * svd.matrixU() * svd.matrixV().transpose());
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void TetrahedronARAPConstraint::addConstraint(
    std::vector<Triplet>& triplets, int& idO) const
{
    // Firstly, we build the PD system
    if (PDSys_Build) {
        idO_ = idO;
        PDSys_Build = false;
    }

    int n = 3;
    for (int i = 0; i < n; ++i) {
        triplets.push_back(Triplet(idO + i, idI_[0],
            -weight_ * (rest_(0, i) + rest_(1, i) + rest_(2, i))));
        triplets.push_back(Triplet(idO + i, idI_[1], weight_ * rest_(0, i)));
        triplets.push_back(Triplet(idO + i, idI_[2], weight_ * rest_(1, i)));
        triplets.push_back(Triplet(idO + i, idI_[3], weight_ * rest_(2, i)));
    }
    idO += n;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE EdgeStrainLimitingConstraint::EdgeStrainLimitingConstraint(
    const std::vector<int>& idI, Scalar weight, const Matrix3X& positions,
    Scalar rangeMin, Scalar rangeMax)
    : Constraint(idI, weight)
    , rangeMin_(rangeMin)
    , rangeMax_(rangeMax)
{
    assert(idI.size() == 2);
    ConstraintType_ = "EdgeStrainLimiting";
    Scalar length = (positions.col(idI_[1]) - positions.col(idI_[0])).norm();
    rest_ = 1. / length;
    weight_ *= std::sqrt(length);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void EdgeStrainLimitingConstraint::project(
    const Matrix3X& positions, Matrix3X& projections, Matrix3X& f_int_nonePD,
    const Matrix3X& oldPositions)
{
    Vector3 edge = positions.col(idI_[1]) - positions.col(idI_[0]);
    Scalar l = edge.norm();
    edge /= l;
    l = clamp(l * rest_, rangeMin_, rangeMax_);
    projections.col(idO_) = weight_ * l * edge;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void EdgeStrainLimitingConstraint::addConstraint(
    std::vector<Triplet>& triplets, int& idO) const
{
    // Firstly, we build the PD system
    if (PDSys_Build) {
        idO_ = idO;
        PDSys_Build = false;
    }

    // The triplets add elements on a Sparse matrix as (row, col, value)
    triplets.push_back(Triplet(idO, idI_[0], -weight_ * rest_));
    triplets.push_back(Triplet(idO, idI_[1], weight_ * rest_));
    idO += 1;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void EdgeStrainLimitingConstraint::setEdgeLength(Scalar length)
{
    rest_ = 1. / length;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE
TriangleStrainLimitingConstraint::TriangleStrainLimitingConstraint(
    const std::vector<int>& idI, Scalar weight, const Matrix3X& positions,
    Scalar rangeMin, Scalar rangeMax)
    : Constraint(idI, weight)
    , rangeMin_(rangeMin)
    , rangeMax_(rangeMax)
{
    assert(idI.size() == 3);
    ConstraintType_ = "TriangleStrainLimiting";
    Matrix32 edges, P;
    edges.col(0) = positions.col(idI_[1]) - positions.col(idI_[0]);
    edges.col(1) = positions.col(idI_[2]) - positions.col(idI_[0]);
    P.col(0) = edges.col(0).normalized();
    P.col(1)
        = (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
    rest_ = (P.transpose() * edges).inverse();
    Scalar A = (P.transpose() * edges).determinant() / 2.;
    weight_ *= std::sqrt(std::abs(A));
}
///////////////////////////////////////////////////////////////////////////////
/* Aq = F & Bp = F_clamped! */
SHAPEOP_INLINE void TriangleStrainLimitingConstraint::project(
    const Matrix3X& positions, Matrix3X& projections, Matrix3X& f_int_nonePD,
    const Matrix3X& oldPositions)
{
    Matrix32 edges, P;
    edges.col(0) = (positions.col(idI_[1]) - positions.col(idI_[0]));
    edges.col(1) = (positions.col(idI_[2]) - positions.col(idI_[0]));
    P.col(0) = edges.col(0).normalized();
    P.col(1)
        = (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
    Matrix22 F = P.transpose() * edges * rest_;
    Eigen::JacobiSVD<Matrix22> svd(
        F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Vector2 S = svd.singularValues();
    // Limits the principal stretches
    S(0) = clamp(S(0), rangeMin_, rangeMax_);
    S(1) = clamp(S(1), rangeMin_, rangeMax_);
    F = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();
    projections.block<3, 2>(0, idO_) = (weight_ * P * F);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void TriangleStrainLimitingConstraint::addConstraint(
    std::vector<Triplet>& triplets, int& idO) const
{
    // Firstly, we build the PD system
    if (PDSys_Build) {
        idO_ = idO;
        PDSys_Build = false;
    }

    int n = 2;
    for (int i = 0; i < n; ++i) {
        triplets.push_back(
            Triplet(idO + i, idI_[0], -weight_ * (rest_(0, i) + rest_(1, i))));
        triplets.push_back(Triplet(idO + i, idI_[1], weight_ * rest_(0, i)));
        triplets.push_back(Triplet(idO + i, idI_[2], weight_ * rest_(1, i)));
    }
    idO += n;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE
TetrahedronStrainLimitingConstraint::TetrahedronStrainLimitingConstraint(
    const std::vector<int>& idI, Scalar weight, const Matrix3X& positions,
    Scalar rangeMin, Scalar rangeMax)
    : Constraint(idI, weight)
    , rangeMin_(rangeMin)
    , rangeMax_(rangeMax)
{
    assert(idI.size() == 4);
    ConstraintType_ = "TetrahedronStrainLimiting";
    Matrix33 edges;
    for (int i = 0; i < 3; ++i)
        edges.col(i) = positions.col(idI_[i + 1]) - positions.col(idI_[0]);
    rest_ = edges.inverse();
    Scalar V = (edges).determinant() / 6.;
    weight_ *= std::sqrt(std::abs(V));
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void TetrahedronStrainLimitingConstraint::project(
    const Matrix3X& positions, Matrix3X& projections, Matrix3X& f_int_nonePD,
    const Matrix3X& oldPositions)
{
    Matrix33 edges;
    for (int i = 0; i < 3; ++i)
        edges.col(i) = positions.col(idI_[i + 1]) - positions.col(idI_[0]);
    Matrix33 F = edges * rest_;
    Eigen::JacobiSVD<Matrix33> svd(
        F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Vector3 S = svd.singularValues();
    S(0) = clamp(S(0), rangeMin_, rangeMax_);
    S(1) = clamp(S(1), rangeMin_, rangeMax_);
    S(2) = clamp(S(2), rangeMin_, rangeMax_);
    if (svd.matrixU().determinant() * svd.matrixV().determinant() < 0.)
        S(2) = -S(2);
    F = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();
    projections.block<3, 3>(0, idO_) = weight_ * F;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void TetrahedronStrainLimitingConstraint::addConstraint(
    std::vector<Triplet>& triplets, int& idO) const
{
    // Firstly, we build the PD system
    if (PDSys_Build) {
        idO_ = idO;
        PDSys_Build = false;
    }

    int n = 3;
    for (int i = 0; i < n; ++i) {
        triplets.push_back(Triplet(idO + i, idI_[0],
            -weight_ * (rest_(0, i) + rest_(1, i) + rest_(2, i))));
        triplets.push_back(Triplet(idO + i, idI_[1], weight_ * rest_(0, i)));
        triplets.push_back(Triplet(idO + i, idI_[2], weight_ * rest_(1, i)));
        triplets.push_back(Triplet(idO + i, idI_[3], weight_ * rest_(2, i)));
    }
    idO += n;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE AreaConstraint::AreaConstraint(const std::vector<int>& idI,
    Scalar weight, const Matrix3X& positions, Scalar rangeMin, Scalar rangeMax)
    : Constraint(idI, weight)
    , rangeMin_(rangeMin)
    , rangeMax_(rangeMax)
{
    assert(idI.size() == 3);
    ConstraintType_ = "Area";
    Matrix32 edges, P;
    edges.col(0) = positions.col(idI_[1]) - positions.col(idI_[0]);
    edges.col(1) = positions.col(idI_[2]) - positions.col(idI_[0]);
    P.col(0) = edges.col(0).normalized();
    P.col(1)
        = (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
    rest_ = (P.transpose() * edges).inverse();
    Scalar A = (P.transpose() * edges).determinant() / 2.;
    weight_ *= std::sqrt(std::abs(A));
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void AreaConstraint::calculateArea(const Matrix3X& positions, Scalar& area)
{
    Matrix32 edges, P;
    edges.col(0) = positions.col(idI_[1]) - positions.col(idI_[0]);
    edges.col(1) = positions.col(idI_[2]) - positions.col(idI_[0]);
    P.col(0) = edges.col(0).normalized();
    P.col(1) = (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
    Scalar A = std::abs((P.transpose() * edges).determinant() / 2.);

    SHAPEOP_OMP_CRITICAL{ area += A; }
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void AreaConstraint::project(const Matrix3X& positions,
    Matrix3X& projections, Matrix3X& f_int_nonePD, const Matrix3X& oldPositions)
{
    Matrix32 edges, P;
    edges.col(0) = (positions.col(idI_[1]) - positions.col(idI_[0]));
    edges.col(1) = (positions.col(idI_[2]) - positions.col(idI_[0]));
    P.col(0) = edges.col(0).normalized();
    P.col(1)
        = (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
    Matrix22 F = P.transpose() * edges * rest_;
    Eigen::JacobiSVD<Matrix22> svd(
        F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Vector2 S = svd.singularValues();
    Vector2 d(0., 0.);
    for (int i = 0; i < SHAPEOP_INNER_ITERATIONS; ++i) {
        // To keep the area constant, you simply set l1*l2 = 1 (li: principal
        // stretches)
        Scalar v = S(0) * S(1);
        Scalar f = v - clamp(v, rangeMin_, rangeMax_);
        Vector2 g(S(1), S(0));
        d = -((f - g.dot(d)) / g.dot(g)) * g;
        S = svd.singularValues() + d;
    }
    F = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();
    projections.block<3, 2>(0, idO_) = (weight_ * P * F);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void AreaConstraint::addConstraint(
    std::vector<Triplet>& triplets, int& idO) const
{
    // Firstly, we build the PD system
    if (PDSys_Build) {
        idO_ = idO;
        PDSys_Build = false;
    }

    int n = 2;
    for (int i = 0; i < n; ++i) {
        triplets.push_back(
            Triplet(idO + i, idI_[0], -weight_ * (rest_(0, i) + rest_(1, i))));
        triplets.push_back(Triplet(idO + i, idI_[1], weight_ * rest_(0, i)));
        triplets.push_back(Triplet(idO + i, idI_[2], weight_ * rest_(1, i)));
    }
    idO += n;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE VolumeConstraint::VolumeConstraint(const std::vector<int>& idI,
    Scalar weight, const Matrix3X& positions, Scalar rangeMin, Scalar rangeMax)
    : Constraint(idI, weight)
    , rangeMin_(rangeMin)
    , rangeMax_(rangeMax)
{
    assert(idI_.size() == 4);
    ConstraintType_ = "Volume";
    Matrix33 edges;
    for (int i = 0; i < 3; ++i)
        edges.col(i) = positions.col(idI_[i + 1]) - positions.col(idI_[0]);
    rest_ = edges.inverse();
    Scalar V = (edges).determinant() / 6.;
    weight_ *= std::sqrt(std::abs(V));
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void VolumeConstraint::calculateVolume(
    const Matrix3X& positions, Scalar& volume)
{
    Matrix33 edges;
    for (int i = 0; i < 3; ++i)
        edges.col(i) = positions.col(idI_[i + 1]) - positions.col(idI_[0]);
    Scalar V = std::abs((edges).determinant() / 6.);

    SHAPEOP_OMP_CRITICAL{ volume += V; }
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void VolumeConstraint::project(const Matrix3X& positions,
    Matrix3X& projections, Matrix3X& f_int_nonePD, const Matrix3X& oldPositions)
{
    Matrix33 edges;
    for (int i = 0; i < 3; ++i)
        edges.col(i) = positions.col(idI_[i + 1]) - positions.col(idI_[0]);
    Matrix33 F = edges * rest_;
    Eigen::JacobiSVD<Matrix33> svd(
        F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Vector3 S = svd.singularValues();
    Vector3 d(0., 0., 0.);
    // For rangeMin_ = rangeMax_ = 1., d stays zero vector and the projection is
    // simply the deformation gradient!
    for (int i = 0; i < SHAPEOP_INNER_ITERATIONS; ++i) {
        // To keep the volume constant, you simply set l1*l2*l3 = 1 (li:
        // principal stretches)
        Scalar v = S(0) * S(1) * S(2);
        Scalar f = v - clamp(v, rangeMin_, rangeMax_);
        Vector3 g(S(1) * S(2), S(0) * S(2), S(0) * S(1));
        d = -((f - g.dot(d)) / g.dot(g)) * g;
        S = svd.singularValues() + d;
    }
    if (svd.matrixU().determinant() * svd.matrixV().determinant() < 0.)
        S(2) = -S(2);
    F = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();
    projections.block<3, 3>(0, idO_) = weight_ * F;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void VolumeConstraint::addConstraint(
    std::vector<Triplet>& triplets, int& idO) const
{
    // Firstly, we build the PD system
    if (PDSys_Build) {
        idO_ = idO;
        PDSys_Build = false;
    }

    int n = 3;
    for (int i = 0; i < n; ++i) {
        triplets.push_back(Triplet(idO + i, idI_[0],
            -weight_ * (rest_(0, i) + rest_(1, i) + rest_(2, i))));
        triplets.push_back(Triplet(idO + i, idI_[1], weight_ * rest_(0, i)));
        triplets.push_back(Triplet(idO + i, idI_[2], weight_ * rest_(1, i)));
        triplets.push_back(Triplet(idO + i, idI_[3], weight_ * rest_(2, i)));
    }
    idO += n;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE BendingConstraint::BendingConstraint(const std::vector<int>& idI,
    Scalar weight, const Matrix3X& positions, Scalar rangeMin, Scalar rangeMax)
    : Constraint(idI, weight)
    , rangeMin_(rangeMin)
    , rangeMax_(rangeMax)
{
    ConstraintType_ = "Bending";
    Matrix3X p(3, idI.size());
    for (int i = 0; i < static_cast<int>(idI_.size()); ++i)
        p.col(i) = positions.col(idI_[i]);
    Scalar l01 = (p.col(0) - p.col(1)).norm();
    Scalar l02 = (p.col(0) - p.col(2)).norm();
    Scalar l12 = (p.col(1) - p.col(2)).norm();
    Scalar r0 = 0.5 * (l01 + l02 + l12);
    Scalar A0 = std::sqrt(r0 * (r0 - l01) * (r0 - l02) * (r0 - l12));
    Scalar l03 = (p.col(0) - p.col(3)).norm();
    Scalar l13 = (p.col(1) - p.col(3)).norm();
    Scalar r1 = 0.5 * (l01 + l03 + l13);
    Scalar A1 = std::sqrt(r1 * (r1 - l01) * (r1 - l03) * (r1 - l13));
    weight_ *= std::sqrt(3.0 / (A0 + A1));
    Scalar cot02 = ((l01 * l01) - (l02 * l02) + (l12 * l12)) / (4.0 * A0);
    Scalar cot12 = ((l01 * l01) + (l02 * l02) - (l12 * l12)) / (4.0 * A0);
    Scalar cot03 = ((l01 * l01) - (l03 * l03) + (l13 * l13)) / (4.0 * A1);
    Scalar cot13 = ((l01 * l01) + (l03 * l03) - (l13 * l13)) / (4.0 * A1);
    w_ = Vector4::Zero();
    w_(0) = cot02 + cot03;
    w_(1) = cot12 + cot13;
    w_(2) = -(cot02 + cot12);
    w_(3) = -(cot03 + cot13);
    n_ = (p * w_).norm();
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void BendingConstraint::project(const Matrix3X& positions,
    Matrix3X& projections, Matrix3X& f_int_nonePD, const Matrix3X& oldPositions)
{
    Vector3 e = Vector3::Zero();
    if (n_ > 1e-6) {
        for (int i = 0; i < static_cast<int>(idI_.size()); ++i)
            e += w_(i) * positions.col(idI_[i]);
        Scalar l = e.norm();
        if (l > 1e-6) {
            e /= l;
            l = n_ * clamp(l / n_, rangeMin_, rangeMax_);
            e *= l;
        }
    }
    projections.col(idO_) = weight_ * e;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void BendingConstraint::addConstraint(
    std::vector<Triplet>& triplets, int& idO) const
{
    // Firstly, we build the PD system
    if (PDSys_Build) {
        idO_ = idO;
        PDSys_Build = false;
    }

    for (int i = 0; i < static_cast<int>(idI_.size()); ++i)
        triplets.push_back(Triplet(idO, idI_[i], weight_ * w_(i)));
    idO += 1;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE ClosenessConstraint::ClosenessConstraint(
    const std::vector<int>& idI, Scalar weight, const Matrix3X& positions)
    : Constraint(idI, weight)
{
    assert(idI.size() == 1);
    ConstraintType_ = "Closeness";
    rest_ = positions.col(idI_[0]);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void ClosenessConstraint::project(const Matrix3X& /*positions*/,
    Matrix3X& projections, Matrix3X& f_int_nonePD, const Matrix3X& oldPositions)
{
    projections.col(idO_) = rest_ * weight_;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void ClosenessConstraint::addConstraint(
    std::vector<Triplet>& triplets, int& idO) const
{
    // Firstly, we build the PD system
    if (PDSys_Build) {
        idO_ = idO;
        PDSys_Build = false;
    }
    
    triplets.push_back(Triplet(idO, idI_[0], weight_));
    idO += 1;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void ClosenessConstraint::setPosition(const Vector3& position)
{
    rest_ = position;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE Vector3 ClosenessConstraint::getPosition() const
{
    return rest_;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
} // namespace npfem
} // namespace plb
///////////////////////////////////////////////////////////////////////////////
#endif // CONSTRAINT_CPP
///////////////////////////////////////////////////////////////////////////////
