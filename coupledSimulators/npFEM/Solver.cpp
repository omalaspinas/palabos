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
#ifndef SOLVER_CPP
#define SOLVER_CPP
///////////////////////////////////////////////////////////////////////////////
#include <iostream>

#include "Solver.h"
#include "LSSolver.h"
#include "Constraint.h"
#include "Force.h"
#include "common.h"
///////////////////////////////////////////////////////////////////////////////
#define THREED
///////////////////////////////////////////////////////////////////////////////
#ifdef NPFEM_SA
#ifdef SHAPEOP_MSVC
// Console Output
#include <windows.h>
#include <iostream>
#endif
#else
#include "core/util.h"
#include "core/plbLogFiles.h"
#endif
///////////////////////////////////////////////////////////////////////////////
#ifdef SHAPEOP_OPENMP
#ifdef SHAPEOP_MSVC
#define SHAPEOP_OMP_PARALLEL __pragma(omp parallel)
#define SHAPEOP_OMP_FOR __pragma(omp for)
#define OMP_CRITICAL __pragma(omp critical)
#else
#define SHAPEOP_OMP_PARALLEL _Pragma("omp parallel")
#define SHAPEOP_OMP_FOR _Pragma("omp for")
#define OMP_CRITICAL _Pragma("omp critical")
#endif
#else
#define SHAPEOP_OMP_PARALLEL
#define SHAPEOP_OMP_FOR
#define OMP_CRITICAL
#endif
///////////////////////////////////////////////////////////////////////////////
#ifndef BENCHMARK
#define WOLFE
#endif // else backtracking
///////////////////////////////////////////////////////////////////////////////
#define MASS_LUMPING_VOLUMETRIC
///////////////////////////////////////////////////////////////////////////////
namespace plb {
namespace npfem {
///////////////////////////////////////////////////////////////////////////////
Solver::Solver()
{
    rotations_ = Matrix33::Identity(3, 3);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE int Solver::addConstraint(const std::shared_ptr<Constraint>& c)
{
    constraints_.push_back(c);
    return static_cast<int>(constraints_.size() - 1);
}
///////////////////////////////////////////////////////////////////////////////
std::vector<std::shared_ptr<Constraint>>& Solver::getConstraints()
{
    return constraints_;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE std::shared_ptr<Constraint>& Solver::getConstraint(int id)
{
    return constraints_[id];
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE int Solver::addForces(const std::shared_ptr<Force>& f)
{
    forces_.push_back(f);
    return static_cast<int>(forces_.size() - 1);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE std::shared_ptr<Force>& Solver::getForce(int id)
{
    return forces_[id];
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver::setPoints(const Matrix3X& p) { points_ = p; }
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver::setVelocities(const Matrix3X& v) { velocities_ = v; }
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE const Matrix3X& Solver::getPoints() const { return points_; }
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver::setTimeStep(Scalar timestep) { delta_ = timestep; }
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE Scalar Solver::getTimeStep() { return delta_; }
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE Matrix3X& Solver::getVelocities()
{
    return velocities_;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver::setConnectivityList(
    const std::vector<std::vector<int>>& connectivity_list)
{
    connectivity_list_ = connectivity_list;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE const std::vector<std::vector<int>>&
    Solver::getConnectivityList() const
{
    return connectivity_list_;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver::set_Palabos_Forces(const Matrix3X& force_matrix)
{
    Palabos_Forces_ = force_matrix;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE const Matrix3X& Solver::get_Palabos_Forces() const
{
    return Palabos_Forces_;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver::set_onSurfaceParticle(
    const std::vector<bool>& onSurfaceParticle)
{
    onSurfaceParticle_ = onSurfaceParticle;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE const std::vector<bool>& Solver::get_onSurfaceParticle() const
{
    return onSurfaceParticle_;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver::shiftPoints(Vector3 vec){
    SHAPEOP_OMP_PARALLEL
    {
        SHAPEOP_OMP_FOR for (int i = 0; i < static_cast<int>(points_.cols());
        ++i)
    {
        points_.col(i) += vec;
    }
    }
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver::rotatePoints(const std::string& axis, const Scalar& theta){
    Matrix33 R;

    if (axis.compare("x") == 0) {
        R.col(0) = Vector3(1., 0., 0.);
        R.col(1) = Vector3(0., std::cos(theta), std::sin(theta));
        R.col(2) = Vector3(0., -std::sin(theta), std::cos(theta));
    }
    else if (axis.compare("y") == 0) {
        R.col(0) = Vector3(std::cos(theta), 0., -std::sin(theta));
        R.col(1) = Vector3(0., 1., 0.);
        R.col(2) = Vector3(std::sin(theta), 0., std::cos(theta));
    }
    else {
        R.col(0) = Vector3(std::cos(theta), std::sin(theta), 0.);
        R.col(1) = Vector3(-std::sin(theta), std::cos(theta), 0.);
        R.col(2) = Vector3(0., 0., 1.);
    }

    SHAPEOP_OMP_PARALLEL
    {
        SHAPEOP_OMP_FOR for (int i = 0; i < static_cast<int>(points_.cols());
        ++i)
    {
        points_.col(i) = R * points_.col(i);
    }
    }

    rotations_ = R*rotations_;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver::calculateArea(Scalar& Area){

    Area = 0.;

    size_t SurfaceMaterialConstraint = 0;
    size_t AreaConstraint = 0;
    for (int i = 0; i < static_cast<int>(constraints_.size()); ++i)
    {
        if (constraints_[i]->get_ConstraintType().compare("SurfaceMaterial") == 0)
            ++SurfaceMaterialConstraint;
        if (constraints_[i]->get_ConstraintType().compare("Area") == 0)
            ++AreaConstraint;
    }

    bool useAreaConstraint = false, useSurfaceMaterialConstraint = false;
    if (SurfaceMaterialConstraint != 0)
        useSurfaceMaterialConstraint = true;
    if (AreaConstraint != 0)
        useAreaConstraint = true;
    if (SurfaceMaterialConstraint != 0 && AreaConstraint != 0)
    {
        useSurfaceMaterialConstraint = true;
        useAreaConstraint = false;
    }
    if (SurfaceMaterialConstraint == 0 && AreaConstraint == 0)
    {
        std::cout << "CANNOT COMPUTE AREA ! \n";
    }


    for (int i = 0; i < static_cast<int>(constraints_.size()); ++i)
    {
        if (useSurfaceMaterialConstraint)
        {
            if (constraints_[i]->get_ConstraintType().compare("SurfaceMaterial") == 0) {
                auto c = std::dynamic_pointer_cast<plb::npfem::SurfaceMaterialConstraint>(getConstraint(i));
                c->calculateArea(points_, Area);
            }
        }
        if (useAreaConstraint)
        {
            if (constraints_[i]->get_ConstraintType().compare("Area") == 0) {
                auto c = std::dynamic_pointer_cast<plb::npfem::AreaConstraint>(getConstraint(i));
                c->calculateArea(points_, Area);
            }
        }
    }
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver::calculateVolume(Scalar& Volume)
{
    Volume = 0.;

    SHAPEOP_OMP_PARALLEL
    {
        SHAPEOP_OMP_FOR for (int i = 0; i < static_cast<int>(constraints_.size()); ++i)
    {
        if (constraints_[i]->get_ConstraintType().compare("Volume") == 0) {
            auto c = std::dynamic_pointer_cast<plb::npfem::VolumeConstraint>(getConstraint(i));
            c->calculateVolume(points_, Volume);
        }
    }
    }
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver::calculateVolumeFromSurfaceAndPrepareConstraint(Scalar& Volume)
{
    Volume = 0.;
    gradC_ = Matrix3X::Zero(3, points_.cols());
    deltaX_ = points_;
    deltaX_.setZero();

    SHAPEOP_OMP_PARALLEL
    {
        SHAPEOP_OMP_FOR for (int i = 0; i < static_cast<int>(constraints_.size()); ++i)
    {
        if (constraints_[i]->get_ConstraintType().compare("SurfaceMaterial") == 0) {
            auto c = std::dynamic_pointer_cast<plb::npfem::SurfaceMaterialConstraint>(getConstraint(i));
            std::vector<int> idI_ = c->get_idI();

            Matrix32 edges, P;
            edges.col(0) = points_.col(idI_[1]) - points_.col(idI_[0]);
            edges.col(1) = points_.col(idI_[2]) - points_.col(idI_[0]);
            P.col(0) = edges.col(0).normalized();
            P.col(1) = (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();

            Scalar Area = std::abs((P.transpose() * edges).determinant() / 2.);
            Vector3 normal = (edges.col(0)).cross(edges.col(1));
            normal.normalize();

            Scalar local_Vol = (1. / 9.) * Area * (points_.col(idI_[0]) + points_.col(idI_[1]) + points_.col(idI_[2])).dot(normal);

            OMP_CRITICAL
            {
                Volume += local_Vol;
            gradC_.col(idI_[0]) += (1. / 3.) * Area * normal;
            gradC_.col(idI_[1]) += (1. / 3.) * Area * normal;
            gradC_.col(idI_[2]) += (1. / 3.) * Area * normal;
            }
        }
    }
    }

        if (applyGlobalVolumeConservation_) {
            Scalar C = Volume - Volume0_;
            deltaX_ = -(C / gradC_.squaredNorm()) * gradC_; // There is a mass weighting, but I consider uniform distribution
        }
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver::mass_lumping()
{
    std::vector<Triplet> triplets;
    for (int i = 0; i < static_cast<int>(constraints_.size()); ++i) {
#ifdef MASS_LUMPING_VOLUMETRIC
        if (constraints_[i]->get_ConstraintType().compare("VolumeMaterial") == 0) {
            auto c = std::dynamic_pointer_cast<plb::npfem::VolumeMaterialConstraint>(getConstraint(i));
            c->mass_lumping(points_, triplets);
        }
#else
        if (constraints_[i]->get_ConstraintType().compare("SurfaceMaterial") == 0) {
            auto c = std::dynamic_pointer_cast<plb::npfem::SurfaceMaterialConstraint>(getConstraint(i));
            c->mass_lumping(points_, triplets);
        }
#endif
    }
    M_.setFromTriplets(triplets.begin(), triplets.end());
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE bool Solver::initialize(Scalar Calpha
    , Scalar Cbeta
    , Scalar timestep
    , Scalar rho
    , bool doModalAnalysis
    , bool applyGlobalVolumeConservation
    , Scalar globalVolumeConservationWeight
    , Scalar externalSolverSpatialResolution
    , Scalar externalSolverTimeResolution)
{
#ifdef NPFEM_SA
#ifdef SHAPEOP_MSVC
    // Console output for dll
    AllocConsole();
    freopen("CONOUT$", "w", stdout);
#endif
#endif
#ifdef DEBUG
    printf(" rho %f  Cbeta %f, Calpha %f, timestep %f, cpu  \n", rho, Cbeta, Calpha, timestep);
#endif

    int n_points = static_cast<int>(points_.cols());
    int n_constraints = static_cast<int>(constraints_.size());

    assert(n_points != 0);
    assert(n_constraints != 0);

    externalSolverSpatialResolution_ = externalSolverSpatialResolution;
    externalSolverTimeResolution_ = externalSolverTimeResolution;

    applyGlobalVolumeConservation_ = applyGlobalVolumeConservation;
    globalVolumeConservationWeight_ = globalVolumeConservationWeight;

    calculateArea(Area0_);
#ifdef THREED
    if (applyGlobalVolumeConservation_)
        calculateVolumeFromSurfaceAndPrepareConstraint(Volume0_);
    else
        calculateVolume(Volume0_);
#endif //THREED

#ifdef NPFEM_SA
    std::cout << "***********************************************************"
        << std::endl;
    std::cout << "CPU VERSION >> Number of Points & Constraints: " << n_points
        << " - " << n_constraints << std::endl;
    std::cout << "Initial Area: " << Area0_ << std::endl;
#ifdef THREED
    std::cout << "Initial Volume: " << Volume0_ << std::endl;
#endif
    std::cout << "***********************************************************"
        << std::endl;
#endif

    // Setup the Projective Dynamics system
    // Build the system like having only PD energies!
    std::vector<Triplet> triplets;
    int idO = 0;
    for (int i = 0; i < n_constraints; ++i) {
        // NONE_PD Energies (may add more : TODO)
        if (constraints_[i]->get_ConstraintType().compare("SurfaceMaterial") == 0)
            continue;
        // PD Energies
        constraints_[i]->addConstraint(triplets, idO);
    }
    projections_.setZero(3, idO);
    SparseMatrix A = SparseMatrix(idO, n_points);
    SparseMatrix At = SparseMatrix(n_points, idO);
    A.setFromTriplets(triplets.begin(), triplets.end());
    At = A.transpose();

    // Quasi-Newton (see Liu_2017 for more details)
    // These matrices correspond to a system that is PD-only
    // Later, we add the nonePD energies & forces, like superposing the 2 systems
    L_ = SparseMatrix(n_points, n_points);
    J_ = SparseMatrix(n_points, idO);
    L_ = At * A;
    J_ = At;
    f_int_nonePD_ = Matrix3X::Zero(3, n_points);

    // Build LAPLACIAN considering the contribution of nonePD energies as well
    triplets.clear();
    idO = 0;
    for (int i = 0; i < n_constraints; ++i)
        constraints_[i]->addConstraint(triplets, idO);
    SparseMatrix B = SparseMatrix(idO, n_points);
    SparseMatrix Bt = SparseMatrix(n_points, idO);
    B.setFromTriplets(triplets.begin(), triplets.end());
    Bt = B.transpose();
    // See Liu_2017 & Thesis for more
    SparseMatrix Laplacian = Bt * B;

    Calpha_ = Calpha; // PBD Damping
    Cbeta_ = Cbeta; // Rayleigh Damping (no mass contribution)
    delta_ = timestep;
    rho_ = rho;
    oldPoints_ = Matrix3X(3, n_points);
    oldPointsHistory_ = 0;
    velocities_ = Matrix3X::Zero(3, n_points);
    momentum_ = Matrix3X::Zero(3, n_points);

    // Mass and Damping Matrices
    M_ = SparseMatrix(n_points, n_points);
    M_tilde_ = SparseMatrix(n_points, n_points);
    M_tilde_inv_ = MatrixXX(n_points, n_points);
    M_star_ = SparseMatrix(n_points, n_points);

    //mass_lumping(); // Go with uniform ditribution
    M_.setIdentity();
#ifdef THREED
    M_ *= (Volume0_ * rho_) / n_points; // mass = Volume * density
#else
    M_ *= 1.;
#endif

    if (doModalAnalysis) {
        // Modal Analysis of the LINEAR System. Kx = lMx <=> Nx = lMx
        Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXX> es(
            (MatrixXX(Laplacian)), (MatrixXX(M_)));

        eigenValues_ = es.eigenvalues();
        eigenVectors_ = es.eigenvectors();

        // ω: the circular frequency
        // λ = ω^2 => ω = sqrt(λ)
        // or in Hertz:
        // natural (circular) frequency ω = 2πf => Natural Period T = 1/f

        naturalPeriods_ = eigenValues_;
        dampingRatios_ = eigenValues_;
        for (int mode = 0; mode < eigenValues_.size(); ++mode) {
            Scalar omega = sqrt(eigenValues_(mode));
            naturalPeriods_(mode) = 2. * acos(-1) / omega;
            dampingRatios_(mode) = 0.5 * Cbeta_ * omega;
#ifndef NPFEM_SA
            global::logfile_nonparallel("NaturalPeriods_and_DampingRatio.log")
                .flushEntry("Mode: " + util::val2str(mode) + ": "
                    + util::val2str(naturalPeriods_(mode)) + " | "
                    + util::val2str(dampingRatios_(mode)));
#endif
        }
    }

    // Damping based on Rayleigh with the "stiffness" part (Laplacian in my case)
    M_tilde_ = M_ + delta_ * (Cbeta_ * Laplacian);
    M_star_ = M_tilde_ / (delta_ * delta_);

    // Inverse of M_tilde
    M_tilde_inv_ = (MatrixXX(M_tilde_)).inverse();

    // Palabos interface
    Palabos_Forces_ = Matrix3X::Zero(3, n_points);
    //bodyID_ = 0;
    Palabos_iT_ = 0;

    // Prefactorize matrix
    solver_ = std::make_shared<plb::npfem::SimplicialLDLTSolver>();
    solver_->initialize(M_star_ + Laplacian); // Hessian Approximation

#ifdef NUM_STABILITY
    MatrixXX N_dense = MatrixXX(M_star_ + Laplacian);
    N_dense_inv_ = N_dense.inverse();
#endif

#ifdef NPFEM_SA
    // If it is called from Palabos, this is set from before
    // For testing purposes ... collisions
    onSurfaceParticle_.clear();
    for (size_t i = 0; i < points_.cols(); ++i)
        onSurfaceParticle_.push_back(true);
#endif

    G_flag_ = 0;
    return solver_->info() == Eigen::Success;
}
///////////////////////////////////////////////////////////////////////////////
std::vector<int> *Solver::get_mesh_graph() {

    if (graph != NULL) {
        //std::cout << "NOT NULL" << std::endl;
        return graph;
    }
    graph = new std::vector<int>[points_.cols()];

    Eigen::SparseMatrix<bool>lookUp = Eigen::SparseMatrix<bool>(points_.cols(), points_.cols());

    for (auto &tri : connectivity_list_) {
        for (int i = 0; i<3; i++) {
            //lookUp.coeffRef(tri[i], tri[(i + 1)%3]) = true;
            lookUp.coeffRef(tri[i], tri[(i + 1) % 3]) = true;
            lookUp.coeffRef(tri[i], tri[(i + 2) % 3]) = true;
        }
    }

    for (int k = 0; k < lookUp.outerSize(); ++k) {
        graph[k].push_back(k);
        for (Eigen::SparseMatrix<bool>::InnerIterator it(lookUp, k); it; ++it) {
            graph[k].push_back(it.row());
        }
    }

    return graph;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE bool Solver::solve(int m
    , unsigned int max_iterations
    , unsigned int max_line_search_loops
    , unsigned int max_attempts_to_solve_stagnation
    , unsigned int convergence_window
    , Scalar tol
    , Scalar gamma
    , Scalar gamma2
    , Scalar collisions_threshold_rep
    , Scalar collisions_weight_rep
    , Scalar collisions_threshold_nonRep
    , Scalar collisions_weight_nonRep
    , Scalar beta_morse
    , Scalar max_u_lb
    , Scalar max_u_lb_fin)
{
#ifdef NPFEM_SA
#ifdef SHAPEOP_MSVC
    // Console output for dll
    AllocConsole();
    freopen("CONOUT$", "w", stdout);
#endif
#endif

    max_u_lb_ = max_u_lb;
    max_u_lb_fin_ = max_u_lb_fin;

    // m = 2, simple (naive) Quasi-Newton & m > 2, L-BFGS
    if (m < 2) m = 15;

    // For convergence issues, move the points at the origin.
    pointsTI_ = points_;
    centroid_ = points_.rowwise().mean();
    points_.colwise() -= centroid_;
    // at the end of the loop, move them back

    // Flatten Colliding Pieces into one suitable data-structure
    size_t numCollidingCandidatesPoints = 0;
    for (size_t i = 0; i < verticesOfCollidingPieces_.size(); ++i)
        numCollidingCandidatesPoints += verticesOfCollidingPieces_[i].cols();

    collidingCandidatesPoints_ = Matrix3X::Zero(3, numCollidingCandidatesPoints);
    normalsCollidingCandidatesPoints_ = Matrix3X::Zero(3, numCollidingCandidatesPoints);

    int shift = 0;

    for (size_t i = 0; i < verticesOfCollidingPieces_.size(); ++i) {
        // Block of size (p,q), starting at (i,j) -> matrix.block(i,j,p,q);
        collidingCandidatesPoints_.block(0, shift, 3, verticesOfCollidingPieces_[i].cols()) = verticesOfCollidingPieces_[i];
        normalsCollidingCandidatesPoints_.block(0, shift, 3, normalsOfCollidingPieces_[i].cols()) = normalsOfCollidingPieces_[i];
        shift += verticesOfCollidingPieces_[i].cols();
    }

    if (numCollidingCandidatesPoints > 0) {
        collidingCandidatesPoints_.colwise() -= centroid_;
        collisionDetection(collisions_threshold_rep, collisions_threshold_nonRep);
        addCollisionConstraints(collisions_threshold_rep, collisions_weight_rep, collisions_threshold_nonRep, collisions_weight_nonRep, beta_morse);
    }

#ifdef NPFEM_SA
    // Convert forces into a dense matrix form
    Palabos_Forces_.setZero(3, points_.cols());
    for (int i = 0; i < static_cast<int>(Palabos_Forces_.cols()); ++i) {
        for (int j = 0; j < static_cast<int>(forces_.size()); ++j) {
            Palabos_Forces_.col(i) += forces_[j]->get(points_, i);
        }
    }
#endif
    
    // Momentum = The new positions in the absence of internal forces
    momentum_ = points_
        + (delta_ * M_tilde_inv_ * M_ * velocities_.transpose()
            + delta_ * delta_ * M_tilde_inv_ * Palabos_Forces_.transpose())
        .transpose();

    // Kharevych variational damping (it is not actively used)
    oldPoints_ = points_;
    if (oldPointsHistory_ == 0) {
        oldPointsEWMA_ = oldPoints_;
        ++oldPointsHistory_;
    }
    else {
        Scalar C_EWMA = 0.8;
        oldPointsEWMA_ = C_EWMA * oldPoints_ + (1. - C_EWMA) * oldPointsEWMA_;
    }

    // First guess
    points_ = momentum_;
    computeInternalForces();

    // counters
    unsigned int it = 0;
    unsigned int avg_line_search_loops = 0;

    // convergence related
    bool stagnation = false;
    unsigned int attempts_to_solve_stagnation = 0;
    Scalar convergence_criterion;
    std::vector<Scalar> objective_history;

    // Since we have the first guess
    objective_history.push_back(evalObjective(points_, momentum_));

    // LBFGS
    std::vector<Matrix3X> iterates, gradients;
    Matrix3X descent_direction = Matrix3X::Zero(3, points_.cols()), q = Matrix3X::Zero(3, points_.cols());

    // Main Optimization Loop to advance body at t+1
    while (true)
    {
        // Store History
        iterates.push_back(points_);
        if (static_cast<int>(iterates.size()) >= m)
            iterates.erase(iterates.begin());

        gradients.push_back((M_star_ * (points_ - momentum_).transpose()
            + L_ * points_.transpose() - J_ * projections_.transpose()
            - f_int_nonePD_.transpose()).transpose());
        if (static_cast<int>(gradients.size()) >= m)
            gradients.erase(gradients.begin());

        if (gradients.back().norm() <= tol)
            break;

        if (static_cast<unsigned int>(objective_history.size()) > convergence_window)
            objective_history.erase(objective_history.begin());

        if (static_cast<unsigned int>(objective_history.size()) >= convergence_window)
        {
            Scalar mean = 0.;
            for (int i = 0; i < static_cast<int>(objective_history.size()); ++i)
                mean += objective_history[i];
            mean /= (Scalar)objective_history.size();

            Scalar stdev = 0.;
            for (int i = 0; i < static_cast<int>(objective_history.size()); ++i)
                stdev += std::pow(objective_history[i] - mean, 2.);
            stdev /= (Scalar)objective_history.size() - 1.;
            stdev = std::pow(stdev, 0.5);

            convergence_criterion = std::abs(stdev / mean);
            if (convergence_criterion <= tol)
                break;
        }

        // LBFGS
        q = -gradients.back();

        std::vector<Matrix3X> s, y;
        for (std::vector<Matrix3X>::iterator iter = iterates.begin();
            iter != iterates.end() - 1; ++iter)
            s.push_back(*(iter + 1) - *iter);

        for (std::vector<Matrix3X>::iterator iter = gradients.begin();
            iter != gradients.end() - 1; ++iter)
            y.push_back(*(iter + 1) - *iter);

        // See Nocedal Chapter 7 for more details
        std::vector<Scalar> rho, alpha;
        for (int i = static_cast<int>(s.size()) - 1; i >= 0; --i) {
            rho.insert(rho.begin(), 1. / (y[i] * s[i].transpose()).trace());
            alpha.insert(
                alpha.begin(), (*rho.begin() * s[i] * q.transpose()).trace());
            q -= *alpha.begin() * y[i];
        }

        SHAPEOP_OMP_PARALLEL
        {
#ifdef NUM_STABILITY
            SHAPEOP_OMP_FOR for (int j = 0; j < 3; ++j) descent_direction.row(j) = (N_dense_inv_ * gradients.back().row(j).transpose()).transpose();
#else
            SHAPEOP_OMP_FOR for (int j = 0; j < 3; ++j) descent_direction.row(j)
            = solver_->solve(q.row(j).transpose(), gradients.back().row(j).transpose()).transpose();
#endif
        }

        Scalar beta;
        for (int i = 0; i < static_cast<int>(s.size()); ++i) {
            beta = rho[i] * (y[i] * descent_direction.transpose()).trace();
            descent_direction += s[i] * (alpha[i] - beta);
        }

        // Line search loop with Wolfe conditions (gamma & gamma2)
        Scalar directional_derivative_prev
            = (gradients.back() * descent_direction.transpose()).trace();
        Scalar alpha_ = 2.;
        bool Armijo_condition = false;
#ifdef WOLFE
        bool curvature_condition = false;
#endif
        unsigned int line_search_loops = 0;
        while (true) {
            alpha_ /= 2.;
            points_ = iterates.back() + alpha_ * descent_direction;

            // Re-evaluate the internal forces (collisions included)
            computeInternalForces();

            // plus one is for stability purposes
            Scalar evalObjective_tmp = evalObjective(points_, momentum_);
            Armijo_condition
                = evalObjective_tmp
                <= (objective_history.back()
                    + gamma * alpha_ * directional_derivative_prev);

#ifdef WOLFE
            // curvature condition
            Matrix3X grad = (M_star_ * (points_ - momentum_).transpose()
                + L_ * points_.transpose()
                - J_ * projections_.transpose() - f_int_nonePD_.transpose())
                .transpose();
            curvature_condition = (grad * descent_direction.transpose()).trace()
                >= gamma2 * directional_derivative_prev;
#endif

            line_search_loops += 1;
            avg_line_search_loops += 1;

#ifdef WOLFE
            if (Armijo_condition && curvature_condition) {
#else
            if (Armijo_condition) {
#endif
                objective_history.push_back(evalObjective_tmp);
                break;
            }

            if (line_search_loops == max_line_search_loops) {
                objective_history.push_back(evalObjective_tmp);
                attempts_to_solve_stagnation++;
                break;
            }
            }

        ++it;

        if (attempts_to_solve_stagnation == max_attempts_to_solve_stagnation) {
            stagnation = true;
            break;
        }

        if (it == max_iterations)
            break;
        }

#ifdef NPFEM_SA
    std::cout << "Loops until convergence: " << it << " | ";
    if (it != 0)
        std::cout << "Average Line Search Loops: "
        << ((Scalar)(avg_line_search_loops)) / ((Scalar)(it))
        << " | attempts_to_solve_stagnation: "
        << attempts_to_solve_stagnation
        << std::endl;
#else
    /*
    if (it != 0)
    logfile_nonparallel(
    "body_" + util::val2str(bodyID_) + "_ShapeOp_Convergence.log")
    .flushEntry("Palabos iT: " + util::val2str(Palabos_iT_)
    + " | Loops until convergence: " + util::val2str(it)
    + " | Average Line Search Loops: "
    + util::val2str(
    ((Scalar)(avg_line_search_loops)) / ((Scalar)(it)))
    + " | attempts_to_solve_stagnation: "
    + util::val2str(attempts_to_solve_stagnation) );
    */
#endif

    // In any case you go through this.
    // If max_attempts or max_iTs reached then this means that the method cannot reduce
    // further the energy. However, the points_ up to this point are a vaild guess, if not the best.
    // Meaning that starting from the momentum (first), we went as far as what we could.
    PBD_Damping();

#ifdef NPFEM_SA
    calculateArea(Area_);
#ifdef THREED
    if (applyGlobalVolumeConservation_)
        calculateVolumeFromSurfaceAndPrepareConstraint(Volume_);
    else
        calculateVolume(Volume_);
#endif // THREED

    Vector3 mean_velocity(0., 0., 0.);
    for (int i = 0; i < static_cast<int>(points_.cols()); ++i)
        mean_velocity += velocities_.col(i);
    mean_velocity /= points_.cols();

    std::cout << "Mean Velocity (micro-m/micro-sec): " << mean_velocity.norm() << std::endl;

    std::cout << "Area Conservation (%): "
        << 100. * std::abs(Area_ - Area0_) / Area0_ << std::endl;

#ifdef THREED
    std::cout << "Volume Conservation (%): "
        << 100. * std::abs(Volume_ - Volume0_) / Volume0_ << std::endl;

    std::cout << "Volume: " << Volume_ << std::endl;
#endif // THREED
#else
    /*
    Vector3 Total_Force_on_Body(0., 0., 0.);
    Vector3 mean_velocity(0., 0., 0.);
    for (int i = 0; i < static_cast<int>(points_.cols()); ++i) {
        Total_Force_on_Body += Palabos_Forces_.col(i);
        mean_velocity += velocities_.col(i);
    }
    mean_velocity /= points_.cols();
    */

    calculateArea(Area_);
    if (applyGlobalVolumeConservation_)
        calculateVolumeFromSurfaceAndPrepareConstraint(Volume_);
    else
        calculateVolume(Volume_);

    /*
    logfile_nonparallel(
    "body_" + util::val2str(bodyID_) + "_ShapeOp_BodyForce.log")
    .flushEntry("Palabos iT: " + util::val2str(Palabos_iT_)
    + " | Total Force on Body (picoNewtons): "
    + util::val2str(Total_Force_on_Body.norm()));

    logfile_nonparallel(
    "body_" + util::val2str(bodyID_) + "_ShapeOp_MeanVelocity.log")
    .flushEntry("Palabos iT: " + util::val2str(Palabos_iT_)
    + " | Mean Velocity (micro-m/micro-sec): "
    + util::val2str(mean_velocity.norm()));

    logfile_nonparallel(
    "body_" + util::val2str(bodyID_) + "_ShapeOp_ConservedQuantities.log")
    .flushEntry("Palabos iT: " + util::val2str(Palabos_iT_)
    + " | Area Conservation (%): "
    + util::val2str(100. * std::abs(Area_ - Area0_) / Area0_)
    + " | Volume Conservation (%): "
    + util::val2str(100. * std::abs(Volume_ - Volume0_) / Volume0_));
    */

    //logfile_nonparallel(
    //    "body_" + util::val2str(bodyID_) + "_Energy.csv")
    //    .flushEntry(util::val2str(evalKinPlusPotEnergy(points_)));
#endif // NPFEM_SA

    //GyrationTensor();

    if (numCollidingCandidatesPoints > 0)
        removeCollisionConstraints();

    // TI: move them back
    points_.colwise() += centroid_;
    collidingCandidatesPoints_.colwise() += centroid_;

    return solver_->info() == Eigen::Success;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE Scalar Solver::evalKinPlusPotEnergy(const Matrix3X& x)
{
    Scalar none_PD_energy = 0.;
    for (int i = 0; i < static_cast<int>(constraints_.size()); ++i)
        none_PD_energy += constraints_[i]->get_E_nonePD();

    Scalar C = 0.5 * (projections_ * projections_.transpose()).trace();

    return 0.5 * (velocities_ * M_tilde_ * velocities_.transpose()).trace()
         + 0.5 * (x * L_ * x.transpose()).trace()
         - (x * J_ * projections_.transpose()).trace() + C
         + none_PD_energy;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE Scalar Solver::evalObjective(const Matrix3X& x, const Matrix3X& y)
{
    Scalar none_PD_energy = 0.;
    for (int i = 0; i < static_cast<int>(constraints_.size()); ++i)
        none_PD_energy += constraints_[i]->get_E_nonePD();

    if (applyGlobalVolumeConservation_)
        none_PD_energy += 0.5 * globalVolumeConservationWeight_ * deltaX_.squaredNorm();

    Scalar C = 0.5 * (projections_ * projections_.transpose()).trace();
    //Scalar C = projections_.squaredNorm() / 2;

#ifdef NUM_STABILITY
    return 0.5 * ((x - y)*M_star_).cwiseProduct((x - y)).sum()
        + 0.5*((x*L_).cwiseProduct(x)).sum()
        - ((x*J_).cwiseProduct(projections_)).sum() + C
        + none_PD_energy;
#else
    return 0.5 * ((x - y) * M_star_ * (x - y).transpose()).trace()
        + 0.5 * (x * L_ * x.transpose()).trace()
        - (x * J_ * projections_.transpose()).trace() + C
        + none_PD_energy;
#endif
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver::GyrationTensor()
{
    G_.setZero();
    for (int m = 0; m < 3; ++m) {
        for (int n = 0; n < 3; ++n) {
            for (int i = 0; i < points_.cols(); ++i) {
                G_(m, n) += (points_.col(i)[m] - Xcm_[m]) *
                    (points_.col(i)[n] - Xcm_[n]);
            }
        }
    }
    G_ /= points_.cols();

    Eigen::SelfAdjointEigenSolver<Matrix33> eigensolver(G_);
    if (eigensolver.info() != Eigen::Success) {
        ;
    }
    else {
        // Ascending order: First is the smallest!
        Vector3 eigenvalues = eigensolver.eigenvalues();
        if (G_flag_ == 0) {
            restShapeEigenvalues_ = eigenvalues;
            ++G_flag_;
        }
#ifdef NPFEM_SA
        std::cout << "The eigenvalues of the Gyration tensor are:\n" << eigensolver.eigenvalues().transpose() << "\n";
#else
        global::logfile_nonparallel(
            "body_" + util::val2str(bodyID_) + "_ShapeOp_GyrationEigen.log")
            .flushEntry(util::val2str(std::abs(eigenvalues[0] - restShapeEigenvalues_[0])));
#endif // NPFEM_SA
    }
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver::PBD_Damping()
{
    velocities_ = (points_ - oldPoints_) / delta_;

    Vector3 Xcm(0., 0., 0.), Vcm(0., 0., 0.);
    Matrix3X tmpX = (M_ * points_.transpose()).transpose();
    Matrix3X tmpV = (M_ * velocities_.transpose()).transpose();
    for (int i = 0; i < static_cast<int>(tmpX.cols()); ++i) {
        Xcm += tmpX.col(i);
        Vcm += tmpV.col(i);
    }
    Scalar TotalMass = (MatrixXX(M_)).trace();
    Xcm /= TotalMass;
    Vcm /= TotalMass;

    Xcm_ = Xcm;

    Matrix3X r = points_;
    r.colwise() -= Xcm;

    Vector3 L(0., 0., 0.);
    for (int i = 0; i < static_cast<int>(tmpX.cols()); ++i) {
        L += (r.col(i)).cross(tmpV.col(i));
    }

    Matrix33 Icm = Matrix33::Zero();
    for (int i = 0; i < static_cast<int>(tmpX.cols()); ++i) {
        // Cross product converted to matrix multiplication
        Matrix33 Rxi = Matrix33::Zero();
        Rxi(0, 1) = -r.col(i)[2];
        Rxi(0, 2) = r.col(i)[1];
        Rxi(1, 0) = r.col(i)[2];
        Rxi(1, 2) = -r.col(i)[0];
        Rxi(2, 0) = -r.col(i)[1];
        Rxi(2, 1) = r.col(i)[0];

        Icm += M_.coeff(i, i) * Rxi * Rxi.transpose();
    }
    Matrix33 Icm_inv = Icm.inverse();

    Vector3 Ocm = Icm_inv * L;

    for (int i = 0; i < static_cast<int>(velocities_.cols()); ++i)
    {
        // See PBD Muller
        Vector3 dvi = Vcm + Ocm.cross(r.col(i)) - velocities_.col(i);
        velocities_.col(i) += Calpha_ * dvi;
    }

#ifndef NPFEM_SA
    // Lattice Units
    Scalar C_vel = externalSolverSpatialResolution_ / externalSolverTimeResolution_;

    Matrix3X lattice_velocities = velocities_ / C_vel;

    size_t velocityLimit = 0;
    for (int i = 0; i < static_cast<int>(velocities_.cols()); ++i)
    {
        if (lattice_velocities.col(i).norm() >= max_u_lb_)
        {
            lattice_velocities.col(i).normalize();
            lattice_velocities.col(i) *= max_u_lb_fin_;

            velocities_.col(i) = lattice_velocities.col(i) * C_vel;

            ++velocityLimit;
        }
    }

#ifdef ENABLE_LOGS
    if (velocityLimit > 0)
    {
        global::logfile_nonparallel("instability.log")
            .flushEntry("At Palabos iT=" + util::val2str(Palabos_iT_) + " | ShapeOp bodyID=" + util::val2str(bodyID_) +
                " | VELOCITY LIMIT ACTIVATED, on " + util::val2str(velocityLimit) + " vertices");
    }
#endif // ENABLE_LOGS

#endif // !NPFEM_SA

    points_ = oldPoints_ + velocities_ * delta_;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver::computeInternalForces()
{
    f_int_nonePD_.setZero(3, points_.cols());
    SHAPEOP_OMP_PARALLEL
    {
#ifdef NUM_STABILITY
        for (int i = static_cast<int>(constraints_.size()) - 1; i >= 0; --i) {
#else
        SHAPEOP_OMP_FOR for (int i = 0; i < static_cast<int>(constraints_.size()); ++i) {
#endif
        constraints_[i]->project(points_, projections_, f_int_nonePD_, oldPointsEWMA_);
        }
    }

#ifdef THREED
    // Spring-Like force
    Scalar Volume;
    if (applyGlobalVolumeConservation_)
    {
        calculateVolumeFromSurfaceAndPrepareConstraint(Volume);
        f_int_nonePD_ += globalVolumeConservationWeight_ * deltaX_;
    }
#endif //THREED
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver::collisionDetection(Scalar collisions_threshold_rep, Scalar collisions_threshold_nonRep)
{
    collidingPointsInfo_.clear();

    // Construct kd-tree - Attention to the TRANSPOSE
    // This tmp_struct is needed for compatibility with different Eigen versions
    // Eigen::Matrix<num_t, Dynamic, Dynamic> mat(nSamples, dim); FROM NANOFLANN
    unsigned int dim = 3;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> tmp_struct(collidingCandidatesPoints_.cols(), dim);
    tmp_struct = collidingCandidatesPoints_.transpose();
    kd_tree extendedPoints_index(tmp_struct, 10 /* default value (optimal 10.-.50) */);
    extendedPoints_index.index->buildIndex();

    // Important note: If L2 norms are used, notice that search radius and all 
    // passed and returned distances are actually squared distances.

    // Check for every point on my Defo Body the NNs and investigate collisions
    for (int i = 0; i < points_.cols(); ++i)
    {
        if (!onSurfaceParticle_[i])
            continue;

        Vector3 bodyPoint = points_.col(i);

        std::vector<Scalar> query_pt(dim);
        query_pt[0] = bodyPoint[0];
        query_pt[1] = bodyPoint[1];
        query_pt[2] = bodyPoint[2];

        // ----------
        // kNN search
        // ----------
        const size_t num_results = 1;
        std::vector<size_t> ret_indexes(num_results);
        std::vector<Scalar> out_dists_sqr(num_results);

        nanoflann::KNNResultSet<Scalar> resultSet(num_results);

        resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
        extendedPoints_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

        Vector3 collidingPoint = collidingCandidatesPoints_.col(ret_indexes[0]);
        Vector3 collidingPointNormal = normalsCollidingCandidatesPoints_.col(ret_indexes[0]);

        Scalar collisions_threshold;
        if (collidingPointNormal.norm() >= 2.0)
            collisions_threshold = collisions_threshold_rep;
        else
            collisions_threshold = collisions_threshold_nonRep;

        if ((collidingPoint - bodyPoint).norm() <= collisions_threshold)
        {
            std::vector<Vector3> container(2);
            
            if (collidingPointNormal.norm() >= 2.0)
                // 0.5 is to decode the info of repulsion
                container[0] = collidingPoint + collisions_threshold * 0.5 * collidingPointNormal; // offset here
            else
                container[0] = collidingPoint + collisions_threshold * collidingPointNormal; // offset here
            
            container[1] = collidingPointNormal;
            collidingPointsInfo_[i] = container;
        }
    }
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver::addCollisionConstraints(
    Scalar collisions_threshold_rep, Scalar collisionsWeight_rep,
    Scalar collisions_threshold_nonRep, Scalar collisionsWeight_nonRep,
    Scalar beta_morse)
{
    std::unordered_map<int, std::vector<Vector3>>::const_iterator it;
    for (it = collidingPointsInfo_.cbegin(); it != collidingPointsInfo_.cend(); ++it)
    {
        std::vector<int> idI;
        idI.push_back(it->first);
        Vector3 collidingPoint = it->second[0];
        Vector3 collidingPointNormal = it->second[1];

        // 1. is a dummy collision weight. Set them below.
        auto c = std::make_shared<CollisionConstraint>(idI, 1., points_);
        
        c->setCollindingPoint(collidingPoint);
        c->setCollindingPointNormal(collidingPointNormal);
        // set the augmented collision energy
        c->setParams(collisions_threshold_rep, collisionsWeight_rep, collisions_threshold_nonRep, collisionsWeight_nonRep, beta_morse);

        addConstraint(c);
    }
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver::removeCollisionConstraints()
{
    std::unordered_map<int, std::vector<Vector3>>::const_iterator it;
    for (it = collidingPointsInfo_.cbegin(); it != collidingPointsInfo_.cend(); ++it)
        constraints_.pop_back();
}
///////////////////////////////////////////////////////////////////////////////
} // namespace npfem
} // namespace plb
///////////////////////////////////////////////////////////////////////////////
#endif // SOLVER_CPP
