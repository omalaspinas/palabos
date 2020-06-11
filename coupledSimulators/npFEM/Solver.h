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
#ifndef SOLVER_H
#define SOLVER_H
///////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <memory>
#include <algorithm>
#include <unordered_map>
#include <math.h>

#include "Types.h"
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* This file contains the main ShapeOp solver.*/
///////////////////////////////////////////////////////////////////////////////
namespace plb {
namespace npfem {
///////////////////////////////////////////////////////////////////////////////
// Forward Declarations
class LSSolver;
class Constraint;
class Force;
///////////////////////////////////////////////////////////////////////////////
class SHAPEOP_API Solver {
public:
    Solver();
    // Copy Constructor
    //Solver(Solver const& solverTmp);

    // Setters
    /* Set the timestep for the dynamics */
    void setTimeStep(Scalar timestep);
    /* Add a constraint to the solver and get back its id */
    int addConstraint(const std::shared_ptr<Constraint>& c);
    /* Add a force to the solver and get back its id */
    int addForces(const std::shared_ptr<Force>& f);
    /* Set the points */
    void setPoints(const Matrix3X& p);
    /* Set the velocities */
    void setVelocities(const Matrix3X& v);
    /* Set the connectivity list of the deformable object (triangle soup) */
    void setConnectivityList(
        const std::vector<std::vector<int>>& connectivity_list);
    /* -- */
    void set_Palabos_Forces(const Matrix3X& force_matrix);
    /* -- */
    void set_onSurfaceParticle(const std::vector<bool>& onSurfaceParticle);

    // Getters
    /* Get the timestep for the dynamics */
    Scalar getTimeStep();
    /* Get a constraint using its id */
    std::vector<std::shared_ptr<Constraint>>& getConstraints();
    std::shared_ptr<Constraint>& getConstraint(int id);
    /* Get a force using its id */
    std::shared_ptr<Force>& getForce(int id);
    /* Get the points */
    const Matrix3X& getPoints() const;
    /* -- */
    Matrix3X& getVelocities();
    /* Get the connectivity list */
    const std::vector<std::vector<int>>& getConnectivityList() const;
    /* -- */
    const Matrix3X& get_Palabos_Forces() const;
    /* -- */
    const std::vector<bool>& get_onSurfaceParticle() const;

	// Initialization
    /* Initialize the ShapeOp linear */
    bool initialize(  Scalar Calpha = 0.
                    , Scalar Cbeta = 0.
                    , Scalar timestep = 1.0
                    , Scalar rho = 1060.
                    , bool doModalAnalysis = false
                    , bool applyGlobalVolumeConservation = true
                    , Scalar globalVolumeConservationWeight = 1.
                    , Scalar externalSolverSpatialResolution = 1.
                    , Scalar externalSolverTimeResolution = 1. );
    std::vector<int> * get_mesh_graph();
    /* -- */
    void mass_lumping();

	// Solve
    /* Solve for one timestep */
    bool solve(  int m = 15
               , unsigned int max_iterations = 1000
               , unsigned int max_line_search_loops = 15
               , unsigned int max_attempts_to_solve_stagnation = 10
               , unsigned int convergence_window = 10
               , Scalar tol = 1e-6
               , Scalar gamma = 1e-4
               , Scalar gamma2 = 0.9
               , Scalar collisions_threshold_rep = 1.0
               , Scalar collisions_weight_rep = 20.
               , Scalar collisions_threshold_nonRep = 1.0
               , Scalar collisions_weight_nonRep = 20.
               , Scalar beta_morse = 40.
               , Scalar max_u_lb = 0.1
               , Scalar max_u_lb_fin = 0.1);
    /* -- */
    void computeInternalForces();
    /* -- */
    void PBD_Damping();
	/* -- */
    Scalar evalKinPlusPotEnergy(const Matrix3X& points);
    /* -- */
    Scalar evalObjective(const Matrix3X& x, const Matrix3X& y);
    /* -- */
    void GyrationTensor();

	// Collisions
	/* -- */
    void collisionDetection(Scalar collisions_threshold_rep, Scalar collisions_threshold_nonRep);
    /* -- */
    void addCollisionConstraints(Scalar collisions_threshold_rep, Scalar collisionsWeight_rep, Scalar collisions_threshold_nonRep, Scalar collisionsWeight_nonRep, Scalar beta_morse);
    /* -- */
    void removeCollisionConstraints();

    /* Manipulate the points */
    void shiftPoints(Vector3 vec);
    void rotatePoints(const std::string& axis, const Scalar& theta);

    // Misc
    int bodyID_, Palabos_iT_, bodyType_;
    Scalar Area0_, Volume0_, Area_, Volume_;
    void calculateArea(Scalar& Area);
    void calculateVolume(Scalar& Volume);
    void calculateVolumeFromSurfaceAndPrepareConstraint(Scalar& Volume);

    // Collisions (TODO: make them private ?)
    // Initial containers that hold the collision candidates
    std::vector<Matrix3X> verticesOfCollidingPieces_;
    std::vector<Matrix3X> normalsOfCollidingPieces_;

    // Modal Analysis
    VectorX eigenValues_;
    VectorX naturalPeriods_;
    VectorX dampingRatios_;
    MatrixXX eigenVectors_;

    //add by joel
    Matrix33 rotations_;
    // Collisions: Convert initial containers into suitable data structures
    // for spatial searching
    Matrix3X collidingCandidatesPoints_, normalsCollidingCandidatesPoints_;

	Scalar Calpha_, Cbeta_;
	Scalar delta_;
	Scalar rho_;

private:
    typedef std::vector<std::shared_ptr<Constraint>> Constraints;
    typedef std::vector<std::shared_ptr<Force>> Forces;
    typedef nanoflann::KDTreeEigenMatrixAdaptor<Eigen::Matrix<Scalar,
                                                Eigen::Dynamic,
                                                Eigen::Dynamic>,
                                                nanoflann::metric_L2>
        kd_tree;

    // Palabos interface
    Matrix3X Palabos_Forces_;
    // Object Connectivity (triangle soup)
    std::vector<std::vector<int>> connectivity_list_;
    // Vector of bools that distinguishes particles that lay on the Surface or
    // inside the Volume
    std::vector<bool> onSurfaceParticle_;

    Matrix3X points_;
    Matrix3X projections_;
    Matrix3X gradC_, deltaX_;
    Constraints constraints_;
    std::shared_ptr<LSSolver> solver_;

    // Translation Invariance (TI)
    Matrix3X pointsTI_;
    Vector3 centroid_;

    SparseMatrix M_, M_tilde_, M_star_;
    MatrixXX M_tilde_inv_;
    Matrix3X oldPoints_, oldPointsEWMA_;
    int oldPointsHistory_;
    Matrix3X velocities_;
    Matrix3X momentum_;
    Forces forces_;

    Scalar globalVolumeConservationWeight_;
    bool applyGlobalVolumeConservation_;
    
    // Gyration Tensor
    Matrix33 G_;
    Vector3 Xcm_, restShapeEigenvalues_;
    int G_flag_;

    // Quasi-Newton & LBFGS
    SparseMatrix L_, J_;
    Matrix3X f_int_nonePD_;

    MatrixXX N_dense_inv_;

    // The collidingPointsInfo stores the ID of the colliding point (key),
    // the colliding point and the colliding point normal
    // vector<>(_1_,_2_)
    std::unordered_map<int, std::vector<Vector3>> collidingPointsInfo_;

    // Avoid Instabilities in the coupling with other solver
    Scalar externalSolverSpatialResolution_;
    Scalar externalSolverTimeResolution_;
    Scalar max_u_lb_, max_u_lb_fin_;

    std::vector<int> *graph = NULL;
};
///////////////////////////////////////////////////////////////////////////////
} // namespace npfem
} // namespace plb
///////////////////////////////////////////////////////////////////////////////
#ifdef SHAPEOP_HEADER_ONLY
#include "Solver.cpp"
#endif
///////////////////////////////////////////////////////////////////////////////
#endif // SOLVER_H
///////////////////////////////////////////////////////////////////////////////
