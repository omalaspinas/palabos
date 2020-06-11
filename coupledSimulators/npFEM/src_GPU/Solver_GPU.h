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
#ifndef SOLVER_GPU_H
#define SOLVER_GPU_H
///////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <memory>
#include <cuda_runtime.h>

#include "Types.h"
#include "Constraint_Flattening.h"
#include "sparse_matrix.h"
#include "common.h"
#include "GPU_data.h"
///////////////////////////////////////////////////////////////////////////////
/** @file
This file contains the main ShapeOp solver.*/
///////////////////////////////////////////////////////////////////////////////
namespace plb {
namespace npfem {
///////////////////////////////////////////////////////////////////////////////
// Forward Declarations
class LSSolver;
class Constraint;
class Force;
typedef std::vector<std::shared_ptr<Constraint> > Constraints;
typedef std::vector<std::shared_ptr<Force> > Forces;
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatrixDyn;
///////////////////////////////////////////////////////////////////////////////
/** \brief ShapeOp Solver. This class implements the main ShapeOp solver based on \cite Bouaziz2012 and \cite Bouaziz2014.*/
class SHAPEOP_API Solver_GPU {
 public:
  Solver_GPU(){};
  Solver_GPU(const MatrixXX &points, 
			 const std::vector< std::vector<int> > &triangles, 
			 const std::vector<bool> &onSurfaceParticle, 
			 std::vector<std::shared_ptr<Constraint>>& constraints, Mesh_info params, float cbeta,
			 const int nb_cells = 1, Scalar phys_timestep = 1, Scalar shapeOp_time_step = 1, int bodyID = 0 , std::vector<int> *graph = NULL);

  void make_gpu_graph(std::vector<int> *graph, int nb);

  ~Solver_GPU();
  void add_collinding_points(ShapeOpScalar  *points, ShapeOpScalar *obstacle_normal, int *nb_per_cell, int nb);
  void free_collinding_points();
  /** \brief Add a constraint to the solver and get back its id.*/
  int addConstraint(const std::shared_ptr<Constraint> &c);
  /** \brief Get a constraint using its id.*/
  std::shared_ptr<Constraint> &getConstraint(int id);
  /** \brief Add a force to the solver and get back its id.*/
  int addForces(const std::shared_ptr<Force> &f);
  /** \brief Get a force using its id.*/
  std::shared_ptr<Force> &getForce(int id);
  /** \brief Set the points.*/
  void setPoints(const Matrix3X &p, const int nb_cells = 1);
  /** \brief Set the timestep for the dynamics.*/
  void setTimeStep(Scalar timestep);
  /** \brief Set the velocity damping for the dynamics.*/
  void setDamping(Scalar damping);
  /** \brief Set the connectivity list of the deformable object (triangle soup).*/
  void setConnectivityList(const std::vector<std::vector<int> > &connectivity_list);
  /** \brief Get the points.*/
  MatrixXXCuda *get_transformedPoints();
  const Matrix3X &getPoints() const;
  /** \brief Get the connectivity list.*/
  const std::vector<std::vector<int> > &getConnectivityList() const;
  /** \brief Initialize the ShapeOp linear system and the different parameters.
  \return true if successfull */

  bool initialize(Scalar timestep = 1.0);
  /** \brief Solve the constraint problem by projecting and merging.
    \return true if successfull */
  void compute_first_centroid();
  void set_initial_positions(const double *centers);
  void reset_position(const double *centers);
  void set_initial_velocity(double x, double y, double z);
  void set_gpu_starting_position(const Matrix3X &points, int cell);
  void set_gpu_starting_velocities(const Matrix3X &vels, int cell);
  ShapeOpScalar *get_centers();
  int get_nb_triangle();
  int get_nb_points();
  int get_nb_cells();

  void send_data_GPU();
  void get_data_from_GPU();

  void rotate_points(const cuda_scalar *matrices);

  double solve( unsigned int max_iterations = ((int)(1e3)), Scalar tol = 1e-6, bool Quasi_Newton = true, Scalar gamma = 0.3, 
              int max_line_search_loops = ((int)(1e3)), int m = 7, Scalar gamma2 = 0.9, Scalar collisions_weight = 500.);

  void make_periodic(float nx,float ny, float nz, float dx, int it);
  void set_Palabos_Forces(const Matrix3X &force_matrix, const int cell_id);
  void shiftPoints(Vector3 vec);
  void flatten_constraints();
  void rotatePoints(const std::string& axis, const Scalar& theta);
  void set_onSurfaceParticle(const std::vector<bool>& onSurfaceParticle);
  const std::vector<bool>& get_onSurfaceParticle() const;
  const MatrixXXCuda& getVelocities() const;
  const Matrix3X& get_Palabos_Forces() const;
  Scalar getTimeStep();
  // Misc
  int bodyID_, Palabos_iT_;

  //copy past by Joel from master cpu version
  std::vector<Matrix3X> verticesOfCollidingPieces_;
  std::vector<Matrix3X> normalsOfCollidingPieces_;
  
  cuda_scalar *force_npd_d_;
  short *triangles_d_;
  short *triangles_h_;
  int  nb_tri_ = 0;
  Mesh_info mesh_info_;

 private:
  //Palabos interface
  Scalar Cbeta_;
  Scalar delta_, solver_step_;
  //Object Connectivity (triangle soup)
  std::vector< std::vector<int> > connectivity_list_;
  std::vector<bool> onSurfaceParticle_;

  //Static
  Matrix3X points_;
  Matrix3X projections_;
  Constraints constraints_;
  std::shared_ptr<LSSolver> solver_;
  SparseMatrix At_;
  SparseMatrix N_;

  MatrixXX N_dense_inv_;
  
  //Dynamic
  bool dynamic_;
  SparseMatrix M_, M_tilde_, M_star_;
  MatrixXX M_tilde_inv_;
  Matrix3X oldPoints_;
  MatrixXXCuda velocities_;
  Matrix3X momentum_;
  Scalar masses_;
  Forces forces_;
  Scalar damping_;

  //Quasi-Newton & LBFGS
  SparseMatrix L_, J_;
  Matrix3X f_int_nonePD_;

  //Collisions
  Matrix3X f_collisions_;
  Scalar E_collisions_;
  std::vector<bool> collision_state_;

  //CUDA
  Scalar *cuda_points;
  Scalar *cuda_forces;
  cudaStream_t stream = 0;
  Constraint_flat *constraints_h_;
  Mesh_data mesh_data_d_, mesh_data_h_;
  Simulation_input simulation_input_d_, simulation_input_h_;
  Simulation_data simulation_data_d_, simulation_data_h_;
  Collision_data collision_data_h_, collision_data_d_;
  cuda_scalar *matrices_d_;
  graph_data graph_h, graph_d;

  MatrixT<Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>  Points_rowMajor_;
  MatrixXXCuda Palabos_Forces_;
  Eigen::Matrix<int, Eigen::Dynamic,    1> ConstraintType_eigen_;
  Eigen::Matrix<int, Eigen::Dynamic,    1> idO_eigen_;
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> rangeMin_eigen_;
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> rangeMax_eigen_;
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Scalar1_eigen_;
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> weight_eigen_;
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> E_nonePD_eigen_;
  Eigen::Matrix<cuda_scalar, Eigen::Dynamic, 1> A_;
  Eigen::Matrix<int, 4,    Eigen::Dynamic, Eigen::RowMajor> idI_eigen_;
  Eigen::Matrix<Scalar, 4, Eigen::Dynamic, Eigen::RowMajor> vectorx_eigen_;
  Eigen::Matrix<Scalar, 4, Eigen::Dynamic, Eigen::RowMajor> matrix22_eigen_;
  Eigen::Matrix<Scalar, 9, Eigen::Dynamic, Eigen::RowMajor> matrix33_eigen_;
};
///////////////////////////////////////////////////////////////////////////////
}
}
///////////////////////////////////////////////////////////////////////////////
#ifdef SHAPEOP_HEADER_ONLY
#include "Solver_GPU.cpp"
#endif
///////////////////////////////////////////////////////////////////////////////
#endif // SOLVER_GPU_H
///////////////////////////////////////////////////////////////////////////////
