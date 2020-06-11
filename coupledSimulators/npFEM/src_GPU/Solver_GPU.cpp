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
#ifndef SOLVER_GPU_CPP
#define SOLVER_GPU_CPP
///////////////////////////////////////////////////////////////////////////////
#include <time.h>
#include <string.h>
///////////////////////////////////////////////////////////////////////////////
#ifdef NPFEM_SA
#ifdef SHAPEOP_MSVC
// Console Output in the case of a dll
#include <windows.h>
#endif
#endif
///////////////////////////////////////////////////////////////////////////////
#include "Solver_GPU.h"
#include "LSSolver.h"
#include "Constraint.h"
#include "Force.h"
#include "device_utilities.h"
#include "sparse_matrix.h"
#include "GPU_data.h"

#include "palabos3D.h"
///////////////////////////////////////////////////////////////////////////////
#ifdef SHAPEOP_OPENMP
#ifdef SHAPEOP_MSVC
#define SHAPEOP_OMP_PARALLEL __pragma(omp parallel)
#define SHAPEOP_OMP_FOR __pragma(omp for)
#else
#define SHAPEOP_OMP_PARALLEL _Pragma("omp parallel")
#define SHAPEOP_OMP_FOR _Pragma("omp for")
#endif
#else
#define SHAPEOP_OMP_PARALLEL
#define SHAPEOP_OMP_FOR
#endif
///////////////////////////////////////////////////////////////////////////////
namespace plb {
namespace npfem {
///////////////////////////////////////////////////////////////////////////////
Solver_GPU::Solver_GPU(const MatrixXX &points,
		const std::vector< std::vector<int> > &triangles,
		const std::vector<bool> &onSurfaceParticle,
		std::vector<std::shared_ptr<Constraint>> &constraints, Mesh_info params, float cbeta,
		const int nb_cells, Scalar phys_timestep, Scalar shapeOp_time_step, int bodyId, std::vector<int> *graph): mesh_info_(params)
{
	solver_step_ = shapeOp_time_step;
	bodyID_ = bodyId;
	setPoints(points, nb_cells);
	constraints_ = constraints;
	setConnectivityList(triangles);
	set_onSurfaceParticle(onSurfaceParticle);
	make_gpu_graph(graph, points.cols());
	Cbeta_ = cbeta;
	initialize(phys_timestep*shapeOp_time_step);
}
/////////////////////////////////////////////////////////////////////////////
void Solver_GPU::make_gpu_graph(std::vector<int> *graph, int nb)
{
	if (graph == NULL) return;

	graph_h.indexs = new int[nb + 1];
	graph_h.indexs[0] = 0;

	for (int i = 1; i <= nb; i++)
		graph_h.indexs[i] = graph_h.indexs[i - 1] + graph[i - 1].size();

	graph_h.nb_edges = graph_h.indexs[nb];
	graph_h.edges = new int[graph_h.nb_edges];

	for (int i = 0; i < nb; i++) {
		for (int j = 0; j < graph[i].size(); j++) {
			graph_h.edges[graph_h.indexs[i] + j] = graph[i][j];
		}
	}
	
    send_graph(&graph_h, &graph_d, nb + 1);
}
////////////////////////////////////////////////////////////////////////////
Solver_GPU::~Solver_GPU()
{
  // TODO : More on this
  /*
  cudaFree(Palabos_Forces_d_), cudaFree(points_d_), cudaFree(projections_d_), cudaFree(oldPoints_d_), cudaFree(velocities_d_), cudaFree(momentum_d_), cudaFree(f_int_nonePD_d_);
  cudaFree(At_dense_d_), cudaFree(N_dense_d_), cudaFree(M_dense_d_), cudaFree(L_dense_d_), cudaFree(J_dense_d_);
  cudaFree(ConstraintType_d_), cudaFree(idO_d_), cudaFree(idI_d_);
  cudaFree(rangeMin_d_), cudaFree(rangeMax_d_), cudaFree(Scalar1_d_), cudaFree(weight_d_), cudaFree(E_nonePD_d_);
  cudaFree(vectorx_d_), cudaFree(matrix22_d_), cudaFree(matrix33_d_);

  delete Palabos_Forces_h_; delete points_h_; delete projections_h_; delete oldPoints_h_; delete velocities_h_; delete momentum_h_; delete f_int_nonePD_h_;
  delete At_dense_h_; delete N_dense_h_; delete M_dense_h_; delete L_dense_h_; delete J_dense_h_;
  delete constraints_h_;
  delete ConstraintType_h_; delete idO_h_; delete idI_h_;
  delete rangeMin_h_; delete rangeMax_h_; delete Scalar1_h_; delete weight_h_; delete E_nonePD_h_;
  delete vectorx_h_; delete matrix22_h_; delete matrix33_h_;
  
  delete cusolverInfo_;
  cusolverDnDestroy(cusolverDnHandle_);
  cublasDestroy(cublasHandle_);
  */
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE int Solver_GPU::addConstraint(const std::shared_ptr<Constraint> &c) 
{
  constraints_.push_back(c);
  return static_cast<int>(constraints_.size() - 1);
}

////////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE std::shared_ptr<Constraint> &Solver_GPU::getConstraint(int id) 
{
	return constraints_[id];
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver_GPU::flatten_constraints() 
{
	constraints_h_ = new Constraint_flat[mesh_info_.n_constraints];

	ConstraintType_eigen_ = Eigen::Matrix<int, Eigen::Dynamic, 1>::Zero(mesh_info_.n_constraints, 1);
	idO_eigen_ = Eigen::Matrix<int, Eigen::Dynamic, 1>::Zero(mesh_info_.n_constraints, 1);
	rangeMin_eigen_ = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Zero(mesh_info_.n_constraints, 1);
	rangeMax_eigen_ = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Zero(mesh_info_.n_constraints, 1);
	Scalar1_eigen_ = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Zero(mesh_info_.n_constraints, 1);
	weight_eigen_ = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Zero(mesh_info_.n_constraints, 1);
	E_nonePD_eigen_ = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Zero(mesh_info_.n_constraints, 1);
	A_ = Eigen::Matrix<cuda_scalar, Eigen::Dynamic, 1>::Zero(mesh_info_.n_constraints, 1);
	idI_eigen_ = Eigen::Matrix<int, 4, Eigen::Dynamic>::Zero(4, mesh_info_.n_constraints);
	vectorx_eigen_ = Eigen::Matrix<Scalar, 4, Eigen::Dynamic>::Zero(4, mesh_info_.n_constraints);
	matrix22_eigen_ = Eigen::Matrix<Scalar, 4, Eigen::Dynamic>::Zero(4, mesh_info_.n_constraints);
	matrix33_eigen_ = Eigen::Matrix<Scalar, 9, Eigen::Dynamic>::Zero(9, mesh_info_.n_constraints);

	SHAPEOP_OMP_PARALLEL
	{
		SHAPEOP_OMP_FOR for (int i = 0; i < mesh_info_.n_constraints; ++i)
		{
			// Common to every constraint
			constraints_h_[i].idO_ = constraints_[i]->get_idO();
			for (int j = 0; j < static_cast<int>(constraints_[i]->get_idI().size()); ++j) constraints_h_[i].idI_[j] = constraints_[i]->get_idI()[j];
			constraints_h_[i].weight_ = constraints_[i]->get_weight();

			if (constraints_[i]->ConstraintType_.compare("Area") == 0)
			{
				constraints_h_[i].ConstraintType_ = 1;
				constraints_h_[i].rangeMin_ = constraints_[i]->getMinRange();
				constraints_h_[i].rangeMax_ = constraints_[i]->getMaxRange();
				Scalar *tmp = constraints_[i]->getMatrix22().data();
				for (int j = 0; j < 4; ++j) constraints_h_[i].matrix22_[j] = tmp[j];
			}
			else if (constraints_[i]->ConstraintType_.compare("Volume") == 0)
			{
				constraints_h_[i].ConstraintType_ = 2;
				constraints_h_[i].rangeMin_ = constraints_[i]->getMinRange();
				constraints_h_[i].rangeMax_ = constraints_[i]->getMaxRange();
				Scalar *tmp = constraints_[i]->getMatrix33().data();
				for (int j = 0; j < 9; ++j) constraints_h_[i].matrix33_[j] = tmp[j];
			}
			else if (constraints_[i]->ConstraintType_.compare("Bending") == 0)
			{
				constraints_h_[i].ConstraintType_ = 3;
				constraints_h_[i].rangeMin_ = constraints_[i]->getMinRange();
				constraints_h_[i].rangeMax_ = constraints_[i]->getMaxRange();
				constraints_h_[i].Scalar1_ = constraints_[i]->getScalar1();
				constraints_h_[i].vectorx_ = constraints_[i]->getVectorX();
				Scalar *tmp = constraints_[i]->getVectorX();
				for (int j = 0; j < 4; ++j) constraints_h_[i].vectorx_[j] = tmp[j];
			}
			else if (constraints_[i]->ConstraintType_.compare("TriangleStrainLimiting") == 0)
			{
				constraints_h_[i].ConstraintType_ = 4;
				constraints_h_[i].rangeMin_ = constraints_[i]->getMinRange();
				constraints_h_[i].rangeMax_ = constraints_[i]->getMaxRange();
				Scalar *tmp = constraints_[i]->getMatrix22().data();
				for (int j = 0; j < 4; ++j) constraints_h_[i].matrix22_[j] = tmp[j];
			}
			else if (constraints_[i]->ConstraintType_.compare("SurfaceMaterial") == 0)
			{
				SurfaceMaterialConstraint *tp = (SurfaceMaterialConstraint*)constraints_[i].get();
				constraints_h_[i].ConstraintType_ = 5;
				constraints_h_[i].rangeMin_ = constraints_[i]->getMinRange();
				constraints_h_[i].rangeMax_ = constraints_[i]->getMaxRange();
				Scalar *tmp = constraints_[i]->getMatrix22().data();
				for (int j = 0; j < 4; ++j) constraints_h_[i].matrix22_[j] = tmp[j];

				A_[i] = tp->A_;
			}
		}
	}

	SHAPEOP_OMP_PARALLEL
	{
		SHAPEOP_OMP_FOR for (int i = 0; i < mesh_info_.n_constraints; ++i)
		{
			ConstraintType_eigen_[i] = constraints_h_[i].ConstraintType_;
			idO_eigen_[i] = constraints_h_[i].idO_;
			rangeMin_eigen_[i] = constraints_h_[i].rangeMin_;
			rangeMax_eigen_[i] = constraints_h_[i].rangeMax_;
			Scalar1_eigen_[i] = constraints_h_[i].Scalar1_;
			weight_eigen_[i] = constraints_h_[i].weight_;
			E_nonePD_eigen_[i] = constraints_h_[i].E_nonePD_;

			idI_eigen_.col(i)      = Eigen::Map<Eigen::Matrix<int   , 4, 1>>(constraints_h_[i].idI_);
			vectorx_eigen_.col(i)  = Eigen::Map<Eigen::Matrix<Scalar, 4, 1>>(constraints_h_[i].vectorx_);
			matrix22_eigen_.col(i) = Eigen::Map<Eigen::Matrix<Scalar, 4, 1>>(constraints_h_[i].matrix22_);
			matrix33_eigen_.col(i) = Eigen::Map<Eigen::Matrix<Scalar, 9, 1>>(constraints_h_[i].matrix33_);
		}
	}

  //printf("passed here 0.5\n?");

  mesh_data_h_.A = A_.data();
  mesh_data_h_.ConstraintType = ConstraintType_eigen_.data();
  mesh_data_h_.idO            = idO_eigen_.data();
  mesh_data_h_.rangeMin       = rangeMin_eigen_.data();
  mesh_data_h_.rangeMax       = rangeMax_eigen_.data();
  mesh_data_h_.Scalar1        = Scalar1_eigen_.data();
  mesh_data_h_.weight         = weight_eigen_.data();
  //E_nonePD_h_       = E_nonePD_eigen_.data();

  mesh_data_h_.idI           = idI_eigen_.data();
  mesh_data_h_.vectorx       = vectorx_eigen_.data();
  mesh_data_h_.matrix22      = matrix22_eigen_.data();
  mesh_data_h_.matrix33      = matrix33_eigen_.data();

  //delete[] constraints_h_;
}

////////////////////////////////////////////////////////////
void Solver_GPU::add_collinding_points(ShapeOpScalar  *points, ShapeOpScalar *obstacle_normal, int *nb_per_cell, int nb)
{
	//carfull this fonction allocate memory on GPU that has to be free later
	collision_data_d_.total_neighbours = nb;
	send_GPU_collinding_points(points, &collision_data_d_.colid_points, obstacle_normal, &collision_data_d_.colid_normals, nb_per_cell, collision_data_d_.nb_neighbours, nb, mesh_info_.nb_cells);
}
//////////////////////////////////////////////////////////////
void Solver_GPU::free_collinding_points()
{
	free_GPU_pointer(collision_data_d_.colid_points);
	free_GPU_pointer(collision_data_d_.colid_normals);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE int Solver_GPU::addForces(const std::shared_ptr<Force> &f)
{
  forces_.push_back(f);
  return static_cast<int>(forces_.size() - 1);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE std::shared_ptr<Force> &Solver_GPU::getForce(int id)
{
  return forces_[id];
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver_GPU::setPoints(const Matrix3X &p, const int nb_cells)
{
  //Vector3 offset;
	
  const int npoints = p.cols();
  points_ = p;
  mesh_info_.nb_cells = nb_cells;
  mesh_info_.n_points = npoints;
  Points_rowMajor_ = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(3*nb_cells, npoints);

  //printf("points pointer %lu %lu \n", Points_rowMajor_.data(), Points_rowMajor_.data() + 3*nb_cells*npoints);

  for (int i = 0; i < 3 * nb_cells; i += 3 ) {
	  //offset << i, 0, 0;
	  Points_rowMajor_.block(i,0, 3, p.cols()) = points_ ;//+ offset*MatrixXX::Ones(1, mesh_info_.n_points);
  }
  //std::cout << Points_rowMajor_.transpose() << std::endl;

  //Points_rowMajor_ = points_;
  simulation_input_h_.points = Points_rowMajor_.data();
  mesh_info_.sum_points = 32;
  while(mesh_info_.sum_points <= mesh_info_.n_points/2) mesh_info_.sum_points *= 2;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE const MatrixXXCuda& Solver_GPU::getVelocities() const
{
	return velocities_;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE const Matrix3X& Solver_GPU::get_Palabos_Forces() const
{
	return Palabos_Forces_;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver_GPU::setTimeStep(Scalar timestep) {
  delta_ = timestep;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE Scalar Solver_GPU::getTimeStep() { return delta_; }
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver_GPU::setDamping(Scalar damping) {
  damping_ = damping;
}
///////////////////////////////////////////////////////////////////////////////
MatrixXXCuda *plb::npfem::Solver_GPU::get_transformedPoints()
{
	return &Points_rowMajor_;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE const Matrix3X &Solver_GPU::getPoints() const {
  return points_;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver_GPU::setConnectivityList(const std::vector<std::vector<int> > &connectivity_list) {

    connectivity_list_ = connectivity_list;
	nb_tri_ = connectivity_list_.size();
	mesh_data_h_.triangles = new short[nb_tri_ * 3];
	mesh_info_.n_triangles = nb_tri_;

	for (int i = 0; i < nb_tri_; i++) {
		mesh_data_h_.triangles[i            ] = connectivity_list_[i][0];
		mesh_data_h_.triangles[i +   nb_tri_] = connectivity_list_[i][1];
		mesh_data_h_.triangles[i + 2*nb_tri_] = connectivity_list_[i][2];
	}
	mesh_info_.sum_triangles = 64;
	while (mesh_info_.sum_triangles <= nb_tri_/2) mesh_info_.sum_triangles *= 2;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE const std::vector<std::vector<int> > &Solver_GPU::getConnectivityList() const{
    return connectivity_list_;
}
void plb::npfem::Solver_GPU::make_periodic(float nx, float ny, float nz, float dx, int it){
	//std::cout << "size " << nx*dx << " " << ny*dx << " " << nz*dx << std::endl;
	device_make_periodic(simulation_input_d_.points, simulation_data_d_.center, nx*dx, ny*dx, nz*dx, mesh_info_.nb_cells, mesh_info_.n_points, it);
}
///////////////////////////////////////////////////////////////////////////////
void Solver_GPU::set_Palabos_Forces(const Matrix3X &force_matrix, const int cell_id){

  if (Palabos_Forces_.cols()>0){
	Palabos_Forces_.block(cell_id*3, 0, 3, force_matrix.cols()) = force_matrix;
  }
  else {
	  Palabos_Forces_ = force_matrix;
  }
}
///////////////////////////////////////////////////////////////////////////////
void Solver_GPU::shiftPoints(Vector3 vec) {
 
  for (int i = 0; i < points_.cols(); ++i) {
	  points_.col(i) += vec;
  }
  //shift_cells_points_g(simulation_input_d_.points, (ShapeOpScalar*)&vec, simulation_data_d_.center, mesh_info_.n_points, mesh_info_.nb_cells);
}
///////////////////////////////////////////////////////////////////////////////
double compute_volume_cpu(double *points, short *triangle, int n_tri, int n_p) {

	Vector3 *grad = new Vector3[n_p];

	for (int i = 0; i < n_p; i++) {
		grad[i] = Vector3(0, 0, 0);
	}
	double volume = 0;

	for (int i = 0; i < n_p; i++) {
		//std::cout << points[i]<<" " << points[i + n_p] <<" "<< points[i + 2 * n_p] << "\n";

	}

	for (int i = 0; i < n_tri; i++) {
		short p0 = triangle[i            ];
		short p1 = triangle[i +     n_tri];
		short p2 = triangle[i + 2 * n_tri];

		Vector3 a(points[p0], points[p0 + n_p], points[p0 + 2 * n_p]);
		Vector3 b(points[p1], points[p1 + n_p], points[p1 + 2 * n_p]);
		Vector3 c(points[p2], points[p2 + n_p], points[p2 + 2 * n_p]);

		//std::cout << p0 <<" " << p1 <<" " << p2 <<"\n";
		//std::cout << a << "\n\n " << b << "\n\n " << c << "\n\n " << "\n\n";

		Vector3 edge0 = b - a;
		Vector3 edge1 = c - a;

		Vector3 n = edge0.cross(edge1);
		double area = n.norm() / 2;
		n = n.normalized();
		//std::cout << n << "\n\n";

		volume += ((a + b + c).dot(n)*area) * 1 / 9;

		grad[p0] += (area*n) / 3;
		grad[p1] += (area*n) / 3;
		grad[p2] += (area*n) / 3;
		//double n[3] = 

	}
	//std::cout << grad[0] <<"\n";

	//std::cout << "volume cpu " << volume << "\n";
	return volume;
}



///////////////////////////////////////////////////////////////////////////////TODO

SHAPEOP_INLINE bool Solver_GPU::initialize(Scalar timestep){
#ifdef NPFEM_SA
#ifdef SHAPEOP_MSVC
    //Console output for dll
    AllocConsole();
    freopen("CONOUT$", "w", stdout);
#endif
#endif

  const int n_points	= mesh_info_.n_points;
  int n_constraints = static_cast<int>(constraints_.size());

  simulation_data_h_.center = new double[3 * mesh_info_.nb_cells];

#ifdef STREAM
  cudaStreamCreate(&stream);
  cudaHostAlloc(&cuda_points, 3*n_points*sizeof(ShapeOpScalar), cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc(&cuda_forces, 3*n_points*sizeof(ShapeOpScalar), cudaHostAllocWriteCombined | cudaHostAllocMapped);
#endif

  assert(n_points != 0);
  assert(n_constraints != 0);

#ifdef NPFEM_SA
  std::cout << "***********************************************************" << std::endl;
  std::cout << "GPU VERSION >> Number of Points & Constraints: " << n_points << " - " << n_constraints << std::endl;
  std::cout << "***********************************************************" << std::endl;
#endif // NPFEM_SA

  std::vector<Triplet> triplets;
  int idO = 0;
  for (int i = 0; i < n_constraints; ++i) {
	  // NONE_PD Energies (may add more : TODO)
	  if (constraints_[i]->get_ConstraintType().compare("SurfaceMaterial") == 0)
		  continue;
	  // PD Energies
	  constraints_[i]->addConstraint(triplets, idO);
  }
  mesh_info_.idO = idO;
  mesh_info_.volume = compute_volume_cpu(simulation_input_h_.points, mesh_data_h_.triangles, mesh_info_.n_triangles, mesh_info_.n_points);
  mesh_info_.n_constraints = n_constraints;
  projections_.setZero(3, idO);
  SparseMatrix A = SparseMatrix(idO, n_points);
  A.setFromTriplets(triplets.begin(), triplets.end());
  At_ = A.transpose();

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
  //printf("Laplacian GPU %f %f %f %f %f \n", Laplacian.coeff(223, 221), Laplacian.coeff(223, 222), Laplacian.coeff(223, 223), Laplacian.coeff(223, 224), Laplacian.coeff(223, 225));

  //Dynamic
  if (velocities_.cols() == 0) {
	  velocities_ = MatrixXXCuda::Zero(3*mesh_info_.nb_cells, n_points);
  }
  delta_ = timestep;
  momentum_ = Matrix3X::Zero(3, n_points);

#ifdef NPFEM_SA
  std::cout << " Calpha" << mesh_info_.calpha << "  Cbeta " << Cbeta_ << "  " << timestep << " rho " << mesh_info_.rho << " delta " << delta_ << std::endl;
#endif

  // Mass and Damping Matrices
  M_ = SparseMatrix(n_points, n_points);
  M_tilde_ = SparseMatrix(n_points, n_points);
  M_tilde_inv_ = MatrixXX(n_points, n_points);
  M_star_ = SparseMatrix(n_points, n_points);

  //mass_lumping(); // Go with uniform ditribution
  M_.setIdentity();
  mesh_info_.mass_total = mesh_info_.volume * mesh_info_.rho;
  M_ *= (mesh_info_.mass_total)/ n_points; // mass = Volume * density

  M_tilde_ = M_ + delta_ * (Cbeta_ * Laplacian);
  M_star_ = M_tilde_ / (delta_ * delta_);
  M_tilde_inv_ = (MatrixXX(M_tilde_)).inverse();

#ifdef NPFEM_SA
  printf("GPU Volume0_ %.17f  rho_ %f n_points %d  number of cells %d\n", mesh_info_.volume, mesh_info_.rho, n_points, mesh_info_.nb_cells);
#endif
  //printf("GPU M_tilde_inv %f %f %f %f %f \n", M_tilde_inv_(223, 221), M_tilde_inv_(223, 222), M_tilde_inv_(223, 223), M_tilde_inv_(223, 224), M_tilde_inv_(223, 225));
  // M_ /= delta_ * delta_;
  N_ = M_star_ + Laplacian;

  //printf("Mtild inv %f %f %f\n", M_tilde_inv_(0,0), M_tilde_inv_(0,1), M_tilde_inv_(0,2));

  //Quasi-Newton
  L_ = SparseMatrix(n_points, n_points);
  J_ = SparseMatrix(At_.rows(), At_.cols());
  L_ = At_ * A;
  J_ = At_;
 
  f_int_nonePD_ = Matrix3X::Zero(3*mesh_info_.nb_cells, n_points);

  //Palabos interface
  Palabos_Forces_ = MatrixXXCuda::Zero(3*mesh_info_.nb_cells,n_points);

  //Collisions
  f_collisions_ = Matrix3X::Zero(3, n_points);
  E_collisions_ = 0.;
  std::vector<bool>().swap(collision_state_);
  for (int i = 0; i < n_points; ++i) collision_state_.push_back(false);
  
  solver_ = std::make_shared<plb::npfem::SimplicialLDLTSolver>();

  //Prefactorize matrix
  solver_->initialize(N_);

  // GPU Initialization
  // Convert Eigen data structures to raw pointer in host memory
  simulation_input_h_.forces_ex = Palabos_Forces_.data();

  //copy past from CPU version by Joel, I dont know what it's needed for
  //bodyID_ = 0;
  Palabos_iT_ = 0;

  simulation_input_h_.velocities    = velocities_.data();

  MatrixDyn N_dense   = MatrixDyn(N_);
  MatrixDyn M_dense_  = MatrixDyn(M_);
  MatrixDyn L_dense_  = MatrixDyn(L_);
  MatrixDyn J_dense_  = MatrixDyn(J_);
  MatrixDyn M_star_dense_ = MatrixDyn(M_star_);
  
 //back to sparse but with our own structure usable on GPU
  mesh_data_h_.L = make_sparse_from_full( L_dense_.data(), L_dense_.rows(), L_dense_.cols());
  mesh_data_h_.J = make_sparse_from_full( J_dense_.data(), J_dense_.rows(), J_dense_.cols());
  mesh_data_h_.M_star = make_sparse_from_full(M_star_dense_.data(), M_star_dense_.rows(), M_star_dense_.cols());

  mesh_data_d_.L = mesh_data_h_.L;
  mesh_data_d_.J = mesh_data_h_.J;
  mesh_data_d_.M_star = mesh_data_h_.M_star;

  mesh_data_h_.M_tild_inv = M_tilde_inv_.data();
  //printf("init pointer l %p %p \n", &mesh_data_d_.L.value, mesh_data_h_.L.value);
  
  VectorX diag = M_dense_.diagonal();
  mesh_data_h_.mass = diag.data();

  N_dense_inv_   = N_dense.inverse();
  mesh_data_h_.N_inv_dense = N_dense_inv_.data();
  //Flatten constraints, quite similar to the .data() of the Eigen Matrices!

  flatten_constraints();
  // After converting all the needed data to raw pointers, we can work without problems on the gpu!
  
  GPU_Init(&mesh_info_,
		   &mesh_data_d_,        &mesh_data_h_,
           &simulation_input_d_, &simulation_input_h_,
           &simulation_data_d_ ,
           &collision_data_d_  , &collision_data_h_, &matrices_d_);

  //return solver_->info() == Eigen::Success;
  return true;
}
//////////////////////////////////////////////////////////////////////////////
void plb::npfem::Solver_GPU::compute_first_centroid() {
	first_center(mesh_info_, simulation_input_d_, simulation_data_d_, simulation_data_d_.center);
}
//////////////////////////////////////////////////////////////////////////////
void plb::npfem::Solver_GPU::set_initial_positions(const double *centers){
	set_cells_initial_position(centers, simulation_data_d_.center, mesh_info_.nb_cells, mesh_info_, mesh_data_d_, simulation_input_d_ , simulation_data_d_);
		
}
void plb::npfem::Solver_GPU::reset_position(const double *centers) {
	//reset_position(const double *center_h, double *center_d, Simulation_input *input_h, Mesh_info info, Mesh_data mesh, Simulation_input *input_d, Simulation_data sim)
	reset_position_d(centers, simulation_data_d_.center, points_.data(), mesh_info_, mesh_data_d_, simulation_input_d_, simulation_data_d_);
}
///////////////////////////////////////////////////////////////////////////////
//TODO make it multible cell
void Solver_GPU::set_initial_velocity(double x, double y, double z){
	Vector3 speed;
	speed << x, y, z;
	printf("init vell nb %d \n", mesh_info_.n_points);
	velocities_ = speed*MatrixXX::Ones(1, mesh_info_.n_points);
	std::cout << velocities_.col(0) << std::endl;
	std::cout << "couille"<< std::endl;
	simulation_input_h_.velocities  = velocities_.data();
	external_forces_from_Host_to_Device(mesh_info_.n_points, simulation_input_h_.velocities, simulation_input_d_.velocities, 0);
	//velocities_.row(1) = Matrix3X::Ones(1, n_points_)*y;
	//velocities_.row(2) = Matrix3X::Ones(1, n_points_)*z;
}
ShapeOpScalar *plb::npfem::Solver_GPU::get_centers(){
	return simulation_data_h_.center;
}

int plb::npfem::Solver_GPU::get_nb_triangle(){
	return mesh_info_.n_triangles;
}
int plb::npfem::Solver_GPU::get_nb_points()
{
	return mesh_info_.n_points;
}
int plb::npfem::Solver_GPU::get_nb_cells(){
	return mesh_info_.nb_cells;
}

//////////////////////////////////////////////////////////////////////////////
void Solver_GPU::send_data_GPU() {
	simulation_input_h_.forces_ex = Palabos_Forces_.data();
	memcpy(cuda_forces, simulation_input_h_.forces_ex, 3 * mesh_info_.n_points * sizeof(ShapeOpScalar));
	external_forces_from_Host_to_Device(mesh_info_.n_points, cuda_forces, simulation_input_d_.forces_ex, stream);
}
///////////////////////////////////////////////////////////////////////////////
void Solver_GPU::get_data_from_GPU() {
	simulation_input_h_.points = Points_rowMajor_.data();
	points_from_Device_to_Host(mesh_info_.n_points, simulation_input_d_.points, simulation_input_h_.points, stream);
	points_.block(0, 0, 3, Points_rowMajor_.cols()) = Points_rowMajor_;
}
//////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver_GPU::rotate_points(const cuda_scalar *matrices) {

	points_time_mat_3X3(matrices_d_,  matrices, simulation_input_d_.points, mesh_info_.nb_cells, mesh_info_.n_points);
}

///////////////////////////////////////////////////////////////////////////////
// solve computes one time step! The algorithms are iterative and so we need some iterations until convergence to the min of the variational problem
SHAPEOP_INLINE double Solver_GPU::solve(unsigned int max_iterations, Scalar tol, bool Quasi_Newton, Scalar gamma, int max_line_search_loops, int m, Scalar gamma2, Scalar collisions_weight) 
{

	clock_t start = clock();

	#ifdef DEBUG
		std::cout << "force_gpu  ex\n" << Palabos_Forces_.col(32) << std::endl;
	#endif

	simulation_input_h_.forces_ex = Palabos_Forces_.data();
	//printf("palabos force %f %f %f  delta_ %f\n", Palabos_Forces_(0, 223), Palabos_Forces_(1, 223), Palabos_Forces_(2, 223), delta_ );
	/*
	std::cout << Palabos_Forces_.col(56) << "\n";

	for (int i = 56; i < 60; i++) {
		printf("%f %f %f \n", Palabos_Forces_h_[i], Palabos_Forces_h_[i + Palabos_Forces_.cols()], Palabos_Forces_h_[i + 2*Palabos_Forces_.cols()]);
	}
	*/
	//compute_volume_cpu(simulation_input_h_.points, mesh_data_h_.triangles, mesh_info_.n_triangles, mesh_info_.n_points);

	#ifdef STREAM
		memcpy(cuda_forces, simulation_input_h_.forces_ex, 3 * mesh_info_.n_points * sizeof(ShapeOpScalar));
		external_forces_from_Host_to_Device(mesh_info_.n_points, cuda_forces, simulation_input_d_.forces_ex, stream);
	#else
		external_forces_from_Host_to_Device(mesh_info_.n_points*mesh_info_.nb_cells, simulation_input_h_.forces_ex, simulation_input_d_.forces_ex, stream);
	#endif	
	//debug_matrix_from_gpu(N_dense_d_, n_points_);
	if (graph_h.nb_edges > 0) {
		median_filter(graph_d, simulation_input_d_.forces_ex, mesh_info_.nb_cells, mesh_info_.n_points);
	}
	
	compute_next_frame_rbc(&mesh_info_, &mesh_data_d_, &simulation_input_d_,
		&simulation_data_d_, &collision_data_d_, delta_, solver_step_, max_iterations, stream);
	
	simulation_input_h_.points = Points_rowMajor_.data();

	#ifdef STREAM
		//points_from_Device_to_Host(mesh_info_.n_points, simulation_input_d_.points, cuda_points, stream);
		//cudaStreamSynchronize(stream);
		//memcpy(simulation_input_h_.points, cuda_points, 3 * mesh_info_.n_points * sizeof(ShapeOpScalar));
	#else
		simulation_input_h_.points = Points_rowMajor_.data();
		simulation_input_h_.velocities = velocities_.data();
		points_from_Device_to_Host(mesh_info_.n_points*mesh_info_.nb_cells, simulation_input_d_.velocities, simulation_input_h_.velocities, stream);
		points_from_Device_to_Host(mesh_info_.n_points*mesh_info_.nb_cells, simulation_input_d_.points, simulation_input_h_.points, stream);
        points_from_Device_to_Host(mesh_info_.nb_cells, simulation_data_d_.center, simulation_data_h_.center, stream);

		//printf("%f %f %f | %d \n\n", simulation_data_h_.center[0], simulation_data_h_.center[1], simulation_data_h_.center[2], bodyID_ );
		
		for(int i = 0; i < 3*mesh_info_.nb_cells; i+=3){
			/*
			std::ofstream ofile;
			ofile.open("./tmp/body_" + plb::util::val2str(bodyID_ + i/3) + "_ComPos.log", std::ofstream::out | std::ofstream::app);
	
			ofile << plb::util::val2str(Palabos_iT_) + "," +
					 plb::util::val2str(simulation_data_h_.center[i]) + "," +
					 plb::util::val2str(simulation_data_h_.center[i + 1]) + "," +
					 plb::util::val2str(simulation_data_h_.center[i + 2]) << std::endl;

			ofile.close();
			*/
			/*
			plb::global::logfile_nonparallel("body_" + plb::util::val2str(bodyID_ + i/3) + "_ComPos.log").flushEntry(plb::util::val2str(Palabos_iT_) + "," +
					plb::util::val2str(simulation_data_h_.center[i    ]) + "," +
					plb::util::val2str(simulation_data_h_.center[i + 1]) + "," +
					plb::util::val2str(simulation_data_h_.center[i + 2]));
			*/
		//points_ += velocities_;
		}
		
		#ifdef DEBUG
			std::cout << "velo GPU: \n" << velocities_.col(102) << std::endl;
		#endif
		//points_.block(0, 0, 3, Points_rowMajor_.cols()) = Points_rowMajor_;
	#endif

	/*
	for (int i = 56; i < 57; i++) {
		std::cout << points_.col(i) << "\n";
		//printf("%f %f %f \n", points_[i], points_[i + Palabos_Forces_.cols()], points_[i + 2 * Palabos_Forces_.cols()]);
	}
	*/
	clock_t end = clock();
	//std::cout << "temps one iter" << (double)(end - start) / CLOCKS_PER_SEC << "\n";
	// DEBUGGING!
	//std::cout << "\n*******************************" << std::endl;
	//std::cout << "Frobenius Norm of projections_: " << projections_.norm() << std::endl;
	//std::cout << "some Projected point: " << proj.block(0, 0, 1, 3) << std::endl;
	//std::cout << "some Projected point: " << projections_.block(0, 0, 1, 3) << std::endl;
	double norm_p = Points_rowMajor_.block(0, 0, 3, Points_rowMajor_.cols()).norm();
	//std::cout << "Frobenius Norm of points_: " << norm_p << std::endl;
	//std::cout << "velocity : " << velocities_.col(0) << std::endl;
	//std::cout << "Frobenius Norm of vel_: " << velocities_.norm() << std::endl;
	//plb::global::logfile_nonparallel("gpu_norm.dat").flushEntry(plb::util::val2str(norm_p));
	//std::cout << "Frobenius Norm of points_: " << Points_rowMajor_.block(3, 0, 3, Points_rowMajor_.cols()).norm() << std::endl;
	//std::cout << "*******************************\n" << std::endl;

	return (double)(end - start);
}
///////////////////////////////////////////////////////////////////////////////
//Rotat point can be done on GPU
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver_GPU::rotatePoints(const std::string& axis, const Scalar& theta)
{
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
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void Solver_GPU::set_onSurfaceParticle(
	const std::vector<bool>& onSurfaceParticle)
{
	onSurfaceParticle_ = onSurfaceParticle;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE const std::vector<bool>& Solver_GPU::get_onSurfaceParticle() const
{
	return onSurfaceParticle_;
}
///////////////////////////////////////////////////////////////////////////////
//TODO tester ca
void Solver_GPU::set_gpu_starting_position(const Matrix3X &points, int cell){
	MatrixXXCuda tp_point = points;
	points_from_Host_to_Device(points.cols(), simulation_input_d_.points, tp_point.data(), cell);
}

void plb::npfem::Solver_GPU::set_gpu_starting_velocities(const Matrix3X & vels, int cell){
	MatrixXXCuda tp_vel = vels;
	points_from_Host_to_Device(vels.cols(), simulation_input_d_.velocities, tp_vel.data(), cell);
}
///////////////////////////////////////////////////////////////////////////////
}
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
#endif // SOLVER_GPU_CPP
///////////////////////////////////////////////////////////////////////////////
