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
#ifndef SHAPEOPWRAPPER_CPP
#define SHAPEOPWRAPPER_CPP

#include "shapeOpWrapper.h"

namespace plb {
namespace npfem {

// Handle file opening error
inline void IsOpen_HandleError(bool open, const char* file, int line)
{
    if (!open)
    {
        std::cout << "File Opening Problem in " << file << " at line " << line << std::endl;
        exit(EXIT_FAILURE);
    }
}
#define ISOPEN_HANDLE_ERROR(err) (IsOpen_HandleError(err, __FILE__, __LINE__))

void setPointsFromCSV(ShapeOp_Solver& s, std::string filename, bool best_fit_T)
{
    char c; // to eat the commas

    double x, y, z;
    std::vector<double> xv, yv, zv;

    int n_points = 0;
    std::ifstream points_file(filename.c_str());
    ISOPEN_HANDLE_ERROR(points_file.is_open());

    while (points_file >> x >> c >> y >> c >> z)
    {
        xv.push_back(x);
        yv.push_back(y);
        zv.push_back(z);
        n_points++;
    }
    points_file.close();

    plb::npfem::Matrix3X q(3, n_points);
    for (int i = 0; i < q.cols(); ++i)
        q.col(i) = plb::npfem::Vector3(xv[i], yv[i], zv[i]);

    if (!best_fit_T)
    {
        // Default Path
        // Set the points with the points that you read previously
        s.setPoints(q);
    }
    else
    {
        // This branch is to avoid setting the initial conditions with deformed bodies
        // such that we avoid energy bursts

        const plb::npfem::Matrix3X& p = s.getPoints();

        ///////////////////////////////////////////////////////////////////

        // See the least-squares fitting using SVD by Sorkine & Rabinovitch

        plb::npfem::Vector3 p_ = p.rowwise().mean();
        plb::npfem::Vector3 q_ = q.rowwise().mean();

        plb::npfem::Matrix3X x = p;
        x.colwise() -= p_;
        plb::npfem::Matrix3X y = q;
        y.colwise() -= q_;

        plb::npfem::Matrix33 S = x * plb::npfem::MatrixXX::Identity(n_points, n_points) * y.transpose();
        Eigen::JacobiSVD<plb::npfem::Matrix33> svd(S, Eigen::ComputeFullU | Eigen::ComputeFullV);
        plb::npfem::Matrix33 U = svd.matrixU();
        plb::npfem::Matrix33 V = svd.matrixV();

        plb::npfem::Vector3 diag(1., 1., (V*U.transpose()).determinant());
        plb::npfem::Matrix33 R = V * diag.asDiagonal() * U.transpose();

        plb::npfem::Vector3 t = q_ - R * p_;

        ///////////////////////////////////////////////////////////////////

		s.setPoints(R*p);
        s.shiftPoints(t);
    }
}

void setVelsFromCSV(ShapeOp_Solver& s, std::string filename)
{
    char c; // to eat the commas

    double x, y, z;
    std::vector<double> xv, yv, zv;

    int n_points = 0;
    std::ifstream vels_file(filename.c_str());
    ISOPEN_HANDLE_ERROR(vels_file.is_open());

    while (vels_file >> x >> c >> y >> c >> z)
    {
        xv.push_back(x);
        yv.push_back(y);
        zv.push_back(z);
        n_points++;
    }
    vels_file.close();

    plb::npfem::Matrix3X q(3, n_points);
    for (int i = 0; i < q.cols(); ++i)
        q.col(i) = plb::npfem::Vector3(xv[i], yv[i], zv[i]);
    
    s.setVelocities(q);
}

void savePointsToCSV(ShapeOp_Solver& s, std::string filename, size_t iT, size_t dt_ShapeOp)
{
    plb::npfem::Matrix3X points = s.getPoints();
    const plb::npfem::Matrix3X& vels = s.getVelocities();

    if (dt_ShapeOp != 1)
    {
        double delta = s.getTimeStep();
        plb::npfem::Matrix3X oldPoints = points - delta * vels;

        double delta_Palabos = delta / ((double)dt_ShapeOp);
        size_t proceedBy = iT % dt_ShapeOp;

        points = oldPoints + vels * delta_Palabos * (double)(proceedBy + 1);
    }

    std::ofstream points_file(filename.c_str());
    ISOPEN_HANDLE_ERROR(points_file.is_open());
    for (int i = 0; i < points.cols(); ++i) {
        points_file << points.col(i)[0] << "," << points.col(i)[1] << ","
                    << points.col(i)[2] << std::endl;
    }
    points_file.close();
}

void saveVelsToCSV(ShapeOp_Solver& s, std::string filename)
{
    const plb::npfem::Matrix3X& velocities = s.getVelocities();

    std::ofstream vels_file(filename.c_str());
    ISOPEN_HANDLE_ERROR(vels_file.is_open());
    for (int i = 0; i < velocities.cols(); ++i) {
        vels_file << velocities.col(i)[0] << "," << velocities.col(i)[1] << ","
            << velocities.col(i)[2] << std::endl;
    }
    vels_file.close();
}

void setConstraintsFromCSV(ShapeOp_Solver& s, std::string filename)
{
    std::ifstream constraints_file(filename.c_str());
    ISOPEN_HANDLE_ERROR(constraints_file.is_open());
    // read every constraint line-by-line
    std::string line;
    while (std::getline(constraints_file, line)) {
        // read details of the constraint
        std::string constraintType;
        std::vector<int> idI;
        plb::npfem::Scalar weight;
        std::vector<plb::npfem::Scalar> scalars;

        std::istringstream ss(line);
        std::string token;
        int field = 0;
        while (std::getline(ss, token, ',')) {
            if (field == 0) {
                constraintType = token;
            } else if (field == 1) {
                std::istringstream ss_tmp(token);
                std::string inds;
                while (std::getline(ss_tmp, inds, ' ')) {
                    idI.push_back(std::stoi(inds));
                }
            } else if (field == 2) {
                weight = std::stof(token);
            } else {
                if (token.length() != 0) {
                    std::istringstream ss_tmp(token);
                    std::string scalar;
                    while (std::getline(ss_tmp, scalar, ' ')) {
                        scalars.push_back(std::stof(scalar));
                    }
                }
            }
            ++field;
        }

        // Add constraint to the solver!
        if (constraintType.compare("VolumeMaterial") == 0) {
            auto c = std::make_shared<plb::npfem::VolumeMaterialConstraint>(
                idI, weight, s.getPoints());
            c->setRangeMin(scalars[0]);
            c->setRangeMax(scalars[1]);
            c->setMiu(scalars[2]);
            c->setLambda(scalars[3]);
            c->setKappa(scalars[4]);
            s.addConstraint(c);
        } else if (constraintType.compare("VolumeDamping") == 0) {
            auto c = std::make_shared<plb::npfem::VolumeDampingConstraint>(
                idI, weight, s.getPoints());
            c->setRangeMin(scalars[0]);
            c->setRangeMax(scalars[1]);
            c->setMiu(scalars[2]);
            c->setLambda(scalars[3]);
            c->setKappa(scalars[4]);
            s.addConstraint(c);
        } else if (constraintType.compare("SurfaceMaterial") == 0) {
            auto c = std::make_shared<plb::npfem::SurfaceMaterialConstraint>(
                idI, weight, s.getPoints());
            c->setRangeMin(scalars[0]);
            c->setRangeMax(scalars[1]);
            c->setMiu(scalars[2]);
            c->setLambda(scalars[3]);
            c->setKappa(scalars[4]);
            s.addConstraint(c);
        } else if (constraintType.compare("Volume") == 0) {
            auto c = std::make_shared<plb::npfem::VolumeConstraint>(
                idI, weight, s.getPoints());
            c->setRangeMin(scalars[0]);
            c->setRangeMax(scalars[1]);
            s.addConstraint(c);
        } else if (constraintType.compare("Area") == 0) {
            auto c = std::make_shared<plb::npfem::AreaConstraint>(
                idI, weight, s.getPoints());
            c->setRangeMin(scalars[0]);
            c->setRangeMax(scalars[1]);
            s.addConstraint(c);
        } else if (constraintType.compare("Bending") == 0) {
            auto c = std::make_shared<plb::npfem::BendingConstraint>(
                idI, weight, s.getPoints());
            c->setRangeMin(scalars[0]);
            c->setRangeMax(scalars[1]);
            s.addConstraint(c);
        } else if (constraintType.compare("EdgeStrainLimiting") == 0) {
            auto c = std::make_shared<plb::npfem::EdgeStrainLimitingConstraint>(
                idI, weight, s.getPoints());
            c->setRangeMin(scalars[0]);
            c->setRangeMax(scalars[1]);
            s.addConstraint(c);
        } else if (constraintType.compare("TriangleStrainLimiting") == 0) {
            auto c
                = std::make_shared<plb::npfem::TriangleStrainLimitingConstraint>(
                    idI, weight, s.getPoints());
            c->setRangeMin(scalars[0]);
            c->setRangeMax(scalars[1]);
            s.addConstraint(c);
        } else {
            std::cout << "Some Constraints are not valid: Check the "
                         "configuration file!"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    constraints_file.close();
}

void setConnectivityListFromCSV(ShapeOp_Solver& s, std::string filename)
{
    std::vector<std::vector<int>> connectivity_csv;
    int n_triangles = 0;
    std::ifstream connectivity_file(filename.c_str());
    ISOPEN_HANDLE_ERROR(connectivity_file.is_open());
    std::string line;
    while (std::getline(connectivity_file, line)) {
        std::vector<int> triangle;
        std::istringstream ss(line);
        std::string token;
        while (std::getline(ss, token, ',')) {
            triangle.push_back(std::stoi(token));
        }
        connectivity_csv.push_back(triangle);
        n_triangles++;
    }
    connectivity_file.close();

    s.setConnectivityList(connectivity_csv);
}

void setForcesFromCSV(ShapeOp_Solver& s, std::string filename)
{
    plb::npfem::Matrix3X forces(3, s.getPoints().cols());

    std::ifstream forces_file(filename.c_str());
    ISOPEN_HANDLE_ERROR(forces_file.is_open());
    // read every force line-by-line
    int i = 0;
    std::string line;
    while (std::getline(forces_file, line)) {
        std::vector<double> force;
        std::istringstream ss(line);
        std::string token;
        while (std::getline(ss, token, ',')) {
            force.push_back(std::stof(token));
        }
        forces.col(i) = plb::npfem::Vector3(force[0], force[1], force[2]);
        ++i;
    }
    forces_file.close();

    addVertexForce(s, forces);
}

void saveForcesToCSV(ShapeOp_Solver& s, std::string filename)
{
    const plb::npfem::Matrix3X& forces = s.get_Palabos_Forces();

    std::ofstream forces_file(filename.c_str());
    ISOPEN_HANDLE_ERROR(forces_file.is_open());
    for (int i = 0; i < forces.cols(); ++i) {
        forces_file << forces.col(i)[0] << "," << forces.col(i)[1] << ","
                    << forces.col(i)[2] << std::endl;
    }
    forces_file.close();
}

void setOnSurfaceParticle(ShapeOp_Solver& s, std::string filename)
{
    std::vector<bool> onSurfaceParticle;

    std::ifstream onSurfaceParticle_file(filename.c_str());
    ISOPEN_HANDLE_ERROR(onSurfaceParticle_file.is_open());

    std::string line;
    while (std::getline(onSurfaceParticle_file, line)) {
        onSurfaceParticle.push_back(std::stoi(line));
    }

    onSurfaceParticle_file.close();

    s.set_onSurfaceParticle(onSurfaceParticle);
}

void addVertexForce(ShapeOp_Solver& s, const plb::npfem::Matrix3X& forces, const int cell_id)
{
    // In the RBC implementation the Force Id coincides with the vertex Id
    // Both ways to edit Vertex forces are equivalent. I keep both for legacy
    // reasons and compatibility with the Grasshopper environment.

    // 1st way to edit the vertex forces
    s.set_Palabos_Forces(forces);
    /*
    // 2nd way to edit the vertex forces
    for (int i = 0; i < forces.cols(); ++i) {
        auto VertexForce
            = std::make_shared<plb::npfem::VertexForce>(forces.col(i), i);
        s.addForces(VertexForce);
    }
	*/
}

void editVertexForce(ShapeOp_Solver& s, const plb::npfem::Matrix3X& forces, const int cell_id)
{
    // In the RBC implementation the Force Id coincides with the vertex Id
    // Both ways to edit Vertex forces are equivalent. I keep both for legacy
    // reasons and compatibility with the Grasshopper environment.

    // 1st way to edit the vertex forces
    s.set_Palabos_Forces(forces);
    /*
    // 2nd way to edit the vertex forces
    for (int i = 0; i < forces.cols(); ++i) {
        auto VertexForce
            = std::dynamic_pointer_cast<plb::npfem::VertexForce>(s.getForce(i));
        // force Id stays unchanged
        VertexForce->setForce(forces.col(i));
    }
	*/
}

}
}

#endif
