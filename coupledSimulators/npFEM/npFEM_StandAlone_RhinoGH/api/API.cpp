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
#ifndef API_CPP
#define API_CPP
///////////////////////////////////////////////////////////////////////////////
#include "API.h"
#include "Solver.h"
#include "Constraint.h"
#include "Force.h"
///////////////////////////////////////////////////////////////////////////////
/* C structure that containts the C++ ShapeOp solver. */
struct ShapeOpSolver {
    /* A shared pointer to a C++ ShapeOp solver. */
    std::shared_ptr<plb::npfem::Solver> s;
};
///////////////////////////////////////////////////////////////////////////////
extern ShapeOpSolver* shapeop_create()
{
    ShapeOpSolver* solver = new ShapeOpSolver;
    solver->s = std::make_shared<plb::npfem::Solver>();
    return solver;
}
extern void shapeop_delete(ShapeOpSolver* op) { delete op; }
extern int shapeop_init(ShapeOpSolver* op)
{
    return static_cast<int>(!op->s->initialize());
}
extern int shapeop_initDynamic(ShapeOpSolver *op, ShapeOpScalar Calpha, ShapeOpScalar Cbeta, ShapeOpScalar timestep
    , ShapeOpScalar rho, bool doModalAnalysis
    , bool applyGlobalVolumeConservation, ShapeOpScalar globalVolumeConservationWeight)
{
    // Last arg for Modal Analysis: Set to false for GH env
    return static_cast<int>(
        !op->s->initialize(Calpha, Cbeta, timestep, rho, false, applyGlobalVolumeConservation, globalVolumeConservationWeight));
}
extern int shapeop_solve(ShapeOpSolver* op, int m, unsigned int max_iterations,
    unsigned int max_line_search_loops,
    unsigned int max_attempts_to_solve_stagnation, unsigned int convergence_window, 
    ShapeOpScalar tol, ShapeOpScalar gamma, ShapeOpScalar gamma2,
    ShapeOpScalar collisions_threshold, ShapeOpScalar collisions_weight)
{
    return static_cast<int>(!op->s->solve(m, max_iterations,
        max_line_search_loops, max_attempts_to_solve_stagnation, convergence_window, tol, gamma,
        gamma2, collisions_threshold, collisions_weight));
}
extern void shapeop_setPoints(
    ShapeOpSolver* op, ShapeOpScalar* points, int nb_points)
{
    Eigen::Map<plb::npfem::Matrix3X> p(points, 3, nb_points);
    op->s->setPoints(p);
}
extern void shapeop_getPoints(
    ShapeOpSolver* op, ShapeOpScalar* points, int nb_points)
{
    Eigen::Map<plb::npfem::Matrix3X> p(points, 3, nb_points);
    p = op->s->getPoints();
}
extern void shapeop_setTimeStep(ShapeOpSolver* op, ShapeOpScalar timestep)
{
    op->s->setTimeStep(timestep);
}
///////////////////////////////////////////////////////////////////////////////
extern int shapeop_addConstraint(ShapeOpSolver* op, const char* constraintType,
    int* ids, int nb_ids, ShapeOpScalar weight)
{

    const std::string ct(constraintType);
    std::vector<int> idv(ids, ids + nb_ids);
    std::shared_ptr<plb::npfem::Constraint> c
        = plb::npfem::Constraint::shapeConstraintFactory(
            ct, idv, weight, op->s->getPoints());
    if (!c) {
        return -1;
    }
    return op->s->addConstraint(c);
}
extern void shapeop_addObstacle(ShapeOpSolver* op, ShapeOpScalar* points,
    ShapeOpScalar* normals, int nb_points)
{
    Eigen::Map<plb::npfem::Matrix3X> p(points, 3, nb_points);
    Eigen::Map<plb::npfem::Matrix3X> n(normals, 3, nb_points);
    
    op->s->verticesOfCollidingPieces_.push_back(p);
    op->s->normalsOfCollidingPieces_.push_back(n);
}
extern shapeop_err shapeop_editConstraint(ShapeOpSolver* op,
    const char* constraintType, int constraint_id, const ShapeOpScalar* scalars,
    int nb_scl)
{
    if (strcmp(constraintType, "SurfaceMaterial") == 0) {
        auto c = std::dynamic_pointer_cast<plb::npfem::SurfaceMaterialConstraint>(
            op->s->getConstraint(constraint_id));
        if (!c) {
            return SO_UNMATCHING_CONSTRAINT_ID;
        }
        if (nb_scl != 5) {
            return SO_INVALID_ARGUMENT_LENGTH;
        }
        c->setRangeMin(scalars[0]);
        c->setRangeMax(scalars[1]);
        c->setMiu(scalars[2]);
        c->setLambda(scalars[3]);
        c->setKappa(scalars[4]);
        return SO_SUCCESS;
    }
    if (strcmp(constraintType, "VolumeMaterial") == 0) {
        auto c = std::dynamic_pointer_cast<plb::npfem::VolumeMaterialConstraint>(
            op->s->getConstraint(constraint_id));
        if (!c) {
            return SO_UNMATCHING_CONSTRAINT_ID;
        }
        if (nb_scl != 5) {
            return SO_INVALID_ARGUMENT_LENGTH;
        }
        c->setRangeMin(scalars[0]);
        c->setRangeMax(scalars[1]);
        c->setMiu(scalars[2]);
        c->setLambda(scalars[3]);
        c->setKappa(scalars[4]);
        return SO_SUCCESS;
    }
    if (strcmp(constraintType, "VolumeDamping") == 0) {
        auto c = std::dynamic_pointer_cast<plb::npfem::VolumeDampingConstraint>(
            op->s->getConstraint(constraint_id));
        if (!c) {
            return SO_UNMATCHING_CONSTRAINT_ID;
        }
        if (nb_scl != 5) {
            return SO_INVALID_ARGUMENT_LENGTH;
        }
        c->setRangeMin(scalars[0]);
        c->setRangeMax(scalars[1]);
        c->setMiu(scalars[2]);
        c->setLambda(scalars[3]);
        c->setKappa(scalars[4]);
        return SO_SUCCESS;
    }
    if (strcmp(constraintType, "TriangleARAP") == 0) {
        auto c = std::dynamic_pointer_cast<plb::npfem::TriangleARAPConstraint>(
            op->s->getConstraint(constraint_id));
        if (!c) {
            return SO_UNMATCHING_CONSTRAINT_ID;
        }
        if (nb_scl != 0) {
            return SO_INVALID_ARGUMENT_LENGTH;
        }
        return SO_SUCCESS;
    }
    if (strcmp(constraintType, "TetrahedronARAP") == 0) {
        auto c = std::dynamic_pointer_cast<plb::npfem::TetrahedronARAPConstraint>(
            op->s->getConstraint(constraint_id));
        if (!c) {
            return SO_UNMATCHING_CONSTRAINT_ID;
        }
        if (nb_scl != 0) {
            return SO_INVALID_ARGUMENT_LENGTH;
        }
        return SO_SUCCESS;
    }
    if (strcmp(constraintType, "EdgeStrainLimiting") == 0) {
        auto c
            = std::dynamic_pointer_cast<plb::npfem::EdgeStrainLimitingConstraint>(
                op->s->getConstraint(constraint_id));
        if (!c) {
            return SO_UNMATCHING_CONSTRAINT_ID;
        }
        if (nb_scl != 2) {
            return SO_INVALID_ARGUMENT_LENGTH;
        }
        c->setRangeMin(scalars[0]);
        c->setRangeMax(scalars[1]);
        return SO_SUCCESS;
    }
    if (strcmp(constraintType, "TriangleStrainLimiting") == 0) {
        auto c = std::
            dynamic_pointer_cast<plb::npfem::TriangleStrainLimitingConstraint>(
                op->s->getConstraint(constraint_id));
        if (!c) {
            return SO_UNMATCHING_CONSTRAINT_ID;
        }
        if (nb_scl != 2) {
            return SO_INVALID_ARGUMENT_LENGTH;
        }
        c->setRangeMin(scalars[0]);
        c->setRangeMax(scalars[1]);
        return SO_SUCCESS;
    }
    if (strcmp(constraintType, "TetrahedronStrainLimiting") == 0) {
        auto c = std::
            dynamic_pointer_cast<plb::npfem::TetrahedronStrainLimitingConstraint>(
                op->s->getConstraint(constraint_id));
        if (!c) {
            return SO_UNMATCHING_CONSTRAINT_ID;
        }
        if (nb_scl != 2) {
            return SO_INVALID_ARGUMENT_LENGTH;
        }
        c->setRangeMin(scalars[0]);
        c->setRangeMax(scalars[1]);
        return SO_SUCCESS;
    }
    if (strcmp(constraintType, "Area") == 0) {
        auto c = std::dynamic_pointer_cast<plb::npfem::AreaConstraint>(
            op->s->getConstraint(constraint_id));
        if (!c) {
            return SO_UNMATCHING_CONSTRAINT_ID;
        }
        if (nb_scl != 2) {
            return SO_INVALID_ARGUMENT_LENGTH;
        }
        c->setRangeMin(scalars[0]);
        c->setRangeMax(scalars[1]);
        return SO_SUCCESS;
    }
    if (strcmp(constraintType, "Volume") == 0) {
        auto c = std::dynamic_pointer_cast<plb::npfem::VolumeConstraint>(
            op->s->getConstraint(constraint_id));
        if (!c) {
            return SO_UNMATCHING_CONSTRAINT_ID;
        }
        if (nb_scl != 2) {
            return SO_INVALID_ARGUMENT_LENGTH;
        }
        c->setRangeMin(scalars[0]);
        c->setRangeMax(scalars[1]);
        return SO_SUCCESS;
    }
    if (strcmp(constraintType, "Bending") == 0) {
        auto c = std::dynamic_pointer_cast<plb::npfem::BendingConstraint>(
            op->s->getConstraint(constraint_id));
        if (!c) {
            return SO_UNMATCHING_CONSTRAINT_ID;
        }
        if (nb_scl != 2) {
            return SO_INVALID_ARGUMENT_LENGTH;
        }
        c->setRangeMin(scalars[0]);
        c->setRangeMax(scalars[1]);
        return SO_SUCCESS;
    }
    if (strcmp(constraintType, "Closeness") == 0) {
        auto c = std::dynamic_pointer_cast<plb::npfem::ClosenessConstraint>(
            op->s->getConstraint(constraint_id));
        if (!c) {
            return SO_UNMATCHING_CONSTRAINT_ID;
        }
        if (nb_scl != 3) {
            return SO_INVALID_ARGUMENT_LENGTH;
        }
        Eigen::Map<const plb::npfem::Vector3> p(scalars, 3, 1);
        c->setPosition(p);
        return SO_SUCCESS;
    }
    return SO_INVALID_CONSTRAINT_TYPE;
}
///////////////////////////////////////////////////////////////////////////////
extern int shapeop_addGravityForce(ShapeOpSolver* op, ShapeOpScalar* force)
{
    Eigen::Map<plb::npfem::Vector3> g(force, 3, 1);
    auto f = std::make_shared<plb::npfem::GravityForce>(g);
    return op->s->addForces(f);
}
extern int shapeop_addVertexForce(
    ShapeOpSolver* op, ShapeOpScalar* force, int id)
{
    Eigen::Map<plb::npfem::Vector3> g(force, 3, 1);
    auto f = std::make_shared<plb::npfem::VertexForce>(g, id);
    return op->s->addForces(f);
}
extern void shapeop_editVertexForce(
    ShapeOpSolver* op, int force_id, ShapeOpScalar* force, int id)
{
    Eigen::Map<plb::npfem::Vector3> g(force, 3, 1);
    auto f = std::dynamic_pointer_cast<plb::npfem::VertexForce>(
        op->s->getForce(force_id));
    f->setId(id);
    f->setForce(g);
}
///////////////////////////////////////////////////////////////////////////////
#endif // API_CPP
///////////////////////////////////////////////////////////////////////////////
