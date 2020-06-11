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
#ifndef API_H
#define API_H
///////////////////////////////////////////////////////////////////////////////
// Boolean vars for C
#include <stdbool.h>

#include "common.h"
///////////////////////////////////////////////////////////////////////////////
/*
* This file implements a C API for the ShapeOp C++ library.
*
* To use the library you need to:
*
* 1) Create the solver with #shapeop_create
*
* 2) Set the vertices with #shapeop_setPoints
*
* 3) Setup the constraints and forces with #shapeop_addConstraint and #shapeop_editConstraint
*
* 4) Initalize the solver
*
* 5) Optimize with #shapeop_solve
*
* 6) Get back the vertices with #shapeop_getPoints
*
* 7) Delete the solver with #shapeop_delete
*/
///////////////////////////////////////////////////////////////////////////////
using namespace plb::npfem;
/* C structure that containts the C++ #ShapeOp::Solver. */
typedef struct ShapeOpSolver ShapeOpSolver;
///////////////////////////////////////////////////////////////////////////////
/* ShapeOp Success and Error type. This list might be extended. To simply test for errors, use =! SO_SUCCESS. */
typedef enum shapeop_err {
  SO_SUCCESS = 0,                 /* ShapeOp Success type indicating that no error happened. */
  SO_INVALID_CONSTRAINT_TYPE = 1, /* ShapeOp Error type indicating an invalid constraint type provided. */
  SO_INVALID_ARGUMENT_LENGTH = 2, /* ShapeOp Error type indicating an invalid length of an array argument. */
  SO_UNMATCHING_CONSTRAINT_ID = 3 /* ShapeOp Error type indicating that the constraint type does not match the id provided. */
} shapeop_err;
///////////////////////////////////////////////////////////////////////////////
#ifdef __cplusplus
extern "C" {
#endif
///////////////////////////////////////////////////////////////////////////////
// Solver
/* Create the ShapeOp solver. For more details see #ShapeOp::Solver. */
SHAPEOP_API ShapeOpSolver *shapeop_create();
/* Delete the ShapeOp solver. For more details see #ShapeOp::Solver. */
SHAPEOP_API void shapeop_delete(ShapeOpSolver *op);
/* Initialize the ShapeOp solver for static geometry processing. For more details see #ShapeOp::Solver. */
SHAPEOP_API int  shapeop_init(ShapeOpSolver *op);
/* Initialize the ShapeOp solver for dynamic geometry processing. For more details see #ShapeOp::Solver. */
SHAPEOP_API int  shapeop_initDynamic(ShapeOpSolver *op, ShapeOpScalar Calpha, ShapeOpScalar Cbeta, ShapeOpScalar timestep
    , ShapeOpScalar rho, bool doModalAnalysis
    , bool applyGlobalVolumeConservation , ShapeOpScalar globalVolumeConservationWeight);
/* Run the optimization. For more details see #ShapeOp::Solver. */
SHAPEOP_API int shapeop_solve(ShapeOpSolver* op, int m,
    unsigned int max_iterations, unsigned int max_line_search_loops,
    unsigned int max_attempts_to_solve_stagnation, unsigned int convergence_window, ShapeOpScalar tol,
    ShapeOpScalar gamma, ShapeOpScalar gamma2,
    ShapeOpScalar collisions_threshold, ShapeOpScalar collisions_weight);
/* Set the vertices to the ShapeOp solver. For more details see #ShapeOp::Solver. */
SHAPEOP_API void shapeop_setPoints(ShapeOpSolver *op, ShapeOpScalar *points, int nb_points);
/* Get the vertices back from the ShapeOp solver. For more details see #ShapeOp::Solver. */
SHAPEOP_API void shapeop_getPoints(ShapeOpSolver *op, ShapeOpScalar *points, int nb_points);
/* Set the timestep of the ShapeOp solver. For more details see #ShapeOp::Solver. */
SHAPEOP_API void shapeop_setTimeStep(ShapeOpSolver *op, ShapeOpScalar timestep);
///////////////////////////////////////////////////////////////////////////////
// Constraints
SHAPEOP_API int shapeop_addConstraint(ShapeOpSolver *op, const char *constraintType, int *ids, int nb_ids, ShapeOpScalar weight);
SHAPEOP_API void shapeop_addObstacle(ShapeOpSolver *op, ShapeOpScalar *points, ShapeOpScalar *normals, int nb_points);
SHAPEOP_API shapeop_err shapeop_editConstraint(ShapeOpSolver *op,
                                               const char *constraintType,
                                               int constraint_id,
                                               const ShapeOpScalar *scalars,
                                               int nb_scl);
///////////////////////////////////////////////////////////////////////////////
// Forces
/* Add a gravity force to the ShapeOp solver. For more details see #ShapeOp::GravityForce.
  \param op The ShapeOp Solver object
  \param force A c-style array of 3 #ShapeOpScalar's specifying the force vector in each dimension
*/
SHAPEOP_API int shapeop_addGravityForce(ShapeOpSolver *op, ShapeOpScalar *force);

/* Add a vertex force to the ShapeOp solver. For more details see #ShapeOp::VertexForce. */
SHAPEOP_API int shapeop_addVertexForce(ShapeOpSolver *op, ShapeOpScalar *force, int id);
/* Edit a vertex force previously added to the ShapeOp solver. For more details see #ShapeOp::VertexForce.*/
SHAPEOP_API void shapeop_editVertexForce(ShapeOpSolver *op, int force_id, ShapeOpScalar *force, int id);
///////////////////////////////////////////////////////////////////////////////
#ifdef __cplusplus
}
#endif
///////////////////////////////////////////////////////////////////////////////
#ifdef SHAPEOP_HEADER_ONLY
#include "API.cpp"
#endif
///////////////////////////////////////////////////////////////////////////////
#endif // API_H
///////////////////////////////////////////////////////////////////////////////
