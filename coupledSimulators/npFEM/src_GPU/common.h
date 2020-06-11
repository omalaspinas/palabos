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
#ifndef COMMON_H
#define COMMON_H
///////////////////////////////////////////////////////////////////////////////
/**
This file is used to define macros and typedefs used by the C++ code and the C
API.*/
///////////////////////////////////////////////////////////////////////////////
/** Defines the scalar type used by the ShapeOp solver (float or double).*/
#ifndef SHAPEOP_SCALAR
// Change both defs
// Depending on the decision the solver may be more unstable
// All the testing is done with double
#define SHAPEOP_SCALAR double // double or float
#define DOUBLE_SHAPEOP // {FLOAT or DOUBLE}_SHAPEOP
#endif

#define MAX_ATTEMPT 7
#define MAX_LINE_SEARCH 5
#define GAMMA 0.0001
#define GAMMA2 0.9
#define TOL 0.00001
#define CONV_WINDOWS 5
#define MEM_SIZE 5

#define BL 0

namespace plb {
namespace npfem {

/** Defines the scalar type used by the ShapeOp solver (float or double).*/
typedef SHAPEOP_SCALAR ShapeOpScalar;
typedef double cuda_scalar;

}
}
///////////////////////////////////////////////////////////////////////////////
#define LU_GPU // Faster than the other 2
//#define CHOLESKY_GPU
//#define QR_GPU //By far slower
///////////////////////////////////////////////////////////////////////////////
/** Defines the API prefix for the current platform.*/
#if defined(_WIN32) || defined(_WIN64)
#pragma warning(disable : 4251) // Disable warnings about templates and std
                                // types exposed in the c++ interface.
#pragma warning(disable : 4244) // conversion from 'double' to
                                // 'plb::npfem::Scalar', possible loss of data
#pragma warning(disable : 4996) // dll window function
#ifdef SHAPEOP_EXPORT
#define SHAPEOP_API __declspec(dllexport)
#else
#define SHAPEOP_API __declspec(dllimport)
#endif
#else
#define SHAPEOP_API
#endif
///////////////////////////////////////////////////////////////////////////////
#ifdef SHAPEOP_HEADER_ONLY
#define SHAPEOP_INLINE inline
#else
#define SHAPEOP_INLINE
#endif
///////////////////////////////////////////////////////////////////////////////
#endif // COMMON_H
///////////////////////////////////////////////////////////////////////////////
