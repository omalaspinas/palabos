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
#ifndef FORCE_H
#define FORCE_H
///////////////////////////////////////////////////////////////////////////////
#include "Types.h"
///////////////////////////////////////////////////////////////////////////////
/*
Forces: keep it like this for legacy reasons and compatibility with Grasshopper
*/
///////////////////////////////////////////////////////////////////////////////
namespace plb {
namespace npfem {
///////////////////////////////////////////////////////////////////////////////
/* 
Base class of any forces. This class defines interface of a ShapeOp force.
*/
class SHAPEOP_API Force {
public:
    virtual ~Force() { ; }
    /* Get force vector.*/
    virtual Vector3 get(const Matrix3X& positions, int id) const = 0;
    virtual int getId() const = 0;
};
///////////////////////////////////////////////////////////////////////////////
/* This class defines a constant force for all vertices. */
class SHAPEOP_API GravityForce : public Force {
public:
    /* Constructor taking the gravity vector as parameter. */
    GravityForce(const Vector3& f);
    virtual ~GravityForce() { ; }
    /* Get gravity vector. */
    virtual Vector3 get(
        const Matrix3X& /*positions*/, int /*id*/) const override final;
    virtual int getId() const override final;

private:
    Vector3 f_;
};
///////////////////////////////////////////////////////////////////////////////
/* This class defines a constant force for a unique vertex. */
class SHAPEOP_API VertexForce : public Force {
public:
    /* Constructor taking the force and the vertex id as parameters. */
    VertexForce(const Vector3& f = Vector3::Zero(), int id = -1);
    virtual ~VertexForce() { ; }
    /* Get force vector. */
    virtual Vector3 get(
        const Matrix3X& /*position*/, int id) const override final;
    /* Set a new vertex id. */
    void setId(int id);
    /* Set a new force.*/
    void setForce(const Vector3& f);
    virtual int getId() const override final;

private:
    Vector3 f_;
    int id_;
};
///////////////////////////////////////////////////////////////////////////////
} // namespace npfem
} // namespace plb
///////////////////////////////////////////////////////////////////////////////
#ifdef SHAPEOP_HEADER_ONLY
#include "Force.cpp"
#endif
///////////////////////////////////////////////////////////////////////////////
#endif // FORCE_H
///////////////////////////////////////////////////////////////////////////////
