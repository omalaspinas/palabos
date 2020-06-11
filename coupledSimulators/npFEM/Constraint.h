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
#ifndef CONSTRAINT_H
#define CONSTRAINT_H
///////////////////////////////////////////////////////////////////////////////
#include <memory>
#include <cmath>

#include "Types.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
///////////////////////////////////////////////////////////////////////////////
namespace plb {
namespace npfem {
///////////////////////////////////////////////////////////////////////////////
class SHAPEOP_API Constraint {
public:
    static std::shared_ptr<Constraint> shapeConstraintFactory(
        const std::string& ConstraintType, const std::vector<int>& idI,
        Scalar weight, const Matrix3X& positions);
    Constraint(const std::vector<int>& idI, Scalar weight);
    virtual ~Constraint() {}

    virtual void project(const Matrix3X& positions, Matrix3X& projections,
        Matrix3X& f_int_nonePD, const Matrix3X& oldPositions)
        = 0;
    virtual void addConstraint(std::vector<Triplet>& triplets, int& idO) const = 0;

    std::size_t nIndices() const { return idI_.size(); }
    Scalar get_E_nonePD() const { return E_nonePD_; }
    std::vector<int> get_idI() const { return idI_; }
    Scalar get_weight() const { return weight_; }
    int get_idO() const { return idO_; }
    std::string get_ConstraintType() { return ConstraintType_; }

    // GPU-related
    virtual Scalar getMinRange() = 0;
    virtual Scalar getMaxRange() = 0;
    virtual Scalar getScalar1() = 0;
    virtual Matrix22 getMatrix22() = 0;
    virtual Matrix33 getMatrix33() = 0;
    virtual Scalar* getVectorX() = 0;

    std::string ConstraintType_;
protected:
    std::vector<int> idI_;
    Scalar weight_;
    mutable int idO_;
    Scalar E_nonePD_;
    mutable bool PDSys_Build;
};
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* Area conservation and membrane material imposition */
class SHAPEOP_API SurfaceMaterialConstraint : public Constraint {
public:
    SurfaceMaterialConstraint(const std::vector<int>& idI, Scalar weight,
        const Matrix3X& positions, Scalar rangeMin = 1.0, Scalar rangeMax = 1.0,
        Scalar miu = 1.0, Scalar lambda = 1.0, Scalar kappa = 1.0);
    virtual ~SurfaceMaterialConstraint() {}

    virtual void project(const Matrix3X& positions, Matrix3X& projections,
        Matrix3X& f_int_nonePD, const Matrix3X& oldPositions) override final;
    virtual void addConstraint(
        std::vector<Triplet>& triplets, int& idO) const override final;

    void setRangeMin(Scalar rMin) { rangeMin_ = rMin; }
    void setRangeMax(Scalar rMax) { rangeMax_ = rMax; }
    void setMiu(Scalar miu) { miu_ = miu; }
    void setLambda(Scalar lambda) { lambda_ = lambda; }
    void setKappa(Scalar kappa) { kappa_ = kappa; }

    void calculateArea(const Matrix3X& positions, Scalar& area);
    // This is used only in the case of a 2d material simulation
    void mass_lumping(
        const Matrix3X& positions, std::vector<Triplet>& triplets);

    virtual Scalar getMinRange() override final { return rangeMin_; }
    virtual Scalar getMaxRange() override final { return rangeMax_; }
    virtual Scalar getScalar1() override final { return 1.; }
    virtual Matrix22 getMatrix22() override final { return rest_; }
    virtual Matrix33 getMatrix33() override final { return Matrix33::Zero(); }
    virtual Scalar* getVectorX() override final { return 0; }

    Scalar A_;
private:
    Matrix22 rest_;
    Scalar rangeMin_;
    Scalar rangeMax_;
    Scalar miu_;
    Scalar lambda_;
    Scalar kappa_;
};
///////////////////////////////////////////////////////////////////////////////
/* Volume conservation and volumetric material imposition */
class SHAPEOP_API VolumeMaterialConstraint : public Constraint {
public:
    VolumeMaterialConstraint(const std::vector<int>& idI, Scalar weight,
        const Matrix3X& positions, Scalar rangeMin = 1.0, Scalar rangeMax = 1.0,
        Scalar miu = 1.0, Scalar lambda = 1.0, Scalar kappa = 1.0);
    virtual ~VolumeMaterialConstraint() {}

    virtual void project(const Matrix3X& positions, Matrix3X& projections,
        Matrix3X& f_int_nonePD, const Matrix3X& oldPositions) override final;
    virtual void addConstraint(
        std::vector<Triplet>& triplets, int& idO) const override final;

    void setRangeMin(Scalar rMin) { rangeMin_ = rMin; }
    void setRangeMax(Scalar rMax) { rangeMax_ = rMax; }
    void setMiu(Scalar miu) { miu_ = miu; }
    void setLambda(Scalar lambda) { lambda_ = lambda; }
    void setKappa(Scalar kappa) { kappa_ = kappa; }

    void calculateVolume(const Matrix3X& positions, Scalar& volume);
    void mass_lumping(const Matrix3X& positions, std::vector<Triplet>& triplets);

    virtual Scalar getMinRange() override final { return rangeMin_; }
    virtual Scalar getMaxRange() override final { return rangeMax_; }
    virtual Scalar getScalar1() override final { return 1.; }
    virtual Matrix22 getMatrix22() override final { return Matrix22::Zero(); }
    virtual Matrix33 getMatrix33() override final { return rest_; }
    virtual Scalar* getVectorX() override final { return 0; }

private:
    Matrix33 rest_;
    Scalar V_;
    Scalar rangeMin_;
    Scalar rangeMax_;
    Scalar miu_;
    Scalar lambda_;
    Scalar kappa_;
};
///////////////////////////////////////////////////////////////////////////////
// I experimented with SurfaceDamping as well, but the solver is VERY unstable
// This is most stable approach, approach the quality of Rayleigh Damping
class SHAPEOP_API VolumeDampingConstraint : public Constraint {
public:
    VolumeDampingConstraint(const std::vector<int>& idI, Scalar weight,
        const Matrix3X& positions, Scalar rangeMin = 1.0, Scalar rangeMax = 1.0,
        Scalar miu = 1.0, Scalar lambda = 1.0, Scalar kappa = 1.0);
    virtual ~VolumeDampingConstraint() {}

    virtual void project(const Matrix3X& positions, Matrix3X& projections,
        Matrix3X& f_int_nonePD, const Matrix3X& oldPositions) override final;
    virtual void addConstraint(
        std::vector<Triplet>& triplets, int& idO) const override final;

    void setRangeMin(Scalar rMin) { rangeMin_ = rMin; }
    void setRangeMax(Scalar rMax) { rangeMax_ = rMax; }
    void setMiu(Scalar miu) { miu_ = miu; }
    void setLambda(Scalar lambda) { lambda_ = lambda; }
    void setKappa(Scalar kappa) { kappa_ = kappa; }

    virtual Scalar getMinRange() override final { return rangeMin_; }
    virtual Scalar getMaxRange() override final { return rangeMax_; }
    virtual Scalar getScalar1() override final { return 1.; }
    virtual Matrix22 getMatrix22() override final { return Matrix22::Zero(); }
    virtual Matrix33 getMatrix33() override final { return rest_; }
    virtual Scalar* getVectorX() override final { return 0; }

private:
    Matrix33 rest_;
    Scalar V_;
    Scalar rangeMin_;
    Scalar rangeMax_;
    Scalar miu_;
    Scalar lambda_;
    Scalar kappa_;
};
///////////////////////////////////////////////////////////////////////////////
// Like Closeness Constraint
class SHAPEOP_API CollisionConstraint : public Constraint {
public:
    CollisionConstraint(
        const std::vector<int>& idI, Scalar weight, const Matrix3X& positions);
    virtual ~CollisionConstraint() {}

    virtual void project(const Matrix3X& positions, Matrix3X& projections,
        Matrix3X& f_int_nonePD, const Matrix3X& oldPositions) override final;
    virtual void addConstraint(
        std::vector<Triplet>& triplets, int& idO) const override final;

    void setCollindingPoint(const Vector3& position);
    void setCollindingPointNormal(const Vector3& normal);
    void setParams(const Scalar& collisions_threshold_rep, const Scalar& collisionsWeight_rep,
        const Scalar& collisions_threshold_nonRep, const Scalar& collisionsWeight_nonRep,
        const Scalar& beta_morse);

    virtual Scalar getMinRange() override final { return 1.; }
    virtual Scalar getMaxRange() override final { return 1.; }
    virtual Scalar getScalar1() override final { return 1.; }
    virtual Matrix22 getMatrix22() override final { return Matrix22::Zero(); }
    virtual Matrix33 getMatrix33() override final { return Matrix33::Zero(); }
    virtual Scalar* getVectorX() override final { return 0; }

private:
    Vector3 collidingPoint_;
    Vector3 collidingPointNormal_;
    Scalar thisWeight_;
    Scalar collisions_threshold_rep_, collisionsWeight_rep_, collisions_threshold_nonRep_, collisionsWeight_nonRep_, beta_morse_;
};
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* Below are the original ShapeOp Energies (Quadratic)                       */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* Project the deformation gradient on SO(3). As-Rigid-as-possible way of
 * deformation. */
/* Energy Density function = Σ(λi-1)^2. */
class SHAPEOP_API TriangleARAPConstraint : public Constraint {
public:
    TriangleARAPConstraint(
        const std::vector<int>& idI, Scalar weight, const Matrix3X& positions);
    virtual ~TriangleARAPConstraint() {}
    virtual void project(const Matrix3X& positions, Matrix3X& projections,
        Matrix3X& f_int_nonePD, const Matrix3X& oldPositions) override final;
    virtual void addConstraint(
        std::vector<Triplet>& triplets, int& idO) const override final;

    virtual Scalar getMinRange() override final { return 1.; }
    virtual Scalar getMaxRange() override final { return 1.; }
    virtual Scalar getScalar1() override final { return 1.; }
    virtual Matrix22 getMatrix22() override final { return rest_; }
    virtual Matrix33 getMatrix33() override final { return Matrix33::Zero(); }
    virtual Scalar* getVectorX() override final { return 0; }

private:
    Matrix22 rest_;
};
///////////////////////////////////////////////////////////////////////////////
class SHAPEOP_API TetrahedronARAPConstraint : public Constraint {
public:
    TetrahedronARAPConstraint(
        const std::vector<int>& idI, Scalar weight, const Matrix3X& positions);
    virtual ~TetrahedronARAPConstraint() {}
    virtual void project(const Matrix3X& positions, Matrix3X& projections,
        Matrix3X& f_int_nonePD, const Matrix3X& oldPositions) override final;
    virtual void addConstraint(
        std::vector<Triplet>& triplets, int& idO) const override final;

    virtual Scalar getMinRange() override final { return 1.; }
    virtual Scalar getMaxRange() override final { return 1.; }
    virtual Scalar getScalar1() override final { return 1.; }
    virtual Matrix22 getMatrix22() override final { return Matrix22::Zero(); }
    virtual Matrix33 getMatrix33() override final { return rest_; }
    virtual Scalar* getVectorX() override final { return 0; }

private:
    Matrix33 rest_;
};
///////////////////////////////////////////////////////////////////////////////
/* Edge strain constraint. Constrains the distance between two points to a
 * range. */
class SHAPEOP_API EdgeStrainLimitingConstraint : public Constraint {
public:
    EdgeStrainLimitingConstraint(const std::vector<int>& idI, Scalar weight,
        const Matrix3X& positions, Scalar rangeMin = 1.0,
        Scalar rangeMax = 1.0);
    virtual ~EdgeStrainLimitingConstraint() {}

    virtual void project(const Matrix3X& positions, Matrix3X& projections,
        Matrix3X& f_int_nonePD, const Matrix3X& oldPositions) override final;
    virtual void addConstraint(
        std::vector<Triplet>& triplets, int& idO) const override final;

    void setEdgeLength(Scalar length);
    void setRangeMin(Scalar rMin) { rangeMin_ = rMin; }
    void setRangeMax(Scalar rMax) { rangeMax_ = rMax; }

    virtual Scalar getMinRange() override final { return rangeMin_; }
    virtual Scalar getMaxRange() override final { return rangeMax_; }
    virtual Scalar getScalar1() override final { return 1.; }
    virtual Matrix22 getMatrix22() override final { return Matrix22::Zero(); }
    virtual Matrix33 getMatrix33() override final { return Matrix33::Zero(); }
    virtual Scalar* getVectorX() override final { return 0; }

private:
    Scalar rest_;
    Scalar rangeMin_;
    Scalar rangeMax_;
};
///////////////////////////////////////////////////////////////////////////////
/* A mesh-independent triangle strain-limiting constraint. */
class SHAPEOP_API TriangleStrainLimitingConstraint : public Constraint {
public:
    TriangleStrainLimitingConstraint(const std::vector<int>& idI, Scalar weight,
        const Matrix3X& positions, Scalar rangeMin = 1.0,
        Scalar rangeMax = 1.0);
    virtual ~TriangleStrainLimitingConstraint() {}

    virtual void project(const Matrix3X& positions, Matrix3X& projections,
        Matrix3X& f_int_nonePD, const Matrix3X& oldPositions) override final;
    virtual void addConstraint(
        std::vector<Triplet>& triplets, int& idO) const override final;

    void setRangeMin(Scalar rMin) { rangeMin_ = rMin; }
    void setRangeMax(Scalar rMax) { rangeMax_ = rMax; }

    virtual Scalar getMinRange() override final { return rangeMin_; }
    virtual Scalar getMaxRange() override final { return rangeMax_; }
    virtual Scalar getScalar1() override final { return 1.; }
    virtual Matrix22 getMatrix22() override final { return rest_; }
    virtual Matrix33 getMatrix33() override final { return Matrix33::Zero(); }
    virtual Scalar* getVectorX() override final { return 0; }

private:
    Matrix22 rest_;
    Scalar rangeMin_;
    Scalar rangeMax_;
};
///////////////////////////////////////////////////////////////////////////////
/* A mesh-independent tetrahedron strain-limiting constraint. */
class SHAPEOP_API TetrahedronStrainLimitingConstraint : public Constraint {
public:
    TetrahedronStrainLimitingConstraint(const std::vector<int>& idI,
        Scalar weight, const Matrix3X& positions, Scalar rangeMin = 1.0,
        Scalar rangeMax = 1.0);
    virtual ~TetrahedronStrainLimitingConstraint() {}

    virtual void project(const Matrix3X& positions, Matrix3X& projections,
        Matrix3X& f_int_nonePD, const Matrix3X& oldPositions) override final;
    virtual void addConstraint(
        std::vector<Triplet>& triplets, int& idO) const override final;

    void setRangeMin(Scalar rMin) { rangeMin_ = rMin; }
    void setRangeMax(Scalar rMax) { rangeMax_ = rMax; }

    virtual Scalar getMinRange() override final { return rangeMin_; }
    virtual Scalar getMaxRange() override final { return rangeMax_; }
    virtual Scalar getScalar1() override final { return 1.; }
    virtual Matrix22 getMatrix22() override final { return Matrix22::Zero(); }
    virtual Matrix33 getMatrix33() override final { return rest_; }
    virtual Scalar* getVectorX() override final { return 0; }

private:
    Matrix33 rest_;
    Scalar rangeMin_;
    Scalar rangeMax_;
};
///////////////////////////////////////////////////////////////////////////////
/* Area constraint. Limits the area of a triangle to a range. */
class SHAPEOP_API AreaConstraint : public Constraint {
public:
    AreaConstraint(const std::vector<int>& idI, Scalar weight,
        const Matrix3X& positions, Scalar rangeMin = 1.0,
        Scalar rangeMax = 1.0);
    virtual ~AreaConstraint() {}

    virtual void project(const Matrix3X& positions, Matrix3X& projections,
        Matrix3X& f_int_nonePD, const Matrix3X& oldPositions) override final;
    virtual void addConstraint(std::vector<Triplet>& triplets, int& idO) const override final;

    void calculateArea(const Matrix3X& positions, Scalar& area);

    void setRangeMin(Scalar rMin) { rangeMin_ = rMin; }
    void setRangeMax(Scalar rMax) { rangeMax_ = rMax; }
    const Matrix22& get_first() const { return rest_; }

    virtual Scalar getMinRange() override final { return rangeMin_; }
    virtual Scalar getMaxRange() override final { return rangeMax_; }
    virtual Scalar getScalar1() override final { return 1.; }
    virtual Matrix22 getMatrix22() override final { return rest_; }
    virtual Matrix33 getMatrix33() override final { return Matrix33::Zero(); }
    virtual Scalar* getVectorX() override final { return 0; }

private:
    Matrix22 rest_;
    Scalar rangeMin_;
    Scalar rangeMax_;
};
///////////////////////////////////////////////////////////////////////////////
/* Volume constraint. Limits the volume of a tetrahedron to a range. */
class SHAPEOP_API VolumeConstraint : public Constraint {
public:
    VolumeConstraint(const std::vector<int>& idI, Scalar weight,
        const Matrix3X& positions, Scalar rangeMin = 1.0,
        Scalar rangeMax = 1.0);
    virtual ~VolumeConstraint() {}

    virtual void project(const Matrix3X& positions, Matrix3X& projections,
        Matrix3X& f_int_nonePD, const Matrix3X& oldPositions) override final;
    virtual void addConstraint(std::vector<Triplet>& triplets, int& idO) const override final;

    void calculateVolume(const Matrix3X& positions, Scalar& volume);

    void setRangeMin(Scalar rMin) { rangeMin_ = rMin; }
    void setRangeMax(Scalar rMax) { rangeMax_ = rMax; }

    virtual Scalar getMinRange() override final { return rangeMin_; }
    virtual Scalar getMaxRange() override final { return rangeMax_; }
    virtual Scalar getScalar1() override final { return 1.; }
    virtual Matrix22 getMatrix22() override final { return Matrix22::Zero(); }
    virtual Matrix33 getMatrix33() override final { return rest_; }
    virtual Scalar* getVectorX() override final { return 0; }

private:
    Matrix33 rest_;
    Scalar rangeMin_;
    Scalar rangeMax_;
};
///////////////////////////////////////////////////////////////////////////////
/* Bending constraint. Limits the bending between two neighboring triangles. */
  /** \brief Constraint constructor. The target bend is set to the bend between the triangles spanned by four vertices in the parameter positions. The parameters rangeMin and rangeMax can be used to specify a target range for the bend [rangeMin*target_bend,rangeMax*target_bend]
    The bending constraint applies to two neighboring triangles sharing an edge
    \param idI are four indices of the vertices of the two triangles ordered as follows:
     <pre>
      |        id2        |
      |       /   \       |
      |     id0---id1     |
      |       \   /       |
      |        id3        |
     </pre>
    \param idI are four indices of the vertices of the tetrahedron
    \param weight The weight of the constraint to be added relative to the other constraints.
    \param positions The positions of all the n vertices stacked in a 3 by n matrix.
    \param rangeMin The factor to determine the minimal bend: rangeMin*target_bend.
    \param rangeMax The factor to determine the maximal bend: rangeMax*target_bend.
  */
class SHAPEOP_API BendingConstraint : public Constraint {
public:
    BendingConstraint(const std::vector<int>& idI, Scalar weight,
        const Matrix3X& positions, Scalar rangeMin = 1.0,
        Scalar rangeMax = 1.0);
    virtual ~BendingConstraint() {}

    virtual void project(const Matrix3X& positions, Matrix3X& projections,
        Matrix3X& f_int_nonePD, const Matrix3X& oldPositions) override final;
    virtual void addConstraint(
        std::vector<Triplet>& triplets, int& idO) const override final;

    void setRangeMin(Scalar rMin) { rangeMin_ = rMin; }
    void setRangeMax(Scalar rMax) { rangeMax_ = rMax; }

    virtual Scalar getMinRange() override final { return rangeMin_; }
    virtual Scalar getMaxRange() override final { return rangeMax_; }
    virtual Scalar getScalar1() override final { return n_; }
    virtual Matrix22 getMatrix22() override final { return Matrix22::Zero(); }
    virtual Matrix33 getMatrix33() override final { return Matrix33::Zero(); }
    virtual Scalar* getVectorX() override final { return w_.data(); }

private:
    VectorX w_;
    Scalar n_;
    Scalar rangeMin_;
    Scalar rangeMax_;
};
///////////////////////////////////////////////////////////////////////////////
/** Closeness constraint. Constrains a vertex to a position in space. Projects
 * onto a given rest position.*/
class SHAPEOP_API ClosenessConstraint : public Constraint {
public:
    ClosenessConstraint(
        const std::vector<int>& idI, Scalar weight, const Matrix3X& positions);
    virtual ~ClosenessConstraint() {}

    virtual void project(const Matrix3X& /*positions*/, Matrix3X& projections,
        Matrix3X& f_int_nonePD, const Matrix3X& oldPositions) override final;
    virtual void addConstraint(
        std::vector<Triplet>& triplets, int& idO) const override final;

    void setPosition(const Vector3& position);
    Vector3 getPosition() const;

    virtual Scalar getMinRange() override final { return 1.; }
    virtual Scalar getMaxRange() override final { return 1.; }
    virtual Scalar getScalar1() override final { return 1.; }
    virtual Matrix22 getMatrix22() override final { return Matrix22::Zero(); }
    virtual Matrix33 getMatrix33() override final { return Matrix33::Zero(); }
    virtual Scalar* getVectorX() override final { return 0; }

private:
    Vector3 rest_;
};
///////////////////////////////////////////////////////////////////////////////
} // namespace npfem
} // namespace plb
///////////////////////////////////////////////////////////////////////////////
#ifdef SHAPEOP_HEADER_ONLY
#include "Constraint.cpp"
#endif
///////////////////////////////////////////////////////////////////////////////
#endif // CONSTRAINT_H
///////////////////////////////////////////////////////////////////////////////
