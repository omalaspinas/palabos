/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2019 FlowKit-Numeca Group Sarl
 * Copyright (C) 2011-2019 University of Geneva
 * 
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
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

/** \file
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- header file.
 */
#ifndef REGULARIZED_BOUNDARY_DYNAMICS_H
#define REGULARIZED_BOUNDARY_DYNAMICS_H

#include "core/globalDefs.h"
#include "boundaryCondition/boundaryDynamics.h"

namespace plb {

/// Regularized velocity boundary dynamics for a straight wall
template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
class RegularizedVelocityBoundaryDynamics :
    public VelocityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>
{
public:
    RegularizedVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics_,
                                        bool automaticPrepareCollision_=true);
    RegularizedVelocityBoundaryDynamics(HierarchicUnserializer& unserializer);

    /// Clone the object, based on its dynamic type
    virtual RegularizedVelocityBoundaryDynamics<T,Descriptor,direction,orientation>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T,Descriptor>& cell) const;
private:
    static int id;
};

/// Regularized velocity boundary dynamics for a straight wall
template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
class RegularizedVelocityConstRhoBoundaryDynamics :
    public VelocityDirichletConstRhoBoundaryDynamics<T,Descriptor,direction,orientation>
{
public:
    RegularizedVelocityConstRhoBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics_,
                                        bool automaticPrepareCollision_=true);
    RegularizedVelocityConstRhoBoundaryDynamics(HierarchicUnserializer& unserializer);

    /// Clone the object, based on its dynamic type
    virtual RegularizedVelocityConstRhoBoundaryDynamics<T,Descriptor,direction,orientation>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T,Descriptor>& cell) const;
private:
    static int id;
};

/// Regularized density Dirichlet boundary dynamics for a straight wall
template<typename T, template<typename U> class Descriptor,
         int direction, int orientation>
class RegularizedDensityBoundaryDynamics : public DensityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>
{
public:
    RegularizedDensityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics_,
                                       bool automaticPrepareCollision_=true);
    RegularizedDensityBoundaryDynamics(HierarchicUnserializer& unserializer);

    /// Clone the object, based on its dynamic type
    virtual RegularizedDensityBoundaryDynamics<T,Descriptor,direction,orientation>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T,Descriptor>& cell) const;
private:
    static int id;
};

}  // namespace plb

#endif  // REGULARIZED_BOUNDARY_DYNAMICS_H
