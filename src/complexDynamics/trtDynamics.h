/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
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

#ifndef TRT_DYNAMICS_H
#define TRT_DYNAMICS_H

#include "core/globalDefs.h"
#include "basicDynamics/isoThermalDynamics.h"


namespace plb {

/// Implementation of the TRT collision step
template<typename T, template<typename U> class Descriptor>
class TRTdynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    // sMinus is the relaxation time associated to the anti-symmetric part.
    // If a value of 0 is provided for sMinus, this parameter will default
    // to a value depending on omega, which guarantees that the wall is
    // exactly half-way between nodes, according to Eq. 11 in Pan et al 2006:
    // https://doi.org/10.1016/j.compfluid.2005.03.008
    // The original formulation of the model is probably in Ginzburg and
    // Adler 1994:
    // https://doi.org/10.1051/jp2:1994123.
    // It corresponds to a "magic parameter" of 3/16.
    TRTdynamics(T omega_, T sMinus_ = T());
    TRTdynamics(HierarchicUnserializer& unserializer);
    
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer& serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer& unserializer);
    /// Clone the object on its dynamic type.
    virtual TRTdynamics<T,Descriptor>* clone() const;
    
    /// Return a unique ID for this class.
    virtual int getId() const;
    
    /* *************** Collision and Equilibrium ************************* */
    
    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);
    
    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(Cell<T,Descriptor>& cell, T rhoBar,
                                 Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat);
    
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
private:
    T sMinus;
    static int id;
};

/// Implementation of incompressible TRT dynamics.
/** This is the TRT equivalent of IncBGKdynamics: the "rho" moment of the
 *  populations appears only as a pressure term in the equilibrium, while
 *  the other terms are multiplied by the constant rho0.
 **/
template<typename T, template<typename U> class Descriptor>
class IncTRTdynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    IncTRTdynamics(T omega_);
    IncTRTdynamics(HierarchicUnserializer& unserializer);
    
    /// Clone the object on its dynamic type.
    virtual IncTRTdynamics<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);
    
    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(Cell<T,Descriptor>& cell, T rhoBar,
                                 Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
                                 
    virtual bool velIsJ() const;

/* *************** Macroscopic variables ***************************** */
    
    /// Velocity is equal to j, not u.
    virtual void computeVelocity( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& u ) const;
private:
    static const T sMinus;
    static int id;
};

}  // namespace plb

#endif  // TRT_DYNAMICS_H

