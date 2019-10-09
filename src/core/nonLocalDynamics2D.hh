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

#ifndef NON_LOCAL_DYNAMICS_2D_HH
#define NON_LOCAL_DYNAMICS_2D_HH

#include "core/globalDefs.h"
#include "core/nonLocalDynamics2D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
NonLocalDynamics2D<T,Descriptor>::NonLocalDynamics2D(Dynamics<T,Descriptor>* baseDynamics_)
    : CompositeDynamics<T,Descriptor>(baseDynamics_, false)
{ }

template<typename T, template<typename U> class Descriptor>
bool NonLocalDynamics2D<T,Descriptor>::isNonLocal() const {
    return true;
}

template<typename T, template<typename U> class Descriptor>
void NonLocalDynamics2D<T,Descriptor>::prepareCollision(Cell<T,Descriptor>& cell) { }



template<typename T, template<typename U> class Descriptor>
NonLocalBoundaryDynamics2D<T,Descriptor>::NonLocalBoundaryDynamics2D(Dynamics<T,Descriptor>* baseDynamics_)
    : NonLocalDynamics2D<T,Descriptor>(baseDynamics_)
{ }

template<typename T, template<typename U> class Descriptor>
bool NonLocalBoundaryDynamics2D<T,Descriptor>::isBoundary() const {
    return true;
}

template<typename T, template<typename U> class Descriptor>
void NonLocalBoundaryDynamics2D<T,Descriptor>::nonLocalAction (
        plint iX, plint iY, BlockLattice2D<T,Descriptor>& lattice )
{
    boundaryCompletion(iX,iY, lattice);
}

}  // namespace plb

#endif  // NON_LOCAL_DYNAMICS_2D_HH

