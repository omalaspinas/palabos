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

/** \file
 * Helper functions for domain initialization -- header file.
 */
#ifndef BC_WRAPPER_2D_HH
#define BC_WRAPPER_2D_HH

#include "plbWrapper/lattice/bcWrapper2D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
OuterBoxBC<T,Descriptor>::OuterBoxBC (
        OnLatticeBoundaryCondition2D<T,Descriptor>* bc_ )
    : bc(bc_)
{ }

template<typename T, template<typename U> class Descriptor>
OuterBoxBC<T,Descriptor>::OuterBoxBC(OuterBoxBC<T,Descriptor> const& rhs)
    : bc (rhs.bc->clone())
{ }

template<typename T, template<typename U> class Descriptor>
OuterBoxBC<T,Descriptor>::~OuterBoxBC()
{ 
    delete bc;
}

template<typename T, template<typename U> class Descriptor>
OuterBoxBC<T,Descriptor>& OuterBoxBC<T,Descriptor>::operator= (
        OuterBoxBC<T,Descriptor> const& rhs )
{
    delete bc;
    bc = rhs.bc->clone();
    return *this;
}

template<typename T, template<typename U> class Descriptor>
OuterBoxBC<T,Descriptor>* OuterBoxBC<T,Descriptor>::clone() const
{
    return new OuterBoxBC<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void OuterBoxBC<T,Descriptor>::setVelocityCondition (
        MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain )
{
    bc->setVelocityConditionOnBlockBoundaries(lattice, domain);
}

template<typename T, template<typename U> class Descriptor>
void OuterBoxBC<T,Descriptor>::setPressureCondition (
        MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain )
{
    bc->setPressureConditionOnBlockBoundaries(lattice, domain);
}

template<typename T, template<typename U> class Descriptor>
OuterBoxBC<T,Descriptor>* generateRegularizedBC() {
    return new OuterBoxBC<T,Descriptor>(createLocalBoundaryCondition2D<T,Descriptor>());
}

}  // namespace plb

#endif  // BC_WRAPPER_2D_HH
