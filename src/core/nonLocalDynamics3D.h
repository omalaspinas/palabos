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

#ifndef NON_LOCAL_DYNAMICS_3D_H
#define NON_LOCAL_DYNAMICS_3D_H

#include "core/globalDefs.h"
#include "core/dynamics.h"
#include "atomicBlock/blockLattice3D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class NonLocalDynamics3D : public CompositeDynamics<T,Descriptor> {
public:
    NonLocalDynamics3D(Dynamics<T,Descriptor>* baseDynamics_);
    virtual bool isNonLocal() const;
    virtual void prepareCollision(Cell<T,Descriptor>& cell);

    virtual void nonLocalAction(plint iX, plint iY, plint iZ, BlockLattice3D<T,Descriptor>& lattice) =0;
};

/*
template<typename T, template<typename U> class Descriptor>
class NonLocalBoundaryDynamics3D : public NonLocalDynamics3D<T,Descriptor> {
public:
    NonLocalBoundaryDynamics3D(Dynamics<T,Descriptor>* baseDynamics_);
    virtual bool isBoundary() const;

    virtual void nonLocalAction(plint iX, plint iY, plint iZ, BlockLattice3D<T,Descriptor>& lattice);
    virtual void boundaryCompletion(plint iX, plint iY, plint iZ, BlockLattice3D<T,Descriptor>& lattice) =0;
};
*/

}  // namespace plb

#endif  // NON_LOCAL_DYNAMICS_3D_H

