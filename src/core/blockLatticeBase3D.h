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
 * Base class for the 2D BlockLattice and MultiBlockLattice -- header file.
 */
#ifndef BLOCK_LATTICE_BASE_3D_H
#define BLOCK_LATTICE_BASE_3D_H

#include "core/globalDefs.h"
#include "core/geometry3D.h"
#include "core/blockStatistics.h"

namespace plb {

template<typename T, template<typename U> class Descriptor> struct Dynamics;
template<typename T, template<typename U> class Descriptor> class Cell;

template<typename T, template<typename U> class Descriptor>
class BlockLatticeBase3D {
public:
    BlockLatticeBase3D();
    virtual ~BlockLatticeBase3D();
    void swap(BlockLatticeBase3D<T,Descriptor>& rhs);
public:
    virtual Cell<T,Descriptor>& get(plint iX, plint iY, plint iZ) =0;
    virtual Cell<T,Descriptor> const& get(plint iX, plint iY, plint iZ) const =0;
    virtual void specifyStatisticsStatus(Box3D domain, bool status) =0;
    virtual void collide(Box3D domain) =0;
    virtual void collide() =0;
    virtual void stream(Box3D domain) =0;
    virtual void stream() =0;
    virtual void collideAndStream(Box3D domain) =0;
    virtual void collideAndStream() =0;
    virtual void incrementTime() =0;
    TimeCounter& getTimeCounter();
    TimeCounter const& getTimeCounter() const;
private:
    TimeCounter timeCounter;
};

}  // namespace plb

#endif

