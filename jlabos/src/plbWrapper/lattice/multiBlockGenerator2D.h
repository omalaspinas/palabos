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
 * Generator functions for all types of multi-blocks, to make them accessible to SWIG;
 * header file.
 */

#ifndef MULTI_BLOCK_GENERATORS_2D_H
#define MULTI_BLOCK_GENERATORS_2D_H

#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiDataField2D.h"


namespace plb {

template<typename T, template<typename U> class Descriptor>
MultiBlockLattice2D<T,Descriptor>*
    generateMultiBlockLattice2D(Box2D const& domain, Dynamics<T,Descriptor> const* backgroundDynamics);

template<typename T1, template<typename U> class Descriptor, typename T2>
MultiNTensorField2D<T2>*
    generateNTensorFieldFromLattice2D (
            MultiBlockLattice2D<T1,Descriptor> const& lattice,
            Box2D const& domain, plint ndim );

}  // namespace plb

#endif  // MULTI_BLOCK_GENERATORS_2D_H
