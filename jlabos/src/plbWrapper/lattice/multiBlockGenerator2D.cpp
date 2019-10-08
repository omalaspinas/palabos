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

#ifdef COMPILE_2D

/** \file
 * Generator functions for all types of multi-blocks, to make them accessible to SWIG;
 * template instantiation.
 */
#include "plbWrapper/lattice/multiBlockGenerator2D.h"
#include "plbWrapper/lattice/multiBlockGenerator2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"


namespace plb {

template
MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D>*
    generateMultiBlockLattice2D (
            Box2D const& domain,
            Dynamics<FLOAT_T,descriptors::DESCRIPTOR_2D> const* backgroundDynamics );

template
MultiNTensorField2D<FLOAT_T>*
    generateNTensorFieldFromLattice2D (
            MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D> const& lattice,
            Box2D const& domain, plint ndim );

template
MultiNTensorField2D<int>*
    generateNTensorFieldFromLattice2D (
            MultiBlockLattice2D<FLOAT_T,descriptors::DESCRIPTOR_2D> const& lattice,
            Box2D const& domain, plint ndim);

}  // namespace plb

#endif  // COMPILE_2D
