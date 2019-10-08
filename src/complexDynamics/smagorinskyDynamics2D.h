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

/* Orestis Malaspinas designed some of the classes and concepts contained
 * in this file. */

#ifndef SMAGORINSKY_DYNAMICS_2D_H
#define SMAGORINSKY_DYNAMICS_2D_H

#include "complexDynamics/smagorinskyDynamics.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
void instantiateStaticSmagorinsky(BlockLattice2D<T,Descriptor>& lattice, Box2D domain);

template<typename T, template<typename U> class Descriptor>
void instantiateStaticSmagorinsky(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain);

template<typename T, template<typename U> class Descriptor, class SmagoFunction>
void instantiateStaticSmagorinsky(BlockLattice2D<T,Descriptor>& lattice, Box2D domain, SmagoFunction smagoFunction, T cSmago0);

template<typename T, template<typename U> class Descriptor, class SmagoFunction>
void instantiateStaticSmagorinsky(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain, SmagoFunction smagoFunction, T cSmago0);

} // namespace plb

#include "complexDynamics/smagorinskyGenerics2D.h"

#endif  // SMAGORINSKY_DYNAMICS_2D_H
