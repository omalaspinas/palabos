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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- template instantiation.
 */

#ifdef COMPILE_3D

#include "dataProcessors/dataAnalysisWrapper3D.h"
#include "dataProcessors/dataAnalysisWrapper3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {

// Atomic-block

template FLOAT_T computeAverageDensity<FLOAT_T, descriptors::DESCRIPTOR_3D> (
                    BlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D>& lattice, Box3D domain );
template FLOAT_T computeAverageDensity<FLOAT_T, descriptors::DESCRIPTOR_3D> (
                    BlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D>& lattice );

template FLOAT_T computeAverageRhoBar<FLOAT_T, descriptors::DESCRIPTOR_3D> (
                    BlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D>& lattice, Box3D domain );
template FLOAT_T computeAverageRhoBar<FLOAT_T, descriptors::DESCRIPTOR_3D> (
                    BlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D>& lattice );

template FLOAT_T computeAverageEnergy<FLOAT_T, descriptors::DESCRIPTOR_3D> (
                    BlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D>& lattice, Box3D domain );
template FLOAT_T computeAverageEnergy<FLOAT_T, descriptors::DESCRIPTOR_3D> (
                    BlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D>& lattice );

// Multi-block

template FLOAT_T computeAverageDensity<FLOAT_T, descriptors::DESCRIPTOR_3D> (
                    MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D>& lattice, Box3D domain );
template FLOAT_T computeAverageDensity<FLOAT_T, descriptors::DESCRIPTOR_3D> (
                    MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D>& lattice );

template FLOAT_T computeAverageRhoBar<FLOAT_T, descriptors::DESCRIPTOR_3D> (
                    MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D>& lattice, Box3D domain );
template FLOAT_T computeAverageRhoBar<FLOAT_T, descriptors::DESCRIPTOR_3D> (
                    MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D>& lattice );

template FLOAT_T computeAverageEnergy<FLOAT_T, descriptors::DESCRIPTOR_3D> (
                    MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D>& lattice, Box3D domain );
template FLOAT_T computeAverageEnergy<FLOAT_T, descriptors::DESCRIPTOR_3D> (
                    MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D>& lattice );

template
void copyPopulations<FLOAT_T,descriptors::DESCRIPTOR_3D> (
        BlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& latticeFrom,
        BlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& latticeTo, Box3D domain );

template
void copyAll<FLOAT_T,descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& latticeFrom,
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& latticeTo, Box3D domain );

template
void copyRegenerate<FLOAT_T,descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& latticeFrom,
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& latticeTo, Box3D domain );
}  // namespace plb

#endif  // COMPILE_3D
