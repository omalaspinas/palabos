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

#include "boundaryCondition/neumannCondition2D.h"
#include "boundaryCondition/neumannCondition2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"

namespace plb {

    template class CopyUnknownPopulationsFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D,  0, -1>;
    template class CopyUnknownPopulationsFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D,  0, +1>;
    template class CopyUnknownPopulationsFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D,  1, -1>;
    template class CopyUnknownPopulationsFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D,  1, +1>;

    template class CopyAllPopulationsFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D,  0,-1>;
    template class CopyAllPopulationsFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D,  0,+1>;
    template class CopyAllPopulationsFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, -1, 0>;
    template class CopyAllPopulationsFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, -1,-1>;
    template class CopyAllPopulationsFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, -1,+1>;
    template class CopyAllPopulationsFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, +1, 0>;
    template class CopyAllPopulationsFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, +1,-1>;
    template class CopyAllPopulationsFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, +1,+1>;

    template class CopyVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D,  0,-1>;
    template class CopyVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D,  0,+1>;
    template class CopyVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, -1, 0>;
    template class CopyVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, -1,-1>;
    template class CopyVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, -1,+1>;
    template class CopyVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, +1, 0>;
    template class CopyVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, +1,-1>;
    template class CopyVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, +1,+1>;

    template class CopyTangentialVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D,  0,-1>;
    template class CopyTangentialVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D,  0,+1>;
    template class CopyTangentialVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, -1, 0>;
    template class CopyTangentialVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, -1,-1>;
    template class CopyTangentialVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, -1,+1>;
    template class CopyTangentialVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, +1, 0>;
    template class CopyTangentialVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, +1,-1>;
    template class CopyTangentialVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, +1,+1>;

    template class CopyNormalVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D,  0,-1>;
    template class CopyNormalVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D,  0,+1>;
    template class CopyNormalVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, -1, 0>;
    template class CopyNormalVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, -1,-1>;
    template class CopyNormalVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, -1,+1>;
    template class CopyNormalVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, +1, 0>;
    template class CopyNormalVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, +1,-1>;
    template class CopyNormalVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, +1,+1>;

    template class CopyDensityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D,  0,-1>;
    template class CopyDensityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D,  0,+1>;
    template class CopyDensityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, -1, 0>;
    template class CopyDensityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, -1,-1>;
    template class CopyDensityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, -1,+1>;
    template class CopyDensityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, +1, 0>;
    template class CopyDensityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, +1,-1>;
    template class CopyDensityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D, +1,+1>;

}

#endif  // COMPILE_2D
