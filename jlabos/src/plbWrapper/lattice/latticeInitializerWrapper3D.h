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
#ifndef SWIG_LATTICE_INITIALIZER_WRAPPER_3D_H
#define SWIG_LATTICE_INITIALIZER_WRAPPER_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "core/dynamics.h"
#include "plbWrapper/lattice/latticeInitializerFunctional3D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
void pypalDefineDynamics( MultiBlockLattice3D<T,Descriptor>& lattice,
                          Box3D domain, Dynamics<T,Descriptor>* dynamics );

template<typename T, template<typename U> class Descriptor>
void maskedDefineDynamics( MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
                           Box3D domain, Dynamics<T,Descriptor>* dynamics );

template<typename T, template<class U> class Descriptor>
void setBoundaryVelocity( MultiBlockLattice3D<T,Descriptor>& lattice,
                          T* velocity, int numDimIs2, Box3D domain );

template<typename T, template<class U> class Descriptor>
void setBoundaryVelocity( MultiBlockLattice3D<T,Descriptor>& lattice,
                          MultiNTensorField3D<T>& velocity,
                          Box3D domain );

template<typename T, template<class U> class Descriptor>
void maskedSetBoundaryVelocity( MultiBlockLattice3D<T,Descriptor>& lattice,
                                MultiNTensorField3D<T>& velocity,
                                MultiNTensorField3D<int>& mask,
                                Box3D domain );

template<typename T, template<class U> class Descriptor>
void initializeAtEquilibrium( MultiBlockLattice3D<T,Descriptor>& lattice,
                              T density, T* velocity, int numDimIs2, Box3D domain );

template<typename T, template<class U> class Descriptor>
void initializeAtEquilibrium( MultiBlockLattice3D<T,Descriptor>& lattice,
                              MultiNTensorField3D<T>& density,
                              MultiNTensorField3D<T>& velocity,
                              Box3D domain );

template<typename T, template<class U> class Descriptor>
void maskedInitializeAtEquilibrium( MultiBlockLattice3D<T,Descriptor>& lattice,
                                    MultiNTensorField3D<T>& density,
                                    MultiNTensorField3D<T>& velocity,
                                    MultiNTensorField3D<int>& mask,
                                    Box3D domain );

template<typename T, template<class U> class Descriptor>
void setExternalVector( MultiBlockLattice3D<T,Descriptor>& lattice,
                        int vectorStartsAt, T* externalVector, int numDimIs2, Box3D domain );

template<typename T, template<class U> class Descriptor>
void setPopulations( MultiBlockLattice3D<T,Descriptor>& lattice,
                     T* populations, int numDimIsQ, Box3D domain );

template<typename T, template<class U> class Descriptor>
void setPopulations( MultiBlockLattice3D<T,Descriptor>& lattice,
                     MultiNTensorField3D<T>& populations,
                     Box3D domain );

template<typename T, template<class U> class Descriptor>
void maskedSetPopulations( MultiBlockLattice3D<T,Descriptor>& lattice,
                           T* populations, int numDimIsQ,
                           MultiNTensorField3D<int>& mask, Box3D domain );

template<typename T, template<class U> class Descriptor>
void maskedSetPopulations( MultiBlockLattice3D<T,Descriptor>& lattice,
                           MultiNTensorField3D<T>& populations,
                           MultiNTensorField3D<int>& mask,
                           Box3D domain );

}  // namespace plb

#endif  // SWIG_LATTICE_INITIALIZER_WRAPPER_3D_H
