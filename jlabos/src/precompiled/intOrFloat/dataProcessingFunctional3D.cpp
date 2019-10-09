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

#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/dataProcessingFunctional3D.hh"

namespace plb {


/* *************** Boxed Data Processor functionals ****************** */

template class BoxProcessingFunctional3D_N<PRECOMP_T>;
template class BoxProcessingFunctional3D_S<PRECOMP_T>;
template class MaskedBoxProcessingFunctional3D_N<PRECOMP_T>;
template class BoxProcessingFunctional3D_NN<PRECOMP_T,int>;
template class BoxProcessingFunctional3D_NN<PRECOMP_T,float>;
template class BoxProcessingFunctional3D_NN<PRECOMP_T,double>;
template class BoxProcessingFunctional3D_SS<PRECOMP_T,int>;
template class BoxProcessingFunctional3D_SS<PRECOMP_T,float>;
template class BoxProcessingFunctional3D_SS<PRECOMP_T,double>;
template class MaskedBoxProcessingFunctional3D_NN<PRECOMP_T,float>;
template class MaskedBoxProcessingFunctional3D_NN<PRECOMP_T,double>;
template class MaskedBoxProcessingFunctional3D_NN<PRECOMP_T,int>;
template class NTensorFieldBoxProcessingFunctional3D<PRECOMP_T>;
template class ScalarFieldBoxProcessingFunctional3D<PRECOMP_T>;
template class MaskedNTensorFieldBoxProcessingFunctional3D<PRECOMP_T>;

/* *************** Bounded Boxed Data Processor functionals ****************** */

template class BoundedBoxProcessingFunctional3D_N<PRECOMP_T>;
template class BoundedMaskedBoxProcessingFunctional3D_N<PRECOMP_T>;
template class BoundedBoxProcessingFunctional3D_NN<PRECOMP_T,PRECOMP_T>;
template class BoundedMaskedBoxProcessingFunctional3D_NN<PRECOMP_T,PRECOMP_T>;
template class BoundedNTensorFieldBoxProcessingFunctional3D<PRECOMP_T>;

}  // namespace plb
