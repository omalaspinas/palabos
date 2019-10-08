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
 * Helper functions for domain initialization -- header file.
 */
#ifndef CO_PROCESSOR_FUNCTIONAL_3D_H
#define CO_PROCESSOR_FUNCTIONAL_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "atomicBlock/atomicContainerBlock3D.h"

namespace plb {

template<typename T>
struct PureDynamics : public ContainerBlockData
{
    PureDynamics() : isPure(true), omega() { }
    PureDynamics(bool isPure_, T omega_) :
        isPure(isPure_), omega(omega_)
    { }
    virtual PureDynamics<T>* clone() const {
        return new PureDynamics<T>(*this);
    }
    bool isPure;
    T omega;
};

template<typename T, template<typename U> class Descriptor>
class IdentifyPureDynamics3D : public BoxProcessingFunctional3D
{
public:
    IdentifyPureDynamics3D(plint dynamicsId_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual IdentifyPureDynamics3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    pluint getMaxChainLength() const;
private:
    plint dynamicsId;
};

}  // namespace plb

#endif  // CO_PROCESSOR_FUNCTIONAL_3D_H

