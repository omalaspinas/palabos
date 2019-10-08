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

#include "atomicBlock/dataProcessorWrapper2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "atomicBlock/atomicBlockOperations2D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** BoxProcessing2D, general case *************************** */

void applyProcessingFunctional(BoxProcessingFunctional2D* functional,
                               Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks)
{
    executeDataProcessor( BoxProcessorGenerator2D(functional, domain),
                          atomicBlocks );
}

void integrateProcessingFunctional(BoxProcessingFunctional2D* functional,
                                   Box2D domain,
                                   std::vector<AtomicBlock2D*> atomicBlocks,
                                   plint level)
{
    addInternalProcessor( BoxProcessorGenerator2D(functional, domain),
                          atomicBlocks, level );
}


/* *************** DotProcessing, general case ***************************** */

void applyProcessingFunctional(DotProcessingFunctional2D* functional,
                               DotList2D const& dotList,
                               std::vector<AtomicBlock2D*> atomicBlocks)
{
    executeDataProcessor( DotProcessorGenerator2D(functional, dotList),
                          atomicBlocks );
}

void integrateProcessingFunctional(DotProcessingFunctional2D* functional,
                                   DotList2D const& dotList,
                                   std::vector<AtomicBlock2D*> atomicBlocks,
                                   plint level)
{
    addInternalProcessor( DotProcessorGenerator2D(functional, dotList),
                          atomicBlocks, level );
}

/* *************** BoundedBoxProcessing2D, general case *************************** */

void applyProcessingFunctional(BoundedBoxProcessingFunctional2D* functional,
                               Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks,
                               plint boundaryWidth )
{
    std::vector<BoxProcessorGenerator2D*> generators;
    functional -> getGenerators(domain, boundaryWidth, generators);
    delete functional;
    for (pluint iGen=0; iGen<generators.size(); ++iGen) {
        executeDataProcessor( *generators[iGen], atomicBlocks );
        delete generators[iGen];
    }
}

void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D* functional,
                                   Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks,
                                   plint boundaryWidth, plint level)
{
    std::vector<BoxProcessorGenerator2D*> generators;
    functional -> getGenerators(domain, boundaryWidth, generators);
    delete functional;
    for (pluint iGen=0; iGen<generators.size(); ++iGen) {
        addInternalProcessor( *generators[iGen], atomicBlocks, level );
        delete generators[iGen];
    }
}

}  // namespace plb
