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
 * Groups all the include files for the 3D multiBlock.
 */

#include "multiBlock/multiBlock3D.h"
#include "multiBlock/coupling3D.h"
#include "multiBlock/group3D.h"
#include "multiBlock/multiContainerBlock3D.h"
#include "multiBlock/multiBlockManagement3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"
#include "multiBlock/serialMultiBlockLattice3D.h"
#include "multiBlock/serialMultiDataField3D.h"
#include "multiBlock/serialBlockCommunicator3D.h"
#include "multiBlock/staticRepartitions3D.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.h"
#include "multiBlock/multiBlockOperations3D.h"
#include "multiBlock/multiBlockSerializer3D.h"
#include "multiBlock/combinedStatistics.h"
#include "multiBlock/multiBlockInfo3D.h"
#include "multiBlock/localMultiBlockInfo3D.h"
#include "multiBlock/nonLocalTransfer3D.h"
#include "multiBlock/nonLocalTransfer2D_3D.h"
#include "multiBlock/multiBlockGenerator3D.h"

