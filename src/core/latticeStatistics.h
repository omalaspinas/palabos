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


#ifndef LATTICE_STATISTICS_H
#define LATTICE_STATISTICS_H

#include "core/globalDefs.h"
#include "core/blockStatistics.h"

namespace plb {

namespace LatticeStatistics {
    enum { avRhoBar=0, avUSqr=1, maxUSqr=0 };
}

inline void gatherStatistics(BlockStatistics& statistics, double rhoBar, double uSqr) {
    statistics.gatherAverage(LatticeStatistics::avRhoBar, rhoBar);
    statistics.gatherAverage(LatticeStatistics::avUSqr, uSqr);
    statistics.gatherMax(LatticeStatistics::maxUSqr, uSqr);
    statistics.incrementStats();
}

} // namespace plb

#endif
