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

#ifndef PLB_DEBUG_H
#define PLB_DEBUG_H

#include <cassert>

#ifdef PLB_DEBUG

    #define PLB_ASSERT( COND )        assert( COND );
    #define PLB_PRECONDITION( COND )  assert( COND );
    #define PLB_POSTCONDITION( COND ) assert( COND );
    #define PLB_STATECHECK( A,B )     assert( (A) == (B) );

#else

    #define PLB_ASSERT( COND )
    #define PLB_PRECONDITION( COND )
    #define PLB_POSTCONDITION( COND )
    #define PLB_STATECHECK( A,B )

#endif  // PLB_DEBUG

namespace plb {

// Programmatically enable core dumps for POSIX systems.
// In a parallel program, this function must be called after plbInit.
void enableCoreDumps();

// Make stdout and stderr unbuffered for better debugging.
void unbufferOutputStdStreams();

}  // namespace plb

#endif  // PLB_DEBUG_H
