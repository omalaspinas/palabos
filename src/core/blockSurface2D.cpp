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
 * Handling the surface area of a block -- implementation.
 */

#include "core/blockSurface2D.h"

namespace plb {

BlockSurface2D::BlockSurface2D(Box2D const& domain_, plint boundaryWidth)
    : domain(domain_),
      bw(boundaryWidth)
{ }

Box2D BlockSurface2D::bulk() const {
    return Box2D( domain.x0+bw,   domain.x1-bw,
                  domain.y0+bw,   domain.y1-bw    );
}

Box2D BlockSurface2D::edge0N() const {
    return Box2D( domain.x0,      domain.x0+bw-1,
                  domain.y0+bw,   domain.y1-bw    );
}

Box2D BlockSurface2D::edge0P() const {
    return Box2D( domain.x1-bw+1, domain.x1,
                  domain.y0+bw,   domain.y1-bw    );
}

Box2D BlockSurface2D::edge1N() const {
    return Box2D( domain.x0+bw,   domain.x1-bw,
                  domain.y0,      domain.y0+bw-1  );
}

Box2D BlockSurface2D::edge1P() const {
    return Box2D( domain.x0+bw,   domain.x1-bw,
                  domain.y1-bw+1, domain.y1       );
}


Box2D BlockSurface2D::cornerNN() const {
    return Box2D( domain.x0,      domain.x0+bw-1,
                  domain.y0,      domain.y0+bw-1  );
}

Box2D BlockSurface2D::cornerPN() const {
    return Box2D( domain.x1-bw+1, domain.x1,
                  domain.y0,      domain.y0+bw-1  );
}

Box2D BlockSurface2D::cornerNP() const {
    return Box2D( domain.x0,      domain.x0+bw-1,
                  domain.y1-bw+1, domain.y1       );
}

Box2D BlockSurface2D::cornerPP() const {
    return Box2D( domain.x1-bw+1, domain.x1,
                  domain.y1-bw+1, domain.y1       );
}

}  // namespace plb
