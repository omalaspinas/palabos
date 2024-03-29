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

/* Code 1.2 in the Palabos tutorial
 */

#include <iomanip>
#include <iostream>

#include "palabos2D.h"
#include "palabos2D.hh"

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR plb::descriptors::D2Q9Descriptor

const plint maxIter = 1000;  // Iterate during 1000 steps.
const plint nx = 600;        // Choice of lattice dimensions.
const plint ny = 600;
const T omega = 1.;  // Choice of the relaxation parameter

T rho0 = 1.;  // All cells have initially density rho ...
// .. except for those inside the disk which have density
//    rho+deltaRho
T deltaRho = 1.e-4;
Array<T, 2> u0(0, 0);

void initializeConstRho(plint, plint, T &rho, Array<T, 2> &u)
{
    u = u0;
    rho = rho0 + deltaRho;
}

void initializeRhoOnDisk(plint iX, plint iY, T &rho, Array<T, 2> &u)
{
    plint radius = nx / 6;
    plint centerX = nx / 3;
    plint centerY = ny / 4;
    u = u0;
    if ((iX - centerX) * (iX - centerX) + (iY - centerY) * (iY - centerY) < radius * radius) {
        rho = rho0 + deltaRho;
    } else {
        rho = rho0;
    }
}

// Initialize the lattice at zero velocity and constant density, except
//   for a slight density excess on a circular sub-domain.
void defineInitialDensityAtCenter(MultiBlockLattice2D<T, DESCRIPTOR> &lattice)
{
    // Initialize constant density everywhere.
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), rho0, u0);

    // And slightly higher density in the central box.
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), initializeRhoOnDisk);

    lattice.initialize();
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice(nx, ny, new BGKdynamics<T, DESCRIPTOR>(omega));

    lattice.periodicity().toggleAll(true);  // Set periodic boundaries.

    defineInitialDensityAtCenter(lattice);

    // Main loop over time iterations.
    for (plint iT = 0; iT < maxIter; ++iT) {
        if (iT % 40 == 0) {  // Write an image every 40th time step.
            pcout << "Writing GIF file at iT=" << iT << endl;
            // Instantiate an image writer with the color map "leeloo".
            ImageWriter<T> imageWriter("leeloo");
            // Write a GIF file with colors rescaled to the range of values
            //   in the matrix
            imageWriter.writeScaledGif(createFileName("u", iT, 6), *computeVelocityNorm(lattice));
        }
        // Execute lattice Boltzmann iteration.
        lattice.collideAndStream();
    }
}
