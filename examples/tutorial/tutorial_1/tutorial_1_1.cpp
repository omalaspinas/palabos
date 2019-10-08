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

/* Code 1.1 in the Palabos tutorial
 */

#include "palabos2D.h"
#include "palabos2D.hh"
#include <iostream>
#include <iomanip>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR plb::descriptors::D2Q9Descriptor

// Initialize the lattice at zero velocity and constant density, except
//   for a slight density excess on a square sub-domain.
void defineInitialDensityAtCenter(MultiBlockLattice2D<T,DESCRIPTOR>& lattice)
{
    // The lattice is of size nx-by-ny
    const plint nx = lattice.getNx();
    const plint ny = lattice.getNy();

    // Create a Box2D which describes the location of cells with a slightly
    //   higher density.
    plint centralSquareRadius = nx/6;
    plint centerX = nx/3;
    plint centerY = ny/4;
    Box2D centralSquare (
            centerX - centralSquareRadius, centerX + centralSquareRadius,
            centerY - centralSquareRadius, centerY + centralSquareRadius );

    // All cells have initially density rho ...
    T rho0 = 1.;
    // .. except for those in the box "centralSquare" which have density
    //    rho+deltaRho
    T deltaRho = 1.e-4;
    Array<T,2> u0((T)0,(T)0);

    // Initialize constant density everywhere.
    initializeAtEquilibrium (
           lattice, lattice.getBoundingBox(), rho0, u0 );

    // And slightly higher density in the central box.
    initializeAtEquilibrium (
           lattice, centralSquare, rho0 + deltaRho, u0 );

    lattice.initialize();
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    const plint maxIter = 1000; // Iterate during 1000 steps.
    const plint nx = 600;       // Choice of lattice dimensions.
    const plint ny = 600;
    const T omega = 1.;        // Choice of the relaxation parameter

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
           nx, ny, new BGKdynamics<T,DESCRIPTOR>(omega) );

    lattice.periodicity().toggleAll(true); // Use periodic boundaries.

    defineInitialDensityAtCenter(lattice);

    // Main loop over time iterations.
    for (plint iT=0; iT<maxIter; ++iT) {
        if (iT%40==0) {  // Write an image every 40th time step.
            pcout << "Writing GIF file at iT=" << iT << endl;
            // Instantiate an image writer with the color map "leeloo".
            ImageWriter<T> imageWriter("leeloo");
            // Write a GIF file with colors rescaled to the range of values
            //   in the matrix
            imageWriter.writeScaledGif (
                    createFileName("u", iT, 6),
                    *computeVelocityNorm(lattice) );
        }
        // Execute lattice Boltzmann iteration.
        lattice.collideAndStream();
    }
}
