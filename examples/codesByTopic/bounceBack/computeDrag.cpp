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

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "cylinder.h"
#include "cylinder.hh"
#include "palabos2D.h"
#include "palabos2D.hh"
#include "poiseuille.h"
#include "poiseuille.hh"

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9Descriptor

/// A functional, used to initialize a pressure boundary to constant density
template <typename T>
class ConstantDensity {
public:
    ConstantDensity(T density_) : density(density_) { }
    T operator()(plint, plint) const
    {
        return density;
    }

private:
    T density;
};

/// A functional, used to instantiate bounce-back nodes at the locations of the cylinder
void defineCylinderGeometry(
    MultiBlockLattice2D<T, DESCRIPTOR> &lattice, IncomprFlowParam<T> const &parameters,
    OnLatticeBoundaryCondition2D<T, DESCRIPTOR> &boundaryCondition, Array<plint, 2> forceIds)
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    Box2D everythingButOutlet(0, nx - 2, 0, ny - 1);
    Box2D outlet(nx - 1, nx - 1, 1, ny - 2);

    // Create Velocity boundary conditions everywhere
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, everythingButOutlet);
    // .. except on right boundary, where we prefer a fixed-pressure condition.
    boundaryCondition.setPressureConditionOnBlockBoundaries(lattice, outlet);

    setBoundaryVelocity(lattice, lattice.getBoundingBox(), PoiseuilleVelocity<T>(parameters));
    setBoundaryDensity(lattice, outlet, ConstantDensity<T>(1.));
    initializeAtEquilibrium(
        lattice, lattice.getBoundingBox(), PoiseuilleVelocityAndDensity<T, DESCRIPTOR>(parameters));

    int cx = nx / 4;
    int cy = ny / 2 + 2;
    int radius = cy / 4;
    // Instead of plain BounceBack, use the dynamics MomentumExchangeBounceBack,
    //   to compute the momentum exchange, and thus, the drag and the lift on
    //   the obstacle, locally during collision.
    defineDynamics(
        lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(cx, cy, radius),
        new MomentumExchangeBounceBack<T, DESCRIPTOR>(forceIds));
    initializeMomentumExchange(
        lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(cx, cy, radius));

    lattice.initialize();
}

void writeGifs(MultiBlockLattice2D<T, DESCRIPTOR> &lattice, plint iter)
{
    const plint imSize = 600;

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(
        createFileName("u", iter, 6), *computeVelocityNorm(lattice), imSize, imSize);
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");

    IncomprFlowParam<T> parameters(
        (T)1e-2,  // uMax
        (T)100.,  // Re
        64,       // N
        5.,       // lx
        1.        // ly
    );
    const T logT = (T)0.01;
#ifndef PLB_REGRESSION
    const T imSave = (T)1.;
    const T maxT = (T)10.1;
#else
    const T maxT = (T)0.5;
#endif

    writeLogFile(parameters, "Poiseuille flow");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice(
        parameters.getNx(), parameters.getNy(),
        new BGKdynamics<T, DESCRIPTOR>(parameters.getOmega()));
    lattice.initialize();

    // The drag and lift acting on the obstacle are computed with help of the
    //   internal statistics object of the lattice. For this purpose, they
    //   need to be registered first, as it is done in the following lines.
    Array<plint, 2> forceIds;
    forceIds[0] = lattice.internalStatSubscription().subscribeSum();
    forceIds[1] = lattice.internalStatSubscription().subscribeSum();

    OnLatticeBoundaryCondition2D<T, DESCRIPTOR> *boundaryCondition =
        createInterpBoundaryCondition2D<T, DESCRIPTOR>();
    // boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    defineCylinderGeometry(lattice, parameters, *boundaryCondition, forceIds);

    // Main loop over time iterations.
    for (plint iT = 0; iT * parameters.getDeltaT() < maxT; ++iT) {
#ifndef PLB_REGRESSION
        if (iT % parameters.nStep(imSave) == 0) {
            pcout << "Saving Gif ..." << endl;
            writeGifs(lattice, iT);
        }
#endif

        if (iT % parameters.nStep(logT) == 0) {
            pcout << "step " << iT << "; lattice time=" << lattice.getTimeCounter().getTime()
                  << "; t=" << iT * parameters.getDeltaT();
        }

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();

        if (iT % parameters.nStep(logT) == 0) {
            pcout << "; av energy=" << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                  << "; av rho=" << getStoredAverageDensity<T>(lattice)
                  << "; drag=" << lattice.getInternalStatistics().getSum(forceIds[0])
                  << "; lift=" << lattice.getInternalStatistics().getSum(forceIds[1]) << endl;
        }
    }

    delete boundaryCondition;
}
