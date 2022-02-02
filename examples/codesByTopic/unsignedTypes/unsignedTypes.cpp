///////////////////////////////////////////////////////////////////////////////
/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 *
 * Contact for Palabos:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * Contact for npFEM:
 * Christos Kotsalos
 * kotsaloscv@gmail.com
 * Computer Science Department
 * University of Geneva
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
///////////////////////////////////////////////////////////////////////////////

/** \file
 *
 * This code demonstrates the use of unsigned data types.
 * These can be used to set flags using bitwise operators
 *
 * An minimal example is given with generic operations functionals
 **/

#include "palabos2D.h"
#include "palabos2D.hh"

using namespace plb;
typedef double T;
typedef unsigned int unsignedType;

// For example purposes, simple functions to set and check for flags are provided
// These can be generalized for any location
// with unsigned long, 64 flag positions are availibale
unsignedType setFlag1(unsignedType fieldValue)
{
    unsigned long loc = 0b1 << 1;
    return (fieldValue | loc);
}

unsignedType setFlag2(unsignedType fieldValue)
{
    unsigned long loc = 0b1 << 2;
    return (fieldValue | loc);
}

bool checkFlag1(unsignedType fieldValue)
{
    unsignedType loc = 0b1 << 1;
    return (fieldValue & loc) != 0u;
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    makeDirectory("tmp", false);

    // Create a scalar field with unsigned data type
    MultiScalarField2D<unsignedType> flagField(100, 100, 0);

    Box2D square1(40, 60, 40, 60);
    Box2D square2(20, 80, 20, 80);

    // set Flags in the specified regions
    apply(setFlag2, flagField, square2);
    apply(setFlag1, flagField, square1);

    // evaluate returns a unqiue_ptr ScalarField with the same data type as the provided ScalarField
    // thus it needs to be converted to boolean.
    std::unique_ptr<MultiScalarField2D<bool>> boolFlag1 =
        copyConvert<unsignedType, bool>(*evaluate(checkFlag1, flagField));

#ifdef PLB_REGRESSION
    pcout << flagField << std::endl;
    pcout << *boolFlag1 << std::endl;
#endif

    // Write images to show that only where the Flag1 was set, the evaluate returns true
    ImageWriter<unsignedType> imageWriter("leeloo.map");
    imageWriter.writeScaledGif("flagField", flagField, 1000, 1000);
    imageWriter.writeScaledGif("flag1", *copyConvert<bool, unsignedType>(*boolFlag1), 1000, 1000);

    // Write VTK to show that only where the Flag1 was set, the evaluate returns true
    VtkImageOutput2D<T> vtkOut("vtkOut", (T)1);
    vtkOut.writeData<unsignedType>(*copyConvert<unsignedType, T>(flagField), "flagField", (T)1);
    vtkOut.writeData<float>(*copyConvert<bool, T>(*boolFlag1), "flag1", (T)1);

    pcout << "Finished running unsignedTypes" << std::endl;
}
