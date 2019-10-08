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

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;

typedef double T;

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./");
    global::IOpolicy().activateParallelIO(false);

    if (argc <4) {
        pcout << "Wrong parameters; the syntax is:" << std::endl
              << (std::string)global::argv(0)
              << " [FLT | DBL | LDBL | INF] outputSTL.stl inputSTL1.stl inputSTL2.stl ..." << std::endl;
        exit(-1);
    }

    std::string precisionStr;
    std::string outFileName;
    std::vector<std::string> inFileNames;
    try {
        global::argv(1).read(precisionStr);
        global::argv(2).read(outFileName);
        for (plint i=3; i<global::argc(); ++i) {
            std::string fileName;
            global::argv(i).read(fileName);
            inFileNames.push_back(fileName);
        }
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters." << std::endl;
        exit(-1);
    }

    Precision precision;
    if (precisionStr == "FLT") {
        precision = FLT;
    } else if (precisionStr == "DBL") {
        precision = DBL;
    } else if (precisionStr == "LDBL") {
        precision = LDBL;
    } else if (precisionStr == "INF") {
        precision = INF;
    } else {
        pcout << "Wrong precision command-line argument." << std::endl;
        exit(-1);
    }

    std::vector<TriangleSet<T>*> inSets(inFileNames.size());
    for (pluint i=0; i<inFileNames.size(); ++i) {
        try {
            inSets[i] = new TriangleSet<T>(inFileNames[i], precision);
        }
        catch (PlbIOException& exception) {
            pcout << "ERROR, could not read STL file " << inFileNames[i]
                  << ": " << exception.what() << std::endl;
            exit(-1);
        }
    }
    TriangleSet<T> outSet;
    outSet.merge(inSets);
    outSet.writeAsciiSTL(outFileName);

    for (pluint i=0; i<inSets.size(); ++i) {
        delete inSets[i];
    }

    return 0;
}

