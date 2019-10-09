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

#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace plb;
using namespace std;


template<typename T>
void plbFileToVtk( std::string fName, std::string identifier, VtkImageOutput3D<T>& vtkOut)
{
    parallelIO::SavedFullMultiBlockSerializer3D *serializer
        = new parallelIO::SavedFullMultiBlockSerializer3D(fName);
    Box3D bbox = serializer->getBoundingBox();
    std::string convertType = serializer->dataType();
    pcout << "Adding the field \"" << identifier << "\" from file " << fName
          << " to the VTK file, with type \"" << convertType << "\""<< std::endl;
    if (convertType=="double") {
        vtkOut.template writeData<double> (
                bbox.getNx(), bbox.getNy(), bbox.getNz(), 
                serializer->getCellDim(), serializer, identifier );
    }
    else if (convertType=="float") {
        vtkOut.template writeData<float> (
                bbox.getNx(), bbox.getNy(), bbox.getNz(), 
                serializer->getCellDim(), serializer, identifier );
    }
    else if (convertType=="int") {
        vtkOut.template writeData<int> (
                bbox.getNx(), bbox.getNy(), bbox.getNz(), 
                serializer->getCellDim(), serializer, identifier );
    }
    else {
        plbIOError("Cannot convert to type "+convertType);
    }
}

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);

    global::directories().setOutputDir("./");
    global::directories().setInputDir("./");

    plint numArgs(global::argc());

    std::vector<std::string> fNames;
    try {
        if (numArgs<3) {
            throw(PlbIOException("Too few arguments"));
        }

        for (plint iArg=1; iArg<numArgs; ++iArg) {
            fNames.push_back(global::argv(iArg));
        }
    }
    catch(PlbIOException const& exception) {
        pcout << exception.what() << std::endl;
        pcout << "Syntax: " << (std::string)global::argv(0)
              << " plb_file_name1 [plb_file_name2, plb_file_name3, ...] output_name"
              << std::endl;
        return -1;
    }

    try {
        double dx = 1.;
        VtkImageOutput3D<double> vtkOut(FileName(fNames.back()).getName(), dx);
        for (plint iFile=0; iFile<(plint)fNames.size()-1; ++iFile) {
            plbFileToVtk(fNames[iFile], FileName(fNames[iFile]).getName(), vtkOut);
        }
    }
    catch(PlbIOException const& exception) {
        pcout << exception.what() << std::endl;
    }
}
