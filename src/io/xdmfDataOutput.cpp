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

#ifdef HDF5

#include "io/xdmfDataOutput.h"

namespace plb {

ParallelXdmfDataWriter3D::ParallelXdmfDataWriter3D(const std::string &fname) :
    xdmf_fname(FileName(fname + ".xdmf").defaultPath(global::directories().getOutputDir())),
    h5_fname(FileName(fname + ".h5").defaultPath(global::directories().getOutputDir()))
{
    if (global::mpi().isMainProcessor()) {
        fhandle = new std::ofstream(xdmf_fname.c_str());
        if (!(*fhandle)) {
            std::cerr << "could not open file " << fname << "\n";
            return;
        }
        (*fhandle) << "<?xml version = \"1.0\"?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        (*fhandle) << "\t<Xdmf Version=\"2.0\">\n";
        (*fhandle) << "\t\t<Domain>\n";
        (*fhandle).flush();
    }
}

ParallelXdmfDataWriter3D::~ParallelXdmfDataWriter3D()
{
    if (global::mpi().isMainProcessor()) {
        (*fhandle) << "\t\t</Domain>\n";
        (*fhandle) << "\t</Xdmf>\n";
        (*fhandle).close();
        delete fhandle;
    }
}

}  // namespace plb

#endif
