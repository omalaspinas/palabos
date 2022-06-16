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

#ifndef HDF_WRAPPER_H
#define HDF_WRAPPER_H

#ifdef HDF5

#include <vector>

#include "core/globalDefs.h"
#include "hdf5.h"

hsize_t openHDFfile(const char *file_name, MPI_Comm comm);
void closeHDFfile(hsize_t fid);

char *readStringHDF5(hsize_t fid);
void writeStringHDF5(hid_t file_id, std::string const &string, int mpi_rank);

// QUESTION: int mpi_rank, MPI_Comm comm are unused. should we remove them?
// for this to work the file has to be already open (the file descriptor is a global var),
// this is an optimization it allows to read data and metadata without having to close and reopen
// the file
std::vector<std::vector<char>> readParallelHDF5(
    hsize_t fid, std::vector<plb::plint> const &my_block_id,
    std::vector<plb::plint> const &i_offset, [[maybe_unused]] int mpi_rank,
    [[maybe_unused]] MPI_Comm comm);

template <typename T>
void writeParallelHDF5(
    const char *file_name, const char *data_set_path, std::vector<plb::plint> const &my_block_id,
    std::vector<plb::plint> const &i_offset, std::vector<std::vector<char>> &data,
    std::string const &metadata, int mpi_rank, MPI_Comm comm, bool opened = false);

#endif
#endif
