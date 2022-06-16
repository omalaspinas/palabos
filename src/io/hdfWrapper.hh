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

#ifndef HDF_WRAPPER_HH
#define HDF_WRAPPER_HH

#ifdef HDF5

#include <cstdlib>

#include "io/hdfWrapper.h"

template <typename T>
void writeParallelHDF5(
    const char *file_name, const char *data_set_path, std::vector<plb::plint> const &my_block_id,
    std::vector<plb::plint> const &i_offset, std::vector<std::vector<char> > &data,
    std::string const &metadata, int mpi_rank, MPI_Comm comm, bool opened)
{
    hid_t h5_type = H5T_C_S1, h5_element_size = sizeof(T);

    if (typeid(T) == typeid(float)) {
        h5_type = H5T_NATIVE_FLOAT;
    }
    if (typeid(T) == typeid(double)) {
        h5_type = H5T_NATIVE_DOUBLE;
    }

    hid_t file_id, dset_id, plist_id; /* file and dataset identifiers */
    hid_t filespace, memspace;        /* file and memory dataspace identifiers */
    herr_t status;
    hsize_t data_piece_size;
    hsize_t offset;
    hsize_t data_set_size = i_offset.back() / h5_element_size;

    if (my_block_id.size() == 0) {
        data_piece_size = 0;
        offset = 0;
    } else if (my_block_id[0] > 0) {
        data_piece_size = i_offset[my_block_id.back()] - i_offset[my_block_id[0] - 1];
        offset = i_offset[my_block_id[0] - 1] / h5_element_size;
    } else {
        data_piece_size = i_offset[my_block_id.back()];
        offset = 0;
    }
    char *data_total = (char *)malloc((int)data_piece_size);
    data_piece_size /= h5_element_size;

    int idx_acc = 0;

    for (int i = 0; i < (int)my_block_id.size(); i++) {
        int block_size;

        if (my_block_id[i] > 0) {
            block_size = i_offset[my_block_id[i]] - i_offset[my_block_id[i] - 1];
        } else {
            block_size = i_offset[0];
        }

        for (int j = 0; j < block_size; j++) {
            data_total[idx_acc + j] = data[i][j];
        }
        idx_acc += block_size;
    }

    MPI_Info info = MPI_INFO_NULL;
    // Create the file.
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);

    if (opened) {
        file_id = H5Fopen(file_name, H5F_ACC_RDWR, plist_id);
    } else {
        file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    }

    H5Pclose(plist_id);

    if (metadata[0]) {
        writeStringHDF5(file_id, metadata, mpi_rank);
    }

    // Create the dataset.
    filespace = H5Screate_simple(1, &data_set_size, NULL);
    memspace = H5Screate_simple(1, &data_piece_size, NULL);

    dset_id = H5Dcreate(
        file_id, data_set_path, h5_type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    // Create the pieces to be written in parallel.
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, NULL, &data_piece_size, NULL);

    // Write data in parallel
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    status = H5Dwrite(dset_id, h5_type, memspace, filespace, plist_id, data_total);

    if (status < 0) {
        printf("hdf5 write didnt work, error code: %d \n", status);
        return;
    }
    free(data_total);
    // Close/release resources.
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);

    H5Pclose(plist_id);
    H5Fclose(file_id);
}

#endif
#endif
