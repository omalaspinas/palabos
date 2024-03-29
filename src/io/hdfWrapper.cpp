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

#include "io/hdfWrapper.h"

#include <cstdlib>

hsize_t openHDFfile(const char *file_name, MPI_Comm comm)
{
    hid_t plist_id;
    MPI_Info info = MPI_INFO_NULL;
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);
    hsize_t fid = H5Fopen(file_name, H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);
    return fid;
}

void closeHDFfile(hsize_t fid)
{
    H5Fclose(fid);
}

char *readStringHDF5(hsize_t fid)
{
    hid_t dset_id;
    herr_t status;
    dset_id = H5Dopen2(fid, "/meta_data", H5P_DEFAULT);
    hid_t dspace = H5Dget_space(dset_id);
    hsize_t dim, maxdims;
    H5Sget_simple_extent_dims(dspace, &dim, &maxdims);

    char *xm_read = (char *)malloc((int)maxdims + 1);
    status = H5Dread(dset_id, H5T_C_S1, H5S_ALL, H5S_ALL, H5P_DEFAULT, xm_read);
    if (status < 0) {
        printf("Reading string with hdf5 didn't work, error code: %d \n", status);
    }
    xm_read[maxdims] = 0;
    H5Dclose(dset_id);
    return xm_read;
}

void writeStringHDF5(hid_t file_id, std::string const &string, int mpi_rank)
{
    hid_t dataset_id, dataspace_id;
    hsize_t nb = string.size();
    herr_t status;
    /* Open an existing file. */
    dataspace_id = H5Screate_simple(1, &nb, NULL);
    dataset_id = H5Dcreate2(
        file_id, "/meta_data", H5T_C_S1, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (mpi_rank == 0) {
        status = H5Dwrite(dataset_id, H5T_C_S1, H5S_ALL, H5S_ALL, H5P_DEFAULT, string.c_str());
        if (status < 0) {
            printf("Writing string with hdf5 didn't work, error code: %d \n", status);
        }
    }
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
}

std::vector<std::vector<char>> readParallelHDF5(
    hsize_t fid, std::vector<plb::plint> const &my_block_id,
    std::vector<plb::plint> const &i_offset)
{
    hid_t dset_id, plist_id;    // file and dataset identifiers //
    hid_t filespace, memspace;  // file and memory dataspace identifiers //
    herr_t status;
    hsize_t data_piece_size;
    hsize_t offset;

    if (my_block_id.size() == 0) {
        data_piece_size = 0;
        offset = 0;
    } else if (my_block_id[0] > 0) {
        data_piece_size = i_offset[my_block_id.back()] - i_offset[my_block_id[0] - 1];
        offset = i_offset[my_block_id[0] - 1];
    } else {
        data_piece_size = i_offset[my_block_id.back()];
        offset = 0;
    }

    // allocate memory for the data to be read////
    char *data_total = (char *)malloc((int)data_piece_size);

    // create dataset
    dset_id = H5Dopen2(fid, "binary_blob", H5P_DEFAULT);

    filespace = H5Dget_space(dset_id);
    memspace = H5Screate_simple(1, &data_piece_size, NULL);

    // create pieces to be read in parallel//
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, NULL, &data_piece_size, NULL);

    // read data in parallel
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    status = H5Dread(dset_id, H5T_C_S1, memspace, filespace, plist_id, data_total);

    if (status < 0) {
        printf("hdf5 read didnt work, error code: %d \n", status);
        std::vector<std::vector<char>> empty;
        return empty;
    }

    // Close/release resources.
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);

    std::vector<std::vector<char>> data(my_block_id.size());

    int idx_acc = 0;
    for (int i = 0; i < (int)my_block_id.size(); i++) {
        int block_size;

        if (my_block_id[i] > 0) {
            block_size = i_offset[my_block_id[i]] - i_offset[my_block_id[i] - 1];
        } else {
            block_size = i_offset[0];
        }
        data[i] = std::vector<char>(block_size);
        for (int j = 0; j < block_size; j++) {
            data[i][j] = data_total[idx_acc + j];
        }
        idx_acc += block_size;
    }
    free(data_total);
    return data;
}

#endif
