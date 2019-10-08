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

#include <dlfcn.h>
#include <mpi.h>
#include <iostream>
#include "utils.h"
char* Mpi::loadMpi(void){
        MPISIZE = 0; 
        MPIRANK = 0;
        //int MPISIZE = 0; 
        //char **MPIRANK = 0;
        //dlopen("libmpi.so.0", RTLD_LAZY);
        dlopen("libmpi.so.0", RTLD_NOW | RTLD_GLOBAL);
        MPI_Init(&MPISIZE, &MPIRANK);
        //int nprocs;
        //MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
        //std::cout << "nb procs :" << nprocs << std::endl;
        return dlerror();
}

int Mpi::finalizeMpi(void){
	dlopen("libmpi.so.0", RTLD_NOW | RTLD_GLOBAL);
	return MPI_Finalize();
}

int Mpi::getRankMpi(void){
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	return rank;
}
