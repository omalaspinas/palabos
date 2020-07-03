if(CI)
###############################################################################

enable_language(CXX)

option(BUILD_HDF5 "compile with hdf5" ON)

if(BUILD_HDF5)
	set(HDF5_PREFER_PARALLEL ON)
	find_package(HDF5 REQUIRED )
	if(HDF5_FOUND)
		if(HDF5_IS_PARALLEL)
			message("HDF5 PARALLEL FOUND")
			add_definitions(-DHDF5)
		else()
			message("You need the parallel version of HDF5 to run this example, the code will compile and run, but without the HFD5 output!")
			option(BUILD_HDF5 "compile with hdf5" OFF)
		endif()
	else()
		message("You don’t have HDF5 library or cmake couldn’t find it, the code will compile and run, but without the HFD5 output!")
		option(BUILD_HDF5 "compile with hdf5" OFF)
	endif()
endif()

set(EXECUTABLE_NAME "aneurysm")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${HDF5_INCLUDE_DIRS})

###############################################################################
set(PROJECT_NAME "aneurysm_hdf5")
add_executable(${PROJECT_NAME} "./${EXECUTABLE_NAME}.cpp")

if(BUILD_HDF5)
	target_link_libraries(${PROJECT_NAME} palabos ${HDF5_LIBRARIES} )
else()
	target_link_libraries(${PROJECT_NAME} palabos )
endif()

if(ENABLE_MPI)
	target_link_libraries(${PROJECT_NAME} ${MPI_CXX_LIBRARIES})
endif()

###############################################################################

# This returns when you compile Palabos and all the examples (root CMakeLists.txt)
return()

###############################################################################
endif()


cmake_minimum_required(VERSION 3.0)

project(aneurysm)
enable_language(CXX)

option(BUILD_HDF5 "compile with hdf5" ON)

if(BUILD_HDF5)
	set(HDF5_PREFER_PARALLEL ON)
	find_package(HDF5 REQUIRED )
	if(HDF5_FOUND)
		if(HDF5_IS_PARALLEL)
			message("HDF5 PARALLEL FOUND")
			add_definitions(-DHDF5)
		else()
			message("You need the parallel version of HDF5 to run this example, the code will compile and run, but without the HFD5 output!")
			option(BUILD_HDF5 "compile with hdf5" OFF)
		endif()
	else()
		message("You don’t have HDF5 library or cmake couldn’t find it, the code will compile and run, but without the HFD5 output!")
		option(BUILD_HDF5 "compile with hdf5" OFF)
	endif()
endif()

set(EXECUTABLE_NAME "aneurysm")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "./")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
    #set(CMAKE_CONFIGURATION_TYPES "Debug;Release")
endif()
message("Generated with config types: ${CMAKE_CONFIGURATION_TYPES}")
message(${CMAKE_BUILD_TYPE})

# Compiler flags
# Append flags: set(CMAKE_XXX_FLAGS "${CMAKE_XXX_FLAGS} ...")
if(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
    message("GCC.")
    set(CMAKE_CXX_FLAGS "-std=c++11  -Wall -Wnon-virtual-dtor -Wno-unused-result")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -DPLB_DEBUG -O0")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
    message("Clang.")
    set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wnon-virtual-dtor")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -DPLB_DEBUG -O0")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    message("MSVC.")
    set(CMAKE_CXX_FLAGS_RELEASE "/Ox /Ot /GS- /GL /DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "/DPLB_DEBUG")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/LTCG /INCREMENTAL:NO /OPT:REF")
else()
    message( FATAL_ERROR "CXX compiler not recognized. CMake will quit." )
endif()

option(ENABLE_MPI "Enable MPI" ON)
if(ENABLE_MPI)
    message("Enabling MPI")
    find_package(MPI REQUIRED)
    if(MPI_CXX_FOUND)
        #set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")
        include_directories(${MPI_CXX_INCLUDE_PATH})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
        add_definitions(-DPLB_MPI_PARALLEL)
    endif()
endif()

if(WIN32)
    option(ENABLE_POSIX "Enable POSIX" OFF)
else()
    option(ENABLE_POSIX "Enable POSIX" ON)
endif()

if(ENABLE_POSIX)
    message("Enabling POSIX")
    add_definitions(-DPLB_USE_POSIX)
endif()

if(APPLE)
    add_definitions(-DPLB_MAC_OS_X)
endif()

if(WIN32 OR CYGWIN)
    add_definitions(-DPLB_WINDOWS)
endif()

###############################################################################
# Palabos Library
###############################################################################
include_directories(${HDF5_INCLUDE_DIRS})
include_directories("../../../src")
include_directories("../../../externalLibraries")
include_directories("../../../externalLibraries/Eigen3")

file(GLOB_RECURSE PALABOS_SRC "../../../src/*.cpp")
file(GLOB_RECURSE EXT_SRC "../../../externalLibraries/tinyxml/*.cpp")

add_library(palabos STATIC ${PALABOS_SRC} ${EXT_SRC})

###############################################################################

add_executable(${EXECUTABLE_NAME} "./${EXECUTABLE_NAME}.cpp")

# Link with the following libraries
if(BUILD_HDF5)
	target_link_libraries(${EXECUTABLE_NAME} palabos ${HDF5_LIBRARIES} )
else()
	target_link_libraries(${EXECUTABLE_NAME} palabos )
endif()
#
if(ENABLE_MPI)
    target_link_libraries(${EXECUTABLE_NAME} ${MPI_CXX_LIBRARIES})
endif()
