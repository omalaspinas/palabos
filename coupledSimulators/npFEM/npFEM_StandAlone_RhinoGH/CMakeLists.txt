cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(npFEM_RhinoGH)

set(DEFAULT ON CACHE INTERNAL "Default options")

enable_language(CXX)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
    #set(CMAKE_CONFIGURATION_TYPES "Debug;Release")
endif()
message("Generated with config types: ${CMAKE_CONFIGURATION_TYPES}")
message(${CMAKE_BUILD_TYPE})

# Compiler flags
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    message("MSVC.")
    add_definitions(-DSHAPEOP_MSVC)
    # For the dll environment
    add_definitions(-DSHAPEOP_EXPORT)
    set(CMAKE_CXX_FLAGS "/bigobj /EHsc")
    set(CMAKE_CXX_FLAGS_RELEASE "/MD /Ox /Ot /fp:fast /GS- /GL /DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "/MTd")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/LTCG /INCREMENTAL:NO /OPT:REF")
else()
    message(FATAL_ERROR "For stand-alone dev, use only MSVC & Grasshopper.")
endif()

option(ENABLE_OPENMP "Enabling OpenMP" ${DEFAULT})
if(ENABLE_OPENMP)
    message("Enabling OpenMP")
    find_package(OpenMP QUIET)
    if(OPENMP_FOUND)
        message("OpenMP found. OpenMP activated in release.")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${OpenMP_EXE_LINKER_FLAGS}")
        add_definitions(-DSHAPEOP_OPENMP)
    else()
        message("OpenMP not found.")
    endif()
endif()

###############################################################################
# Palabos Library
###############################################################################

file(TO_CMAKE_PATH $ENV{PALABOS_ROOT} PALABOS_ROOT)

###############################################################################
# npFEM Library
###############################################################################

# Development of npFEM (stand-alone/SA version)
add_definitions(-DNPFEM_SA)

# Find Eigen
include_directories("${PALABOS_ROOT}/externalLibraries")
include_directories("${PALABOS_ROOT}/externalLibraries/Eigen3")

# Find NanoFlann
include_directories("${PALABOS_ROOT}/externalLibraries/nanoflann")

include_directories("${PALABOS_ROOT}/coupledSimulators/npFEM")
include_directories("${PALABOS_ROOT}/coupledSimulators/npFEM/npFEM_StandAlone_RhinoGH/api")

file(GLOB NPFEM_H "${PALABOS_ROOT}/coupledSimulators/npFEM/*.h")
file(GLOB NPFEM_CPP "${PALABOS_ROOT}/coupledSimulators/npFEM/*.cpp")

file(GLOB API_HEADERS "${PALABOS_ROOT}/coupledSimulators/npFEM/npFEM_StandAlone_RhinoGH/api/*.h")
file(GLOB API_SOURCES "${PALABOS_ROOT}/coupledSimulators/npFEM/npFEM_StandAlone_RhinoGH/api/*.cpp")

# By adding the NPFEM_H, we manage to show the headers in a separate folder 
# in the project tree (quite handy for development)

# Build NPFEM Library (For stand-alone go with dynamic lib)
set(HEADERS ${NPFEM_H} ${API_HEADERS})
set(SOURCES ${NPFEM_CPP} ${API_SOURCES})
add_library(npFEM SHARED ${SOURCES} ${HEADERS})

###############################################################################