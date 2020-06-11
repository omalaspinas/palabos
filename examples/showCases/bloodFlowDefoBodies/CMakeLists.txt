if(CI)
###############################################################################

enable_language(CXX)

set(EXECUTABLE_NAME "bloodFlowDefoBodies")
set(EXECUTABLE_NAME_GPU "bloodFlowDefoBodies_gpu")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    message("MSVC.")
	add_definitions(/DMSVC)
    add_definitions(/DSHAPEOP_MSVC)
    add_definitions(/DSHAPEOP_EXPORT)
endif()

###############################################################################
# npFEM Library (Compile npFEM separately)
###############################################################################

include_directories("../../../coupledSimulators")
include_directories("../../../coupledSimulators/npFEM")

file(GLOB NPFEM_H "../../../coupledSimulators/npFEM/*.h")
file(GLOB NPFEM_CPP "../../../coupledSimulators/npFEM/*.cpp")

add_library(npFEM STATIC ${NPFEM_CPP} ${NPFEM_H})

###############################################################################

file(GLOB DEPS_H "./*.h")
file(GLOB PARAMS_1 "./*.xml")
file(GLOB PARAMS_2 "./*.pos")

add_executable(${EXECUTABLE_NAME} "./${EXECUTABLE_NAME}.cpp" ${DEPS_H} ${PARAMS_1} ${PARAMS_2})
target_link_libraries(${EXECUTABLE_NAME} palabos npFEM)
if(ENABLE_MPI)
	target_link_libraries(${EXECUTABLE_NAME} ${MPI_CXX_LIBRARIES})
endif()

###############################################################################

# This returns when you compile Palabos and all the examples (root CMakeLists.txt)
return()

###############################################################################
endif()

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(bloodFlowDefoBodies)

enable_language(CXX)

option(ENABLE_CUDA "Enabling CUDA" OFF)
if(ENABLE_CUDA)
  enable_language(CUDA)
endif()

set(EXECUTABLE_NAME "bloodFlowDefoBodies")
set(EXECUTABLE_NAME_GPU "bloodFlowDefoBodies_gpu")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")

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
    if(CYGWIN)
        set(CMAKE_CXX_FLAGS "-Wa,-mbig-obj -std=c++11 -Wall -Wnon-virtual-dtor")
    else()
        set(CMAKE_CXX_FLAGS "-std=c++11 -w")
        #set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wnon-virtual-dtor")
    endif()
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -DPLB_DEBUG -O0")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
    message("Clang.")
    set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wnon-virtual-dtor")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -DPLB_DEBUG -O0")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    message("MSVC.")
	add_definitions(/DMSVC)
    add_definitions(/DSHAPEOP_MSVC)
    add_definitions(/DSHAPEOP_EXPORT)
    set(CMAKE_CXX_FLAGS "/bigobj /EHsc")
    set(CMAKE_CXX_FLAGS_RELEASE "/MD /Ox /Ot /fp:fast /GS- /GL /DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "/MTd /DPLB_DEBUG")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/LTCG /INCREMENTAL:NO /OPT:REF")
else()
    message( FATAL_ERROR "CXX compiler not recognized. CMake will quit." )
endif()

option(ENABLE_MPI "Enable MPI" ON)
if(ENABLE_MPI)
    message("Enabling MPI")
    find_package(MPI REQUIRED)
    if(MPI_CXX_FOUND)
        include_directories(${MPI_CXX_INCLUDE_DIRS})
        include_directories(SYSTEM ${MPI_INCLUDE_PATH})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
        add_definitions(-DPLB_MPI_PARALLEL)
    endif()
endif()

option(ENABLE_OPENMP "Enabling OpenMP" OFF)
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

option(ENABLE_LOGS "Enabling Log Files" OFF)
if (ENABLE_LOGS)
    message("Attention: Log files of huge size are generated.")
    add_definitions(-DENABLE_LOGS)
endif()

if(ENABLE_CUDA)
    message("Enabling CUDA")

    find_package(CUDA REQUIRED)
    if(NOT CUDA_FOUND)
        message(FATAL_ERROR "CUDA not found.")
    else()
        message("CUDA found.")
    endif()

    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_60")
    if(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
        set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG --use_fast_math -lineinfo -g")
        set(CMAKE_CUDA_FLAGS_DEBUG "-O0 -G -g")
    endif()
endif()

###############################################################################
# Palabos Library
###############################################################################

#file(TO_CMAKE_PATH $ENV{PALABOS_ROOT} PALABOS_ROOT)

include_directories("../../../src")
include_directories("../../../externalLibraries")
include_directories("../../../externalLibraries/Eigen3")
add_definitions(-DPLB_USE_EIGEN)

file(GLOB_RECURSE PALABOS_SRC "../../../src/*.cpp")
file(GLOB_RECURSE EXT_SRC     "../../../externalLibraries/tinyxml/*.cpp")

add_library(palabos STATIC ${PALABOS_SRC} ${EXT_SRC})

###############################################################################
# npFEM Library (Compile npFEM separately)
###############################################################################

include_directories("../../../coupledSimulators")
include_directories("../../../coupledSimulators/npFEM")
if(ENABLE_CUDA)
    include_directories("${CUDA_TOOLKIT_ROOT_DIR}")
    include_directories("${CUDA_TOOLKIT_ROOT_DIR}/include/")
    include_directories("../../../coupledSimulators/npFEM/src_GPU")
endif()

file(GLOB NPFEM_H "../../../coupledSimulators/npFEM/*.h")
file(GLOB NPFEM_CPP "../../../coupledSimulators/npFEM/*.cpp")
if(ENABLE_CUDA)
    # cuh files
    file(GLOB NPFEM_CUDA_H "../../../coupledSimulators/npFEM/src_GPU/*.h")
    # cu files
    file(GLOB NPFEM_CUDA_SOURCES_CPP "../../../coupledSimulators/npFEM/src_GPU/*.cpp")
    file(GLOB NPFEM_CUDA_SOURCES_CU  "../../../coupledSimulators/npFEM/src_GPU/*.cu")
    set(NPFEM_CUDA_SOURCES ${NPFEM_CUDA_SOURCES_CPP} ${NPFEM_CUDA_SOURCES_CU})  
endif()

# By adding the NPFEM_H, we manage to show the headers in a separate folder 
# in the project tree (quite handy for development)

if(NOT ENABLE_CUDA)
    add_library(npFEM STATIC ${NPFEM_CPP} ${NPFEM_H})
else()
    set(HEADERS ${NPFEM_H} ${NPFEM_CUDA_H})
    set(SOURCES ${NPFEM_CPP} ${NPFEM_CUDA_SOURCES})
    add_library(npFEM STATIC ${HEADERS} ${SOURCES})
    target_link_libraries(npFEM ${CUDA_LINK_LIBRARIES_KEYWORD})
endif()

###############################################################################

file(GLOB DEPS_H "./*.h")
file(GLOB PARAMS_1 "./*.xml")
file(GLOB PARAMS_2 "./*.pos")

add_executable(${EXECUTABLE_NAME} "./${EXECUTABLE_NAME}.cpp" ${DEPS_H} ${PARAMS_1} ${PARAMS_2})
target_link_libraries(${EXECUTABLE_NAME} palabos npFEM)
if(ENABLE_MPI)
	target_link_libraries(${EXECUTABLE_NAME} ${MPI_CXX_LIBRARIES})
endif()

# GPU Version of the code
if(ENABLE_CUDA)
    add_executable(${EXECUTABLE_NAME_GPU} EXCLUDE_FROM_ALL "./${EXECUTABLE_NAME}.cpp" ${DEPS_H} ${PARAMS_1} ${PARAMS_2})
    target_compile_options(${EXECUTABLE_NAME_GPU} PUBLIC -DNPFEM_CUDA=1)
    set_property(TARGET ${EXECUTABLE_NAME_GPU} PROPERTY CUDA_SEPARABLE_COMPILATION ON)
    target_link_libraries(${EXECUTABLE_NAME_GPU} palabos npFEM)
    if(ENABLE_MPI)
        target_link_libraries(${EXECUTABLE_NAME_GPU} ${MPI_CXX_LIBRARIES})
    endif()

    # >> make gpu, to compile the GPU version
    add_custom_target(gpu)
    add_dependencies(gpu ${EXECUTABLE_NAME_GPU})
endif()