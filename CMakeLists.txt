cmake_minimum_required(VERSION 2.8)
# # # This CMakeLists.txt file can be used to test the compilation of palabos using CMake.
# # # It can be also used as a template to build your own application and serve as "project builder" for IDE for palabos
# # # developers.

if (${CMAKE_VERSION} VERSION_LESS "3.0.0")
    project(palabos) # change to the name of your project if you want to compile your application
    set(CMAKE_LANGUAGES_COMPILER CXX)
    find_package(MPI REQUIRED)
else ()
    project(palabos LANGUAGES CXX) # change to the name of your project if you woant to compile your application
endif ()

set(DEFAULT ON CACHE INTERNAL "Default options")

# # # If you want to compile your application uncomment the next line
# set(EXECUTABLE_NAME "my_executable")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()
message("Generated with config types: ${CMAKE_CONFIGURATION_TYPES}")
message(${CMAKE_BUILD_TYPE})

set(PALABOS_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
# # # If you want to compile your application, or this file is outside the palabos root folder
# # # add to PATH the palabos folder, comment set(PALABOS_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
# # # and uncomment the next line
#file(TO_CMAKE_PATH $ENV{PALABOS_ROOT} PALABOS_ROOT)

# Compiler flags
# Append flags: set(CMAKE_XXX_FLAGS "${CMAKE_XXX_FLAGS} ...")
if (${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
    message("GCC.")
    if (CYGWIN)
        set(CMAKE_CXX_FLAGS "-Wa,-mbig-obj -std=c++17 -Wall -Wnon-virtual-dtor")
    else ()
        set(CMAKE_CXX_FLAGS "-std=c++17 -Wall -Wnon-virtual-dtor")
    endif ()
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -DPLB_DEBUG -O0")
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
    message("Clang.")
    set(CMAKE_CXX_FLAGS "-std=c++17 -Wall -Wnon-virtual-dtor")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -DPLB_DEBUG -O0")
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    message("MSVC.")
    add_definitions(/DMSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "/MT /Ox /Ot /GS- /GL /DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "/MTd /DPLB_DEBUG")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/LTCG /INCREMENTAL:NO /OPT:REF")
else ()
    message(FATAL_ERROR "CXX compiler not recognized. CMake will quit.")
endif ()

option(ENABLE_MPI "Enable MPI" ${DEFAULT})
if (ENABLE_MPI)
    message("Enabling MPI")
    find_package(MPI REQUIRED)
    if (MPI_CXX_FOUND)
        set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")
        include_directories(${MPI_CXX_INCLUDE_DIRS})
        add_definitions(-DPLB_MPI_PARALLEL)
    endif ()
endif ()

#option(ENABLE_OPENMP "Enabling OpenMP" ${DEFAULT})
#if (ENABLE_OPENMP)
#    message("Enabling OpenMP")
#    find_package(OpenMP QUIET)
#    if (OPENMP_FOUND)
#        message("OpenMP found. OpenMP activated in release.")
#        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
#        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${OpenMP_EXE_LINKER_FLAGS}")
#    else ()
#        message("OpenMP not found.")
#    endif ()
#endif ()

if (WIN32)
    option(ENABLE_POSIX "Enable POSIX" OFF)
else ()
    option(ENABLE_POSIX "Enable POSIX" ${DEFAULT})
endif ()
if (ENABLE_POSIX)
    message("Enabling POSIX")
    add_definitions(-DPLB_USE_POSIX)
endif ()

if (APPLE)
    add_definitions(-DPLB_MAC_OS_X)
endif ()

if (WIN32 OR CYGWIN)
    add_definitions(-DPLB_WINDOWS)
endif ()

include_directories("${PALABOS_ROOT}/src")
include_directories("${PALABOS_ROOT}/externalLibraries")
include_directories("${PALABOS_ROOT}/externalLibraries/Eigen3")

file(GLOB_RECURSE PALABOS_SRC "${PALABOS_ROOT}/src/*.cpp")
file(GLOB_RECURSE EXT_SRC "${PALABOS_ROOT}/externalLibraries/tinyxml/*.cpp")

add_library(palabos STATIC ${PALABOS_SRC} ${EXT_SRC})

# # # EXAMPLE OF HOW TO COMPILE YOUR OWN APPLICATION:
#file(GLOB DEPS_H "mypath/*.h")
#file(GLOB PARAMS "mypathxml/*.xml")
#add_executable(${EXECUTABLE_NAME} "./${EXECUTABLE_NAME}.cpp" ${DEPS_H} ${PARAMS})
#
# # # Link with the compiled libraries
#target_link_libraries(${EXECUTABLE_NAME} palabos)
#if(ENABLE_MPI)
#    target_link_libraries(${EXECUTABLE_NAME} ${MPI_CXX_LIBRARIES})
#endif()