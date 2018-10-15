#!/bin/bash
# You can use your interpreter of choice (bash, sh, zsh, ...)

# For any question please contact:
#   - kratos@listas.cimne.upc.edu
#   - croig@cimne.upc.edu

# Or visit our forums at
#   - http://kratos-wiki.cimne.upc.edu/forum/index.php

########################################################################################################################
#### CLEANUP                                                                                                        ####
########################################################################################################################
# Clean up the terminal
clear

# Clean up the temporal files from other compilations
rm -rf CMakeCache.txt 2> /dev/null
rm -rf *.cmake 2> /dev/null
rm -rf CMakeFiles 2> /dev/null

########################################################################################################################
#### CONFIGURATION                                                                                                  ####
########################################################################################################################
# Configuration parameters

# Kratos root
#    Indicate the source directory of kratos. Should be the location that containts the cmake_build folder
#    with this script. By default it points to .. which should work if you run the command from the cmake_build dir
#    - Example:
#        KRATOS_ROOT=/home/youruser/kratos
# --------------------------------------------------------------------------------------------------------------
KRATOS_ROOT=".."

# CMake
#    Indicates the cmake binary you want to use. In ubuntu systems you can simply call the default one
#    - Example
#        CMAKE_BIN=cmake
# --------------------------------------------------------------------------------------------------------------
CMAKE_BIN=cmake

# Compilers
#    Indicate where your compiler can be found. In ubuntu you can simply call the default one
#    - Example:
#        C_COMPILER=gcc
#        CXX_COMPILER=g++
# --------------------------------------------------------------------------------------------------------------
C_COMPILER=gcc-6
CXX_COMPILER=g++-6

# Build type
#    Indicate the build type. Possible values are "Release", "RelWithDebInfo" or "Debug"
#    - Example:
#        BUILD_TYPE="RelWithDebInfo"
# --------------------------------------------------------------------------------------------------------------
BUILD_TYPE="Custom"

# Flags for performance
#    Indicate the compiler performance related flags. Please notice that OX flags are added by cmake
#    depending on the Build type selected.
#    - Example:
#        C_PERF_FLAGS="-msse3 -fopenmp"
#        CXX_PERF_FLAGS="-msse3 -fopenmp"
# --------------------------------------------------------------------------------------------------------------
C_PERF_FLAGS="-O0 -fopenmp"
CXX_PERF_FLAGS="-O0 -fopenmp"

# Flags for warnings control
#    Indicate the warning related flags. Wall is enabled by default. This is the proper place to ignore
#    specific warnings
#    - Example:
#        C_WARN_FLAGS="-Wall"
#        CXX_WARN_FLAGS="-Wall"
# --------------------------------------------------------------------------------------------------------------
C_WARN_FLAGS="-Wall"
CXX_WARN_FLAGS="-Wall -Werror=suggest-override -Werror=unused-variable \
                -Werror=misleading-indentation -Werror=return-type \
                -Werror=sign-compare -Werror=unused-but-set-variable"

# Other flags
#    Indicate any other flag you want to add here
#    specific warnings
#    - Example:
#        C_CUSTOM_FLAGS="-fmessage-length=20"
#        CXX_CUSTOM_FLAGS="-fmessage-length=20"
# --------------------------------------------------------------------------------------------------------------
C_CUSTOM_FLAGS=""
CXX_CUSTOM_FLAGS=""

CMAKE_LIBS=(
  # Boost
  #    Indicate your boost root directory in case you don't want to use the system default
  #    - Example:
  #        -DBOOST_ROOT="${KRATOS_ROOT}/external_libraries/boost_1_61_0"
  # --------------------------------------------------------------------------------------------------------------
  # Should be compiling against the system boost
  # -DBOOST_ROOT="${KRATOS_ROOT}/external_libraries/boost_1_61_0"

  # Python
  #    Indicate your python binary dir in case you don't want to use the system default or you
  #    have multiple versions and you want to select one in particular
  #    - Example (for ubuntu 14.04):
  #        -DPYTHON_EXECUTABLE="/usr/bin/python
  # --------------------------------------------------------------------------------------------------------------
  -DPYTHON_EXECUTABLE="/usr/bin/python3.4"
)

########################################################################################################################
#### APPLICATIONS                                                                                                   ####
########################################################################################################################
# List of applications to compile. Set to ON/OFF

CMAKE_APPLICATION=(
  -DEXTERNAL_SOLVERS_APPLICATION=ON
  -DSTRUCTURAL_MECHANICS_APPLICATION=ON
  -DFLUID_DYNAMICS_APPLICATION=ON
  -DDEM_APPLICATION=ON
  -DSWIMMING_DEM_APPLICATION=OFF
)

########################################################################################################################
#### ADVANCED                                                                                                       ####
########################################################################################################################

 ## You shouldn't touch anything below here unless you know what you are doing

########################################################################################################################
#### COMPILER                                                                                                       ####
########################################################################################################################
# Build block for cmake

CMAKE_BUILD=(
   # CMake C compiler
  -DCMAKE_C_COMPILER=${C_COMPILER}
  -DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} ${C_PERF_FLAGS} ${C_WARN_FLAGS} ${C_CUSTOM_FLAGS}"

  # CMake C++ compiler
  # Please DO NOT REMOVE THE "-std=c++11" FLAG.
  -DCMAKE_CXX_COMPILER=${CXX_COMPILER}
  -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -std=c++11 ${CXX_PERF_FLAGS} ${CXX_WARN_FLAGS} ${CXX_CUSTOM_FLAGS}"

  # Build type
  # NOTE: This is better commented for travis since we don't want to use
  # a default configuration. (-O0 is prefered here)
  # -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"

  # Install info
  -DCMAKE_INSTALL_RPATH="${KRATOS_ROOT}/libs"
  -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE

  # Python version to build pybind modules
  -DPYBIND11_PYTHON_VERSION="3.4"
)

########################################################################################################################
#### EXTRA                                                                                                          ####
########################################################################################################################
# Set additional arguments ( for example trillinos etc....)

CMAKE_EXTRA=(
  # Skip problems
  -DCMAKE_SKIP_RPATH=ON

  # Runkratos
  -DINSTALL_EMBEDDED_PYTHON=ON

  # Metis
  -DMETIS_APPLICATION=OFF
  -DUSE_METIS_5=ON
  -DMETIS_ROOT_DIR="/home/roigcarlo/CompiledLibs/metis-5.1.0"

  # Trillinos
  -DTRILINOS_APPLICATION=OFF
  -DTRILINOS_ROOT="/home/youruser/compiled_libraries/trilinos-10.2.0"

  -DUSE_COTIRE=ON
)

########################################################################################################################
#### COMPILE                                                                                                        ####
########################################################################################################################

# This line issues the configuration command with cmake that will generate the Makefile
${CMAKE_BIN} ${KRATOS_ROOT} "${CMAKE_BUILD[@]}" "${CMAKE_LIBS[@]}" "${CMAKE_APPLICATION[@]}" "${CMAKE_EXTRA[@]}"
