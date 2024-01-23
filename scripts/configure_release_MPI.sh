#!/bin/bash
# Please do not modify this script

# You can use your interpreter of choice (bash, sh, zsh, ...)

# For any question please contact with us in:
#   - https://github.com/KratosMultiphysics/Kratos

# Optional parameters:
# You can find a list will all the compiation options in INSTALL.md or here:
#   - https://github.com/KratosMultiphysics/Kratos/wiki/Compilation-options

# Function to add apps
add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

# Set compiler
# export CC=gcc
# export CXX=g++
# export MPI_C=mpicc
# export MPI_CXX=mpicxx

# export CMAKE_CXX_COMPILER=mpicxx
# export CMAKE_C_COMPILER=mpicc
# export CMAKE_Fortran_COMPILER=mpif90
# export Trilinos_ENABLE_CXX11=ON
# export MPI_BASE_DIR="/usr/lib/x86_64-linux-gnu/openmpi"
# export MPI_INCLUDE_DIRS="${MPI_ROOT}/include"
# export MPI_INCLUDE="/usr/include/x86_64-linux-gnu/openmpi"
# export MPI_LIB="/usr/lib/x86_64-linux-gnu/openmpi/lib"
# export MPI_BIN="/usr/bin"

# Set variables
export KRATOS_SOURCE="${KRATOS_SOURCE:-"$( cd "$(dirname "$0")" ; pwd -P )"/..}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

# Set basic configuration
export KRATOS_BUILD_TYPE=${KRATOS_BUILD_TYPE:-"Release"}
export PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE:-"/usr/bin/python3"}
export BOOST_ROOT=${BOOST_ROOT:-"/usr/include/boost"}

# export TRILINOS_INCLUDE_DIR="/usr/include/trilinos"
# export TRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu"
# export TRILINOS_LIBRARY_PREFIX="trilinos_"

#export METIS_DIR="/usr/include"

# Set applications to compile
export KRATOS_APPLICATIONS=
#add_app ${KRATOS_APP_DIR}/ExternalSolversApplication
add_app ${KRATOS_APP_DIR}/LinearSolversApplication
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/ConstitutiveLawsApplication
add_app ${KRATOS_APP_DIR}/HDF5Application
#add_app ${KRATOS_APP_DIR}/DropletDynamicsApplication
#add_app ${KRATOS_APP_DIR}/ConvectionDiffusionApplication
add_app ${KRATOS_APP_DIR}/MetisApplication
add_app ${KRATOS_APP_DIR}/TrilinosApplication
#add_app ${KRATOS_APP_DIR}/MultilevelMonteCarloApplication
#add_app ${KRATOS_APP_DIR}/MeshingApplication
#add_app ${KRATOS_APP_DIR}/MappingApplication
#add_app ${KRATOS_APP_DIR}/ExaquteSandboxApplication
#add_app ${KRATOS_APP_DIR}/StatisticsApplication
add_app ${KRATOS_APP_DIR}/ContactStructuralMechanicsApplication

# Clean
clear
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

echo "Kratos build type is ${KRATOS_BUILD_TYPE}"

# Configure
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" \
${KRATOS_CMAKE_OPTIONS_FLAGS} \
-DUSE_MPI=ON \
-DKRATOS_BUILD_TESTING=ON \
-DPYBIND11_PYTHON_VERSION="3.8" \
-DCMAKE_CXX_FLAGS="${KRATOS_CMAKE_CXX_FLAGS} -std=c++11 -O0 -fopenmp -Wall" \
-DTRILINOS_INCLUDE_DIR="/usr/include/trilinos" \
-DTRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu" \
-DTRILINOS_LIBRARY_PREFIX="trilinos_" #\
# -DMETIS_DIR="/usr/include" \
# -DBLAS_LIBRARIES="/usr/lib/x86_64-linux-gnu/blas/libblas.so.3" \
# -DLAPACK_LIBRARIES="/usr/lib/x86_64-linux-gnu/liblapack.so.3" \
# -DUSE_COTIRE=ON \
# -DINCLUDE_MMG=ON

# Buid
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j8 #$(nproc)
