#!/bin/bash
# Please do not modify this script

# You can use your interpreter of choice (bash, sh, zsh, ...)

# For any question please contact with us in:
#   - https://github.com/KratosMultiphysics/Kratos

# Optional parameters:
# You can find a list with all the compilation options in INSTALL.md or here:
#   - https://github.com/KratosMultiphysics/Kratos/wiki/Compilation-options

# Function to add apps
add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

# Set compiler
export CC=${CC:-gcc}
export CXX=${CXX:-g++}

# Set variables
export KRATOS_SOURCE="${KRATOS_SOURCE:-"$( cd "$(dirname "$0")" ; pwd -P )"/..}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
# export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

# Set basic configuration
export KRATOS_BUILD_TYPE=${KRATOS_BUILD_TYPE:-"Release"}
export PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE:-"/usr/bin/python3"}

# Set applications to compile
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/HDF5Application
add_app ${KRATOS_APP_DIR}/MetisApplication

# Clean
clear
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

# Configure
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" \
-DUSE_MPI=ON                                                       \
-DUSE_EIGEN_MKL=OFF                                                 \
-DTRILINOS_INCLUDE_DIR="/usr/include/trilinos"                      \
-DTRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu"                  \
-DTRILINOS_LIBRARY_PREFIX="trilinos_"                               \
-DKRATOS_GENERATE_PYTHON_STUBS=ON

# Build
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j$(nproc)