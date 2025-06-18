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
export CC=gcc
export CXX=g++

# Set variables
export KRATOS_SOURCE="${KRATOS_SOURCE:-"$( cd "$(dirname "$0")" ; pwd -P )"/..}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
# export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

# Set basic configuration
export KRATOS_BUILD_TYPE=${KRATOS_BUILD_TYPE:-"Release"}
export PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE:-"/usr/bin/python3"}

#export KRATOS_BUILD_TYPE=${KRATOS_BUILD_TYPE:-"Debug"}
#export PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE:-"/usr/bin/python3-dbg"}



# Set applications to compile
export KRATOS_APPLICATIONS=
#add_app ${KRATOS_APP_DIR}/LinearSolversApplication
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication
add_app ${KRATOS_APP_DIR}/ConvectionDiffusionApplication
add_app ${KRATOS_APP_DIR}/MeshingApplication
add_app ${KRATOS_APP_DIR}/PfemMeltingApplication

# Clean
clear
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

# Configure
#cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" -DUSE_MPI=OFF -DUSE_EIGEN_MKL=OFF

# Configure
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" -DUSE_MPI=OFF -DINSTALL_RUNKRATOS=OFF -DKRATOS_BUILD_TESTING=OFF -DUSE_TRIANGLE_NONFREE_TPL=ON -DUSE_TETGEN_NONFREE_TPL=ON -DUSE_TETGEN_NONFREE_TPL_URL=https://github.com/PFEM/tetgen-1.5.0/archive/kratos.zip -DUSE_TRIANGLE_NONFREE_TPL_URL=https://github.com/PFEM/triangle-1.6.0/

# Buid
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j2

