#!/bin/bash
# You can use your interpreter of choice (bash, sh, zsh, ...)

# For any question please contact with us in:
#   - https://github.com/KratosMultiphysics/Kratos

# Set compiler
export CC=gcc
export CXX=g++

# Set variables
export KRATOS_SOURCE="${KRATOS_SOURCE:-${PWD}}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"

# Set applications to compile
export KRATOS_BUILD_TYPE=${KRATOS_BUILD_TYPE:="Release"}
export KRATOS_APPLICATIONS="${KRATOS_APP_DIR}/ExternalSolversApplication;"
# export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}${KRATOS_APP_DIR}/StructuralMechanicsApplication;"
# export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}${KRATOS_APP_DIR}/FluidDynamicsApplication;"
# export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}${KRATOS_APP_DIR}/MeshMovingApplication;"
# export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}${KRATOS_APP_DIR}/DEMApplication;"
# export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}${KRATOS_APP_DIR}/CSharpWrapperApplication;"
export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}${KRATOS_APP_DIR}/MetisApplication;"
export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}${KRATOS_APP_DIR}/TrilinosApplication;"

# Clean
clear
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

echo "Kratos build type is ${KRATOS_BUILD_TYPE}"

# Configure
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" \
-DMPI_NEEDED=ON \
-DTRILINOS_APPLICATION=ON \
-DUSE_COTIRE=ON

# Buid
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target all_unity    -- -j4
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install/fast -- -j4
