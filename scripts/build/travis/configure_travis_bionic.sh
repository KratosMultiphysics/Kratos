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

export KRATOS_BUILD_TYPE=${KRATOS_BUILD_TYPE:="Custom"}
export KRATOS_APPLICATIONS="${KRATOS_APP_DIR}/ExternalSolversApplication;"
export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}${KRATOS_APP_DIR}/StructuralMechanicsApplication;"
export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}${KRATOS_APP_DIR}/FluidDynamicsApplication;"
export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}${KRATOS_APP_DIR}/MeshMovingApplication;"
export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}${KRATOS_APP_DIR}/DEMApplication;"
export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}${KRATOS_APP_DIR}/CSharpWrapperApplication;"
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
-DPYTHON_EXECUTABLE="/usr/bin/python3.6" -DMPI_NEEDED=ON \
-DPYBIND11_PYTHON_VERSION="3.6" \
-DCMAKE_CXX_FLAGS="-O0 -fopenmp -Wall \
-Wignored-qualifiers -Werror=ignored-qualifiers -Werror=suggest-override -Werror=unused-variable \
-Werror=misleading-indentation -Werror=return-type \
-Werror=sign-compare -Werror=unused-but-set-variable \
-Werror=unused-local-typedefs -Werror=reorder -Werror=maybe-uninitialized" \
-DINSTALL_PYTHON_USING_LINKS=ON \
-DTRILINOS_APPLICATION=ON \
-DTRILINOS_INCLUDE_DIR="/usr/include/trilinos" \
-DTRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu" \
-DTRILINOS_LIBRARY_PREFIX="trilinos_" \
-DUSE_COTIRE=ON

# Buid
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target all_unity    -- -j1
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install/fast -- -j1

