#!/bin/bash
# You can use your interpreter of choice (bash, sh, zsh, ...)

# For any question please contact with us in:
#   - https://github.com/KratosMultiphysics/Kratos

# Optional parameters:
# You can find a list will all the compiation options in INSTALL.md or here:
#   - https://github.com/KratosMultiphysics/Kratos/wiki/Compilation-options

add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

# Set compiler
export CC=gcc-6
export CXX=g++-6

# Set variables
export KRATOS_SOURCE="${KRATOS_SOURCE:-${PWD}}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"

# Set build type
export KRATOS_BUILD_TYPE=${KRATOS_BUILD_TYPE:="Custom"}
export PYTHON_EXECUTABLE="/usr/bin/python3.4"
export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

# Set applications to compile
add_app ${KRATOS_APP_DIR}/ExternalSolversApplication;
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication;
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication;
add_app ${KRATOS_APP_DIR}/MeshMovingApplication;
add_app ${KRATOS_APP_DIR}/DEMApplication;
add_app ${KRATOS_APP_DIR}/CSharpWrapperApplication;
add_app ${KRATOS_APP_DIR}/MetisApplication;
add_app ${KRATOS_APP_DIR}/TrilinosApplication;

# Clean
clear
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

echo "Kratos build type is ${KRATOS_BUILD_TYPE}"

# Configure
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" \
-DUSE_MPI=ON \
-DPYBIND11_PYTHON_VERSION="3.4" \
-DCMAKE_CXX_FLAGS="-O0 -fopenmp -Wall \
-Wignored-qualifiers -Werror=ignored-qualifiers -Werror=suggest-override -Werror=unused-variable \
-Werror=misleading-indentation -Werror=return-type \
-Werror=sign-compare -Werror=unused-but-set-variable \
-Werror=unused-local-typedefs -Werror=reorder -Werror=maybe-uninitialized" \
-DTRILINOS_INCLUDE_DIR="/usr/include/trilinos" \
-DTRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu" \
-DTRILINOS_LIBRARY_PREFIX="trilinos_" \
-DUSE_COTIRE=ON

# Buid
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target all_unity    -- -j1
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install/fast -- -j1
