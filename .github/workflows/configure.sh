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

# Set variables
export KRATOS_SOURCE="${KRATOS_SOURCE:-${PWD}}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
export PYTHON_EXECUTABLE="/usr/bin/python3.6"
export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

# Set applications to compile
add_app ${KRATOS_APP_DIR}/ExternalSolversApplication;
add_app ${KRATOS_APP_DIR}/EigenSolversApplication;
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication;
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication;
add_app ${KRATOS_APP_DIR}/MeshMovingApplication;
add_app ${KRATOS_APP_DIR}/DEMApplication;
add_app ${KRATOS_APP_DIR}/CSharpWrapperApplication;
add_app ${KRATOS_APP_DIR}/MetisApplication;
add_app ${KRATOS_APP_DIR}/TrilinosApplication;
add_app ${KRATOS_APP_DIR}/ShapeOptimizationApplication;
add_app ${KRATOS_APP_DIR}/CoSimulationApplication;
add_app ${KRATOS_APP_DIR}/CableNetApplication;
add_app ${KRATOS_APP_DIR}/RANSApplication;
add_app ${KRATOS_APP_DIR}/MappingApplication;
add_app ${KRATOS_APP_DIR}/FSIApplication;
add_app ${KRATOS_APP_DIR}/MeshingApplication;
add_app ${KRATOS_APP_DIR}/CompressiblePotentialFlowApplication;
add_app ${KRATOS_APP_DIR}/HDF5Application;
add_app ${KRATOS_APP_DIR}/ContactStructuralMechanicsApplication;

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
-DPYBIND11_PYTHON_VERSION="3.6" \
-DCMAKE_CXX_FLAGS="${KRATOS_CMAKE_CXX_FLAGS} -std=c++11 -O0 -fopenmp -Wall" \
-DTRILINOS_INCLUDE_DIR="/usr/include/trilinos" \
-DTRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu" \
-DTRILINOS_LIBRARY_PREFIX="trilinos_" \
-DUSE_COTIRE=ON \
-DINCLUDE_MMG=ON                                    \
-DMMG_INCLUDE_DIR="/usr/local/include/mmg/"         \
-DMMG2D_INCLUDE_DIR="/usr/local/include/mmg/mmg2d/" \
-DMMG3D_INCLUDE_DIR="/usr/local/include/mmg/mmg3d/" \
-DMMGS_INCLUDE_DIR="/usr/local/include/mmg/mmgs/"   \
-DMMG_LIBRARY="/usr/local/lib/libmmg.so"            \
-DMMG2D_LIBRARY="/usr/local/lib/libmmg2d.so"        \
-DMMG3D_LIBRARY="/usr/local/lib/libmmg3d.so"        \
-DMMGS_LIBRARY="/usr/local/lib/libmmgs.so"

# Buid
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target all_unity    -- -j1
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install/fast -- -j1

