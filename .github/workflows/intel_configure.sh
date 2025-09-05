#!/bin/bash
# You can use your interpreter of choice (bash, sh, zsh, ...)

# For any question please contact with us in:
#   - https://github.com/KratosMultiphysics/Kratos

# Optional parameters:
# You can find a list with all the compilation options in INSTALL.md or here:
#   - https://github.com/KratosMultiphysics/Kratos/wiki/Compilation-options

add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

# Set variables
export KRATOS_SOURCE="${KRATOS_SOURCE:-${PWD}}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

# Set applications to compile
# add_app ${KRATOS_APP_DIR}/ConvectionDiffusionApplication;
# add_app ${KRATOS_APP_DIR}/LinearSolversApplication;
# add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication;
# add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication;
# add_app ${KRATOS_APP_DIR}/FluidDynamicsBiomedicalApplication;
# add_app ${KRATOS_APP_DIR}/MeshMovingApplication;
# add_app ${KRATOS_APP_DIR}/DEMApplication;
# add_app ${KRATOS_APP_DIR}/CSharpWrapperApplication;
# add_app ${KRATOS_APP_DIR}/MetisApplication;
# add_app ${KRATOS_APP_DIR}/TrilinosApplication;
# add_app ${KRATOS_APP_DIR}/ShapeOptimizationApplication;
# add_app ${KRATOS_APP_DIR}/CoSimulationApplication;
# add_app ${KRATOS_APP_DIR}/CableNetApplication;
# add_app ${KRATOS_APP_DIR}/RANSApplication;
# add_app ${KRATOS_APP_DIR}/MappingApplication;
# add_app ${KRATOS_APP_DIR}/FSIApplication;
# add_app ${KRATOS_APP_DIR}/MeshingApplication;
# add_app ${KRATOS_APP_DIR}/CompressiblePotentialFlowApplication;
# add_app ${KRATOS_APP_DIR}/HDF5Application;
# add_app ${KRATOS_APP_DIR}/ContactStructuralMechanicsApplication;
# add_app ${KRATOS_APP_DIR}/IgaApplication;
# add_app ${KRATOS_APP_DIR}/MPMApplication;
# add_app ${KRATOS_APP_DIR}/ChimeraApplication;
# add_app ${KRATOS_APP_DIR}/MultilevelMonteCarloApplication;
# add_app ${KRATOS_APP_DIR}/StatisticsApplication;
# add_app ${KRATOS_APP_DIR}/SwimmingDEMApplication;
# add_app ${KRATOS_APP_DIR}/ConstitutiveLawsApplication;

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
-DPYTHON_EXECUTABLE="/usr/bin/python3.10" \
-DCMAKE_CXX_FLAGS="${KRATOS_CMAKE_CXX_FLAGS} -O3 -Wall -Werror-all -diag-disable 1478 -diag-disable 1786" \
-DTRILINOS_INCLUDE_DIR="/usr/include/trilinos" \
-DTRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu" \
-DTRILINOS_LIBRARY_PREFIX="trilinos_" \
-DCMAKE_UNITY_BUILD=ON

# Build
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j2
