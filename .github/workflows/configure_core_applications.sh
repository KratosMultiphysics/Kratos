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
export PYTHON_EXECUTABLE="/usr/bin/python3.10"
export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

# Set applications to compile
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication;
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication;
add_app ${KRATOS_APP_DIR}/ConstitutiveLawsApplication;
add_app ${KRATOS_APP_DIR}/GeoMechanicsApplication;
add_app ${KRATOS_APP_DIR}/ConvectionDiffusionApplication;
add_app ${KRATOS_APP_DIR}/RomApplication;

# Clean
clear
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

echo "Kratos build type is ${KRATOS_BUILD_TYPE}"

# Configure
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" \
${KRATOS_CMAKE_OPTIONS_FLAGS}                                       \
-DUSE_MPI=ON                                                        \
-DBOOST_ROOT="/workspace/boost/boost_1_87_0"                        \
-DBoost_NO_SYSTEM_PATHS=ON                                          \
-DEXCLUDE_KRATOS_CORE=ON                                            \
-DEXCLUDE_AUTOMATIC_DEPENDENCIES=ON                                 \
-DREMOVE_INSTALL_DIRECTORIES=OFF                                    \
-DCMAKE_CXX_FLAGS="${KRATOS_CMAKE_CXX_FLAGS} -O0 -Wall"             \
-DCMAKE_POLICY_VERSION_MINIMUM=3.5                                  \
-DKRATOS_USE_FUTURE=ON                                              \
-DKRATOS_USE_LEGACY=OFF                                             \
-DTRILINOS_INCLUDE_DIR="/usr/include/trilinos"                      \
-DTRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu"                  \
-DTRILINOS_LIBRARY_PREFIX="trilinos_"                               \
-DUSE_EIGEN_SUITESPARSE:BOOL=ON                                     \
-DCMAKE_UNITY_BUILD=ON

# Build
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j${KRATOS_CI_CORES}
