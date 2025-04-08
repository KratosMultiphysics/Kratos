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
add_app ${KRATOS_APP_DIR}/HDF5Application;
add_app ${KRATOS_APP_DIR}/MedApplication;
add_app ${KRATOS_APP_DIR}/MappingApplication;
add_app ${KRATOS_APP_DIR}/MeshMovingApplication;
add_app ${KRATOS_APP_DIR}/MeshingApplication;
add_app ${KRATOS_APP_DIR}/StatisticsApplication;


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
-DEXCLUDE_KRATOS_CORE=ON                                            \
-DEXCLUDE_AUTOMATIC_DEPENDENCIES=ON                                 \
-DREMOVE_INSTALL_DIRECTORIES=OFF                                    \
-DCMAKE_CXX_FLAGS="${KRATOS_CMAKE_CXX_FLAGS} -O0 -Wall"             \
-DCMAKE_POLICY_VERSION_MINIMUM=3.5                                  \
-DTRILINOS_INCLUDE_DIR="/usr/include/trilinos"                      \
-DTRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu"                  \
-DTRILINOS_LIBRARY_PREFIX="trilinos_"                               \
-DCMAKE_UNITY_BUILD=ON                                              \
-DINCLUDE_MMG=ON                                                    \
-DCMAKE_C_COMPILER_LAUNCHER=sccache                                 \
-DCMAKE_CXX_COMPILER_LAUNCHER=sccache

# Build
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j${KRATOS_CI_CORES}
