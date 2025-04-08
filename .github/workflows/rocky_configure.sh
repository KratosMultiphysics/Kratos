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
export PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE:-"/usr/bin/python3"}
export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

# Boost setup
BOOST_VERSION="1.86.0"
BOOST_DIR_NAME="boost_1_86_0"
BOOST_TARBALL_URL="https://archives.boost.io/release/${BOOST_VERSION}/source/${BOOST_DIR_NAME}.tar.gz"
BOOST_DOWNLOAD_DIR="${KRATOS_SOURCE}/external_libraries"
BOOST_EXTRACT_DIR="${BOOST_DOWNLOAD_DIR}/${BOOST_DIR_NAME}"

# Download and extract Boost if not already done
mkdir -p "${BOOST_DOWNLOAD_DIR}"
if [ ! -d "${BOOST_EXTRACT_DIR}" ]; then
    echo "Downloading Boost ${BOOST_VERSION}..."
    curl -L "${BOOST_TARBALL_URL}" -o "${BOOST_DOWNLOAD_DIR}/${BOOST_DIR_NAME}.tar.gz"
    echo "Extracting Boost..."
    tar -xzf "${BOOST_DOWNLOAD_DIR}/${BOOST_DIR_NAME}.tar.gz" -C "${BOOST_DOWNLOAD_DIR}"
fi

# Define BOOST_ROOT
export BOOST_ROOT="${BOOST_EXTRACT_DIR}"
echo "BOOST_ROOT is set to ${BOOST_ROOT}"

# Add applications
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication;
add_app ${KRATOS_APP_DIR}/MappingApplication;
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication;
add_app ${KRATOS_APP_DIR}/MeshingApplication;
add_app ${KRATOS_APP_DIR}/LinearSolversApplication;
add_app ${KRATOS_APP_DIR}/ConstitutiveLawsApplication;
add_app ${KRATOS_APP_DIR}/CoSimulationApplication;

# Clean
clear
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

echo "Kratos build type is ${KRATOS_BUILD_TYPE}"

# Configure
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" \
${KRATOS_CMAKE_OPTIONS_FLAGS}                                       \
-DUSE_MPI=OFF                                                       \
-DPYBIND11_PYTHON_VERSION="3.8"                                     \
-DBOOST_ROOT="${BOOST_ROOT}"                                        \
-DCMAKE_CXX_FLAGS="${KRATOS_CMAKE_CXX_FLAGS} -O0 -Wall"             \
-DCMAKE_UNITY_BUILD=ON                                              \
-DCMAKE_C_COMPILER_LAUNCHER=sccache                                 \
-DCMAKE_CXX_COMPILER_LAUNCHER=sccache

# Build
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j2
