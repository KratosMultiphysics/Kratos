 
# Function to add apps
add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

# Set compiler
export CC=gcc
export CXX=g++

export CMAKE_C_COMPILER=/usr/bin/gcc
export MAKE_CXX_COMPILER=/usr/bin/g++
export CMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -std=c++11 "

# Set variables
export KRATOS_SOURCE="${KRATOS_SOURCE:-"$( cd "$(dirname "$0")" ; pwd -P )"/..}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
export KRATOS_INSTALL_PYTHON_USING_LINKS=ON
export INSTALL_EMBEDDED_PYTHON=ON

# Set basic configuration
export KRATOS_BUILD_TYPE="Release"
export PYTHON_LIBRARY="/usr/lib/python3.6/config-3.6m-x86_64-linux-gnu/libpython3.6m.so"
export PYTHON_EXECUTABLE="/usr/bin/python3"
export BOOST_ROOT="/home/elaf/boost_1_59_0"
export PYTHON_INCLUDE_DIR="/usr/include/python3.6"

# Set applications to compile
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/ExternalSolversApplication
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication
add_app ${KRATOS_APP_DIR}/ConvectionDiffusionApplication
add_app ${KRATOS_APP_DIR}/MeshingApplication
add_app ${KRATOS_APP_DIR}/FSIApplication
add_app ${KRATOS_APP_DIR}/ULFapplication

# Clean
clear
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

# Configure
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" -DUSE_MPI=OFF -DINSTALL_RUNKRATOS=ON -DKRATOS_BUILD_TESTING=ON -DUSE_TETGEN_NONFREE_TPL=OFF -DTETGEN_INCLUDE_DIR="/home/elaf/tetgen1.4.3" -DTETGEN_LIBRARY_DIR="/home/elaf/tetgen1.4.3"

# Buid
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j4

