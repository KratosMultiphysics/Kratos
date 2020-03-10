#!/bin/bash

# Function to add apps
add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

# Set compiler
export CC=gcc
export CXX=g++

# Set variables
export KRATOS_SOURCE=${KRATOS_ROOT}
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
# export KRATOS_INSTALL_PYTHON_USING_LINKS=ON


export KRATOS_BUILD_TYPE="Release"
export PYTHON_EXECUTABLE=$1

# Set applications to compile
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication
add_app ${KRATOS_APP_DIR}/DEMApplication
add_app ${KRATOS_APP_DIR}/ContactStructuralMechanicsApplication

# Clean
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

${CMAKE} -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" \
-DUSE_MPI=OFF                                                          \
-DCMAKE_C_COMPILER=/opt/rh/devtoolset-8/root/usr/bin/gcc               \
-DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-8/root/usr/bin/g++             \
-DCMAKE_CXX_FLAGS="-msse3 -std=c++11 "                                 \
-DCMAKE_C_FLAGS="-msse3"                                               \
-DBOOST_ROOT="/workspace/boost/boost_1_71_0"                           \
-DLAPACK_LIBRARIES="/usr/lib64/liblapack.so.3"                         \
-DBLAS_LIBRARIES="/usr/lib64/libblas.so.3"                             \
-DINSTALL_RUNKRATOS=OFF
