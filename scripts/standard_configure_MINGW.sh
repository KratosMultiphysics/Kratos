#!/bin/bash
# Please do not modify this script

# You can use your interpreter of choice (bash, sh, zsh, ...)

# For any question please contact with us in:
#   - https://github.com/KratosMultiphysics/Kratos

# Optional parameters:
# You can find a list will all the compiation options in INSTALL.md or here:
#   - https://github.com/KratosMultiphysics/Kratos/wiki/Compilation-options

# Function to add apps
add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

# Set compiler
export CC=gcc
export CXX=g++

# Set variables
export KRATOS_SOURCE="kratos_folder/Kratos"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
# export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

# Set basic configuration
export KRATOS_BUILD_TYPE=${KRATOS_BUILD_TYPE:-"Release"}
export PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE:-"location_python/python.exe"}

# Set applications to compile
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/ExternalSolversApplication
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication
add_app ${KRATOS_APP_DIR}/ConvectionDiffusionApplication

# Configure
cmake ..                                                                                            \
-G "MinGW Makefiles"                                                                                \
-DWIN32=TRUE                                                                                        \
-DCMAKE_INSTALL_PREFIX="${KRATOS_SOURCE}/bin/${KRATOS_BUILD_TYPE}"                                  \
-DCMAKE_BUILD_TYPE="${KRATOS_BUILD_TYPE}"                                                           \
-DCMAKE_EXE_LINKER_FLAGS="-s"                                                                       \
-DCMAKE_SHARED_LINKER_FLAGS="-s"                                                                    \
-DLAPACK_LIBRARIES="location_msys64/msys64/mingw64/lib/liblapack.dll.a"                             \
-DBLAS_LIBRARIES="location_msys64/msys64/mingw64/lib/libblas.dll.a"                                 \
-H"${KRATOS_SOURCE}"                                                                                \
-B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}"                                                            \
-DUSE_MPI=OFF                                                                                       \

# Buid
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j$(nproc)
