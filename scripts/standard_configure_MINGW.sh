#!/bin/bash
# Please do not modify this script

# You can use your interpreter of choice (bash, sh, zsh, ...)

# For any question please contact with us in:
#   - https://github.com/KratosMultiphysics/Kratos

# Optional parameters:
# You can find a list with all the compilation options in INSTALL.md or here:
#   - https://github.com/KratosMultiphysics/Kratos/wiki/Compilation-options

# Function to add apps
add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

# Set compiler #NOTE: Currently only GCC is supported, linking error on recent versions of Clang/LLVM, see https://github.com/llvm/llvm-project/issues/53433
export CC=${CC:-gcc}
export CXX=${CXX:-g++}

# Set variables
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
# export KRATOS_INSTALL_PYTHON_USING_LINKS=ON
export KRATOS_SHARED_MEMORY_PARALLELIZATION=${KRATOS_SHARED_MEMORY_PARALLELIZATION:-"OpenMP"}

# Set basic configuration
export PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE:-"C:/Windows/py.exe"}

# Set applications to compile
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/LinearSolversApplication
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication

# Configure
cmake ..                                                                                            \
-G "MinGW Makefiles"                                                                                \
-DWIN32=TRUE                                                                                        \
-DCMAKE_INSTALL_PREFIX="${KRATOS_SOURCE}/bin/${KRATOS_BUILD_TYPE}"                                  \
-DCMAKE_BUILD_TYPE="${KRATOS_BUILD_TYPE}"                                                           \
-DCMAKE_EXE_LINKER_FLAGS="-s"                                                                       \
-DCMAKE_SHARED_LINKER_FLAGS="-s"                                                                    \
-DCMAKE_POLICY_VERSION_MINIMUM=3.5                                                                  \
-H"${KRATOS_SOURCE}"                                                                                \
-B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}"                                                            \
-DUSE_MPI=OFF                                                                                       \
-DKRATOS_SHARED_MEMORY_PARALLELIZATION="${KRATOS_SHARED_MEMORY_PARALLELIZATION}"                    \
-DUSE_EIGEN_MKL=OFF

# Build
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j$(nproc)