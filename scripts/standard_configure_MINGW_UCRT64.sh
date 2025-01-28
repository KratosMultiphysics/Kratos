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

# Set compiler
export CC=${CC:-clang}
export CXX=${CXX:-clang++}

# Set variables
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
# export KRATOS_INSTALL_PYTHON_USING_LINKS=ON
export KRATOS_SHARED_MEMORY_PARALLELIZATION=${KRATOS_SHARED_MEMORY_PARALLELIZATION:-"OpenMP"}

# Set basic configuration
export PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE:-"C:/Windows/py.exe"}

export CMAKE_EXPORT_COMPILE_COMMANDS=${CMAKE_EXPORT_COMPILE_COMMANDS:-"OFF"}
#export KRATOS_CMAKE_CXX_FLAGS="-Wno-maybe-uninitialized -Wignored-qualifiers"

# Set applications to compile
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/LinearSolversApplication
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication

# Set compilation tool
export COMPILATION_TOOL=${COMPILATION_TOOL:-"Ninja"}

# Determine the folder where the compiler is located
export MSYS_BIN_FOLDER=$(dirname $(which ${CXX}))

# Set CMake strip
export CMAKE_STRIP=${CMAKE_STRIP:-"${MSYS_BIN_FOLDER}/strip.exe"}

# Configure
cmake ..                                                                                            \
-G "${COMPILATION_TOOL}"                                                                            \
-DCMAKE_INSTALL_PREFIX="${KRATOS_SOURCE}/bin/${KRATOS_BUILD_TYPE}"                                  \
-DCMAKE_BUILD_TYPE="${KRATOS_BUILD_TYPE}"                                                           \
-H"${KRATOS_SOURCE}"                                                                                \
-B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}"                                                            \
-DCMAKE_STRIP="${CMAKE_STRIP}"                                                                      \
-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}                                                                    \
-DCMAKE_CXX_FLAGS="${KRATOS_CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS}"                                    \
-DCMAKE_EXPORT_COMPILE_COMMANDS=${CMAKE_EXPORT_COMPILE_COMMANDS}                                    \
-DKRATOS_BUILD_TESTING=ON                                                                           \
-DKRATOS_SHARED_MEMORY_PARALLELIZATION="${KRATOS_SHARED_MEMORY_PARALLELIZATION}"                    \
-DUSE_EIGEN_MKL=OFF

# Build
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j$(nproc)