#!/bin/bash

# Darwin/macOS wheel build configuration script
# Supports both Docker CI environments and local Homebrew development
# Usage: ./configure.sh <python_executable> [install_prefix]

# Function to add apps
add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

# Set variables
export KRATOS_SOURCE=${KRATOS_ROOT:-.}
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
# export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

export KRATOS_BUILD_TYPE="Release"
export PYTHON_EXECUTABLE=$1
export KRATOS_INSTALL_PREFIX=${2:-"/workspace/kratos/Kratos/bin/Release/Python-314"}

# Detect environment: Docker CI or local Homebrew development
if [ -d "/opt/homebrew/opt" ]; then
    # Local macOS Homebrew development environment
    LLVM_PATH="/opt/homebrew/opt/llvm"
    BOOST_ROOT="/opt/homebrew/opt/boost"
    COMPILER_CXX="${LLVM_PATH}/bin/clang++"
    COMPILER_C="${LLVM_PATH}/bin/clang"
else
    # Docker CI environment (fallback)
    LLVM_PATH="/workspace/llvm"
    BOOST_ROOT="/workspace/boost/boost_1_87_0"
    COMPILER_CXX="clang++"
    COMPILER_C="clang"
fi

# Determine CPU architecture and set appropriate flags
ARCH=$(uname -m)
if [ "$ARCH" = "arm64" ] || [ "$ARCH" = "aarch64" ]; then
    # Apple Silicon - no SSE flags needed
    CPU_FLAGS=""
else
    # Intel Mac - can use SSE3
    CPU_FLAGS="-msse3"
fi

# Set applications to compile (common subset for wheels)
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication
add_app ${KRATOS_APP_DIR}/DEMApplication
add_app ${KRATOS_APP_DIR}/ContactStructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/MPMApplication;
add_app ${KRATOS_APP_DIR}/ConvectionDiffusionApplication;
add_app ${KRATOS_APP_DIR}/DamApplication;
add_app ${KRATOS_APP_DIR}/PoromechanicsApplication;
add_app ${KRATOS_APP_DIR}/FSIApplication;
add_app ${KRATOS_APP_DIR}/SwimmingDEMApplication;
add_app ${KRATOS_APP_DIR}/LinearSolversApplication;
add_app ${KRATOS_APP_DIR}/ConstitutiveLawsApplication;
add_app ${KRATOS_APP_DIR}/MeshingApplication;
add_app ${KRATOS_APP_DIR}/MetisApplication;
add_app ${KRATOS_APP_DIR}/DemStructuresCouplingApplication;
add_app ${KRATOS_APP_DIR}/MeshMovingApplication;
add_app ${KRATOS_APP_DIR}/CSharpWrapperApplication;
add_app ${KRATOS_APP_DIR}/ShapeOptimizationApplication;
add_app ${KRATOS_APP_DIR}/CoSimulationApplication;
add_app ${KRATOS_APP_DIR}/CableNetApplication;
add_app ${KRATOS_APP_DIR}/RANSApplication;
add_app ${KRATOS_APP_DIR}/MappingApplication;
add_app ${KRATOS_APP_DIR}/CompressiblePotentialFlowApplication;
add_app ${KRATOS_APP_DIR}/IgaApplication;
add_app ${KRATOS_APP_DIR}/ChimeraApplication;
add_app ${KRATOS_APP_DIR}/StatisticsApplication;
add_app ${KRATOS_APP_DIR}/RomApplication;
add_app ${KRATOS_APP_DIR}/ShallowWaterApplication;
add_app ${KRATOS_APP_DIR}/OptimizationApplication;
add_app ${KRATOS_APP_DIR}/GeoMechanicsApplication;
add_app ${KRATOS_APP_DIR}/SystemIdentificationApplication;

# Clean
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

# CMake configuration for macOS
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}"    \
-DCMAKE_POLICY_VERSION_MINIMUM=3.5                                     \
-DKRATOS_USE_FUTURE=ON                                                 \
-DCMAKE_INSTALL_PREFIX="${KRATOS_INSTALL_PREFIX}"                      \
-DUSE_TRIANGLE_NONFREE_TPL=ON                                          \
-DMAKE_TRILINOS_OPTIONAL=ON                                            \
-DCMAKE_C_COMPILER="${COMPILER_C}"                                     \
-DCMAKE_CXX_COMPILER="${COMPILER_CXX}"                                 \
-DCMAKE_CXX_FLAGS="${CPU_FLAGS} -std=c++17"                            \
-DCMAKE_C_FLAGS="${CPU_FLAGS}"                                         \
-DBOOST_ROOT="${BOOST_ROOT}"                                           \
-DINCLUDE_MMG=ON                                                       \
-DKRATOS_BUILD_TESTING=OFF                                             \
-DKRATOS_GENERATE_PYTHON_STUBS=ON