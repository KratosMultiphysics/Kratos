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

# Set variables
export KRATOS_SOURCE="${KRATOS_SOURCE:-"$( cd "$(dirname "$0")" ; pwd -P )"/..}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

# Set basic configuration
export KRATOS_BUILD_TYPE=Release
export BOOST_ROOT=/Users/kratos/Tools/boost_1_87_0
export PYTHON_EXECUTABLE=python3

# Set applications to compile
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
add_app ${KRATOS_APP_DIR}/LinearSolversApplication;
add_app ${KRATOS_APP_DIR}/ConstitutiveLawsApplication;
add_app ${KRATOS_APP_DIR}/MeshingApplication;
add_app ${KRATOS_APP_DIR}/DemStructuresCouplingApplication;
add_app ${KRATOS_APP_DIR}/MeshMovingApplication;
add_app ${KRATOS_APP_DIR}/ShapeOptimizationApplication;
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
clear
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

export LIBCXX_ENABLE_VENDOR_AVAILABILITY_ANNOTATIONS=ON

# Configure
/Applications/CMake.app/Contents/bin/cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}"	\
 -DCMAKE_EXPORT_COMPILE_COMMANDS=ON 										                                \
 -DCMAKE_INSTALL_PREFIX=$3                                                                                  \
 -DCMAKE_C_COMPILER=/opt/homebrew/opt/gcc/bin/gcc-15 								                        \
 -DCMAKE_CXX_COMPILER=/opt/homebrew/opt/gcc/bin/g++-15 							                            \
 -DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -I/opt/homebrew/include"            						                \
 -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -I/opt/homebrew/include"         					                \
 -DOpenMP_C_FLAGS="-fopenmp" 											                                    \
 -DOpenMP_CXX_FLAGS="-fopenmp" 											                                    \
 -DOpenMP_C_LIB_NAMES="gomp" 											                                    \
 -DOpenMP_CXX_LIB_NAMES="gomp" 											                                    \
 -DOpenMP_gomp_LIBRARY=/opt/homebrew/lib/gcc/15/libgomp.dylib  							                    \
 -DCMAKE_OSX_SYSROOT=$(xcrun --show-sdk-path)  									                            \
 -DUSE_EIGEN_MKL=OFF                                                                                        \
 -DKRATOS_BUILD_TESTING=OFF 											                                    \
 -DKRATOS_USE_FUTURE=ON 											                                        \
 -DKRATOS_GENERATE_PYTHON_STUBS=ON

# Build
/Applications/CMake.app/Contents/bin/cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install -- -j$(sysctl -n hw.physicalcpu)