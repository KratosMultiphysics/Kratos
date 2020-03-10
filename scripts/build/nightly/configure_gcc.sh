#!/bin/bash
# You can use your interpreter of choice (bash, sh, zsh, ...)

# For any question please contact with us in:
#   - https://github.com/KratosMultiphysics/Kratos

# Optional parameters:
# You can find a list will all the compiation options in INSTALL.md or here:
#   - https://github.com/KratosMultiphysics/Kratos/wiki/Compilation-options

add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

# Set compiler
export CC=gcc
export CXX=g++

# Set variables
export KRATOS_SOURCE="${KRATOS_SOURCE:-${PWD}}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"

# Set build type
export KRATOS_BUILD_TYPE=${KRATOS_BUILD_TYPE:="Custom"}

# Set applications to compile
add_app ${KRATOS_APP_DIR}/ExternalSolversApplication
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/ContactStructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/FluidDynamicsApplication
add_app ${KRATOS_APP_DIR}/DEMApplication
add_app ${KRATOS_APP_DIR}/FSIDynamicsApplication
add_app ${KRATOS_APP_DIR}/MeshingApplication
add_app ${KRATOS_APP_DIR}/MeshMovingApplication
add_app ${KRATOS_APP_DIR}/MappingApplication
add_app ${KRATOS_APP_DIR}/ConvectionDifussionApplication
add_app ${KRATOS_APP_DIR}/HDF5Application
add_app ${KRATOS_APP_DIR}/EigenSolversApplication
add_app ${KRATOS_APP_DIR}/ShapeOptimizationApplication
add_app ${KRATOS_APP_DIR}/IgaApplication
add_app ${KRATOS_APP_DIR}/DamApplication
add_app ${KRATOS_APP_DIR}/ParticleMechanicsApplication
add_app ${KRATOS_APP_DIR}/PoromechanicsApplication
add_app ${KRATOS_APP_DIR}/CompressiblePotentialFlowApplication
add_app ${KRATOS_APP_DIR}/CableNetApplication
add_app ${KRATOS_APP_DIR}/SolidMechanicsApplication
add_app ${KRATOS_APP_DIR}/ConstitutiveModelsApplication
add_app ${KRATOS_APP_DIR}/CoSimulationApplication
add_app ${KRATOS_APP_DIR}/CSharpWrapperApplication

# Clean
clear
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

# Configure
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}"                     \
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3 -fopenmp"                                      \
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -std=c++11 -fopenmp"                       \
-DPYTHON_EXECUTABLE="/usr/bin/python3.5"                                                \
-DBOOST_DIR="${HOME}/CompiledLibs/boost_1_59_0"                                         \
-DINCLUDE_FEAST=ON                                                                      \
-DINCLUDE_MMG=ON                                                                        \
-DMMG_INCLUDE_DIR="${HOME}/MMGPrecompiled-master/include"                               \
-DMMG2D_INCLUDE_DIR="${HOME}/MMGPrecompiled-master/include/mmg/mmg2d/"                  \
-DMMG3D_INCLUDE_DIR="${HOME}/MMGPrecompiled-master/include/mmg/mmg3d/"                  \
-DMMGS_INCLUDE_DIR="${HOME}/MMGPrecompiled-master/include/mmg/mmgs/"                    \
-DMMG_LIBRARY="${HOME}/MMGPrecompiled-master/lib/libmmg.a"                              \
-DMMG2D_LIBRARY="${HOME}/MMGPrecompiled-master/lib/libmmg2d.a"                          \
-DMMG3D_LIBRARY="${HOME}/MMGPrecompiled-master/lib/libmmg3d.a"                          \
-DMMGS_LIBRARY="${HOME}/MMGPrecompiled-master/lib/libmmgs.a"                            \
-DUSE_EIGEN_MKL=OFF                                                                     \
-DEIGEN_ROOT="${HOME}/eigen/"                                                           \
-DANURBS_ROOT="${HOME}/ANurbs"

# Buid
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target all_unity -- -j1
cmake --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" --target install/fast -- -j1


