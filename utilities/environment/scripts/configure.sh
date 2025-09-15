#!/bin/bash
# Please do not modify this script

# You can use your interpreter of choice (bash, sh, zsh, ...)

# For any question please contact with us in:
#   - https://github.com/KratosMultiphysics/Kratos

# Optional parameters:
# You can find a list will all the compiation options in INSTALL.md or here:
#   - https://github.com/KratosMultiphysics/Kratos/wiki/Compilation-options

source /opt/intel/oneapi/setvars.sh intel64

# Function to add apps
add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

# assumes the python environment is activated.
kratos_install_dir=$(python -c 'import site; print(site.getsitepackages()[0])')
python_venv_name=$(which python | rev | cut -d"/" -f 3 | rev)

# Set compiler
export compiler_type=$(echo $python_venv_name | rev | cut -d"_" -f2 | rev)
export KRATOS_BUILD_TYPE=$(echo $python_venv_name | rev | cut -d"_" -f1 | rev)
export KRATOS_CPP_CONFIG_NAME="${compiler_type}_${KRATOS_BUILD_TYPE}"

current_kratos_name="<KRATOS_NAME>"
current_kratos_name_from_venv=$(echo $python_venv_name | rev | cut -c$((${#KRATOS_CPP_CONFIG_NAME}+2))- | rev)
if [ ! "${current_kratos_name}" = "${current_kratos_name_from_venv}" ]; then
    echo "Mismatching kratos ($current_kratos_name) and python env($python_venv_name) name."
    exit
fi



case $compiler_type in
    "gcc")
        export CC=gcc
        export CXX=g++
        ;;
    "clang")
        export CC=clang
        export CXX=clang++
        ;;
    "intel")
        export CC=icx
        export CXX=icpx
        ;;
    *)
        echo "-- Unsupported compiler type provided in environment variable KRATOS_CPP_CONFIG_NAME. [ python_venv_name = \"$python_venv_name\" ]. Followings are the accepted compiler types:"
        printf "\tgcc\n"
        printf "\tclang\n"
        printf "\tintel\n"
        exit
        ;;
esac

echo "===== Compiling kratos using \"$compiler_type\" in \"$KRATOS_BUILD_TYPE\" with \"$python_venv_name\" python environemnt"

# Set variables
export KRATOS_SOURCE="${KRATOS_SOURCE:-"$( cd "$(dirname "$0")" ; pwd -P )"/..}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
export KRATOS_INSTALL_PYTHON_USING_LINKS=ON

# Set basic configuration
# export KRATOS_BUILD_TYPE=${KRATOS_BUILD_TYPE:-${KRATOS_BUILD_TYPE}}
export PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE:-"python"}

# Set applications to compile
export KRATOS_APPLICATIONS=
add_app ${KRATOS_APP_DIR}/HDF5Application
add_app ${KRATOS_APP_DIR}/MetisApplication
add_app ${KRATOS_APP_DIR}/TrilinosApplication
add_app ${KRATOS_APP_DIR}/LinearSolversApplication
add_app ${KRATOS_APP_DIR}/StructuralMechanicsApplication
add_app ${KRATOS_APP_DIR}/ShapeOptimizationApplication
add_app ${KRATOS_APP_DIR}/OptimizationApplication
add_app ${KRATOS_APP_DIR}/MeshMovingApplication

# Clean
clear
rm -rf "${KRATOS_BUILD}/${KRATOS_CPP_CONFIG_NAME}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_CPP_CONFIG_NAME}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_CPP_CONFIG_NAME}/CMakeFiles"

# Configure
cmake -H"${KRATOS_SOURCE}" -B"${KRATOS_BUILD}/${KRATOS_CPP_CONFIG_NAME}"   \
-DCMAKE_INSTALL_PREFIX="${kratos_install_dir}"                             \
-DCMAKE_EXPORT_COMPILE_COMMANDS=ON                                         \
-DUSE_MPI=ON                                                               \
-DUSE_EIGEN_MKL=ON                                                         \
-DCMAKE_CXX_COMPILER_LAUNCHER=ccache                                       \
-DKRATOS_GENERATE_PYTHON_STUBS=ON                                          \

# Buid
cmake --build "${KRATOS_BUILD}/${KRATOS_CPP_CONFIG_NAME}" --target install -j10

ln -sf ${KRATOS_BUILD}/${KRATOS_CPP_CONFIG_NAME}/compile_commands.json ${KRATOS_BUILD}/compile_commands.json
rm ${KRATOS_BUILD}/test
echo "-- Linking ${kratos_install_dir}/test with ${KRATOS_BUILD}/test"
ln -sf ${kratos_install_dir}/test ${KRATOS_BUILD}/test
