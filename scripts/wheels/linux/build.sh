#!/bin/bash
PYTHONS=("38" "39" "310" "311" "312" "313")
export KRATOS_VERSION="10.2.3"

BASE_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
export KRATOS_ROOT="/workspace/kratos/Kratos"
WHEEL_ROOT="/workspace/wheel"
WHEEL_OUT="/data_swap_guest"
CORE_LIB_DIR="/workspace/coreLibs"

# Created the wheel building directory.
setup_wheel_dir () {
    cd $KRATOS_ROOT
    mkdir $WHEEL_ROOT
    cp scripts/wheels/setup.py ${WHEEL_ROOT}/setup.py
    mkdir ${WHEEL_ROOT}/KratosMultiphysics
    mkdir ${WHEEL_ROOT}/KratosMultiphysics/.libs        # This dir is necessary mainly to store the shared libs in windows
}

# Creates the wheel for KratosCore.
build_core_wheel () {
    setup_wheel_dir
    cd $KRATOS_ROOT

    PREFIX_LOCATION=$1

    mkdir ${WHEEL_ROOT}/KratosMultiphysics

    cp ${PREFIX_LOCATION}/KratosMultiphysics/*       ${WHEEL_ROOT}/KratosMultiphysics
    cp ${KRATOS_ROOT}/kratos/KratosMultiphysics.json ${WHEEL_ROOT}/wheel.json

    cd $WHEEL_ROOT

    $PYTHON_LOCATION setup.py bdist_wheel

    cd ${WHEEL_ROOT}/dist

    auditwheel repair *.whl

    mkdir $CORE_LIB_DIR
    unzip -j wheelhouse/KratosMultiphysics* 'KratosMultiphysics.libs/*' -d $CORE_LIB_DIR

    cp wheelhouse/* ${WHEEL_OUT}/
    cd
    rm -r $WHEEL_ROOT
}

# Creates the wheel for each KratosApplication.
build_application_wheel () {
    setup_wheel_dir

    cp ${KRATOS_ROOT}/applications/${1}/${1}.json ${WHEEL_ROOT}/wheel.json
    cd $WHEEL_ROOT

    $PYTHON_LOCATION setup.py bdist_wheel

    auditwheel repair dist/*.whl

    optimize_wheel

    cp ${WHEEL_ROOT}/wheelhouse/* ${WHEEL_OUT}/

    cd
    rm -r $WHEEL_ROOT
}

# Kreates the wheel bundle.
build_kratos_all_wheel () {
    setup_wheel_dir
    cp ${KRATOS_ROOT}/kratos/KratosMultiphysics-all.json ${WHEEL_ROOT}/wheel.json
    cp ${KRATOS_ROOT}/scripts/wheels/linux/setup_kratos_all.py ${WHEEL_ROOT}/setup.py
    cd ${WHEEL_ROOT}
    $PYTHON_LOCATION setup.py bdist_wheel
    cp dist/* ${WHEEL_OUT}/

    cd
    rm -r $WHEEL_ROOT
}

# Removes duplicated libraries from existing wheels.
optimize_wheel(){

    cd ${WHEEL_ROOT}/wheelhouse
    ARCHIVE_NAME=$(ls .)
    mkdir tmp
    unzip ${ARCHIVE_NAME} -d tmp
    rm $ARCHIVE_NAME

    echo "Begin exclude list for ${APPNAME}"
    echo "List of core libs to be excluded:"
    echo "$(ls ${CORE_LIB_DIR})"

    # Clean excluded libraries
    for LIBRARY in $(ls tmp/Kratos${APPNAME}.libs)
    do
        if [ -f "${CORE_LIB_DIR}/${LIBRARY}" ] || grep $(echo $LIBRARY | cut -f1 -d"-") "${WHEEL_ROOT}/excluded.txt" ; then
            echo "-- Removing ${LIBRARY} - already present in dependent wheel."

            rm tmp/Kratos${APPNAME}.libs/${LIBRARY}     # Try to remove from app dir

            sed -i "/${LIBRARY}/d" tmp/*.dist-info/RECORD
        else
            echo "-- Keeping ${LIBRARY}"
        fi
    done

    # Alson clean the possible copies done in the setup.py
    rm tmp/KratosMultiphysics/.libs/libKratos*

    cd tmp
    zip -r ../${ARCHIVE_NAME} ./*
    cd ..
    rm -r tmp
}

# Buils the KratosXCore components for the kernel and applications
build_core () {
	cd $KRATOS_ROOT

	PYTHON_LOCATION=$1
    PREFIX_LOCATION=$2

	cp /workspace/kratos/Kratos/scripts/wheels/linux/configure.sh ./configure.sh
	chmod +x configure.sh
	./configure.sh $PYTHON_LOCATION $PREFIX_LOCATION

    cmake --build "${KRATOS_ROOT}/build/Release" --target KratosKernel -- -j$(nproc)
}

# Buils the KratosXInterface components for the kernel and applications given an specific version of python
build_interface () {
    cd $KRATOS_ROOT

	PYTHON_LOCATION=$1
    PREFIX_LOCATION=$2

	cp /workspace/kratos/Kratos/scripts/wheels/linux/configure.sh ./configure.sh
	chmod +x configure.sh
	./configure.sh $PYTHON_LOCATION $PREFIX_LOCATION

    cmake --build "${KRATOS_ROOT}/build/Release" --target KratosPythonInterface -- -j$(nproc)
    cmake --build "${KRATOS_ROOT}/build/Release" --target install -- -j$(nproc)
}

# Core can be build independently of the python version.
# Install path should be useless here.
echo "Starting core build"
build_core python3.8 ${KRATOS_ROOT}/bin/core
echo "Finished core build"

for PYTHON_VERSION in  "${PYTHONS[@]}"
do
    PYTHON_TMP=$(ls /opt/python | grep $PYTHON_VERSION | cut -d "-" -f 2)
    export PYTHON=${PYTHON_TMP#cp}
    echo "Starting build for python${PYTHON_VERSION}"

	PYTHON_LOCATION=/opt/python/$(ls /opt/python | grep $PYTHON_VERSION)/bin/python
    PREFIX_LOCATION=$KRATOS_ROOT/bin/Release/python_$PYTHON

    $PYTHON_LOCATION -m pip install mypy

    build_interface $PYTHON_LOCATION $PREFIX_LOCATION

	cd $KRATOS_ROOT
	export HASH=$(git show -s --format=%h) # Used in version number
	export LD_LIBRARY_PATH=${PREFIX_LOCATION}/libs:$BASE_LD_LIBRARY_PATH
	echo $LD_LIBRARY_PATH

    echo "Building Core Wheel"
    build_core_wheel $PREFIX_LOCATION

    echo "Building App Wheels"
    for APPLICATION in $(ls -d ${PREFIX_LOCATION}/applications/*)
    do
        APPNAME=$(basename "$APPLICATION")
        echo "Building ${APPNAME} Wheel"
        build_application_wheel $APPNAME
    done

    echo "Building Bundle Wheel"
    build_kratos_all_wheel $PREFIX_LOCATION

	echo finished build for python${PYTHON_VERSION}

	export LD_LIBRARY_PATH=$BASE_LD_LIBRARY_PATH

done

