#!/bin/bash
PYTHONS=("cp36" "cp37" "cp38" "cp39" "cp310")
export KRATOS_VERSION="9.0.0"

BASE_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
export KRATOS_ROOT="/workspace/kratos/Kratos"
WHEEL_ROOT="/workspace/wheel"
WHEEL_OUT="/data_swap_guest"
CORE_LIB_DIR="/workspace/coreLibs"

setup_wheel_dir () {
    cd $KRATOS_ROOT
    mkdir $WHEEL_ROOT
    cp scripts/wheels/setup.py ${WHEEL_ROOT}/setup.py
    mkdir ${WHEEL_ROOT}/KratosMultiphysics
    mkdir ${WHEEL_ROOT}/KratosMultiphysics/.libs
}

build_core_wheel () {
    setup_wheel_dir
    cd $KRATOS_ROOT

    PREFIX_LOCATION=$1
    
    cp ${PREFIX_LOCATION}/KratosMultiphysics/* ${WHEEL_ROOT}/KratosMultiphysics
    cp scripts/wheels/linux/KratosMultiphysics.json ${WHEEL_ROOT}/wheel.json
    cp scripts/wheels/__init__.py ${WHEEL_ROOT}/KratosMultiphysics/__init__.py
    
    cd $WHEEL_ROOT
    $PYTHON_LOCATION setup.py bdist_wheel
    cd ${WHEEL_ROOT}/dist
    auditwheel repair *.whl
    
    mkdir $CORE_LIB_DIR
    unzip -j wheelhouse/KratosMultiphysics* 'KratosMultiphysics/.libs/*' -d $CORE_LIB_DIR
    
    cp wheelhouse/* ${WHEEL_OUT}/
    cd
    rm -r $WHEEL_ROOT
}

build_application_wheel () {
    setup_wheel_dir
    cp ${KRATOS_ROOT}/scripts/wheels/linux/applications/${1} ${WHEEL_ROOT}/wheel.json
    cd $WHEEL_ROOT
    $PYTHON_LOCATION setup.py bdist_wheel
    cd dist
    auditwheel repair *.whl
    
    optimize_wheel
    cp ${WHEEL_ROOT}/dist/wheelhouse/* ${WHEEL_OUT}/
    
    cd
    rm -r $WHEEL_ROOT
}

build_kratos_all_wheel () {
    setup_wheel_dir
    cp ${KRATOS_ROOT}/scripts/wheels/linux/KratosMultiphysics-all.json ${WHEEL_ROOT}/wheel.json
    cp ${KRATOS_ROOT}/scripts/wheels/linux/setup_kratos_all.py ${WHEEL_ROOT}/setup.py
    cd ${WHEEL_ROOT}
    $PYTHON_LOCATION setup.py bdist_wheel
    cp dist/* ${WHEEL_OUT}/

    cd
    rm -r $WHEEL_ROOT
}

optimize_wheel(){
    cd ${WHEEL_ROOT}/dist/wheelhouse
    ARCHIVE_NAME=$(ls .)
    mkdir tmp
    unzip ${ARCHIVE_NAME} -d tmp
    rm $ARCHIVE_NAME
    
    for LIBRARY in $(ls tmp/KratosMultiphysics/.libs)
    do
        if [ -f "${CORE_LIB_DIR}/${LIBRARY}" ] || grep -Fxq $(echo $LIBRARY | cut -f1 -d"-") "${WHEEL_ROOT}/excluded.txt" ; then
            echo "removing ${LIBRARY} - already present in dependent wheel."
            rm tmp/KratosMultiphysics/.libs/${LIBRARY}
            sed -i "/${LIBRARY}/d" tmp/*.dist-info/RECORD
        fi
    done
    cd tmp
    zip -r ../${ARCHIVE_NAME} ./*
    cd ..
    rm -r tmp
}

build_core () {
	cd $KRATOS_ROOT
	# git clean -ffxd

	PYTHON_LOCATION=$1
    PREFIX_LOCATION=$2

	cp /workspace/kratos/Kratos/scripts/wheels/linux/configure.sh ./configure.sh
	chmod +x configure.sh
	./configure.sh $PYTHON_LOCATION $PREFIX_LOCATION

    cmake --build "${KRATOS_ROOT}/build/Release" --target KratosKernel -- -j$(nproc)
}

build_interface () {
    cd $KRATOS_ROOT
	# git clean -ffxd

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
echo starting core build
build_core python3.6 ${KRATOS_ROOT}/bin/core
echo finished core build

for PYTHON_VERSION in  "${PYTHONS[@]}"
do
    PYTHON_TMP=$(ls /opt/python | grep $PYTHON_VERSION | cut -d "-" -f 2)
    export PYTHON=${PYTHON_TMP#cp}
    echo starting build for python${PYTHON_VERSION}

	PYTHON_LOCATION=/opt/python/$(ls /opt/python | grep $PYTHON_VERSION)/bin/python
    PREFIX_LOCATION=$KRATOS_ROOT/bin/Release/python_$PYTHON

    build_interface $PYTHON_LOCATION $PREFIX_LOCATION
	
	cd $KRATOS_ROOT
	export HASH=$(git show -s --format=%h) # Used in version number
	export LD_LIBRARY_PATH=${PREFIX_LOCATION}/libs:$BASE_LD_LIBRARY_PATH
	echo $LD_LIBRARY_PATH

    build_core_wheel $PREFIX_LOCATION

    for APPLICATION in $(ls ${KRATOS_ROOT}/scripts/wheels/linux/applications)
    do
        build_application_wheel $APPLICATION
    done

    build_kratos_all_wheel $PREFIX_LOCATION

	echo finished build for python${PYTHON_VERSION}

	export LD_LIBRARY_PATH=$BASE_LD_LIBRARY_PATH

done

