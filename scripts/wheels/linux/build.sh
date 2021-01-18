#!/bin/bash
PYTHONS=("35" "36" "37" "38")
export KRATOS_VERSION="7.0.3"

BASE_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
export KRATOS_ROOT="/workspace/kratos/Kratos"
WHEEL_ROOT="/workspace/wheel"
WHEEL_OUT="/workspace/out"
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
    
    cp bin/Release/KratosMultiphysics/* ${WHEEL_ROOT}/KratosMultiphysics
    cp scripts/wheels/linux/KratosMultiphysics.json ${WHEEL_ROOT}/wheel.json
    cp scripts/wheels/__init__.py ${WHEEL_ROOT}/KratosMultiphysics/__init__.py
    
    cd $WHEEL_ROOT
    $PYTHON_LOCATION setup.py bdist_wheel
    cd ${WHEEL_ROOT}/dist
    auditwheel repair *.whl
    
    mkdir $CORE_LIB_DIR
    unzip -j wheelhouse/KratosMultiphysics* 'KratosMultiphysics/.libs/*' -d $CORE_LIB_DIR
    
    cp wheelhouse/* $WHEEL_OUT
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

build_application_wheel () {
    setup_wheel_dir
    cp ${KRATOS_ROOT}/scripts/wheels/linux/applications/${1} ${WHEEL_ROOT}/wheel.json
    cd $WHEEL_ROOT
    $PYTHON_LOCATION setup.py bdist_wheel
    cd dist
    auditwheel repair *.whl
    
    optimize_wheel
    cp ${WHEEL_ROOT}/dist/wheelhouse/* ${WHEEL_OUT}
    
    cd
    rm -r $WHEEL_ROOT
}

build_kratos_all_wheel () {
    setup_wheel_dir
    cp ${KRATOS_ROOT}/scripts/wheels/linux/KratosMultiphysics-all.json ${WHEEL_ROOT}/wheel.json
    cp ${KRATOS_ROOT}/scripts/wheels/linux/setup_kratos_all.py ${WHEEL_ROOT}/setup.py
    cd ${WHEEL_ROOT}
    $PYTHON_LOCATION setup.py bdist_wheel
    cp dist/* ${WHEEL_OUT}

    cd
    rm -r $WHEEL_ROOT
}

build () {
	cd $KRATOS_ROOT
	git clean -ffxd

	PYTHON_LOCATION=$1

	cp /workspace/kratos/Kratos/scripts/wheels/linux/configure.sh ./configure.sh
	chmod +x configure.sh
	./configure.sh $PYTHON_LOCATION

	cmake --build "${KRATOS_ROOT}/build/Release" --target install -- -j$2
}


for PYTHON_VERSION in  "${PYTHONS[@]}"
do
    PYTHON_TMP=$(ls /opt/python | grep $PYTHON_VERSION | cut -d "-" -f 2)
    export PYTHON=${PYTHON_TMP#cp}
    echo starting build for python${PYTHON}
	PYTHON_LOCATION=/opt/python/$(ls /opt/python | grep $PYTHON_VERSION)/bin/python
    build $PYTHON_LOCATION $1
	
	
	cd $KRATOS_ROOT
	export HASH=$(git show -s --format=%h) #used in version number
	export LD_LIBRARY_PATH=${KRATOS_ROOT}/bin/Release/libs:$BASE_LD_LIBRARY_PATH
	echo $LD_LIBRARY_PATH

    build_core_wheel

    for APPLICATION in $(ls ${KRATOS_ROOT}/scripts/wheels/linux/applications)
    do
        build_application_wheel $APPLICATION
    done

    build_kratos_all_wheel

	echo finished build for python${PYTHON_VERSION}

	export LD_LIBRARY_PATH=$BASE_LD_LIBRARY_PATH

done

