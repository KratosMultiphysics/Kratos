#!/bin/bash

BASE_LD_LIBRARY_PATH=$LD_LIBRARY_PATH

PYTHONS=("35" "36" "37")

for PYTHON in  "${PYTHONS[@]}"
do
	echo starting build for python${PYTHON}

	cd /workspace/kratos/Kratos
	git clean -ffxd

	PYTHON_LOCATION=/opt/python/cp${PYTHON}-cp${PYTHON}m/bin/python

	cd cmake_build
	cp /workspace/kratos/Kratos/scripts/wheels/linux/configure.sh ./configure.sh
	chmod +x configure.sh
	./configure.sh $PYTHON_LOCATION

	make -j$1
	make install

	cd /workspace/kratos/Kratos
	export HASH=$(git show -s --format=%h) #used in version number
	export LD_LIBRARY_PATH=$(pwd)/libs:$BASE_LD_LIBRARY_PATH
	echo $LD_LIBRARY_PATH

	mkdir /workspace/wheel
	cp -r /workspace/kratos/Kratos/KratosMultiphysics /workspace/wheel
	mkdir /workspace/wheel/KratosMultiphysics/.libs
	cp /workspace/kratos/Kratos/libs/Kratos* /workspace/wheel/KratosMultiphysics/.libs

	cp /workspace/kratos/Kratos/scripts/wheels/linux/setup.py /workspace/wheel/setup.py
	cp /workspace/kratos/Kratos/scripts/wheels/linux/README.md /workspace/wheel/README.md
	cp /workspace/kratos/Kratos/scripts/wheels/linux/__init__.py /workspace/wheel/KratosMultiphysics/__init__.py


	cd /workspace/wheel

	$PYTHON_LOCATION setup.py bdist_wheel
	cd dist
	ls
	auditwheel repair *.whl
	mv wheelhouse/*manylinux2010_x86_64.whl /workspace/out

	echo finished build for python${PYTHON}
	echo cleaning up

	#cleanup
	cd /workspace
	rm -rf /workspace/wheel
	export LD_LIBRARY_PATH=$BASE_LD_LIBRARY_PATH

done

