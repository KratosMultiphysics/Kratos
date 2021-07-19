#!/bin/bash
export KRATOS_VERSION="9.0.0.dev1"

cpus=$1
pythons=$2
useCotire=$3

BASE_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
export KRATOS_ROOT="/workspace/kratos/Kratos"
wheelRoot="/workspace/wheel"
wheelOut="/workspace/out"
coreLibDir="/workspace/coreLibs"

setup_wheel_dir() {
  cd $KRATOS_ROOT
  mkdir $wheelRoot
  cp scripts/wheels/setup.py ${wheelRoot}/setup.py
  mkdir ${wheelRoot}/KratosMultiphysics
  mkdir ${wheelRoot}/KratosMultiphysics/.libs
}

build_core_wheel() {
  setup_wheel_dir
  cd $KRATOS_ROOT

  cp bin/Release/KratosMultiphysics/* ${wheelRoot}/KratosMultiphysics
  cp scripts/wheels/linux/KratosMultiphysics.json ${wheelRoot}/wheel.json
  /bin/cp ${KRATOS_ROOT}/scripts/wheels/__init__.py ${wheelRoot}/KratosMultiphysics/__init__.py

  cd $wheelRoot
  $pythonLocation setup.py bdist_wheel
  cd ${wheelRoot}/dist
  auditwheel repair *.whl

  mkdir $coreLibDir
  unzip -j wheelhouse/KratosMultiphysics* 'KratosMultiphysics/.libs/*' -d $coreLibDir

  cp wheelhouse/* $wheelOut
  cd
  rm -r $wheelRoot
}

optimize_wheel() {
  cd ${wheelRoot}/dist/wheelhouse
  archiveName=$(ls .)
  mkdir tmp
  unzip "$archiveName" -d tmp
  rm "$archiveName"

  for library in $(ls tmp/KratosMultiphysics/.libs); do
    if [ -f "${coreLibDir}/${library}" ] || grep -Fxq $(echo "$library" | cut -f1 -d"-") "${wheelRoot}/excluded.txt"; then
      echo "removing ${library} - already present in dependent wheel."
      rm tmp/KratosMultiphysics/.libs/"$library"
      sed -i "/${library}/d" tmp/*.dist-info/RECORD
    fi
  done
  cd tmp
  zip -r ../"$archiveName" ./*
  cd ..
  rm -r tmp
}

build_application_wheel() {
  setup_wheel_dir
  cp ${KRATOS_ROOT}/scripts/wheels/linux/applications/"$1" ${wheelRoot}/wheel.json
  cd $wheelRoot
  $pythonLocation setup.py bdist_wheel
  cd dist
  auditwheel repair *.whl

  optimize_wheel
  cp ${wheelRoot}/dist/wheelhouse/* ${wheelOut}

  cd
  rm -r $wheelRoot
}

build_kratos_all_wheel() {
  setup_wheel_dir
  cp ${KRATOS_ROOT}/scripts/wheels/linux/KratosMultiphysics-all.json ${wheelRoot}/wheel.json
  cp ${KRATOS_ROOT}/scripts/wheels/linux/setup_kratos_all.py ${wheelRoot}/setup.py
  cd ${wheelRoot}
  $pythonLocation setup.py bdist_wheel
  cp dist/* ${wheelOut}

  cd
  rm -r $wheelRoot
}

build() {
  cd $KRATOS_ROOT
  git clean -ffxd

  pythonLocation=$1

  cp /workspace/kratos/Kratos/scripts/wheels/linux/configure.sh ./configure.sh
  chmod +x configure.sh
  ./configure.sh "$pythonLocation" "$useCotire"

  if [[ $useCotire == "ON" ]];
  then
    cmake --build "${KRATOS_ROOT}/build/Release" --target all_unity -- -j1 && \
    cmake --build "${KRATOS_ROOT}/build/Release" --target zlibstatic -- -j1 && \
    cmake --build "${KRATOS_ROOT}/build/Release" --target install/fast -- -j1
  else
    cmake --build "${KRATOS_ROOT}/build/Release" --target install -- -j"$cpus"
  fi

}

for pythonVersion in "${pythons[@]}"; do
  PYTHON_TMP=$(ls /opt/python | grep "$pythonVersion" | cut -d "-" -f 2)
  export PYTHON=${PYTHON_TMP#cp}
  echo starting build for python"$PYTHON"
  pythonLocation=/opt/python/$(ls /opt/python | grep "$pythonVersion")/bin/python
  build "$pythonLocation"

  cd $KRATOS_ROOT
  export LD_LIBRARY_PATH=${KRATOS_ROOT}/bin/Release/libs:$BASE_LD_LIBRARY_PATH
  echo "$LD_LIBRARY_PATH"

  build_core_wheel

  for application in $(ls ${KRATOS_ROOT}/scripts/wheels/linux/applications); do
    if [[ $application != _* ]];
    then
      echo - starting build for "$application"
      build_application_wheel "$application"
    fi

  done

  build_kratos_all_wheel

  echo finished build for python"$pythonVersion"

  export LD_LIBRARY_PATH=$BASE_LD_LIBRARY_PATH

done
