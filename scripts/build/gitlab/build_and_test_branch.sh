#!/bin/bash

set -e

echo Loading cache
cp -R /buildcache/* ./cmake_build
cp scripts/build/gitlab/configure_gitlab_trusty.sh cmake_build/configure.sh

export PYTHONPATH=${PYTHONPATH}:${CI_PROJECT_DIR}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CI_PROJECT_DIR}/libs

cd cmake_build
CCACHE_COMPILERCHECK=content
bash configure.sh
make -j8
make -j8 runkratos
make install/fast
echo Build complete

echo Running tests
cd ../kratos/python_scripts
python3 run_tests.py -l small -c python3
