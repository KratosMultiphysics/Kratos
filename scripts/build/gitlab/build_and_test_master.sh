#!/bin/bash

set -e

cp scripts/build/gitlab/configure_gitlab_trusty.sh cmake_build/configure.sh

export PYTHONPATH=${PYTHONPATH}:${CI_PROJECT_DIR}/bin/Linux/Custom
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CI_PROJECT_DIR}/bin/Linux/Custom/libs

cd cmake_build
CMAKE_BUILD_TYPE=Custom
CCACHE_COMPILERCHECK=content
bash configure.sh
make -j8
make -j8 runkratos
make install/fast
echo Build complete

echo Running tests
cd ../kratos/python_scripts
python3 run_tests.py -l small -c python3
