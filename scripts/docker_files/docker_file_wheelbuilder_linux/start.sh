#!/bin/bash

BRANCH=ci/updating-wheel-builder

cd /workspace/kratos
git clone --depth 1 --single-branch -b $BRANCH https://github.com/KratosMultiphysics/Kratos.git

cd /workspace/kratos/Kratos/scripts/wheels/linux/
chmod +x build.sh
./build.sh