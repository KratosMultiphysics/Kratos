#!/usr/bin/env bash

path_to_test_folder=$(pwd)

runcompss \
    --lang=python \
    --python_interpreter=python3 \
    --pythonpath=$path_to_test_folder/../xmc/classDefs_solverWrapper/problemDefs_KratosMultiphysics/poisson_2d \
    ./test_xmcAlgorithm.py
