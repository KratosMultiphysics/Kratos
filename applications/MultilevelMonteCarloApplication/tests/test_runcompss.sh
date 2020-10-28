#!/bin/sh

script_dir=$(dirname "$(readlink -f "$0")")
cd "$script_dir"

runcompss --lang=python --python_interpreter=python3 \
    --pythonpath="${script_dir}/../xmc/classDefs_solverWrapper/problemDefs_KratosMultiphysics/poisson_2d" \
    test_xmcAlgorithm.py

runcompss --lang=python --python_interpreter=python3 test_momentEstimator_deterministic.py

runcompss --lang=python --python_interpreter=python3 parallel-only_tests.py
