#!/bin/sh

script_dir=$(dirname "$(readlink -f "$0")")
cd "$script_dir"

./test_runcompss_xmcAlgorithm.sh

runcompss --lang=python --python_interpreter=python3 test_momentEstimator_deterministic.py

runcompss --lang=python --python_interpreter=python3 parallel-only_tests.py

./test_runcompss_xmcAlgorithm_mpi.sh
