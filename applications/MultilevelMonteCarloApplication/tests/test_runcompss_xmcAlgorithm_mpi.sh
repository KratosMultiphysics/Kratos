#!/usr/bin/env bash

path_to_test_folder=$(pwd)

path_to_replace="RectangularCylinder2D_81nodes"
new_path="$path_to_test_folder/test_xmcAlgorithm_mpi/problem_settings/RectangularCylinder2D_81nodes"

# set absolute path in Kratos parameters
sed -i "s|$path_to_replace|$new_path|g" "test_xmcAlgorithm_mpi/problem_settings/ProjectParametersRectangularCylinder2D_Fractional_MPI.json"

export OMP_NUM_THREADS=1
export computing_units_mlmc_execute_0=1;
export computing_procs_mlmc_execute_0=2;
runcompss \
    --lang=python \
    --cpu_affinity="disabled" \
    --python_interpreter=python3 \
    --pythonpath=$path_to_test_folder/test_xmcAlgorithm_mpi \
    ./test_xmcAlgorithm_mpi.py

# revert change in Kratos parameters
sed -i "s|$new_path|$path_to_replace|g" "test_xmcAlgorithm_mpi/problem_settings/ProjectParametersRectangularCylinder2D_Fractional_MPI.json"
