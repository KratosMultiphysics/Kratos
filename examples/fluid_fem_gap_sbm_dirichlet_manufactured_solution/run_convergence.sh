#!/bin/bash
set -euo pipefail

KRATOS_EXAMPLE_DIRECTORY="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
KRATOS_REPOSITORY_ROOT="$(cd "${KRATOS_EXAMPLE_DIRECTORY}/../.." && pwd -P)"
KRATOS_PYTHON_EXECUTABLE="${PYTHON_EXECUTABLE:-/Library/Frameworks/Python.framework/Versions/3.10/bin/python3}"

export OMP_NUM_THREADS=4
export OMP_THREAD_LIMIT=4
export PYTHONPATH="${KRATOS_REPOSITORY_ROOT}/bin/Release${PYTHONPATH:+:${PYTHONPATH}}"
export DYLD_LIBRARY_PATH="${KRATOS_REPOSITORY_ROOT}/bin/Release/libs${DYLD_LIBRARY_PATH:+:${DYLD_LIBRARY_PATH}}"

exec "${KRATOS_PYTHON_EXECUTABLE}" "${KRATOS_EXAMPLE_DIRECTORY}/run_convergence.py" "$@"
