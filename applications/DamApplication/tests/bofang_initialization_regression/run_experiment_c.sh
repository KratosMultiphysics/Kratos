#!/usr/bin/env bash
# Runner for variant C (historical Kratos at ece5cfe...) in the pre_13472 worktree.
set -euo pipefail

WORKTREE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
EXPERIMENT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="${EXPERIMENT_DIR}/case"
INSTRU_DIR="${EXPERIMENT_DIR}/instrumentation"
OUT_DIR="${EXPERIMENT_DIR}/results/pre_13472"
KRATOS_BIN="${WORKTREE_DIR}/bin/Release"

MONITORED_NODES="1,4,7,10,13,16,17,18"

export LD_LIBRARY_PATH="${KRATOS_BIN}/libs:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${KRATOS_BIN}:${CASE_DIR}:${EXPERIMENT_DIR}:${PYTHONPATH:-}"

mkdir -p "${OUT_DIR}"
(cd "${WORKTREE_DIR}" && git rev-parse HEAD) > "${OUT_DIR}/HEAD.txt"

cd "${CASE_DIR}"
python3 "${INSTRU_DIR}/instrumented_dam_analysis.py" \
    ProjectParameters.json "${OUT_DIR}" "C" "${MONITORED_NODES}" legacy

echo "Variant C completed. Results in ${OUT_DIR}"
