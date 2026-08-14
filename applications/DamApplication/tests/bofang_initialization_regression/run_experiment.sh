#!/usr/bin/env bash
# Runner for the Bofang initialization-regression experiment.
#
# Variants:
#   A : current master (no code change)
#   B : current master + Bofang process reverted to ExecuteInitialize()
#   B2: current master + Bofang process reverted to ExecuteInitialize()
#       + dam_analysis.py pre-solver-initialize process call reverted to
#       ExecuteInitialize()  (faithful reproduction of the pre-#13472 lifecycle)
#   C : historical Kratos at ece5cfe... (run separately from its own worktree)
#
# Usage:
#   run_experiment.sh <variant> <kratos_bin_dir> <patch_dir>
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
CASE_DIR="${SCRIPT_DIR}/case"
INSTRU_DIR="${SCRIPT_DIR}/instrumentation"
RESULTS_DIR="${SCRIPT_DIR}/results"

VARIANT="$1"
KRATOS_BIN="${2:-${REPO_ROOT}/bin/Release}"
BUILD_DIR="${3:-${REPO_ROOT}/build/master}"
PATCH_DIR="${4:-/tmp/opencode/bofang_patches}"

MONITORED_NODES="1,4,7,10,13,16,17,18"

case "$VARIANT" in
  A)
    OUT_DIR="${RESULTS_DIR}/current_master"
    LEGACY_FLAG=""
    PATCH_FILE=""
    ;;
  B)
    OUT_DIR="${RESULTS_DIR}/current_master_legacy_bofang"
    LEGACY_FLAG=""
    PATCH_FILE="${PATCH_DIR}/bofang_legacy_process_only.patch"
    ;;
  B2)
    OUT_DIR="${RESULTS_DIR}/current_master_legacy_bofang_faithful"
    LEGACY_FLAG="legacy"
    PATCH_FILE="${PATCH_DIR}/bofang_legacy_faithful.patch"
    ;;
  *)
    echo "Unknown variant: ${VARIANT}"
    exit 1
    ;;
esac

if [ -n "${PATCH_FILE}" ]; then
    (cd "${REPO_ROOT}" && git apply "${PATCH_FILE}")
    echo "Applied patch ${PATCH_FILE}"
    # rebuild the DamApplication pybind module (the Bofang process is header-only
    # and is instantiated by the python interface)
    cmake --build "${BUILD_DIR}" --target KratosDamApplication install -- -j8
fi

export LD_LIBRARY_PATH="${KRATOS_BIN}/libs:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${KRATOS_BIN}:${CASE_DIR}:${SCRIPT_DIR}:${PYTHONPATH:-}"

mkdir -p "${OUT_DIR}"

# record build info
(cd "${REPO_ROOT}" && git rev-parse HEAD) > "${OUT_DIR}/HEAD.txt" 2>/dev/null || true
(cd "${REPO_ROOT}" && git diff --stat) > "${OUT_DIR}/git_diff.txt" 2>/dev/null || true

cd "${CASE_DIR}"
python3 "${INSTRU_DIR}/instrumented_dam_analysis.py" \
    ProjectParameters.json "${OUT_DIR}" "${VARIANT}" "${MONITORED_NODES}" ${LEGACY_FLAG}

# restore the tree for B variants
if [ -n "${PATCH_FILE}" ]; then
    (cd "${REPO_ROOT}" && git apply -R "${PATCH_FILE}")
    echo "Reverted patch ${PATCH_FILE}"
    # rebuild the pybind module back to the pristine (A) state
    cmake --build "${BUILD_DIR}" --target KratosDamApplication install -- -j8
fi

echo "Variant ${VARIANT} completed. Results in ${OUT_DIR}"
