#!/usr/bin/env bash
# =============================================================================
# KaHIP Build Configuration Script
# =============================================================================
# Comprehensive build script for the KaHIP graph partitioning library.
# Supports serial (KaFFPa), parallel (ParHIP/kaffpaE), Python bindings,
# ILP solvers, 64-bit edge mode, and installation to custom prefix.
#
# Usage:
#   ./configure.sh [OPTIONS]
#
# Run from the build/ directory or any directory — it auto-locates the source.
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Locate source root (parent of this script's directory)
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${SCRIPT_DIR}"

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
BUILD_TYPE="Release"
NCORES=""
VERBOSE=0
DRY_RUN=0

# Feature flags (mirrors CMake options)
NOMPI=0
PARHIP=1
DETERMINISTIC_PARHIP=0
BUILDPYTHONMODULE=0
USE_ILP=0
USE_TCMALLOC=0
MODE_64BIT=0
NONATIVEOPTIMIZATIONS=0
OPTIMIZED_OUTPUT=0
DEPLOY=1
DEPLOY_DIR="${SOURCE_DIR}/deploy"
GENERATOR=""
EXTRA_CMAKE_ARGS=()

# ---------------------------------------------------------------------------
# Color output helpers
# ---------------------------------------------------------------------------
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
BLUE='\033[0;34m'; CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'

info()    { echo -e "${BLUE}[INFO]${RESET}  $*"; }
ok()      { echo -e "${GREEN}[OK]${RESET}    $*"; }
warn()    { echo -e "${YELLOW}[WARN]${RESET}  $*"; }
error()   { echo -e "${RED}[ERROR]${RESET} $*" >&2; }
section() { echo -e "\n${BOLD}${CYAN}>>> $*${RESET}"; }

# ---------------------------------------------------------------------------
# Usage / help
# ---------------------------------------------------------------------------
usage() {
cat <<EOF
${BOLD}KaHIP Build Configuration Script${RESET}

${BOLD}USAGE${RESET}
    $(basename "$0") [OPTIONS]

${BOLD}BUILD OPTIONS${RESET}
    --build-type TYPE       CMake build type: Release|Debug|RelWithDebInfo|MinSizeRel
                            Default: Release
    --jobs N                Number of parallel compile jobs (default: auto-detect)
    --generator NAME        CMake generator, e.g. "Ninja" or "Unix Makefiles"
    --no-deploy             Skip copying artifacts to deploy/ directory
    --deploy-dir PATH       Destination for deploy artifacts (default: <src>/deploy)
    --dry-run               Print cmake/make commands without executing them
    --verbose               Enable verbose make output (make VERBOSE=1)

${BOLD}FEATURE FLAGS${RESET}
    --no-mpi                Disable all MPI-dependent targets (ParHIP, kaffpaE)
    --no-parhip             Build MPI targets but skip ParHIP specifically
    --deterministic-parhip  Enforce deterministic output in ParHIP (may reduce quality)
    --python                Build Python (pybind11) bindings
    --ilp                   Build ILP improver and exact solver (requires Gurobi)
    --tcmalloc              Link against tcmalloc for faster memory allocation
    --64bit                 Enable 64-bit edge index mode (kahip_idx = int64_t)
    --no-native             Disable -march=native optimizations (for portable builds)
    --optimized-output      Enable KAFFPAOUTPUT optimized output mode

${BOLD}EXTRA CMAKE ARGUMENTS${RESET}
    Any unrecognized arguments are passed verbatim to cmake.
    Examples:
        $(basename "$0") -DGUROBI_HOME=/opt/gurobi
        $(basename "$0") -DCMAKE_CXX_COMPILER=g++-13

${BOLD}EXAMPLES${RESET}
    # Standard release build (with MPI/ParHIP if available)
    $(basename "$0")

    # Debug build without MPI
    $(basename "$0") --build-type Debug --no-mpi

    # Full-featured build with Python bindings and 64-bit edges
    $(basename "$0") --python --64bit --tcmalloc

    # Portable build (no native CPU optimizations), install to /usr/local
    $(basename "$0") --no-native --prefix /usr/local

    # Build with Ninja generator and 8 parallel jobs
    $(basename "$0") --generator Ninja --jobs 8

    # ILP support (requires Gurobi, set GUROBI_HOME)
    $(basename "$0") --ilp -DGUROBI_HOME=/opt/gurobi1000/linux64

    # Deterministic parallel partitioning
    $(basename "$0") --deterministic-parhip

    # CI-friendly: no native opts, no MPI, just the core library
    $(basename "$0") --no-mpi --no-native --no-deploy

EOF
}

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --build-type)       BUILD_TYPE="$2";            shift 2 ;;
        --jobs)             NCORES="$2";                shift 2 ;;
        --generator)        GENERATOR="$2";             shift 2 ;;
        --no-deploy)        DEPLOY=0;                   shift   ;;
        --deploy-dir)       DEPLOY_DIR="$2";            shift 2 ;;
        --dry-run)          DRY_RUN=1;                  shift   ;;
        --verbose)          VERBOSE=1;                  shift   ;;
        --no-mpi)           NOMPI=1;                    shift   ;;
        --no-parhip)        PARHIP=0;                   shift   ;;
        --deterministic-parhip) DETERMINISTIC_PARHIP=1; shift  ;;
        --python)           BUILDPYTHONMODULE=1;        shift   ;;
        --ilp)              USE_ILP=1;                  shift   ;;
        --tcmalloc)         USE_TCMALLOC=1;             shift   ;;
        --64bit)            MODE_64BIT=1;               shift   ;;
        --no-native)        NONATIVEOPTIMIZATIONS=1;    shift   ;;
        --optimized-output) OPTIMIZED_OUTPUT=1;         shift   ;;
        -h|--help)          usage; exit 0 ;;
        -D*)                EXTRA_CMAKE_ARGS+=("$1");   shift   ;;
        *)
            error "Unknown option: $1"
            echo "Run '$(basename "$0") --help' for usage."
            exit 1
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Auto-detect number of cores
# ---------------------------------------------------------------------------
if [[ -z "$NCORES" ]]; then
    case "$(uname -s)" in
        Darwin) NCORES=$(sysctl -n hw.logicalcpu 2>/dev/null || echo 4) ;;
        Linux)  NCORES=$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4) ;;
        *)      NCORES=4 ;;
    esac
fi

# ---------------------------------------------------------------------------
# Dependency checks
# ---------------------------------------------------------------------------
section "Checking dependencies"

check_cmd() {
    local cmd="$1" label="${2:-$1}"
    if command -v "$cmd" &>/dev/null; then
        ok "$label found: $(command -v "$cmd")"
        return 0
    else
        warn "$label not found"
        return 1
    fi
}

check_cmd cmake  "CMake"
check_cmd make   "Make"   || check_cmd ninja "Ninja" || warn "No make/ninja found — build may fail"
check_cmd ccache "ccache" && info "ccache detected, will be used by CMake automatically"

HAS_MPI=0
if [[ $NOMPI -eq 0 ]]; then
    if check_cmd mpicc "MPI C compiler" && check_cmd mpicxx "MPI C++ compiler"; then
        HAS_MPI=1
        MPI_VERSION=$(mpicc --version 2>&1 | head -1 || echo "unknown")
        info "MPI version: ${MPI_VERSION}"
    else
        warn "MPI not found. MPI-dependent targets (ParHIP, kaffpaE) will be skipped."
        warn "To suppress this warning, pass --no-mpi explicitly."
        HAS_MPI=0
    fi
fi

if [[ $BUILDPYTHONMODULE -eq 1 ]]; then
    check_cmd python3 "Python 3" || { error "Python 3 required for --python"; exit 1; }
    PYTHON_VERSION=$(python3 --version 2>&1)
    info "Python: $PYTHON_VERSION"
    if ! python3 -c "import pybind11" 2>/dev/null; then
        warn "pybind11 not found in Python path. Attempting to locate via pip..."
        if python3 -m pip show pybind11 &>/dev/null; then
            PYBIND11_DIR=$(python3 -m pip show pybind11 | grep Location | cut -d ' ' -f 2)/pybind11/share/cmake/pybind11
            export pybind11_DIR="$PYBIND11_DIR"
            info "pybind11 cmake dir: $PYBIND11_DIR"
        else
            error "pybind11 not installed. Run: pip install pybind11"
            exit 1
        fi
    fi
fi

if [[ $USE_ILP -eq 1 ]]; then
    if [[ -z "${GUROBI_HOME:-}" ]]; then
        warn "GUROBI_HOME is not set. Gurobi detection may fail."
        warn "Set GUROBI_HOME to your Gurobi installation directory."
    else
        info "Gurobi home: $GUROBI_HOME"
    fi
fi

if [[ $USE_TCMALLOC -eq 1 ]]; then
    if ldconfig -p 2>/dev/null | grep -q libtcmalloc || find /usr/lib /usr/local/lib -name "libtcmalloc*" 2>/dev/null | grep -q .; then
        ok "tcmalloc found"
    else
        warn "tcmalloc not found. Linking may succeed if CMake finds it elsewhere."
    fi
fi

# Check for optional Metis (enables fast_node_ordering)
if ldconfig -p 2>/dev/null | grep -q libmetis || find /usr/lib /usr/local/lib -name "libmetis*" 2>/dev/null | grep -q .; then
    ok "METIS found — fast_node_ordering will be built"
else
    info "METIS not found — fast_node_ordering will be skipped (this is optional)"
fi

# CMake version check
CMAKE_VERSION=$(cmake --version | head -1 | awk '{print $3}')
info "CMake version: $CMAKE_VERSION"
CMAKE_MAJOR=$(echo "$CMAKE_VERSION" | cut -d. -f1)
CMAKE_MINOR=$(echo "$CMAKE_VERSION" | cut -d. -f2)
if [[ $CMAKE_MAJOR -lt 3 ]] || [[ $CMAKE_MAJOR -eq 3 && $CMAKE_MINOR -lt 10 ]]; then
    error "CMake >= 3.10 required (found $CMAKE_VERSION)"
    exit 1
fi

# ---------------------------------------------------------------------------
# Build cmake argument list
# ---------------------------------------------------------------------------
section "Configuring CMake arguments"

CMAKE_ARGS=(
    "-DCMAKE_BUILD_TYPE=${BUILD_TYPE}"
    "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON"
)
[[ -n "$GENERATOR" ]]         && CMAKE_ARGS+=("-G${GENERATOR}")
[[ $NOMPI -eq 1 ]]            && CMAKE_ARGS+=("-DNOMPI=ON")
[[ $PARHIP -eq 0 ]]           && CMAKE_ARGS+=("-DPARHIP=OFF")
[[ $DETERMINISTIC_PARHIP -eq 1 ]] && CMAKE_ARGS+=("-DDETERMINISTIC_PARHIP=ON")
[[ $BUILDPYTHONMODULE -eq 1 ]] && CMAKE_ARGS+=("-DBUILDPYTHONMODULE=ON")
[[ $USE_ILP -eq 1 ]]          && CMAKE_ARGS+=("-DUSE_ILP=ON")
[[ $USE_TCMALLOC -eq 1 ]]     && CMAKE_ARGS+=("-DUSE_TCMALLOC=ON")
[[ $MODE_64BIT -eq 1 ]]       && CMAKE_ARGS+=("-D64BITMODE=ON")
[[ $NONATIVEOPTIMIZATIONS -eq 1 ]] && CMAKE_ARGS+=("-DNONATIVEOPTIMIZATIONS=ON")
[[ $OPTIMIZED_OUTPUT -eq 1 ]] && CMAKE_ARGS+=("-DOPTIMIZED_OUTPUT=ON")

# Gurobi hint
[[ -n "${GUROBI_HOME:-}" ]] && CMAKE_ARGS+=("-DGUROBI_HOME=${GUROBI_HOME}")

# pybind11 dir if we found it above
[[ -n "${pybind11_DIR:-}" ]] && CMAKE_ARGS+=("-Dpybind11_DIR=${pybind11_DIR}")

# Pass through any extra -D args
CMAKE_ARGS+=("${EXTRA_CMAKE_ARGS[@]}")

# Print summary
echo ""
info "Source directory : $SOURCE_DIR"
info "Build directory  : $BUILD_DIR"
info "Build type       : $BUILD_TYPE"
info "Parallel jobs    : $NCORES"
[[ $DEPLOY -eq 1 ]] && info "Deploy directory : $DEPLOY_DIR"
echo ""
info "CMake options:"
for arg in "${CMAKE_ARGS[@]}"; do
    echo "    $arg"
done

# ---------------------------------------------------------------------------
# Helper: run or print command
# ---------------------------------------------------------------------------
run() {
    if [[ $DRY_RUN -eq 1 ]]; then
        echo "[DRY-RUN] $*"
    else
        "$@"
    fi
}

# ---------------------------------------------------------------------------
# Configure
# ---------------------------------------------------------------------------
section "Running cmake configure"

CMAKE_CMD=(cmake "${CMAKE_ARGS[@]}" "$SOURCE_DIR")
info "Command: ${CMAKE_CMD[*]}"
echo ""

cd "$BUILD_DIR"
run "${CMAKE_CMD[@]}"

# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------
section "Building KaHIP (${NCORES} jobs)"

if [[ $VERBOSE -eq 1 ]]; then
    run cmake --build . --config "$BUILD_TYPE" --parallel "$NCORES" -- VERBOSE=1
else
    run cmake --build . --config "$BUILD_TYPE" --parallel "$NCORES"
fi

ok "Build completed successfully."

# ---------------------------------------------------------------------------
# Deploy artifacts
# ---------------------------------------------------------------------------
if [[ $DEPLOY -eq 1 && $DRY_RUN -eq 0 ]]; then
    section "Deploying artifacts to $DEPLOY_DIR"

    rm -rf "$DEPLOY_DIR"
    mkdir -p "$DEPLOY_DIR"

    # ---- Serial executables ----
    echo ""
    info "Copying serial executables..."
    SERIAL_BINS=(
        kaffpa
        global_multisection
        evaluator
        edge_evaluator
        graphchecker
        node_separator
        node_ordering
        partition_to_vertex_separator
        label_propagation
        edge_partitioning
        interface_test
    )
    for bin in "${SERIAL_BINS[@]}"; do
        if [[ -f "${BUILD_DIR}/${bin}" ]]; then
            cp "${BUILD_DIR}/${bin}" "${DEPLOY_DIR}/"
            ok "  ${bin}"
        fi
    done

    # Optional serial executables (may not be built depending on options)
    OPTIONAL_BINS=(kaffpaE fast_node_ordering ilp_improve ilp_exact)
    for bin in "${OPTIONAL_BINS[@]}"; do
        if [[ -f "${BUILD_DIR}/${bin}" ]]; then
            cp "${BUILD_DIR}/${bin}" "${DEPLOY_DIR}/"
            ok "  ${bin} (optional)"
        fi
    done

    # ---- Libraries ----
    echo ""
    info "Copying libraries..."

    # Shared library
    for lib in "${BUILD_DIR}"/libkahip.so* "${BUILD_DIR}"/libkahip.dylib "${BUILD_DIR}"/kahip.dll; do
        [[ -f "$lib" ]] && cp "$lib" "${DEPLOY_DIR}/" && ok "  $(basename "$lib")"
    done

    # Static library
    if [[ -f "${BUILD_DIR}/libkahip_static.a" ]]; then
        cp "${BUILD_DIR}/libkahip_static.a" "${DEPLOY_DIR}/libkahip.a"
        ok "  libkahip.a (static)"
    fi

    # pkg-config
    if [[ -f "${BUILD_DIR}/kahip.pc" ]]; then
        cp "${BUILD_DIR}/kahip.pc" "${DEPLOY_DIR}/"
        ok "  kahip.pc"
    fi

    # ---- Headers ----
    echo ""
    info "Copying headers..."
    cp "${SOURCE_DIR}/interface/kaHIP_interface.h" "${DEPLOY_DIR}/"
    ok "  kaHIP_interface.h"

    # ---- ParHIP parallel artifacts ----
    echo ""
    PARHIP_SRC_DIR="${BUILD_DIR}/parallel/parallel_src"
    if [[ -d "$PARHIP_SRC_DIR" ]]; then
        info "Copying ParHIP parallel artifacts..."
        mkdir -p "${DEPLOY_DIR}/parallel"

        # ParHIP executables
        [[ -f "${PARHIP_SRC_DIR}/parhip" ]] && cp "${PARHIP_SRC_DIR}/parhip" "${DEPLOY_DIR}/parhip" && ok "  parhip"
        for f in "${PARHIP_SRC_DIR}"/toolbox*; do
            [[ -f "$f" ]] && cp "$f" "${DEPLOY_DIR}/toolbox" && ok "  toolbox"
        done
        for f in "${PARHIP_SRC_DIR}"/graph2binary*; do
            [[ -f "$f" && ! "$f" == *_external* ]] && cp "$f" "${DEPLOY_DIR}/graph2binary" && ok "  graph2binary"
        done
        for f in "${PARHIP_SRC_DIR}"/graph2binary_external*; do
            [[ -f "$f" ]] && cp "$f" "${DEPLOY_DIR}/graph2binary_external" && ok "  graph2binary_external"
        done
        for f in "${PARHIP_SRC_DIR}"/dsp*; do
            [[ -f "$f" ]] && cp "$f" "${DEPLOY_DIR}/distributed_edge_partitioning" && ok "  distributed_edge_partitioning"
        done
        for f in "${PARHIP_SRC_DIR}"/readbgf*; do
            [[ -f "$f" ]] && cp "$f" "${DEPLOY_DIR}/readbgf" && ok "  readbgf"
        done

        # ParHIP library — keep the libparhip_interface* name, FindKaHIP.cmake
        # searches for NAMES parhip_interface / parhip_interface_static.
        for lib in "${PARHIP_SRC_DIR}"/libparhip_interface*.a "${PARHIP_SRC_DIR}"/libparhip_interface*.so*; do
            if [[ -f "$lib" ]]; then
                cp "$lib" "${DEPLOY_DIR}/" && ok "  $(basename "$lib")"
            fi
        done

        # Modified KaHIP (used internally by ParHIP)
        MODIFIED_KAHIP_DIR="${BUILD_DIR}/parallel/modified_kahip"
        if [[ -d "$MODIFIED_KAHIP_DIR" ]]; then
            for lib in "${MODIFIED_KAHIP_DIR}"/lib*.a; do
                [[ -f "$lib" ]] && cp "$lib" "${DEPLOY_DIR}/parallel/libkahip.a" && ok "  parallel/libkahip.a"
            done
        fi

        # ParHIP header
        cp "${SOURCE_DIR}/parallel/parallel_src/interface/parhip_interface.h" "${DEPLOY_DIR}/"
        ok "  parhip_interface.h"

        # ParHIP pkg-config
        for pc in "${PARHIP_SRC_DIR}"/*.pc; do
            [[ -f "$pc" ]] && cp "$pc" "${DEPLOY_DIR}/" && ok "  $(basename "$pc")"
        done
    else
        info "ParHIP was not built — skipping parallel artifacts"
    fi

    # ---- Python module ----
    if [[ $BUILDPYTHONMODULE -eq 1 ]]; then
        echo ""
        info "Copying Python module artifacts..."
        for f in "${BUILD_DIR}"/kahip*.so "${BUILD_DIR}"/kahip*.pyd; do
            [[ -f "$f" ]] && cp "$f" "${DEPLOY_DIR}/" && ok "  $(basename "$f")"
        done
        cp "${SOURCE_DIR}/misc/pymodule/callkahipfrompython.py" "${DEPLOY_DIR}/" 2>/dev/null && ok "  callkahipfrompython.py"
    fi

    # ---- compile_commands.json (for IDE/clangd support) ----
    if [[ -f "${BUILD_DIR}/compile_commands.json" ]]; then
        cp "${BUILD_DIR}/compile_commands.json" "${DEPLOY_DIR}/"
        ok "  compile_commands.json"
    fi

    echo ""
    ok "Deployed artifacts:"
    echo "==============================================="
    (cd "${DEPLOY_DIR}" && ls -1 --color=auto 2>/dev/null || ls -1)
    echo "==============================================="
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
section "Build Summary"

echo ""
printf "  %-30s %s\n" "Source:"          "$SOURCE_DIR"
printf "  %-30s %s\n" "Build:"           "$BUILD_DIR"
printf "  %-30s %s\n" "Build type:"      "$BUILD_TYPE"
printf "  %-30s %s\n" "64-bit edges:"    "$([ $MODE_64BIT -eq 1 ] && echo yes || echo no)"
printf "  %-30s %s\n" "MPI/ParHIP:"      "$([ $NOMPI -eq 1 ] && echo disabled || ([ $HAS_MPI -eq 1 ] && echo enabled || echo 'not found'))"
printf "  %-30s %s\n" "Python bindings:" "$([ $BUILDPYTHONMODULE -eq 1 ] && echo yes || echo no)"
printf "  %-30s %s\n" "ILP (Gurobi):"    "$([ $USE_ILP -eq 1 ] && echo yes || echo no)"
printf "  %-30s %s\n" "tcmalloc:"        "$([ $USE_TCMALLOC -eq 1 ] && echo yes || echo no)"
printf "  %-30s %s\n" "Native opts:"     "$([ $NONATIVEOPTIMIZATIONS -eq 1 ] && echo disabled || echo enabled)"
[[ $DEPLOY -eq 1 && $DRY_RUN -eq 0 ]] && printf "  %-30s %s\n" "Deployed to:" "$DEPLOY_DIR"
echo ""

ok "KaHIP build complete."
