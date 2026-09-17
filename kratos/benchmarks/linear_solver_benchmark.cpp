//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <fstream>
#include <filesystem>
#include <iostream>

// External includes
#include <benchmark/benchmark.h>

// Project includes
#include "includes/kratos_application.h"
#include "includes/kratos_parameters.h"
#include "includes/matrix_market_interface.h"
#include "factories/linear_solver_factory.h"
#include "spaces/ublas_space.h"

namespace Kratos
{

// Type definitions
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType  = UblasSpace<double, Matrix, Vector>;
using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;
using SparseMatrixType = typename SparseSpaceType::MatrixType;
using VectorType       = typename SparseSpaceType::VectorType;

/// Global file path variables (set via command-line arguments; defaults provided)
static std::string g_lhs_file = "A.mm";
static std::string g_rhs_file = "b.mm.rhs";
static std::string g_solver_config_file = "solver_settings.json";

/// Default AMGCL solver settings used as fallback when no config file is found
static constexpr const char* g_default_solver_settings = R"({
    "solver_type"                    : "AMGCL",
    "preconditioner_type"            : "amg",
    "smoother_type"                  : "ilu0",
    "krylov_type"                    : "gmres",
    "coarsening_type"                : "aggregation",
    "max_iteration"                  : 100,
    "gmres_krylov_space_dimension"   : 100,
    "verbosity"                      : 0,
    "tolerance"                      : 1e-6,
    "scaling"                        : false,
    "block_size"                     : 1,
    "use_block_matrices_if_possible" : true
})";

/**
 * @brief Builds a small synthetic diagonally-dominant SPD system of size @p N.
 * @details Used as fallback when no matrix/vector files are provided.
 * @param N System size
 * @param rA Output sparse matrix
 * @param rB Output RHS vector
 */
void BuildSyntheticSystem(std::size_t N, SparseMatrixType& rA, VectorType& rB)
{
    rA.resize(N, N, false);
    rA.clear();
    rB.resize(N, false);
    // Tridiagonal: -1  2  -1  (scaled so diagonal dominates)
    for (std::size_t i = 0; i < N; ++i) {
        if (i > 0)     rA(i, i - 1) = -1.0;
        rA(i, i)       = 4.0;
        if (i < N - 1) rA(i, i + 1) = -1.0;
        rB[i] = 1.0;
    }
}

/**
 * @brief Reads the LHS matrix from a Matrix Market file, or builds a synthetic system if the file is absent.
 * @param rFileName The path to the .mm file (may be empty or non-existent for synthetic fallback)
 * @param rA The sparse matrix to be filled
 * @param rB RHS vector (only populated when using the synthetic fallback)
 * @return true if a real file was loaded, false if the synthetic fallback was used
 */
bool ReadLHS(const std::string& rFileName, SparseMatrixType& rA, VectorType& rB)
{
    if (!rFileName.empty() && std::filesystem::exists(rFileName)) {
        const bool success = ReadMatrixMarketMatrix(rFileName.c_str(), rA);
        KRATOS_ERROR_IF_NOT(success) << "Failed to read LHS matrix from: " << rFileName << std::endl;
        return true;
    }
    std::cout << "[WARNING] LHS file not found (\"" << rFileName << "\"). Using synthetic 1000x1000 tridiagonal system.\n";
    BuildSyntheticSystem(1000, rA, rB);
    return false;
}

/**
 * @brief Reads the RHS vector from a Matrix Market file, or uses the already-populated @p rB from the synthetic system.
 * @param rFileName The path to the .mm file (may be empty or non-existent for synthetic fallback)
 * @param rB The vector to be filled
 * @param rUseSynthetic If true, @p rB is already filled and this call is a no-op
 */
void ReadRHS(const std::string& rFileName, VectorType& rB, bool rUseSynthetic)
{
    if (rUseSynthetic) return;  // synthetic RHS already built alongside LHS
    KRATOS_ERROR_IF(!std::filesystem::exists(rFileName)) << "RHS file does not exist: " << rFileName << std::endl;
    const bool success = ReadMatrixMarketVector(rFileName.c_str(), rB);
    KRATOS_ERROR_IF_NOT(success) << "Failed to read RHS vector from: " << rFileName << std::endl;
}

/**
 * @brief Creates a linear solver from a JSON configuration file, or from built-in AMGCL defaults if absent.
 * @param rFileName The path to the .json file
 * @return A shared pointer to the created linear solver
 */
LinearSolverType::Pointer CreateSolverFromJSON(const std::string& rFileName)
{
    Parameters solver_params;
    if (!rFileName.empty() && std::filesystem::exists(rFileName)) {
        std::ifstream config_file(rFileName);
        KRATOS_ERROR_IF_NOT(config_file.is_open()) << "Could not open config file: " << rFileName << std::endl;
        solver_params = Parameters(config_file);
    } else {
        std::cout << "[WARNING] Config file not found (\"" << rFileName << "\"). Using built-in AMGCL defaults.\n";
        solver_params = Parameters(g_default_solver_settings);
    }

    using LinearSolverFactoryType = LinearSolverFactory<SparseSpaceType, LocalSpaceType>;
    return LinearSolverFactoryType().Create(solver_params);
}

/**
 * @brief Benchmark for the full Solve method (Initialize + SolveStep + Finalize + Clear)
 */
static void BM_LinearSolverSolve(benchmark::State& state)
{
    // Register Kratos application (needed for solver factories)
    KratosApplication application("KratosApplication");
    application.Register();

    // Read LHS and RHS (fallback to synthetic system if files absent)
    SparseMatrixType A;
    VectorType b;
    const bool synthetic = !ReadLHS(g_lhs_file, A, b);
    ReadRHS(g_rhs_file, b, synthetic);

    // Create solver (fallback to AMGCL defaults if config file absent)
    auto p_solver = CreateSolverFromJSON(g_solver_config_file);

    const auto system_size = SparseSpaceType::Size(b);
    KRATOS_ERROR_IF(SparseSpaceType::Size1(A) != system_size || SparseSpaceType::Size2(A) != system_size)
        << "Dimension mismatch: A is " << SparseSpaceType::Size1(A) << "x" << SparseSpaceType::Size2(A)
        << " but b has size " << system_size << std::endl;

    for (auto _ : state) {
        state.PauseTiming();
        // Create a fresh copy of the matrix and RHS for each iteration
        // (some solvers modify the matrix during factorization)
        SparseMatrixType A_copy(A);
        VectorType b_copy(b);
        VectorType x = ZeroVector(system_size);
        state.ResumeTiming();

        p_solver->Solve(A_copy, x, b_copy);
    }

    // Report custom counters
    state.counters["SystemSize"] = static_cast<double>(system_size);
    state.counters["NNZ"] = static_cast<double>(A.nnz());
}

/**
 * @brief Benchmark for the InitializeSolutionStep method only (e.g., factorization)
 */
static void BM_LinearSolverInitializeSolutionStep(benchmark::State& state)
{
    // Register Kratos application
    KratosApplication application("KratosApplication");
    application.Register();

    // Read LHS and RHS (fallback to synthetic system if files absent)
    SparseMatrixType A;
    VectorType b;
    const bool synthetic = !ReadLHS(g_lhs_file, A, b);
    ReadRHS(g_rhs_file, b, synthetic);

    // Create solver (fallback to AMGCL defaults if config file absent)
    auto p_solver = CreateSolverFromJSON(g_solver_config_file);

    const auto system_size = SparseSpaceType::Size(b);
    VectorType x = ZeroVector(system_size);

    // Initialize once
    p_solver->Initialize(A, x, b);

    for (auto _ : state) {
        state.PauseTiming();
        SparseMatrixType A_copy(A);
        VectorType b_copy(b);
        VectorType x_copy = ZeroVector(system_size);
        state.ResumeTiming();

        p_solver->InitializeSolutionStep(A_copy, x_copy, b_copy);
    }

    state.counters["SystemSize"] = static_cast<double>(system_size);
    state.counters["NNZ"] = static_cast<double>(A.nnz());
}

/**
 * @brief Benchmark for the PerformSolutionStep method only (solve after factorization)
 */
static void BM_LinearSolverPerformSolutionStep(benchmark::State& state)
{
    // Register Kratos application
    KratosApplication application("KratosApplication");
    application.Register();

    // Read LHS and RHS (fallback to synthetic system if files absent)
    SparseMatrixType A;
    VectorType b;
    const bool synthetic = !ReadLHS(g_lhs_file, A, b);
    ReadRHS(g_rhs_file, b, synthetic);

    // Create solver (fallback to AMGCL defaults if config file absent)
    auto p_solver = CreateSolverFromJSON(g_solver_config_file);

    const auto system_size = SparseSpaceType::Size(b);
    VectorType x = ZeroVector(system_size);

    // Initialize and prepare the solver (factorize once)
    p_solver->Initialize(A, x, b);
    p_solver->InitializeSolutionStep(A, x, b);

    for (auto _ : state) {
        state.PauseTiming();
        VectorType b_copy(b);
        VectorType x_copy = ZeroVector(system_size);
        state.ResumeTiming();

        p_solver->PerformSolutionStep(A, x_copy, b_copy);
    }

    state.counters["SystemSize"] = static_cast<double>(system_size);
    state.counters["NNZ"] = static_cast<double>(A.nnz());
}

// Register benchmarks
BENCHMARK(BM_LinearSolverSolve)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_LinearSolverInitializeSolutionStep)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_LinearSolverPerformSolutionStep)->Unit(benchmark::kMillisecond);

}  // namespace Kratos

/**
 * @brief Custom main to parse extra command-line arguments for file paths
 * @details Usage:
 *   ./linear_solver_benchmark --lhs=matrix.mm --rhs=vector.mm --config=solver.json [google benchmark flags]
 */
int main(int argc, char** argv)
{
    // Parse custom arguments before passing to Google Benchmark
    std::vector<char*> remaining_args;
    remaining_args.push_back(argv[0]);

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg.rfind("--lhs=", 0) == 0) {
            Kratos::g_lhs_file = arg.substr(6);
        } else if (arg.rfind("--rhs=", 0) == 0) {
            Kratos::g_rhs_file = arg.substr(6);
        } else if (arg.rfind("--config=", 0) == 0) {
            Kratos::g_solver_config_file = arg.substr(9);
        } else {
            remaining_args.push_back(argv[i]);
        }
    }

    // Print usage/help if explicitly requested
    bool show_help = false;
    for (auto* arg : remaining_args) {
        const std::string s(arg);
        if (s == "--help" || s == "-h") { show_help = true; break; }
    }
    if (show_help) {
        std::cout << "Usage: " << argv[0] << " [--lhs=<matrix.mm>] [--rhs=<vector.mm>] [--config=<solver.json>] [benchmark flags]\n"
                  << "\nOptional arguments (defaults shown):\n"
                  << "  --lhs=<file.mm>       Path to the LHS sparse matrix in Matrix Market format  [A.mm]\n"
                  << "  --rhs=<file.mm>       Path to the RHS vector in Matrix Market format          [b.mm]\n"
                  << "  --config=<file.json>  Path to the JSON file with linear solver configuration  [solver_settings.json]\n"
                  << "\nExample solver configuration (solver_settings.json):\n"
                  << "  {\n"
                  << "      \"solver_type\": \"amgcl\"\n"
                  << "  }\n"
                  << "\nAll Google Benchmark flags are also supported (e.g., --benchmark_format=json).\n";
        return 0;
    }

    // Validate that the required files are specified (they may be empty only if the user explicitly cleared them)
    if (Kratos::g_lhs_file.empty() || Kratos::g_rhs_file.empty() || Kratos::g_solver_config_file.empty()) {
        std::cerr << "Error: one or more required file paths are empty.\n"
                  << "Run with --help for usage information.\n";
        return 1;
    }

    // Print configuration
    std::cout << "Linear Solver Benchmark Configuration:\n"
              << "  LHS matrix:      " << Kratos::g_lhs_file << "\n"
              << "  RHS vector:      " << Kratos::g_rhs_file << "\n"
              << "  Solver config:   " << Kratos::g_solver_config_file << "\n\n";

    // Pass remaining args to Google Benchmark
    int new_argc = static_cast<int>(remaining_args.size());
    benchmark::Initialize(&new_argc, remaining_args.data());
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
