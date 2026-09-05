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
#include <cmath>

// External includes
#include <benchmark/benchmark.h>

// Project includes
#include "spaces/ublas_space.h"
#include "spaces/eigen_space.h"
#include "linear_solvers/cg_solver.h"
#include "utilities/parallel_utilities.h"

// Side-by-side performance comparison of the sparse linear-algebra backends
// (uBLAS vs Eigen). Both space implementations are always compiled, so a
// single binary benchmarks both under identical compiler flags — the
// KRATOS_LINEAR_ALGEBRA_BACKEND option only selects which one the default
// space aliases point to, and this file names both spaces explicitly.
//
// Every benchmark is registered twice through BENCHMARK_TEMPLATE, once per
// space, on identical data. The system matrix is a synthetic banded
// diagonally-dominant (SPD) matrix mimicking a FEM stencil, built by writing
// the CSR arrays directly (exactly as the builder-and-solvers do).
//
// The parallel behavior follows the usual Kratos shared-memory settings: run
// with OMP_NUM_THREADS=1 for serial numbers and higher values for the
// threaded ones.

namespace Kratos
{

namespace
{

using UblasSparse = TUblasSparseSpace<double>;
using EigenSparse = TEigenSparseSpace<double>;

/// Each sparse family is paired with the dense space of its own backend, so
/// neither benchmark leg mixes uBLAS and Eigen types.
template <class TSpaceType>
struct PairedDenseSpace
{
    using Type = TUblasDenseSpace<double>;
};

template <>
struct PairedDenseSpace<EigenSparse>
{
    using Type = TEigenDenseSpace<double>;
};

constexpr std::size_t BandHalfWidth = 4; // 9 entries per interior row, FEM-stencil-like

template <class TSpaceType>
struct BandedSystem
{
    typename TSpaceType::MatrixType A;
    typename TSpaceType::VectorType x;
    typename TSpaceType::VectorType y;
};

/// Builds a banded, diagonally dominant matrix by writing the CSR arrays
/// directly (the same access pattern the block builder-and-solver uses).
template <class TSpaceType>
BandedSystem<TSpaceType> MakeBandedSystem(const std::size_t Size)
{
    BandedSystem<TSpaceType> system;

    // Count the nonzeros
    std::size_t nnz = 0;
    for (std::size_t i = 0; i < Size; ++i) {
        const std::size_t begin = (i < BandHalfWidth) ? 0 : i - BandHalfWidth;
        const std::size_t end = (i + BandHalfWidth + 1 > Size) ? Size : i + BandHalfWidth + 1;
        nnz += end - begin;
    }

    system.A = typename TSpaceType::MatrixType(Size, Size, nnz);

    auto* row_indices = system.A.index1_data().begin();
    auto* col_indices = system.A.index2_data().begin();
    auto* values = system.A.value_data().begin();

    std::size_t counter = 0;
    row_indices[0] = 0;
    for (std::size_t i = 0; i < Size; ++i) {
        const std::size_t begin = (i < BandHalfWidth) ? 0 : i - BandHalfWidth;
        const std::size_t end = (i + BandHalfWidth + 1 > Size) ? Size : i + BandHalfWidth + 1;
        for (std::size_t j = begin; j < end; ++j) {
            col_indices[counter] = j;
            values[counter] = (i == j) ? 2.0 * BandHalfWidth + 1.0 : -0.5;
            ++counter;
        }
        row_indices[i + 1] = counter;
    }
    system.A.set_filled(Size + 1, nnz);

    system.x = typename TSpaceType::VectorType(Size);
    system.y = typename TSpaceType::VectorType(Size);
    for (std::size_t i = 0; i < Size; ++i) {
        system.x[i] = 1.0 + 0.001 * static_cast<double>(i % 17);
        system.y[i] = 0.0;
    }

    return system;
}

} // namespace

// --- Sparse matrix kernels -------------------------------------------------

template <class TSpaceType>
static void BM_SpMV(benchmark::State& rState)
{
    auto system = MakeBandedSystem<TSpaceType>(rState.range(0));
    for (auto _ : rState) {
        TSpaceType::Mult(system.A, system.x, system.y);
        benchmark::DoNotOptimize(system.y[0]);
    }
}

template <class TSpaceType>
static void BM_TransposeSpMV(benchmark::State& rState)
{
    auto system = MakeBandedSystem<TSpaceType>(rState.range(0));
    for (auto _ : rState) {
        TSpaceType::TransposeMult(system.A, system.x, system.y);
        benchmark::DoNotOptimize(system.y[0]);
    }
}

template <class TSpaceType>
static void BM_MatrixFrobeniusNorm(benchmark::State& rState)
{
    auto system = MakeBandedSystem<TSpaceType>(rState.range(0));
    for (auto _ : rState) {
        benchmark::DoNotOptimize(TSpaceType::TwoNorm(system.A));
    }
}

template <class TSpaceType>
static void BM_SetToZeroMatrix(benchmark::State& rState)
{
    auto system = MakeBandedSystem<TSpaceType>(rState.range(0));
    for (auto _ : rState) {
        TSpaceType::SetToZero(system.A);
        benchmark::DoNotOptimize(system.A.value_data().begin());
    }
}

// Graph construction the way the block builder-and-solver performs it:
// (rows, cols, nnz) construction plus direct CSR array filling.
template <class TSpaceType>
static void BM_CSRConstruction(benchmark::State& rState)
{
    for (auto _ : rState) {
        auto system = MakeBandedSystem<TSpaceType>(rState.range(0));
        benchmark::DoNotOptimize(system.A.value_data().begin());
    }
}

// --- Vector kernels ----------------------------------------------------------

template <class TSpaceType>
static void BM_Dot(benchmark::State& rState)
{
    auto system = MakeBandedSystem<TSpaceType>(rState.range(0));
    TSpaceType::Mult(system.A, system.x, system.y);
    for (auto _ : rState) {
        benchmark::DoNotOptimize(TSpaceType::Dot(system.x, system.y));
    }
}

template <class TSpaceType>
static void BM_VectorTwoNorm(benchmark::State& rState)
{
    auto system = MakeBandedSystem<TSpaceType>(rState.range(0));
    for (auto _ : rState) {
        benchmark::DoNotOptimize(TSpaceType::TwoNorm(system.x));
    }
}

template <class TSpaceType>
static void BM_ScaleAndAdd(benchmark::State& rState)
{
    auto system = MakeBandedSystem<TSpaceType>(rState.range(0));
    typename TSpaceType::VectorType z(rState.range(0));
    for (auto _ : rState) {
        TSpaceType::ScaleAndAdd(1.5, system.x, -0.5, system.y, z); // z = 1.5 x - 0.5 y
        benchmark::DoNotOptimize(z[0]);
    }
}

template <class TSpaceType>
static void BM_UnaliasedAdd(benchmark::State& rState)
{
    auto system = MakeBandedSystem<TSpaceType>(rState.range(0));
    for (auto _ : rState) {
        TSpaceType::UnaliasedAdd(system.y, 0.001, system.x); // y += 0.001 x
        benchmark::DoNotOptimize(system.y[0]);
    }
}

// --- End-to-end iterative solve ---------------------------------------------
// CG chains SpMV, Dot and the vector updates through the space, so it is a
// representative aggregate of the backend performance. The matrix and the
// starting point are identical for both backends, so the iteration counts are
// identical too.

template <class TSpaceType>
static void BM_CGSolve(benchmark::State& rState)
{
    auto system = MakeBandedSystem<TSpaceType>(rState.range(0));

    typename TSpaceType::VectorType b(rState.range(0));
    TSpaceType::Mult(system.A, system.x, b); // manufactured solution: x

    for (auto _ : rState) {
        rState.PauseTiming();
        typename TSpaceType::VectorType solution(rState.range(0));
        TSpaceType::SetToZero(solution);
        typename TSpaceType::VectorType rhs = b;
        CGSolver<TSpaceType, typename PairedDenseSpace<TSpaceType>::Type> solver(1e-10, 500);
        rState.ResumeTiming();

        solver.Solve(system.A, solution, rhs);
        benchmark::DoNotOptimize(solution[0]);
    }
}

// --- Expression-level operations ---------------------------------------------
// The ublas idiom (prod, inner_prod, trans, sum, noalias, ...) spelled
// IDENTICALLY for both families: for the uBLAS types the calls resolve to the
// boost::numeric::ublas expression templates injected by ublas_interface.h,
// for the Eigen types to the compat operations of eigen_compat_operations.h.

struct UblasExprFamily
{
    using MatrixType = Matrix; // boost::numeric::ublas::matrix<double>
    using VectorType = Vector; // boost::numeric::ublas::vector<double>
};

struct EigenExprFamily
{
    using MatrixType = EigenMatrix<double>;
    using VectorType = EigenVector<double>;
};

namespace
{

template <class TFamily>
typename TFamily::MatrixType MakeDenseMatrix(const std::size_t Size)
{
    typename TFamily::MatrixType matrix(Size, Size);
    for (std::size_t i = 0; i < Size; ++i)
        for (std::size_t j = 0; j < Size; ++j)
            matrix(i, j) = 1.0 / (1.0 + static_cast<double>(i + j));
    return matrix;
}

template <class TFamily>
typename TFamily::VectorType MakeDenseVector(const std::size_t Size)
{
    typename TFamily::VectorType vector(Size);
    for (std::size_t i = 0; i < Size; ++i)
        vector[i] = 1.0 + 0.001 * static_cast<double>(i % 17);
    return vector;
}

} // namespace

template <class TFamily>
static void BM_ExprDenseProdMM(benchmark::State& rState)
{
    auto A = MakeDenseMatrix<TFamily>(rState.range(0));
    auto B = MakeDenseMatrix<TFamily>(rState.range(0));
    typename TFamily::MatrixType C(rState.range(0), rState.range(0));
    for (auto _ : rState) {
        noalias(C) = prod(A, B);
        benchmark::DoNotOptimize(C(0, 0));
    }
}

template <class TFamily>
static void BM_ExprDenseProdMV(benchmark::State& rState)
{
    auto A = MakeDenseMatrix<TFamily>(rState.range(0));
    auto x = MakeDenseVector<TFamily>(rState.range(0));
    typename TFamily::VectorType y(rState.range(0));
    for (auto _ : rState) {
        noalias(y) = prod(A, x);
        benchmark::DoNotOptimize(y[0]);
    }
}

template <class TFamily>
static void BM_ExprTransProdMM(benchmark::State& rState)
{
    auto A = MakeDenseMatrix<TFamily>(rState.range(0));
    auto B = MakeDenseMatrix<TFamily>(rState.range(0));
    typename TFamily::MatrixType C(rState.range(0), rState.range(0));
    for (auto _ : rState) {
        noalias(C) = prod(trans(A), B);
        benchmark::DoNotOptimize(C(0, 0));
    }
}

template <class TFamily>
static void BM_ExprTrans(benchmark::State& rState)
{
    auto A = MakeDenseMatrix<TFamily>(rState.range(0));
    typename TFamily::MatrixType C(rState.range(0), rState.range(0));
    for (auto _ : rState) {
        noalias(C) = trans(A);
        benchmark::DoNotOptimize(C(0, 0));
    }
}

template <class TFamily>
static void BM_ExprOuterProd(benchmark::State& rState)
{
    auto x = MakeDenseVector<TFamily>(rState.range(0));
    auto y = MakeDenseVector<TFamily>(rState.range(0));
    typename TFamily::MatrixType C(rState.range(0), rState.range(0));
    for (auto _ : rState) {
        noalias(C) = outer_prod(x, y);
        benchmark::DoNotOptimize(C(0, 0));
    }
}

template <class TFamily>
static void BM_ExprInnerProd(benchmark::State& rState)
{
    auto x = MakeDenseVector<TFamily>(rState.range(0));
    auto y = MakeDenseVector<TFamily>(rState.range(0));
    for (auto _ : rState) {
        benchmark::DoNotOptimize(inner_prod(x, y));
    }
}

template <class TFamily>
static void BM_ExprSum(benchmark::State& rState)
{
    auto x = MakeDenseVector<TFamily>(rState.range(0));
    for (auto _ : rState) {
        benchmark::DoNotOptimize(sum(x));
    }
}

template <class TFamily>
static void BM_ExprNorm2(benchmark::State& rState)
{
    auto x = MakeDenseVector<TFamily>(rState.range(0));
    for (auto _ : rState) {
        benchmark::DoNotOptimize(norm_2(x));
    }
}

template <class TFamily>
static void BM_ExprNorm1(benchmark::State& rState)
{
    auto x = MakeDenseVector<TFamily>(rState.range(0));
    for (auto _ : rState) {
        benchmark::DoNotOptimize(norm_1(x));
    }
}

template <class TFamily>
static void BM_ExprNormInf(benchmark::State& rState)
{
    auto x = MakeDenseVector<TFamily>(rState.range(0));
    for (auto _ : rState) {
        benchmark::DoNotOptimize(norm_inf(x));
    }
}

template <class TFamily>
static void BM_ExprAxpy(benchmark::State& rState)
{
    auto x = MakeDenseVector<TFamily>(rState.range(0));
    auto y = MakeDenseVector<TFamily>(rState.range(0));
    typename TFamily::VectorType z(rState.range(0));
    for (auto _ : rState) {
        noalias(z) = 2.0 * x + y;
        benchmark::DoNotOptimize(z[0]);
    }
}

// Sparse product spelled as an expression, noalias(y) = prod(A, x)
// (the space Mult is benchmarked separately above)
template <class TSpaceType>
static void BM_ExprSparseProdMV(benchmark::State& rState)
{
    auto system = MakeBandedSystem<TSpaceType>(rState.range(0));
    for (auto _ : rState) {
        noalias(system.y) = prod(system.A, system.x);
        benchmark::DoNotOptimize(system.y[0]);
    }
}

// --- Registration -------------------------------------------------------------

#define KRATOS_REGISTER_BACKEND_BENCHMARK(name)                                       \
    BENCHMARK_TEMPLATE(name, UblasSparse)->Name(#name "/ublas")->Arg(1<<14)->Arg(1<<20); \
    BENCHMARK_TEMPLATE(name, EigenSparse)->Name(#name "/eigen")->Arg(1<<14)->Arg(1<<20);

KRATOS_REGISTER_BACKEND_BENCHMARK(BM_SpMV)
KRATOS_REGISTER_BACKEND_BENCHMARK(BM_TransposeSpMV)
KRATOS_REGISTER_BACKEND_BENCHMARK(BM_MatrixFrobeniusNorm)
KRATOS_REGISTER_BACKEND_BENCHMARK(BM_SetToZeroMatrix)
KRATOS_REGISTER_BACKEND_BENCHMARK(BM_CSRConstruction)
KRATOS_REGISTER_BACKEND_BENCHMARK(BM_Dot)
KRATOS_REGISTER_BACKEND_BENCHMARK(BM_VectorTwoNorm)
KRATOS_REGISTER_BACKEND_BENCHMARK(BM_ScaleAndAdd)
KRATOS_REGISTER_BACKEND_BENCHMARK(BM_UnaliasedAdd)

BENCHMARK_TEMPLATE(BM_CGSolve, UblasSparse)->Name("BM_CGSolve/ublas")->Arg(1<<14)->Arg(1<<18)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(BM_CGSolve, EigenSparse)->Name("BM_CGSolve/eigen")->Arg(1<<14)->Arg(1<<18)->Unit(benchmark::kMillisecond);

#undef KRATOS_REGISTER_BACKEND_BENCHMARK

// Expression-level operations: dense matrix ops at element-local (16) and
// mid (128/512) sizes, vector ops at 4k and 1M entries
#define KRATOS_REGISTER_EXPR_MATRIX_BENCHMARK(name, ...)                                      \
    BENCHMARK_TEMPLATE(name, UblasExprFamily)->Name(#name "/ublas")->__VA_ARGS__;             \
    BENCHMARK_TEMPLATE(name, EigenExprFamily)->Name(#name "/eigen")->__VA_ARGS__;

KRATOS_REGISTER_EXPR_MATRIX_BENCHMARK(BM_ExprDenseProdMM, Arg(16)->Arg(128))
KRATOS_REGISTER_EXPR_MATRIX_BENCHMARK(BM_ExprTransProdMM, Arg(16)->Arg(128))
KRATOS_REGISTER_EXPR_MATRIX_BENCHMARK(BM_ExprDenseProdMV, Arg(16)->Arg(512))
KRATOS_REGISTER_EXPR_MATRIX_BENCHMARK(BM_ExprTrans, Arg(16)->Arg(512))
KRATOS_REGISTER_EXPR_MATRIX_BENCHMARK(BM_ExprOuterProd, Arg(16)->Arg(512))
KRATOS_REGISTER_EXPR_MATRIX_BENCHMARK(BM_ExprInnerProd, Arg(1<<12)->Arg(1<<20))
KRATOS_REGISTER_EXPR_MATRIX_BENCHMARK(BM_ExprSum, Arg(1<<12)->Arg(1<<20))
KRATOS_REGISTER_EXPR_MATRIX_BENCHMARK(BM_ExprNorm2, Arg(1<<12)->Arg(1<<20))
KRATOS_REGISTER_EXPR_MATRIX_BENCHMARK(BM_ExprNorm1, Arg(1<<12)->Arg(1<<20))
KRATOS_REGISTER_EXPR_MATRIX_BENCHMARK(BM_ExprNormInf, Arg(1<<12)->Arg(1<<20))
KRATOS_REGISTER_EXPR_MATRIX_BENCHMARK(BM_ExprAxpy, Arg(1<<12)->Arg(1<<20))

#undef KRATOS_REGISTER_EXPR_MATRIX_BENCHMARK

BENCHMARK_TEMPLATE(BM_ExprSparseProdMV, UblasSparse)->Name("BM_ExprSparseProdMV/ublas")->Arg(1<<14)->Arg(1<<17);
BENCHMARK_TEMPLATE(BM_ExprSparseProdMV, EigenSparse)->Name("BM_ExprSparseProdMV/eigen")->Arg(1<<14)->Arg(1<<17);

// --- Fixed-size (element-local) dense operations ------------------------------
// Side-by-side comparison of the backend-selected fixed-size types (array_1d,
// BoundedMatrix, BoundedVector — Eigen-backed when KRATOS_LINEAR_ALGEBRA_BACKEND
// is "eigen") against the boost::numeric::ublas bounded types named explicitly.
// Under the uBLAS backend both families resolve to the same implementations, so
// equal timings there are the expected sanity baseline; under the Eigen backend
// the "/backend" numbers measure the Eigen-backed types.
// The kernels are the idioms element CalculateLocalSystem code is made of.

template <std::size_t TSize>
struct BackendFixedFamily
{
    using VectorType = array_1d<double, TSize>;
    using MatrixType = BoundedMatrix<double, TSize, TSize>;
};

template <std::size_t TSize>
struct UblasFixedFamily
{
    using VectorType = boost::numeric::ublas::bounded_vector<double, TSize>;
    using MatrixType = boost::numeric::ublas::bounded_matrix<double, TSize, TSize>;
};

namespace
{

template <class TFamily>
typename TFamily::MatrixType MakeFixedMatrix()
{
    typename TFamily::MatrixType matrix;
    for (std::size_t i = 0; i < matrix.size1(); ++i)
        for (std::size_t j = 0; j < matrix.size2(); ++j)
            matrix(i, j) = (i == j) ? 2.0 + static_cast<double>(i) : 1.0 / (1.0 + static_cast<double>(i + j));
    return matrix;
}

template <class TFamily>
typename TFamily::VectorType MakeFixedVector()
{
    typename TFamily::VectorType vector;
    for (std::size_t i = 0; i < vector.size(); ++i)
        vector[i] = 1.0 + 0.001 * static_cast<double>(i % 17);
    return vector;
}

} // namespace

template <class TFamily>
static void BM_FixedProdMV(benchmark::State& rState)
{
    const auto R = MakeFixedMatrix<TFamily>();
    const auto x = MakeFixedVector<TFamily>();
    typename TFamily::VectorType y;
    for (auto _ : rState) {
        noalias(y) = prod(R, x);
        benchmark::DoNotOptimize(y[0]);
    }
}

template <class TFamily>
static void BM_FixedProdMM(benchmark::State& rState)
{
    const auto A = MakeFixedMatrix<TFamily>();
    const auto B = MakeFixedMatrix<TFamily>();
    typename TFamily::MatrixType C;
    for (auto _ : rState) {
        noalias(C) = prod(A, B);
        benchmark::DoNotOptimize(C(0, 0));
    }
}

template <class TFamily>
static void BM_FixedTransProdMM(benchmark::State& rState)
{
    const auto A = MakeFixedMatrix<TFamily>();
    const auto B = MakeFixedMatrix<TFamily>();
    typename TFamily::MatrixType C;
    for (auto _ : rState) {
        noalias(C) = prod(trans(A), B);
        benchmark::DoNotOptimize(C(0, 0));
    }
}

/// The local-axes rotation kernel of the coordinate transformation utilities:
/// two rotations plus an accumulation, chained through a temporary.
template <class TFamily>
static void BM_FixedRotationKernel(benchmark::State& rState)
{
    const auto R = MakeFixedMatrix<TFamily>();
    const auto x = MakeFixedVector<TFamily>();
    typename TFamily::VectorType tmp;
    typename TFamily::VectorType result = MakeFixedVector<TFamily>();
    for (auto _ : rState) {
        noalias(tmp) = prod(R, x);
        noalias(result) = prod(trans(R), tmp);
        benchmark::DoNotOptimize(result[0]);
    }
}

template <class TFamily>
static void BM_FixedAxpy(benchmark::State& rState)
{
    const auto x = MakeFixedVector<TFamily>();
    const auto y = MakeFixedVector<TFamily>();
    typename TFamily::VectorType z;
    for (auto _ : rState) {
        noalias(z) = 2.0 * x + y;
        benchmark::DoNotOptimize(z[0]);
    }
}

template <class TFamily>
static void BM_FixedInnerProd(benchmark::State& rState)
{
    const auto x = MakeFixedVector<TFamily>();
    const auto y = MakeFixedVector<TFamily>();
    for (auto _ : rState) {
        benchmark::DoNotOptimize(inner_prod(x, y));
    }
}

template <class TFamily>
static void BM_FixedNorm2(benchmark::State& rState)
{
    const auto x = MakeFixedVector<TFamily>();
    for (auto _ : rState) {
        benchmark::DoNotOptimize(norm_2(x));
    }
}

template <class TFamily>
static void BM_FixedOuterProd(benchmark::State& rState)
{
    const auto x = MakeFixedVector<TFamily>();
    const auto y = MakeFixedVector<TFamily>();
    typename TFamily::MatrixType C;
    for (auto _ : rState) {
        noalias(C) = outer_prod(x, y);
        benchmark::DoNotOptimize(C(0, 0));
    }
}

#define KRATOS_REGISTER_FIXED_BENCHMARK(name)                                                     \
    BENCHMARK_TEMPLATE(name, BackendFixedFamily<3>)->Name(#name "/backend/3");                    \
    BENCHMARK_TEMPLATE(name, UblasFixedFamily<3>)->Name(#name "/ublas/3");                        \
    BENCHMARK_TEMPLATE(name, BackendFixedFamily<9>)->Name(#name "/backend/9");                    \
    BENCHMARK_TEMPLATE(name, UblasFixedFamily<9>)->Name(#name "/ublas/9");

KRATOS_REGISTER_FIXED_BENCHMARK(BM_FixedProdMV)
KRATOS_REGISTER_FIXED_BENCHMARK(BM_FixedProdMM)
KRATOS_REGISTER_FIXED_BENCHMARK(BM_FixedTransProdMM)
KRATOS_REGISTER_FIXED_BENCHMARK(BM_FixedRotationKernel)
KRATOS_REGISTER_FIXED_BENCHMARK(BM_FixedAxpy)
KRATOS_REGISTER_FIXED_BENCHMARK(BM_FixedInnerProd)
KRATOS_REGISTER_FIXED_BENCHMARK(BM_FixedNorm2)
KRATOS_REGISTER_FIXED_BENCHMARK(BM_FixedOuterProd)

#undef KRATOS_REGISTER_FIXED_BENCHMARK

} // namespace Kratos

BENCHMARK_MAIN();
