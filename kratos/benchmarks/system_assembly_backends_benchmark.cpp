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
#include <map>
#include <memory>

// External includes
#include <benchmark/benchmark.h>

// Project includes
#include "containers/model.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_2d_4.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "linear_solvers/linear_solver.h"
#include "processes/structured_mesh_generator_process.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "spaces/default_spaces.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

// Full-mesh system build/assembly benchmark of the linear-algebra backend
// selected at configure time (uBLAS or Eigen): where
// linear_algebra_backends_benchmark.cpp measures the individual operations,
// this file measures something closer to a real case — the sparse graph
// construction and the scattered CSR assembly of a complete FE mesh through
// the actual ResidualBasedBlockBuilderAndSolver, with a real element
// (DistanceCalculationElementSimplex, Poisson step) providing the local
// contributions. The two backends are mutually exclusive in one build, so
// comparing them means running the binary of each build; the benchmark names
// carry the backend.
//
// The parallel behavior follows the usual Kratos shared-memory settings: run
// with OMP_NUM_THREADS=1 for serial numbers and higher values for the
// threaded ones (the assembly scatter is parallel with atomic updates).

namespace Kratos
{

namespace
{

#ifdef KRATOS_USE_EIGEN_BACKEND
#define KRATOS_BACKEND_NAME "eigen"
#else
#define KRATOS_BACKEND_NAME "ublas"
#endif

using DefaultSparse = TDefaultSparseSpace<double>;
using DenseSpace = TDefaultDenseSpace<double>;

/// Builds (once per mesh size) a structured simplicial
/// mesh of DistanceCalculationElementSimplex elements on the unit
/// square/cube, with the DISTANCE DOF added and smooth signed nodal DISTANCE
/// values, ready for a FRACTIONAL_STEP = 1 (Poisson) build.
ModelPart& GetBenchmarkModelPart(const std::size_t Dimension, const std::size_t NumberOfDivisions)
{
    // Deliberately leaked, like the model cache below: the application's
    // teardown must not run after the ParallelEnvironment singleton is gone.
    static const bool registered = [] {
        (new KratosApplication("KratosApplication"))->Register();
        return true;
    }();
    (void)registered;

    // Deliberately leaked: a ModelPart's teardown talks to the
    // ParallelEnvironment singleton, which is destroyed before function-local
    // statics at exit.
    static std::map<std::pair<std::size_t, std::size_t>, Model*> model_cache;

    auto& p_model = model_cache[{Dimension, NumberOfDivisions}];
    if (p_model != nullptr) {
        return p_model->GetModelPart("AssemblyBenchmark");
    }

    p_model = new Model();
    auto& r_model_part = p_model->CreateModelPart("AssemblyBenchmark");
    r_model_part.SetBufferSize(1);
    r_model_part.AddNodalSolutionStepVariable(DISTANCE);

    Parameters mesher_parameters(R"({
        "number_of_divisions"        : 0,
        "element_name"               : "",
        "condition_name"             : "LineCondition2D2N",
        "create_skin_sub_model_part" : false
    })");
    mesher_parameters["number_of_divisions"].SetInt(NumberOfDivisions);

    if (Dimension == 2) {
        mesher_parameters["element_name"].SetString("DistanceCalculationElementSimplex2D3N");
        auto p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
        auto p_point_2 = Kratos::make_intrusive<Node>(2, 0.0, 1.0, 0.0);
        auto p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
        auto p_point_4 = Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0);
        Quadrilateral2D4<Node> domain_geometry(p_point_1, p_point_2, p_point_3, p_point_4);
        StructuredMeshGeneratorProcess(domain_geometry, r_model_part, mesher_parameters).Execute();
    } else {
        mesher_parameters["element_name"].SetString("DistanceCalculationElementSimplex3D4N");
        mesher_parameters["condition_name"].SetString("SurfaceCondition3D3N");
        auto p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
        auto p_point_2 = Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0);
        auto p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
        auto p_point_4 = Kratos::make_intrusive<Node>(4, 0.0, 1.0, 0.0);
        auto p_point_5 = Kratos::make_intrusive<Node>(5, 0.0, 0.0, 1.0);
        auto p_point_6 = Kratos::make_intrusive<Node>(6, 1.0, 0.0, 1.0);
        auto p_point_7 = Kratos::make_intrusive<Node>(7, 1.0, 1.0, 1.0);
        auto p_point_8 = Kratos::make_intrusive<Node>(8, 0.0, 1.0, 1.0);
        Hexahedra3D8<Node> domain_geometry(p_point_1, p_point_2, p_point_3, p_point_4,
                                           p_point_5, p_point_6, p_point_7, p_point_8);
        StructuredMeshGeneratorProcess(domain_geometry, r_model_part, mesher_parameters).Execute();
    }

    VariableUtils().AddDof(DISTANCE, r_model_part);

    // A smooth signed level-set-like field, so the element's source term takes
    // both signs across the mesh.
    block_for_each(r_model_part.Nodes(), [](Node& rNode) {
        rNode.FastGetSolutionStepValue(DISTANCE) =
            std::sin(4.0 * rNode.X()) * std::cos(4.0 * rNode.Y()) * std::cos(4.0 * rNode.Z()) + 0.1;
    });

    // Poisson step of the variational redistancing element.
    r_model_part.GetProcessInfo().SetValue(FRACTIONAL_STEP, 1);
    r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, static_cast<int>(Dimension));

    return r_model_part;
}

/// The per-leg builder machinery: a block builder-and-solver and static
/// scheme instantiated over the given sparse space (dense space fixed to the
/// compiled backend's, see the file note), with the DOF set and equation
/// numbering prepared.
template <class TSparseSpace>
struct AssemblyFixture
{
    using LinearSolverType = LinearSolver<TSparseSpace, DenseSpace>;
    using BuilderAndSolverType = ResidualBasedBlockBuilderAndSolver<TSparseSpace, DenseSpace, LinearSolverType>;
    using SchemeType = ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, DenseSpace>;

    explicit AssemblyFixture(ModelPart& rModelPart)
        : mrModelPart(rModelPart),
          mpScheme(std::make_shared<SchemeType>())
    {
        mBuilderAndSolver.SetDofSetIsInitializedFlag(false);
        mBuilderAndSolver.SetUpDofSet(mpScheme, mrModelPart);
        mBuilderAndSolver.SetUpSystem(mrModelPart);
    }

    /// Allocates and graph-constructs the system, then runs the standard
    /// initialization sequence so Build/BuildRHS can be called repeatedly.
    void InitializeSystem()
    {
        mBuilderAndSolver.ResizeAndInitializeVectors(mpScheme, mpA, mpDx, mpb, mrModelPart);
        mBuilderAndSolver.InitializeSolutionStep(mrModelPart, *mpA, *mpDx, *mpb);
        mpScheme->Initialize(mrModelPart);
        mpScheme->InitializeSolutionStep(mrModelPart, *mpA, *mpDx, *mpb);
        mpScheme->InitializeNonLinIteration(mrModelPart, *mpA, *mpDx, *mpb);
    }

    ModelPart& mrModelPart;
    typename SchemeType::Pointer mpScheme;
    BuilderAndSolverType mBuilderAndSolver{};
    typename TSparseSpace::MatrixPointerType mpA;
    typename TSparseSpace::VectorPointerType mpDx;
    typename TSparseSpace::VectorPointerType mpb;
};

} // namespace

// --- Sparse graph construction ----------------------------------------------
// DOF set and equation numbering are prepared once; each iteration allocates
// fresh system containers and runs ConstructMatrixStructure through
// ResizeAndInitializeVectors — the first-solve graph-building cost.
template <class TSparseSpace>
void BM_ConstructSystemStructure(benchmark::State& rState)
{
    auto& r_model_part = GetBenchmarkModelPart(rState.range(0), rState.range(1));
    AssemblyFixture<TSparseSpace> fixture(r_model_part);

    for (auto _ : rState) {
        typename TSparseSpace::MatrixPointerType p_A;
        typename TSparseSpace::VectorPointerType p_Dx;
        typename TSparseSpace::VectorPointerType p_b;
        fixture.mBuilderAndSolver.ResizeAndInitializeVectors(fixture.mpScheme, p_A, p_Dx, p_b, r_model_part);
        benchmark::DoNotOptimize(p_A->nnz());
    }
}

// --- Full LHS + RHS build ----------------------------------------------------
// The realistic per-nonlinear-iteration cost: zero the system (as the
// strategies do) and assemble every element's local system into the CSR
// matrix and the RHS vector.
template <class TSparseSpace>
void BM_BuildLHSAndRHS(benchmark::State& rState)
{
    auto& r_model_part = GetBenchmarkModelPart(rState.range(0), rState.range(1));
    AssemblyFixture<TSparseSpace> fixture(r_model_part);
    fixture.InitializeSystem();

    for (auto _ : rState) {
        TSparseSpace::SetToZero(*fixture.mpA);
        TSparseSpace::SetToZero(*fixture.mpb);
        fixture.mBuilderAndSolver.Build(fixture.mpScheme, r_model_part, *fixture.mpA, *fixture.mpb);
        benchmark::DoNotOptimize(TSparseSpace::TwoNorm(*fixture.mpb));
    }
}

// Mesh sizes: {dimension, number_of_divisions}. 2D: ~33k / ~526k triangles
// (~17k / ~263k DOFs); 3D: ~25k / ~384k tetrahedra (~5k / ~69k DOFs).
#define KRATOS_REGISTER_ASSEMBLY_BENCHMARK(name)                    \
    BENCHMARK_TEMPLATE(name, DefaultSparse)                         \
        ->Name(#name "/" KRATOS_BACKEND_NAME)                       \
        ->Args({2, 128})->Args({2, 512})->Args({3, 16})->Args({3, 40}) \
        ->Unit(benchmark::kMillisecond);

KRATOS_REGISTER_ASSEMBLY_BENCHMARK(BM_ConstructSystemStructure)
KRATOS_REGISTER_ASSEMBLY_BENCHMARK(BM_BuildLHSAndRHS)

#undef KRATOS_REGISTER_ASSEMBLY_BENCHMARK
#undef KRATOS_BACKEND_NAME

} // namespace Kratos

BENCHMARK_MAIN();
