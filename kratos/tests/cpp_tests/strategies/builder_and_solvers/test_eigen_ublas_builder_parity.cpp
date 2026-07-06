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
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"

/* Utility includes */
#include "containers/model.h"
#include "spaces/ublas_space.h"
#include "spaces/eigen_space.h"

/* Element include */
#include "geometries/line_2d_2.h"
#include "tests/test_utilities/test_bar_element.h"

/* Linear solvers */
#include "linear_solvers/reorderer.h"
#include "linear_solvers/direct_solver.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"

/* The most basic scheme (static) */
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

/* The builder and solver */
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

namespace Kratos::Testing
{

namespace EigenUblasParityInternals
{

// NOTE: this test compiles the block builder-and-solver with BOTH the uBLAS
// and the Eigen sparse spaces unconditionally, independently of the
// KRATOS_LINEAR_ALGEBRA_BACKEND setting: it is the permanent compile and
// behavior parity check between the two backends.

/// Generates the same truss structure used by test_builder_and_solver.cpp
static inline void SetUpTestModelPart(ModelPart& rModelPart, const bool WithConstraint)
{
    using GeometryType = Geometry<Node>;

    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
    rModelPart.AddNodalSolutionStepVariable(REACTION);
    rModelPart.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

    Node::Pointer pnode1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    Node::Pointer pnode2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    Node::Pointer pnode3 = rModelPart.CreateNewNode(3, 2.0, 0.0, 0.0);

    auto p_prop = rModelPart.CreateNewProperties(1, 0);
    p_prop->SetValue(YOUNG_MODULUS, 206900000000.0);
    p_prop->SetValue(NODAL_AREA, 0.01);

    GeometryType::Pointer pgeom1 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode1, pnode2})});
    rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>(1, pgeom1, p_prop));
    GeometryType::Pointer pgeom2 = Kratos::make_shared<Line2D2<Node>>(PointerVector<Node>{std::vector<Node::Pointer>({pnode2, pnode3})});
    rModelPart.AddElement(Kratos::make_intrusive<TestBarElement>(2, pgeom2, p_prop));

    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X, REACTION_X);
        r_node.AddDof(DISPLACEMENT_Y, REACTION_Y);
        r_node.AddDof(DISPLACEMENT_Z, REACTION_Z);
    }

    const auto& r_process_info = rModelPart.GetProcessInfo();
    for (auto& r_elem : rModelPart.Elements()) {
        r_elem.Initialize(r_process_info);
        r_elem.InitializeSolutionStep(r_process_info);
    }

    for (auto& r_node : rModelPart.Nodes()) {
        (r_node.FastGetSolutionStepValue(DISPLACEMENT)).clear();
        (r_node.FastGetSolutionStepValue(DISPLACEMENT, 1)).clear();
        (r_node.FastGetSolutionStepValue(DISPLACEMENT, 2)).clear();
    }

    for (auto& r_node : rModelPart.Nodes()) {
        r_node.Fix(DISPLACEMENT_Y);
        r_node.Fix(DISPLACEMENT_Z);
    }
    pnode1->Fix(DISPLACEMENT_X);

    if (WithConstraint) {
        rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, *pnode2, DISPLACEMENT_X, *pnode3, DISPLACEMENT_X, 1.0, 0.0);
    }
}

/// Builds the system with the block builder-and-solver for the given sparse
/// space and returns (A, b, Dx after solving with skyline LU).
template<class TSparseSpaceType>
struct SystemResults
{
    typename TSparseSpaceType::MatrixPointerType pA;
    typename TSparseSpaceType::VectorPointerType pb;
    typename TSparseSpaceType::VectorPointerType pDx;
};

template<class TSparseSpaceType>
SystemResults<TSparseSpaceType> BuildAndSolveSystem(
    Model& rModel,
    const bool WithConstraint
    )
{
    using LocalSpaceType = TUblasDenseSpace<double>;
    using LinearSolverType = LinearSolver<TSparseSpaceType, LocalSpaceType>;
    using ReordererType = Reorderer<TSparseSpaceType, LocalSpaceType>;
    using SkylineLUFactorizationSolverType = SkylineLUFactorizationSolver<TSparseSpaceType, LocalSpaceType, ReordererType>;
    using SchemeType = ResidualBasedIncrementalUpdateStaticScheme<TSparseSpaceType, LocalSpaceType>;
    using BuilderAndSolverType = ResidualBasedBlockBuilderAndSolver<TSparseSpaceType, LocalSpaceType, LinearSolverType>;

    ModelPart& r_model_part = rModel.CreateModelPart("Main", 3);
    SetUpTestModelPart(r_model_part, WithConstraint);

    auto p_scheme = Kratos::make_shared<SchemeType>();
    auto p_solver = typename LinearSolverType::Pointer(new SkylineLUFactorizationSolverType());
    auto p_builder_and_solver = Kratos::make_shared<BuilderAndSolverType>(p_solver);

    p_builder_and_solver->SetDofSetIsInitializedFlag(false);
    p_builder_and_solver->Clear();
    p_scheme->Clear();

    SystemResults<TSparseSpaceType> results;

    p_builder_and_solver->SetUpDofSet(p_scheme, r_model_part);
    p_builder_and_solver->SetUpSystem(r_model_part);
    p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, results.pA, results.pDx, results.pb, r_model_part);

    auto& rA = *results.pA;
    auto& rDx = *results.pDx;
    auto& rb = *results.pb;

    p_builder_and_solver->InitializeSolutionStep(r_model_part, rA, rDx, rb);
    p_scheme->InitializeSolutionStep(r_model_part, rA, rDx, rb);
    p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);

    p_builder_and_solver->Build(p_scheme, r_model_part, rA, rb);
    if (r_model_part.MasterSlaveConstraints().size() != 0) {
        p_builder_and_solver->ApplyConstraints(p_scheme, r_model_part, rA, rb);
    }
    p_builder_and_solver->ApplyDirichletConditions(p_scheme, r_model_part, rA, rDx, rb);

    p_builder_and_solver->SystemSolve(rA, rDx, rb);

    return results;
}

template<class TEigenResults, class TUblasResults>
void CompareSystems(TEigenResults& rEigen, TUblasResults& rUblas)
{
    constexpr double tolerance = 1e-9;

    auto& r_eigen_A = *rEigen.pA;
    auto& r_ublas_A = *rUblas.pA;

    // Same sizes and same sparsity pattern
    KRATOS_EXPECT_EQ(r_eigen_A.size1(), r_ublas_A.size1());
    KRATOS_EXPECT_EQ(r_eigen_A.size2(), r_ublas_A.size2());
    KRATOS_EXPECT_EQ(r_eigen_A.nnz(), static_cast<std::size_t>(r_ublas_A.nnz()));

    const auto eigen_row_indices = r_eigen_A.index1_data();
    const auto eigen_col_indices = r_eigen_A.index2_data();
    const auto eigen_values = r_eigen_A.value_data();
    const auto& ublas_row_indices = r_ublas_A.index1_data();
    const auto& ublas_col_indices = r_ublas_A.index2_data();
    const auto& ublas_values = r_ublas_A.value_data();

    for (std::size_t i = 0; i < r_eigen_A.size1() + 1; ++i) {
        KRATOS_EXPECT_EQ(static_cast<std::size_t>(eigen_row_indices[i]), ublas_row_indices[i]);
    }
    const std::size_t nnz = r_eigen_A.nnz();
    for (std::size_t k = 0; k < nnz; ++k) {
        KRATOS_EXPECT_EQ(static_cast<std::size_t>(eigen_col_indices[k]), ublas_col_indices[k]);
        KRATOS_EXPECT_NEAR(eigen_values[k], ublas_values[k], tolerance * (1.0 + std::abs(ublas_values[k])));
    }

    // Same RHS and same solution
    const auto& r_eigen_b = *rEigen.pb;
    const auto& r_ublas_b = *rUblas.pb;
    const auto& r_eigen_Dx = *rEigen.pDx;
    const auto& r_ublas_Dx = *rUblas.pDx;
    KRATOS_EXPECT_EQ(static_cast<std::size_t>(r_eigen_b.size()), static_cast<std::size_t>(r_ublas_b.size()));
    for (std::size_t i = 0; i < static_cast<std::size_t>(r_ublas_b.size()); ++i) {
        KRATOS_EXPECT_NEAR(r_eigen_b[i], r_ublas_b[i], tolerance * (1.0 + std::abs(r_ublas_b[i])));
        KRATOS_EXPECT_NEAR(r_eigen_Dx[i], r_ublas_Dx[i], tolerance * (1.0 + std::abs(r_ublas_Dx[i])));
    }
}

} // namespace EigenUblasParityInternals

KRATOS_TEST_CASE_IN_SUITE(EigenUblasBlockBuilderAndSolverParity, KratosCoreFastSuite)
{
    Model eigen_model, ublas_model;
    auto eigen_results = EigenUblasParityInternals::BuildAndSolveSystem<TEigenSparseSpace<double>>(eigen_model, false);
    auto ublas_results = EigenUblasParityInternals::BuildAndSolveSystem<TUblasSparseSpace<double>>(ublas_model, false);
    EigenUblasParityInternals::CompareSystems(eigen_results, ublas_results);
}

KRATOS_TEST_CASE_IN_SUITE(EigenUblasBlockBuilderAndSolverParityWithConstraints, KratosCoreFastSuite)
{
    Model eigen_model, ublas_model;
    auto eigen_results = EigenUblasParityInternals::BuildAndSolveSystem<TEigenSparseSpace<double>>(eigen_model, true);
    auto ublas_results = EigenUblasParityInternals::BuildAndSolveSystem<TUblasSparseSpace<double>>(ublas_model, true);
    EigenUblasParityInternals::CompareSystems(eigen_results, ublas_results);
}

} // namespace Kratos::Testing
