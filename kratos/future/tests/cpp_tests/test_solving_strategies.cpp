//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes
#include <limits>

// External includes

// Project includes
#include "constraints/linear_master_slave_constraint.h"
#include "containers/model.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "testing/testing.h"

#ifdef KRATOS_USE_FUTURE
#include "future/linear_solvers/amgcl_solver.h"
#include "future/linear_solvers/skyline_lu_factorization_solver.h"
#include "future/solving_strategies/schemes/static_scheme.h"
#include "future/solving_strategies/strategies/linear_strategy.h"
#endif

namespace Kratos::Testing
{

namespace
{

class TestElement : public Element
{
public:

    using IndexType = Element::IndexType;

    using DofsArrayType = Element::DofsArrayType;

    using DofsVectorType = Element::DofsVectorType;

    using EquationIdVectorType = Element::EquationIdVectorType;

    TestElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {}

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
        CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    }

    void CalculateLeftHandSide(
        MatrixType &rLeftHandSideMatrix,
        const ProcessInfo &rCurrentProcessInfo) override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        if (rLeftHandSideMatrix.size1() != n_nodes || rLeftHandSideMatrix.size2() != n_nodes) {
            rLeftHandSideMatrix.resize(n_nodes, n_nodes, false);
        }
        const double aux = 1.0 / r_geom.Length();
        rLeftHandSideMatrix(0,0) = aux;
        rLeftHandSideMatrix(0,1) = -aux;
        rLeftHandSideMatrix(1,0) = -aux;
        rLeftHandSideMatrix(1,1) = aux;
    }

    void CalculateRightHandSide(
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo) override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        if (rRightHandSideVector.size() != n_nodes) {
            rRightHandSideVector.resize(n_nodes);
        }
        const double f = 1.0;
        const double h = r_geom.Length();
        const double aux_f = f * h / 2.0;
        const double aux_phi =  (r_geom[0].GetSolutionStepValue(DISTANCE) - r_geom[1].GetSolutionStepValue(DISTANCE)) / h;
        rRightHandSideVector(0) = aux_f - aux_phi;
        rRightHandSideVector(1) = aux_f + aux_phi;
    }

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        if (rResult.size() != n_nodes) {
            rResult.resize(n_nodes);
        }
        for (unsigned int i = 0; i < n_nodes; ++i) {
            rResult[i] = r_geom[i].GetDof(DISTANCE).EquationId();
        }
    }

    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        if (rElementalDofList.size() != n_nodes) {
            rElementalDofList.resize(n_nodes);
        }
        for (unsigned int i = 0; i < n_nodes; i++) {
            rElementalDofList[i] = r_geom[i].pGetDof(DISTANCE);
        }
    }

};

static void SetUpTestSchemesModelPart(
    const std::size_t NumElems,
    const double ElemSize,
    ModelPart& rModelPart)
{
    const int domain_size = 2;
    const std::size_t buffer_size = 3;
    rModelPart.SetBufferSize(buffer_size);
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, domain_size);

    rModelPart.AddNodalSolutionStepVariable(DISTANCE);

    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    for (std::size_t i_elem = 0; i_elem < NumElems; ++i_elem) {
        auto p_node_a = rModelPart.pGetNode(i_elem + 1);
        auto p_node_b = rModelPart.CreateNewNode(i_elem + 2, (i_elem + 1)*ElemSize, 0.0, 0.0);
        auto p_geom = Kratos::make_shared<Line2D2<Node>>(p_node_a, p_node_b);
        rModelPart.AddElement(Kratos::make_intrusive<TestElement>(i_elem + 1, p_geom));
    }

    rModelPart.GetNodalSolutionStepVariablesList().AddDof(&DISTANCE);
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(DISTANCE);
    }

    for (auto& r_elem : rModelPart.Elements()) {
        r_elem.Set(ACTIVE, true);
    }
}

}

KRATOS_TEST_CASE_IN_SUITE(StaticScheme, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 2;
    const double elem_size = 1.0;
    SetUpTestSchemesModelPart(num_elems, elem_size, r_test_model_part);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "build_type" : "block"
        }
    })");
    auto p_scheme = Kratos::make_unique<Future::StaticScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, scheme_settings);

    // Create the DOF set and set the global ids
    // Note that in a standard case this happens at the strategy level
    ModelPart::DofsArrayType dof_set;
    p_scheme->SetUpDofArray(dof_set);
    p_scheme->SetUpSystemIds(dof_set);

    // Set up the matrix graph and arrays
    // Note that in a standard case this happens at the strategy level
    auto p_lhs = Kratos::make_shared<CsrMatrix<>>();
    auto p_dx = Kratos::make_shared<SystemVector<>>();
    auto p_rhs = Kratos::make_shared<SystemVector<>>();
    p_scheme->ResizeAndInitializeVectors(dof_set, p_lhs, p_dx, p_rhs);

    // Call the build
    p_scheme->Build(*p_lhs, *p_rhs);

    // Check resultant matrices
    const double tol = 1.0e-12;
    std::vector<double> expected_rhs = {0.5,1.0,0.5};
    BoundedMatrix<double,3,3> expected_lhs;
    expected_lhs(0,0) = 1.0; expected_lhs(0,1) = -1.0; expected_lhs(0,2) = 0.0;
    expected_lhs(1,0) = -1.0; expected_lhs(1,1) = 2.0; expected_lhs(1,2) = -1.0;
    expected_lhs(2,0) = 0.0; expected_lhs(2,1) = -1.0; expected_lhs(2,2) = 1.0;
    KRATOS_CHECK_VECTOR_NEAR((*p_rhs), expected_rhs, tol); // Note that as there are not non-zero entries in the sparse vector we can use the standard macro
    for (unsigned int i = 0; i < p_lhs->size1(); ++i) {
        for (unsigned int j = 0; j < p_lhs->size2(); ++j) {
            const double expected_val = expected_lhs(i,j);
            if (std::abs(expected_val) > tol) {
                KRATOS_CHECK_NEAR(p_lhs->operator()(i,j), expected_val, tol); // Note that we check if the expected value is non-zero as this is a CSR matrix
            }
        }
    }
#else
    true;
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategy, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 2;
    const double elem_size = 1.0;
    SetUpTestSchemesModelPart(num_elems, elem_size, r_test_model_part);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "build_type" : "block"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<CsrMatrix<>, SystemVector<>>;
    using LinearSolverType = Future::LinearSolver<CsrMatrix<>, SystemVector<>>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Apply Dirichlet BCs
    auto p_node_1 = r_test_model_part.pGetNode(1);
    p_node_1->FastGetSolutionStepValue(DISTANCE, 0) = 1.0;
    p_node_1->Fix(DISTANCE);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->InitializeSolutionStep();
    p_strategy->Predict();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();
    p_strategy->Clear();

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 2.5, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 3.0, 1.0e-12);
#else
    true;
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyWithPeriodicityConstraint, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 3;
    const double elem_size = 1.0;
    SetUpTestSchemesModelPart(num_elems, elem_size, r_test_model_part);

    // Create a periodicity constraint with jump
    const double jump = 1.0;
    const double weight = 1.0;
    auto& r_master_node = r_test_model_part.GetNode(1);
    auto& r_slave_node = r_test_model_part.GetNode(3);
    auto p_const_1 = r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, r_master_node, DISTANCE, r_slave_node, DISTANCE, weight, jump);
    p_const_1->Set(ACTIVE, true);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "build_type" : "block"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<CsrMatrix<>, SystemVector<>>;
    using LinearSolverType = Future::LinearSolver<CsrMatrix<>, SystemVector<>>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Apply Dirichlet BCs
    auto p_node_1 = r_test_model_part.pGetNode(1);
    p_node_1->FastGetSolutionStepValue(DISTANCE, 0) = 1.0;
    p_node_1->Fix(DISTANCE);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->InitializeSolutionStep();
    p_strategy->Predict();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();
    p_strategy->Clear();

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 1.5, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 0.0, 1.0e-12);
#else
    true;
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyWithMultipleDofsConstraints, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 3;
    const double elem_size = 1.0;
    SetUpTestSchemesModelPart(num_elems, elem_size, r_test_model_part);

    // Create a constraint involving multiple master DOFs
    std::vector<typename Dof<double>::Pointer> slave_dofs(1);
    std::vector<typename Dof<double>::Pointer> master_dofs(2);
    slave_dofs[0] = r_test_model_part.GetNode(4).pGetDof(DISTANCE);
    master_dofs[0] = r_test_model_part.GetNode(1).pGetDof(DISTANCE);
    master_dofs[1] = r_test_model_part.GetNode(2).pGetDof(DISTANCE);

    Vector jump_vector(1);
    jump_vector[0] = 1.0;

    Matrix relation_matrix(1, 2);
    relation_matrix(0,0) = 0.5;
    relation_matrix(0,1) = 0.5;

    auto p_const_1 = r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, master_dofs, slave_dofs, relation_matrix, jump_vector);
    p_const_1->Set(ACTIVE, true);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "build_type" : "block"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<CsrMatrix<>, SystemVector<>>;
    using LinearSolverType = Future::LinearSolver<CsrMatrix<>, SystemVector<>>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Apply Dirichlet BCs
    auto p_node_1 = r_test_model_part.pGetNode(1);
    p_node_1->FastGetSolutionStepValue(DISTANCE, 0) = 1.0;
    p_node_1->Fix(DISTANCE);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->InitializeSolutionStep();
    p_strategy->Predict();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();
    p_strategy->Clear();

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 2.888888888888888, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 3.166666666666666, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(4).FastGetSolutionStepValue(DISTANCE), 0.0, 1.0e-12);
#else
    true;
#endif
}

}  // namespace Kratos::Testing.

