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
#include "containers/model.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "testing/testing.h"
#include "solving_strategies/schemes/new_scheme.h"

namespace Kratos::Testing
{

// typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
// typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

// typedef ModelPart::DofsArrayType DofsArrayType;

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
        const double aux = f * h / 2.0;
        rRightHandSideVector(0) = aux;
        rRightHandSideVector(1) = aux;
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

static void SetUpTestSchemesModelPart(ModelPart& rModelPart)
{
    const int domain_size = 2;
    const std::size_t buffer_size = 3;
    rModelPart.SetBufferSize(buffer_size);
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, domain_size);

    rModelPart.AddNodalSolutionStepVariable(DISTANCE);

    auto p_node_1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_node_3 = rModelPart.CreateNewNode(3, 2.0, 0.0, 0.0);

    auto p_geom_1 = Kratos::make_shared<Line2D2<Node>>(p_node_1, p_node_2);
    auto p_geom_2 = Kratos::make_shared<Line2D2<Node>>(p_node_2, p_node_3);

    rModelPart.AddElement(Kratos::make_intrusive<TestElement>(1, p_geom_1));
    rModelPart.AddElement(Kratos::make_intrusive<TestElement>(2, p_geom_2));

    rModelPart.GetNodalSolutionStepVariablesList().AddDof(&DISTANCE);
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(DISTANCE);
    }

    for (auto& r_elem : rModelPart.Elements()) {
        r_elem.Set(ACTIVE, true);
    }
}

}

KRATOS_TEST_CASE_IN_SUITE(NewScheme, KratosCoreFastSuite)
{
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    SetUpTestSchemesModelPart(r_test_model_part);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_type" : "block",
        "echo_level" : 2
    })");
    auto p_scheme = Kratos::make_unique<NewScheme<>>(scheme_settings);

    // Create the DOF set and set the global ids
    p_scheme->SetUpDofArray(r_test_model_part);
    p_scheme->SetUpSystem(r_test_model_part);

    // Set up the matrix graph and arrays
    // Note that in a standard case this happens at the strategy level
    SparseGraph<> sparse_graph;
    AssemblyUtilities::SetUpSparseGraph(r_test_model_part, sparse_graph);
    auto p_dx = Kratos::make_shared<NewScheme<>::SystemVectorType>(sparse_graph);
    auto p_rhs = Kratos::make_shared<NewScheme<>::SystemVectorType>(sparse_graph);
    auto p_lhs = Kratos::make_shared<NewScheme<>::SystemMatrixType>(sparse_graph);
    sparse_graph.Clear();

    // // Initialize the global system arrays
    // p_scheme->ResizeAndInitializeVectors(r_test_model_part, p_lhs, p_dx, p_rhs);

    // Call the build
    p_scheme->Build(r_test_model_part, *p_lhs, *p_rhs);

    // Check resultant matrices
    const double tol = 1.0e-12;
    std::vector<double> expected_rhs = {0.5,1.0,0.5};
    BoundedMatrix<double,3,3> expected_lhs;
    expected_lhs(0,0) = 2.0; expected_lhs(0,1) = -1.0; expected_lhs(0,2) = 0.0;
    expected_lhs(1,0) = -1.0; expected_lhs(1,1) = 2.0; expected_lhs(1,2) = -1.0;
    expected_lhs(2,0) = 0.0; expected_lhs(2,1) = -1.0; expected_lhs(2,2) = 2.0;
    KRATOS_CHECK_VECTOR_NEAR((*p_rhs), expected_rhs, tol);
}

}  // namespace Kratos::Testing.

