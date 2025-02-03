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
    const int domain_size = 3;
    const std::size_t buffer_size = 3;
    rModelPart.SetBufferSize(buffer_size);
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, domain_size);

    rModelPart.AddNodalSolutionStepVariable(DISTANCE);

    auto p_node_1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_node_3 = rModelPart.CreateNewNode(3, 2.0, 0.0, 0.0);

    auto p_geom_1 = Kratos::make_shared<Line2D2<Node>>(p_node_1, p_node_2);
    auto p_geom_2 = Kratos::make_shared<Line2D2<Node>>(p_node_2, p_node_3);

    auto p_elem_1 = Kratos::make_shared<TestElement>(1, p_geom_1);
    auto p_elem_2 = Kratos::make_shared<TestElement>(2, p_geom_2);

    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(DISTANCE);
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
        "build_type" : "block"
    })");
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    auto p_scheme = Kratos::make_unique<NewScheme<SparseSpaceType, LocalSpaceType>>(scheme_settings);

    // Create the DOF set
    p_scheme->SetUpDofArray(r_test_model_part);

}

}  // namespace Kratos::Testing.

