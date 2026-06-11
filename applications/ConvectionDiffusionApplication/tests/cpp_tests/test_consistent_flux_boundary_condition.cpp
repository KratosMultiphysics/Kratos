// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/global_pointer_variables.h"
#include "testing/testing.h"

// Application includes
#include "convection_diffusion_application.h"
#include "../test_utilities/convection_diffusion_testing_utilities.h"

namespace Kratos::Testing
{

namespace
{

void AddTemperatureDofs(ModelPart& rModelPart)
{
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(TEMPERATURE);
    }
}

void SetUniformConductivity(ModelPart& rModelPart, const double Value)
{
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.FastGetSolutionStepValue(CONDUCTIVITY) = Value;
    }
}

void SetParentElementToCondition(Condition& rCondition, Element& rElement)
{
    rCondition.SetValue(NEIGHBOUR_ELEMENTS, GlobalPointersVector<Element>({&rElement}));
}

void CheckVector(const Vector& rVector, const std::vector<double>& rExpected, const double Tolerance)
{
    KRATOS_EXPECT_EQ(rVector.size(), rExpected.size());
    for (IndexType i = 0; i < rVector.size(); ++i) {
        KRATOS_EXPECT_NEAR(rVector[i], rExpected[i], Tolerance);
    }
}

void CheckMatrix(const Matrix& rMatrix, const std::vector<double>& rExpected, const double Tolerance)
{
    KRATOS_EXPECT_EQ(rMatrix.size1() * rMatrix.size2(), rExpected.size());
    for (IndexType i = 0; i < rMatrix.size1(); ++i) {
        for (IndexType j = 0; j < rMatrix.size2(); ++j) {
            KRATOS_EXPECT_NEAR(rMatrix(i, j), rExpected[i * rMatrix.size2() + j], Tolerance);
        }
    }
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(ConsistentFluxBoundaryCondition2D2NTriangleParent, KratosConvectionDiffusionFastSuite)
{
    Model model;
    auto& r_test_model_part = model.CreateModelPart("TestModelPart");
    ConvectionDiffusionTestingUtilities::SetEntityUnitTestModelPart(r_test_model_part);

    r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    AddTemperatureDofs(r_test_model_part);
    SetUniformConductivity(r_test_model_part, 2.0);

    r_test_model_part.GetNode(1).FastGetSolutionStepValue(TEMPERATURE) = 1.0;
    r_test_model_part.GetNode(2).FastGetSolutionStepValue(TEMPERATURE) = 2.0;
    r_test_model_part.GetNode(3).FastGetSolutionStepValue(TEMPERATURE) = 3.0;

    auto p_properties = r_test_model_part.pGetProperties(0);
    r_test_model_part.CreateNewElement("LaplacianElement2D3N", 1, std::vector<ModelPart::IndexType>{1, 2, 3}, p_properties);
    r_test_model_part.CreateNewCondition("ConsistentFluxBoundaryCondition2D2N", 1, std::vector<ModelPart::IndexType>{1, 2}, p_properties);

    auto p_condition = r_test_model_part.pGetCondition(1);
    auto p_element = r_test_model_part.pGetElement(1);
    SetParentElementToCondition(*p_condition, *p_element);

    KRATOS_EXPECT_EQ(p_condition->Check(r_test_model_part.GetProcessInfo()), 0);

    Matrix lhs;
    Vector rhs;
    p_condition->CalculateLocalSystem(lhs, rhs, r_test_model_part.GetProcessInfo());

    // Sign convention: the consistent-flux BC restores the weak-form boundary
    // term  - \int_\Gamma N_i * D * (n . grad c)  (n outward), so its LHS is
    // -K_b with K_b(i,j) = \int_\Gamma N_i * D * (n . grad N_j), and its
    // RHS = -LHS*c = +K_b*c. See the OutwardDiffusiveFluxSign test below for an
    // implementation-independent check of this sign.
    CheckMatrix(lhs, {
        -1.0, 0.0, 1.0,
        -1.0, 0.0, 1.0,
        0.0, 0.0, 0.0}, 1.0e-12);
    CheckVector(rhs, {-2.0, -2.0, 0.0}, 1.0e-12);
}

KRATOS_TEST_CASE_IN_SUITE(ConsistentFluxBoundaryCondition2D2NQuadrilateralParent, KratosConvectionDiffusionFastSuite)
{
    Model model;
    auto& r_test_model_part = model.CreateModelPart("TestModelPart");
    ConvectionDiffusionTestingUtilities::SetEntityUnitTestModelPart(r_test_model_part);

    r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_test_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    r_test_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
    AddTemperatureDofs(r_test_model_part);
    SetUniformConductivity(r_test_model_part, 2.0);

    r_test_model_part.GetNode(1).FastGetSolutionStepValue(TEMPERATURE) = 1.0;
    r_test_model_part.GetNode(2).FastGetSolutionStepValue(TEMPERATURE) = 2.0;
    r_test_model_part.GetNode(3).FastGetSolutionStepValue(TEMPERATURE) = 3.0;
    r_test_model_part.GetNode(4).FastGetSolutionStepValue(TEMPERATURE) = 4.0;

    auto p_properties = r_test_model_part.pGetProperties(0);
    r_test_model_part.CreateNewElement("EulerianConvDiff2D4N", 1, std::vector<ModelPart::IndexType>{1, 2, 3, 4}, p_properties);
    r_test_model_part.CreateNewCondition("ConsistentFluxBoundaryCondition2D2N", 1, std::vector<ModelPart::IndexType>{1, 2}, p_properties);

    auto p_condition = r_test_model_part.pGetCondition(1);
    auto p_element = r_test_model_part.pGetElement(1);
    SetParentElementToCondition(*p_condition, *p_element);

    KRATOS_EXPECT_EQ(p_condition->Check(r_test_model_part.GetProcessInfo()), 0);

    Matrix lhs;
    Vector rhs;
    p_condition->CalculateLocalSystem(lhs, rhs, r_test_model_part.GetProcessInfo());

    CheckMatrix(lhs, {
        -2.0 / 3.0, -1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0,
        -1.0 / 3.0, -2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0}, 1.0e-12);
    CheckVector(rhs, {-7.0 / 3.0, -5.0 / 3.0, 0.0, 0.0}, 1.0e-12);
}

KRATOS_TEST_CASE_IN_SUITE(ConsistentFluxBoundaryCondition3D3NTetrahedronParent, KratosConvectionDiffusionFastSuite)
{
    Model model;
    auto& r_test_model_part = model.CreateModelPart("TestModelPart");
    ConvectionDiffusionTestingUtilities::SetEntityUnitTestModelPart(r_test_model_part);

    r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    r_test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
    AddTemperatureDofs(r_test_model_part);
    SetUniformConductivity(r_test_model_part, 2.0);

    r_test_model_part.GetNode(1).FastGetSolutionStepValue(TEMPERATURE) = 1.0;
    r_test_model_part.GetNode(2).FastGetSolutionStepValue(TEMPERATURE) = 2.0;
    r_test_model_part.GetNode(3).FastGetSolutionStepValue(TEMPERATURE) = 3.0;
    r_test_model_part.GetNode(4).FastGetSolutionStepValue(TEMPERATURE) = 4.0;

    auto p_properties = r_test_model_part.pGetProperties(0);
    r_test_model_part.CreateNewElement("LaplacianElement3D4N", 1, std::vector<ModelPart::IndexType>{1, 2, 3, 4}, p_properties);
    r_test_model_part.CreateNewCondition("ConsistentFluxBoundaryCondition3D3N", 1, std::vector<ModelPart::IndexType>{1, 2, 3}, p_properties);

    auto p_condition = r_test_model_part.pGetCondition(1);
    auto p_element = r_test_model_part.pGetElement(1);
    SetParentElementToCondition(*p_condition, *p_element);

    KRATOS_EXPECT_EQ(p_condition->Check(r_test_model_part.GetProcessInfo()), 0);

    Matrix lhs;
    Vector rhs;
    p_condition->CalculateLocalSystem(lhs, rhs, r_test_model_part.GetProcessInfo());

    Vector rhs_only;
    p_condition->CalculateRightHandSide(rhs_only, r_test_model_part.GetProcessInfo());

    CheckMatrix(lhs, {
        -1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0,
        -1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0,
        -1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0,
        0.0, 0.0, 0.0, 0.0}, 1.0e-12);
    CheckVector(rhs, {-1.0, -1.0, -1.0, 0.0}, 1.0e-12);
    CheckVector(rhs_only, {-1.0, -1.0, -1.0, 0.0}, 1.0e-12);
}

// Implementation-independent sign check: the consistent-flux BC must add to the
// RHS the TRUE outward diffusive flux of the current field projected on the
// shape functions, i.e. RHS_i = +\int_\Gamma N_i * D * (n . grad c) dGamma with
// n the OUTWARD normal. We impose a linear field with a known gradient and
// compare against the closed-form result, which does not depend on how K_b is
// assembled internally and therefore pins the overall sign.
KRATOS_TEST_CASE_IN_SUITE(ConsistentFluxBoundaryCondition2D2NOutwardDiffusiveFluxSign, KratosConvectionDiffusionFastSuite)
{
    Model model;
    auto& r_test_model_part = model.CreateModelPart("TestModelPart");
    ConvectionDiffusionTestingUtilities::SetEntityUnitTestModelPart(r_test_model_part);

    // Reference triangle; the condition is the edge {1,2} lying on y = 0.
    // The parent centroid is at (1/3, 1/3), so the outward normal is (0, -1, 0).
    r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    AddTemperatureDofs(r_test_model_part);

    const double conductivity = 2.0;
    SetUniformConductivity(r_test_model_part, conductivity);

    // Linear field c = 3*x + 5*y -> grad c = (3, 5, 0), so n . grad c = -5.
    r_test_model_part.GetNode(1).FastGetSolutionStepValue(TEMPERATURE) = 0.0; // c(0,0)
    r_test_model_part.GetNode(2).FastGetSolutionStepValue(TEMPERATURE) = 3.0; // c(1,0)
    r_test_model_part.GetNode(3).FastGetSolutionStepValue(TEMPERATURE) = 5.0; // c(0,1)

    auto p_properties = r_test_model_part.pGetProperties(0);
    r_test_model_part.CreateNewElement("LaplacianElement2D3N", 1, std::vector<ModelPart::IndexType>{1, 2, 3}, p_properties);
    r_test_model_part.CreateNewCondition("ConsistentFluxBoundaryCondition2D2N", 1, std::vector<ModelPart::IndexType>{1, 2}, p_properties);

    auto p_condition = r_test_model_part.pGetCondition(1);
    auto p_element = r_test_model_part.pGetElement(1);
    SetParentElementToCondition(*p_condition, *p_element);

    KRATOS_EXPECT_EQ(p_condition->Check(r_test_model_part.GetProcessInfo()), 0);

    Vector rhs;
    p_condition->CalculateRightHandSide(rhs, r_test_model_part.GetProcessInfo());

    // RHS_i = D * (n . grad c) * \int_edge N_i. The edge {1,2} has length 1, so
    // \int N_1 = \int N_2 = 1/2 and \int N_3 = 0 (N_3 vanishes on y = 0).
    // RHS = 2 * (-5) * {1/2, 1/2, 0} = {-5, -5, 0}.
    const double n_dot_grad_c = -5.0;
    const double expected = conductivity * n_dot_grad_c * 0.5;
    CheckVector(rhs, {expected, expected, 0.0}, 1.0e-12);
}

} // namespace Kratos::Testing