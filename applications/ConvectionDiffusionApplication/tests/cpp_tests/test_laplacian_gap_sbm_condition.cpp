// BSD License, Kratos default license: kratos/license.txt
#include <cmath>
#include "testing/testing.h"
#include "convection_diffusion_application.h"
#include "containers/model.h"
#include "includes/convection_diffusion_settings.h"
#include "integration/line_gauss_legendre_integration_points.h"
#include "convection_diffusion_application_variables.h"

namespace Kratos::Testing {
namespace {
Condition::Pointer CreateLaplacianGapCondition(ModelPart& rPart, bool Connector)
{
    rPart.AddNodalSolutionStepVariable(TEMPERATURE);
    auto p_settings = Kratos::make_shared<ConvectionDiffusionSettings>();
    p_settings->SetUnknownVariable(TEMPERATURE);
    rPart.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, p_settings);
    auto p_properties = rPart.CreateNewProperties(1);
    p_properties->SetValue(PENALTY_COEFFICIENT, 40.0);
    rPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    rPart.CreateNewNode(4, 0.0, 1.0, 0.0);
    rPart.CreateNewNode(5, 2.0, 0.0, 0.0);
    rPart.CreateNewNode(6, 2.0, 1.0, 0.0);
    for (auto& r_node : rPart.Nodes()) r_node.AddDof(TEMPERATURE);
    auto p_first = rPart.CreateNewElement("Element2D4N", 1, {1, 2, 3, 4}, p_properties);
    auto p_second = rPart.CreateNewElement("Element2D4N", 2, {2, 5, 6, 3}, p_properties);
    auto p_condition = rPart.CreateNewCondition(Connector ? "LaplacianGapConnectorCondition2D" : "LaplacianGapDirichletCondition2D",
        1, Connector ? std::vector<ModelPart::IndexType>{1, 2, 3, 4, 5, 6} :
                       std::vector<ModelPart::IndexType>{1, 2, 3, 4}, p_properties);
    GlobalPointersVector<Element> parents;
    parents.push_back(GlobalPointer<Element>(p_first.get()));
    if (Connector) parents.push_back(GlobalPointer<Element>(p_second.get()));
    p_condition->SetValue(NEIGHBOUR_ELEMENTS, parents);
    Matrix segment = ZeroMatrix(2, 3);
    segment(0, 0) = segment(1, 0) = 1.2;
    segment(1, 1) = 1.0;
    p_condition->SetValue(Connector ? SURROGATE_BOUNDARY_FACE_COORDINATES : GAP_INTERFACE_COORDINATES, segment);
    if (!Connector) p_condition->SetValue(GAP_DIRICHLET_VALUES, Vector(3, 2.0));
    return p_condition;
}
}

KRATOS_TEST_CASE_IN_SUITE(LaplacianGapDirichletConstantAndLinear, KratosConvectionDiffusionFastSuite)
{
    Model model;
    auto& r_part = model.CreateModelPart("Scalar");
    auto p_condition = CreateLaplacianGapCondition(r_part, false);
    for (auto& r_node : r_part.Nodes()) r_node.FastGetSolutionStepValue(TEMPERATURE) = 2.0;
    Matrix lhs;
    Vector rhs;
    p_condition->CalculateLocalSystem(lhs, rhs, r_part.GetProcessInfo());
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        KRATOS_EXPECT_NEAR(rhs[i], 0.0, 1e-11);
        for (std::size_t j = 0; j < rhs.size(); ++j) KRATOS_EXPECT_NEAR(lhs(i,j), lhs(j,i), 1e-12);
    }
    const auto& r_rule = LineGaussLegendreIntegrationPoints3::IntegrationPoints();
    Vector values(r_rule.size());
    for (std::size_t g = 0; g < r_rule.size(); ++g) values[g] = 1.0 + 1.2 + 2.0 * 0.5 * (1.0 + r_rule[g].X());
    p_condition->SetValue(GAP_DIRICHLET_VALUES, values);
    for (auto& r_node : r_part.Nodes()) r_node.FastGetSolutionStepValue(TEMPERATURE) = 1.0 + r_node.X() + 2.0*r_node.Y();
    p_condition->CalculateLocalSystem(lhs, rhs, r_part.GetProcessInfo());
    double sum = 0.0;
    for (double value : rhs) sum += value;
    KRATOS_EXPECT_NEAR(sum, 1.0, 1e-11);
}

KRATOS_TEST_CASE_IN_SUITE(LaplacianGapConnectorSymmetryAndConsistency, KratosConvectionDiffusionFastSuite)
{
    Model model;
    auto& r_part = model.CreateModelPart("Scalar");
    auto p_condition = CreateLaplacianGapCondition(r_part, true);
    Matrix lhs;
    Vector rhs;
    for (auto& r_node : r_part.Nodes()) r_node.FastGetSolutionStepValue(TEMPERATURE) = 2.0;
    p_condition->CalculateLocalSystem(lhs, rhs, r_part.GetProcessInfo());
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        KRATOS_EXPECT_NEAR(rhs[i], 0.0, 1e-11);
        for (std::size_t j = 0; j < rhs.size(); ++j) KRATOS_EXPECT_NEAR(lhs(i,j), lhs(j,i), 1e-12);
    }
    for (auto& r_node : r_part.Nodes()) r_node.FastGetSolutionStepValue(TEMPERATURE) = r_node.X() + 2.0*r_node.Y();
    p_condition->CalculateLocalSystem(lhs, rhs, r_part.GetProcessInfo());
    double sum = 0.0;
    for (double value : rhs) sum += value;
    KRATOS_EXPECT_NEAR(sum, 0.0, 1e-11);
    KRATOS_EXPECT_GT(norm_2(rhs), 1e-8);
}

KRATOS_TEST_CASE_IN_SUITE(LaplacianGapConnectorSmallSegment, KratosConvectionDiffusionFastSuite)
{
    Model model;
    auto& r_part = model.CreateModelPart("Scalar");
    auto p_condition = CreateLaplacianGapCondition(r_part, true);
    Matrix segment = ZeroMatrix(2, 3);
    segment(0, 0) = segment(1, 0) = 1.0;
    segment(0, 1) = 0.5;
    segment(1, 1) = 0.5 + 1e-8;
    p_condition->SetValue(SURROGATE_BOUNDARY_FACE_COORDINATES, segment);
    for (auto& r_node : r_part.Nodes()) r_node.FastGetSolutionStepValue(TEMPERATURE) = 1.0;
    Matrix lhs;
    Vector rhs;
    p_condition->CalculateLocalSystem(lhs, rhs, r_part.GetProcessInfo());
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        KRATOS_EXPECT_NEAR(rhs[i], 0.0, 1e-12);
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            KRATOS_EXPECT_TRUE(std::isfinite(lhs(i,j)));
            KRATOS_EXPECT_LT(std::abs(lhs(i,j)), 1e-5);
        }
    }
}
}
