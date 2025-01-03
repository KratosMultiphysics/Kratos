// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "containers/model.h"
#include "custom_processes/geo_extrapolate_integration_point_values_to_nodes_process.h"
#include "geo_mechanics_application_variables.h"
#include "geometries/quadrilateral_2d_4.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

class StubElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(StubElement);

    StubElement(IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry) {}

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override
    {
        rOutput = mIntegrationDoubleValues;
    }

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override
    {
        rOutput = mIntegrationMatrixValues;
    }

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                      std::vector<array_1d<double, 3>>&    rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        rOutput = mIntegrationArrayValues;
    }

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override
    {
        rOutput = mIntegrationVectorValues;
    }

    using Element::CalculateOnIntegrationPoints;

    IntegrationMethod GetIntegrationMethod() const override
    {
        return GeometryData::IntegrationMethod::GI_GAUSS_2;
    }

    std::vector<double>              mIntegrationDoubleValues = {};
    std::vector<Matrix>              mIntegrationMatrixValues = {};
    std::vector<Vector>              mIntegrationVectorValues = {};
    std::vector<array_1d<double, 3>> mIntegrationArrayValues  = {};
};

ModelPart& CreateModelPartWithTwoStubElements(Model& model)
{
    auto& model_part = model.CreateModelPart("MainModelPart");
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(HYDRAULIC_HEAD);
    model_part.AddNodalSolutionStepVariable(CAUCHY_STRESS_TENSOR);
    model_part.AddNodalSolutionStepVariable(CAUCHY_STRESS_VECTOR);
    model_part.AddNodalSolutionStepVariable(FLUID_FLUX_VECTOR);
    auto node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto node_2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto node_3 = model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    auto node_4 = model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
    auto node_5 = model_part.CreateNewNode(5, 2.0, 0.0, 0.0);
    auto node_6 = model_part.CreateNewNode(6, 2.0, 1.0, 0.0);

    auto geometry_1 = std::make_shared<Quadrilateral2D4<Node>>(node_1, node_2, node_3, node_4);
    auto element_1  = Kratos::make_intrusive<StubElement>(1, geometry_1);
    model_part.AddElement(element_1);

    auto geometry_2 = std::make_shared<Quadrilateral2D4<Node>>(node_2, node_5, node_6, node_3);
    auto element_2  = Kratos::make_intrusive<StubElement>(2, geometry_2);
    model_part.AddElement(element_2);

    return model_part;
}

KRATOS_TEST_CASE_IN_SUITE(TestExtrapolationProcess_InitializesNodalArea, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    //   This test uses the following two-element system.
    //   4------3------6
    //   |  El1 |  El2 |
    //   1------2------5
    //
    Model model;
    auto& model_part = CreateModelPartWithTwoStubElements(model);

    auto parameters = Parameters(R"(
    {
        "model_part_name"            : "MainModelPart",
        "echo_level"                 : 0,
        "list_of_variables"          : ["HYDRAULIC_HEAD"]
    })");

    GeoExtrapolateIntegrationPointValuesToNodesProcess process(model_part, parameters);

    process.ExecuteBeforeSolutionLoop();

    std::vector<double> expected_values = {1.0, 2.0, 2.0, 1.0, 1.0, 1.0};
    std::vector<double> actual_values;
    std::transform(model_part.Nodes().begin(), model_part.Nodes().end(), std::back_inserter(actual_values),
                   [](const auto& node) { return node.GetValue(NODAL_AREA); });

    KRATOS_EXPECT_VECTOR_NEAR(actual_values, expected_values, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(TestExtrapolationProcess_ExtrapolatesCorrectlyForConstantField,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    //   This test uses the following two-element system.
    //   4------3------6
    //   |  El1 |  El2 |
    //   1------2------5
    //

    Model model;
    auto& model_part = CreateModelPartWithTwoStubElements(model);

    dynamic_cast<StubElement&>(model_part.Elements()[1]).mIntegrationDoubleValues = {1.0, 1.0, 1.0, 1.0};
    dynamic_cast<StubElement&>(model_part.Elements()[2]).mIntegrationDoubleValues = {1.0, 1.0, 1.0, 1.0};

    auto parameters = Parameters(R"(
    {
        "model_part_name"            : "MainModelPart",
        "echo_level"                 : 0,
        "list_of_variables"          : ["HYDRAULIC_HEAD"]
    })");

    GeoExtrapolateIntegrationPointValuesToNodesProcess process(model_part, parameters);
    process.Execute();

    for (auto& node : model_part.Nodes()) {
        KRATOS_EXPECT_NEAR(node.FastGetSolutionStepValue(HYDRAULIC_HEAD), 1.0, 1e-6);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TestExtrapolationProcess_ExtrapolatesCorrectlyForTwoConstantButDifferentFields,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    //   This test uses the following two-element system.
    //   4------3------6
    //   |  El1 |  El2 |
    //   1------2------5
    //
    Model model;
    auto& model_part = CreateModelPartWithTwoStubElements(model);

    dynamic_cast<StubElement&>(model_part.Elements()[1]).mIntegrationDoubleValues = {1.0, 1.0, 1.0, 1.0};
    dynamic_cast<StubElement&>(model_part.Elements()[2]).mIntegrationDoubleValues = {2.0, 2.0, 2.0, 2.0};

    auto parameters = Parameters(R"(
    {
        "model_part_name"            : "MainModelPart",
        "echo_level"                 : 0,
        "list_of_variables"          : ["HYDRAULIC_HEAD"]
    })");

    GeoExtrapolateIntegrationPointValuesToNodesProcess process(model_part, parameters);
    process.Execute();

    std::vector<double> expected_values = {1.0, 1.5, 1.5, 1.0, 2.0, 2.0};
    std::vector<double> actual_values;
    std::transform(model_part.Nodes().begin(), model_part.Nodes().end(), std::back_inserter(actual_values),
                   [](const auto& node) { return node.FastGetSolutionStepValue(HYDRAULIC_HEAD); });

    KRATOS_EXPECT_VECTOR_NEAR(actual_values, expected_values, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(TestExtrapolationProcess_ExtrapolatesCorrectlyForLinearFields,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    //   This test uses the following two-element system.
    //   4------3------6
    //   |  El1 |  El2 |
    //   1------2------5
    //
    Model model;
    auto& model_part = CreateModelPartWithTwoStubElements(model);

    // Linear field in x between -1 and 1
    dynamic_cast<StubElement&>(model_part.Elements()[1]).mIntegrationDoubleValues = {
        -0.57735, 0.57735, 0.57735, -0.57735};

    // Linear field in y between -1 and 1
    dynamic_cast<StubElement&>(model_part.Elements()[2]).mIntegrationDoubleValues = {
        -0.57735, -0.57735, 0.57735, 0.57735};

    auto                                               parameters = Parameters(R"(
     {
         "model_part_name"            : "MainModelPart",
         "echo_level"                 : 0,
         "list_of_variables"          : ["HYDRAULIC_HEAD"]
     })");
    GeoExtrapolateIntegrationPointValuesToNodesProcess process(model_part, parameters);
    process.Execute();

    std::vector<double> expected_values = {-1, 0, 1, -1, -1, 1};
    std::vector<double> actual_values;
    std::transform(model_part.Nodes().begin(), model_part.Nodes().end(), std::back_inserter(actual_values),
                   [](const auto& node) { return node.FastGetSolutionStepValue(HYDRAULIC_HEAD); });

    KRATOS_EXPECT_VECTOR_NEAR(actual_values, expected_values, 1e-6)
}

KRATOS_TEST_CASE_IN_SUITE(TestExtrapolationProcess_ExtrapolatesMatrixCorrectlyForLinearFields,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    //   This test uses the following two-element system.
    //   4------3------6
    //   |  El1 |  El2 |
    //   1------2------5
    //
    Model model;
    auto& model_part = CreateModelPartWithTwoStubElements(model);

    // Linear field in x between -1 and 1
    dynamic_cast<StubElement&>(model_part.Elements()[1]).mIntegrationMatrixValues = {
        ScalarMatrix(3, 3, -0.57735), ScalarMatrix(3, 3, 0.57735), ScalarMatrix(3, 3, 0.57735),
        ScalarMatrix(3, 3, -0.57735)};

    // Linear field in y between -1 and 1
    dynamic_cast<StubElement&>(model_part.Elements()[2]).mIntegrationMatrixValues = {
        ScalarMatrix(3, 3, -0.57735), ScalarMatrix(3, 3, -0.57735), ScalarMatrix(3, 3, 0.57735),
        ScalarMatrix(3, 3, 0.57735)};

    auto parameters = Parameters(R"(
     {
         "model_part_name"            : "MainModelPart",
         "echo_level"                 : 0,
         "list_of_variables"          : ["CAUCHY_STRESS_TENSOR"]
     })");

    GeoExtrapolateIntegrationPointValuesToNodesProcess process(model_part, parameters);
    process.Execute();

    std::vector<Matrix> expected_values = {ScalarMatrix(3, 3, -1), ScalarMatrix(3, 3, 0),
                                           ScalarMatrix(3, 3, 1),  ScalarMatrix(3, 3, -1),
                                           ScalarMatrix(3, 3, -1), ScalarMatrix(3, 3, 1)};
    std::vector<Matrix> actual_values;
    std::transform(
        model_part.Nodes().begin(), model_part.Nodes().end(), std::back_inserter(actual_values),
        [](const auto& node) { return node.FastGetSolutionStepValue(CAUCHY_STRESS_TENSOR); });

    for (std::size_t i = 0; i < expected_values.size(); ++i) {
        KRATOS_EXPECT_MATRIX_NEAR(actual_values[i], expected_values[i], 1e-6)
    }
}

KRATOS_TEST_CASE_IN_SUITE(TestExtrapolationProcess_ExtrapolatesVectorCorrectlyForLinearFields,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    //   This test uses the following two-element system.
    //   4------3------6
    //   |  El1 |  El2 |
    //   1------2------5
    //
    Model model;
    auto& model_part = CreateModelPartWithTwoStubElements(model);

    // Linear field in x between -1 and 1
    dynamic_cast<StubElement&>(model_part.Elements()[1]).mIntegrationVectorValues = {
        ScalarVector(6, -0.57735), ScalarVector(6, 0.57735), ScalarVector(6, 0.57735),
        ScalarVector(6, -0.57735)};

    // Linear field in y between -1 and 1
    dynamic_cast<StubElement&>(model_part.Elements()[2]).mIntegrationVectorValues = {
        ScalarVector(6, -0.57735), ScalarVector(6, -0.57735), ScalarVector(6, 0.57735), ScalarVector(6, 0.57735)};

    auto parameters = Parameters(R"(
     {
         "model_part_name"            : "MainModelPart",
         "echo_level"                 : 0,
         "list_of_variables"          : ["CAUCHY_STRESS_VECTOR"]
     })");

    GeoExtrapolateIntegrationPointValuesToNodesProcess process(model_part, parameters);
    process.Execute();

    std::vector<Vector> expected_values = {ScalarVector(6, -1), ScalarVector(6, 0),
                                           ScalarVector(6, 1),  ScalarVector(6, -1),
                                           ScalarVector(6, -1), ScalarVector(6, 1)};
    std::vector<Vector> actual_values;
    std::transform(
        model_part.Nodes().begin(), model_part.Nodes().end(), std::back_inserter(actual_values),
        [](const auto& node) { return node.FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR); });

    for (std::size_t i = 0; i < expected_values.size(); ++i) {
        KRATOS_EXPECT_VECTOR_NEAR(actual_values[i], expected_values[i], 1e-6)
    }
}

KRATOS_TEST_CASE_IN_SUITE(TestExtrapolationProcess_ExtrapolatesArrayCorrectlyForLinearFields,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    //   This test uses the following two-element system.
    //   4------3------6
    //   |  El1 |  El2 |
    //   1------2------5
    //
    Model model;
    auto& model_part = CreateModelPartWithTwoStubElements(model);

    // Linear field in x between -1 and 1
    dynamic_cast<StubElement&>(model_part.Elements()[1]).mIntegrationArrayValues = {
        ScalarVector(3, -0.57735), ScalarVector(3, 0.57735), ScalarVector(3, 0.57735),
        ScalarVector(3, -0.57735)};

    // Linear field in y between -1 and 1
    dynamic_cast<StubElement&>(model_part.Elements()[2]).mIntegrationArrayValues = {
        ScalarVector(3, -0.57735), ScalarVector(3, -0.57735), ScalarVector(3, 0.57735), ScalarVector(3, 0.57735)};

    auto parameters = Parameters(R"(
     {
         "model_part_name"            : "MainModelPart",
         "echo_level"                 : 0,
         "list_of_variables"          : ["FLUID_FLUX_VECTOR"]
     })");

    GeoExtrapolateIntegrationPointValuesToNodesProcess process(model_part, parameters);
    process.Execute();

    std::vector<Vector>              expected_values = {ScalarVector(3, -1), ScalarVector(3, 0),
                                                        ScalarVector(3, 1),  ScalarVector(3, -1),
                                                        ScalarVector(3, -1), ScalarVector(3, 1)};
    std::vector<array_1d<double, 3>> actual_values;
    std::transform(model_part.Nodes().begin(), model_part.Nodes().end(), std::back_inserter(actual_values),
                   [](const auto& node) { return node.FastGetSolutionStepValue(FLUID_FLUX_VECTOR); });

    for (std::size_t i = 0; i < expected_values.size(); ++i) {
        KRATOS_EXPECT_VECTOR_NEAR(actual_values[i], expected_values[i], 1e-6)
    }
}

} // namespace Kratos::Testing
