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
#include "test_setup_utilities/model_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <numbers>
#include <string>

using namespace std::numbers;
using namespace std::string_literals;

namespace Kratos::Testing
{

class StubElementForNodalExtrapolationTest : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(StubElementForNodalExtrapolationTest);

    StubElementForNodalExtrapolationTest(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, std::move(pGeometry))
    {
    }

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
    auto element_1  = Kratos::make_intrusive<StubElementForNodalExtrapolationTest>(1, geometry_1);
    model_part.AddElement(element_1);

    auto geometry_2 = std::make_shared<Quadrilateral2D4<Node>>(node_2, node_5, node_6, node_3);
    auto element_2  = Kratos::make_intrusive<StubElementForNodalExtrapolationTest>(2, geometry_2);
    model_part.AddElement(element_2);

    return model_part;
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

    dynamic_cast<StubElementForNodalExtrapolationTest&>(model_part.Elements()[1]).mIntegrationDoubleValues = {
        1.0, 1.0, 1.0, 1.0};
    dynamic_cast<StubElementForNodalExtrapolationTest&>(model_part.Elements()[2]).mIntegrationDoubleValues = {
        1.0, 1.0, 1.0, 1.0};

    const auto parameters = Parameters(R"(
    {
        "model_part_name"            : "MainModelPart",
        "echo_level"                 : 0,
        "list_of_variables"          : ["HYDRAULIC_HEAD"]
    })");

    GeoExtrapolateIntegrationPointValuesToNodesProcess process(model, parameters);
    process.ExecuteBeforeSolutionLoop();
    process.ExecuteFinalizeSolutionStep();

    for (auto& node : model_part.Nodes()) {
        KRATOS_EXPECT_NEAR(node.FastGetSolutionStepValue(HYDRAULIC_HEAD), 1.0, Defaults::absolute_tolerance);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TestExtrapolationProcess_ExtrapolatesCorrectlyForConstantFieldWithInactiveElement,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    //   This test uses the following two-element system, where the second element is inactive.
    //   4------3------6
    //   |  El1 |  El2 |
    //   1------2------5
    //

    Model model;
    auto& model_part = CreateModelPartWithTwoStubElements(model);

    dynamic_cast<StubElementForNodalExtrapolationTest&>(model_part.Elements()[1]).mIntegrationDoubleValues = {
        1.0, 1.0, 1.0, 1.0};
    dynamic_cast<StubElementForNodalExtrapolationTest&>(model_part.Elements()[2]).mIntegrationDoubleValues = {
        1.0, 1.0, 1.0, 1.0};
    model_part.Elements()[2].Set(ACTIVE, false);

    const auto parameters = Parameters(R"(
    {
        "model_part_name"            : "MainModelPart",
        "echo_level"                 : 0,
        "list_of_variables"          : ["HYDRAULIC_HEAD"]
    })");

    GeoExtrapolateIntegrationPointValuesToNodesProcess process(model, parameters);
    process.ExecuteBeforeSolutionLoop();
    process.ExecuteFinalizeSolutionStep();

    for (auto& node : model_part.Elements()[1].GetGeometry()) {
        KRATOS_EXPECT_NEAR(node.FastGetSolutionStepValue(HYDRAULIC_HEAD), 1.0, Defaults::absolute_tolerance);
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

    dynamic_cast<StubElementForNodalExtrapolationTest&>(model_part.Elements()[1]).mIntegrationDoubleValues = {
        1.0, 1.0, 1.0, 1.0};
    dynamic_cast<StubElementForNodalExtrapolationTest&>(model_part.Elements()[2]).mIntegrationDoubleValues = {
        2.0, 2.0, 2.0, 2.0};

    const auto parameters = Parameters(R"(
    {
        "model_part_name"            : "MainModelPart",
        "echo_level"                 : 0,
        "list_of_variables"          : ["HYDRAULIC_HEAD"]
    })");

    GeoExtrapolateIntegrationPointValuesToNodesProcess process(model, parameters);
    process.ExecuteBeforeSolutionLoop();
    process.ExecuteFinalizeSolutionStep();

    std::vector<double> expected_values = {1.0, 1.5, 1.5, 1.0, 2.0, 2.0};
    std::vector<double> actual_values;
    actual_values.reserve(model_part.Nodes().size());
    std::transform(model_part.Nodes().begin(), model_part.Nodes().end(), std::back_inserter(actual_values),
                   [](const auto& node) { return node.FastGetSolutionStepValue(HYDRAULIC_HEAD); });

    KRATOS_EXPECT_VECTOR_NEAR(actual_values, expected_values, Defaults::absolute_tolerance)
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
    dynamic_cast<StubElementForNodalExtrapolationTest&>(model_part.Elements()[1]).mIntegrationDoubleValues = {
        -inv_sqrt3, inv_sqrt3, inv_sqrt3, -inv_sqrt3};

    // Linear field in y between -1 and 1
    dynamic_cast<StubElementForNodalExtrapolationTest&>(model_part.Elements()[2]).mIntegrationDoubleValues = {
        -inv_sqrt3, -inv_sqrt3, inv_sqrt3, inv_sqrt3};

    const auto                                         parameters = Parameters(R"(
     {
         "model_part_name"            : "MainModelPart",
         "echo_level"                 : 0,
         "list_of_variables"          : ["HYDRAULIC_HEAD"]
     })");
    GeoExtrapolateIntegrationPointValuesToNodesProcess process(model, parameters);
    process.ExecuteBeforeSolutionLoop();
    process.ExecuteFinalizeSolutionStep();

    std::vector<double> expected_values = {-1, 0, 1, -1, -1, 1};
    std::vector<double> actual_values;
    actual_values.reserve(model_part.Nodes().size());
    std::transform(model_part.Nodes().begin(), model_part.Nodes().end(), std::back_inserter(actual_values),
                   [](const auto& node) { return node.FastGetSolutionStepValue(HYDRAULIC_HEAD); });

    KRATOS_EXPECT_VECTOR_NEAR(actual_values, expected_values, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(TestExtrapolationProcess_ExtrapolatesCorrectlyForLinearFields_EvenIfUnrelatedEmptyModelPartsAreSupplied,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    //   This test uses the following two-element system.
    //   4------3------6
    //   |  El1 |  El2 |
    //   1------2------5
    //
    Model model;
    auto& model_part = CreateModelPartWithTwoStubElements(model);
    model.CreateModelPart("foo"s);
    model.CreateModelPart("bar"s);

    // Linear field in x between -1 and 1
    dynamic_cast<StubElementForNodalExtrapolationTest&>(model_part.Elements()[1]).mIntegrationDoubleValues = {
        -inv_sqrt3, inv_sqrt3, inv_sqrt3, -inv_sqrt3};

    // Linear field in y between -1 and 1
    dynamic_cast<StubElementForNodalExtrapolationTest&>(model_part.Elements()[2]).mIntegrationDoubleValues = {
        -inv_sqrt3, -inv_sqrt3, inv_sqrt3, inv_sqrt3};

    const auto                                         parameters = Parameters(R"(
     {
         "model_part_name_list" : ["MainModelPart", "foo", "bar"],
         "echo_level"           : 0,
         "list_of_variables"    : ["HYDRAULIC_HEAD"]
     })");
    GeoExtrapolateIntegrationPointValuesToNodesProcess process(model, parameters);
    process.ExecuteBeforeSolutionLoop();
    process.ExecuteFinalizeSolutionStep();

    std::vector<double> expected_values = {-1, 0, 1, -1, -1, 1};
    std::vector<double> actual_values;
    actual_values.reserve(model_part.Nodes().size());
    std::transform(model_part.Nodes().begin(), model_part.Nodes().end(), std::back_inserter(actual_values),
                   [](const auto& node) { return node.FastGetSolutionStepValue(HYDRAULIC_HEAD); });

    KRATOS_EXPECT_VECTOR_NEAR(actual_values, expected_values, Defaults::absolute_tolerance)
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
    dynamic_cast<StubElementForNodalExtrapolationTest&>(model_part.Elements()[1]).mIntegrationMatrixValues = {
        ScalarMatrix(3, 3, -inv_sqrt3), ScalarMatrix(3, 3, inv_sqrt3),
        ScalarMatrix(3, 3, inv_sqrt3), ScalarMatrix(3, 3, -inv_sqrt3)};

    // Linear field in y between -1 and 1
    dynamic_cast<StubElementForNodalExtrapolationTest&>(model_part.Elements()[2]).mIntegrationMatrixValues = {
        ScalarMatrix(3, 3, -inv_sqrt3), ScalarMatrix(3, 3, -inv_sqrt3),
        ScalarMatrix(3, 3, inv_sqrt3), ScalarMatrix(3, 3, inv_sqrt3)};

    const auto parameters = Parameters(R"(
     {
         "model_part_name"            : "MainModelPart",
         "echo_level"                 : 0,
         "list_of_variables"          : ["CAUCHY_STRESS_TENSOR"]
     })");

    GeoExtrapolateIntegrationPointValuesToNodesProcess process(model, parameters);
    process.ExecuteBeforeSolutionLoop();
    process.ExecuteFinalizeSolutionStep();

    std::vector<Matrix> expected_values = {ScalarMatrix(3, 3, -1), ScalarMatrix(3, 3, 0),
                                           ScalarMatrix(3, 3, 1),  ScalarMatrix(3, 3, -1),
                                           ScalarMatrix(3, 3, -1), ScalarMatrix(3, 3, 1)};
    std::vector<Matrix> actual_values;
    actual_values.reserve(model_part.Nodes().size());
    std::transform(
        model_part.Nodes().begin(), model_part.Nodes().end(), std::back_inserter(actual_values),
        [](const auto& node) { return node.FastGetSolutionStepValue(CAUCHY_STRESS_TENSOR); });

    for (std::size_t i = 0; i < expected_values.size(); ++i) {
        KRATOS_EXPECT_MATRIX_NEAR(actual_values[i], expected_values[i], Defaults::absolute_tolerance)
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
    dynamic_cast<StubElementForNodalExtrapolationTest&>(model_part.Elements()[1]).mIntegrationVectorValues = {
        ScalarVector(6, -inv_sqrt3), ScalarVector(6, inv_sqrt3), ScalarVector(6, inv_sqrt3),
        ScalarVector(6, -inv_sqrt3)};

    // Linear field in y between -1 and 1
    dynamic_cast<StubElementForNodalExtrapolationTest&>(model_part.Elements()[2]).mIntegrationVectorValues = {
        ScalarVector(6, -inv_sqrt3), ScalarVector(6, -inv_sqrt3), ScalarVector(6, inv_sqrt3),
        ScalarVector(6, inv_sqrt3)};

    const auto parameters = Parameters(R"(
     {
         "model_part_name"            : "MainModelPart",
         "echo_level"                 : 0,
         "list_of_variables"          : ["CAUCHY_STRESS_VECTOR"]
     })");

    GeoExtrapolateIntegrationPointValuesToNodesProcess process(model, parameters);
    process.ExecuteBeforeSolutionLoop();
    process.ExecuteFinalizeSolutionStep();

    std::vector<Vector> expected_values = {ScalarVector(6, -1), ScalarVector(6, 0),
                                           ScalarVector(6, 1),  ScalarVector(6, -1),
                                           ScalarVector(6, -1), ScalarVector(6, 1)};
    std::vector<Vector> actual_values;
    actual_values.reserve(model_part.Nodes().size());
    std::transform(
        model_part.Nodes().begin(), model_part.Nodes().end(), std::back_inserter(actual_values),
        [](const auto& node) { return node.FastGetSolutionStepValue(CAUCHY_STRESS_VECTOR); });

    for (std::size_t i = 0; i < expected_values.size(); ++i) {
        KRATOS_EXPECT_VECTOR_NEAR(actual_values[i], expected_values[i], Defaults::absolute_tolerance)
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
    dynamic_cast<StubElementForNodalExtrapolationTest&>(model_part.Elements()[1]).mIntegrationArrayValues = {
        ScalarVector(3, -inv_sqrt3), ScalarVector(3, inv_sqrt3), ScalarVector(3, inv_sqrt3),
        ScalarVector(3, -inv_sqrt3)};

    // Linear field in y between -1 and 1
    dynamic_cast<StubElementForNodalExtrapolationTest&>(model_part.Elements()[2]).mIntegrationArrayValues = {
        ScalarVector(3, -inv_sqrt3), ScalarVector(3, -inv_sqrt3), ScalarVector(3, inv_sqrt3),
        ScalarVector(3, inv_sqrt3)};

    const auto parameters = Parameters(R"(
     {
         "model_part_name_list"       : ["MainModelPart"],
         "echo_level"                 : 0,
         "list_of_variables"          : ["FLUID_FLUX_VECTOR"]
     })");

    GeoExtrapolateIntegrationPointValuesToNodesProcess process(model, parameters);
    process.ExecuteBeforeSolutionLoop();
    process.ExecuteFinalizeSolutionStep();

    std::vector<Vector>              expected_values = {ScalarVector(3, -1), ScalarVector(3, 0),
                                                        ScalarVector(3, 1),  ScalarVector(3, -1),
                                                        ScalarVector(3, -1), ScalarVector(3, 1)};
    std::vector<array_1d<double, 3>> actual_values;
    actual_values.reserve(model_part.Nodes().size());
    std::transform(model_part.Nodes().begin(), model_part.Nodes().end(), std::back_inserter(actual_values),
                   [](const auto& node) { return node.FastGetSolutionStepValue(FLUID_FLUX_VECTOR); });

    for (std::size_t i = 0; i < expected_values.size(); ++i) {
        KRATOS_EXPECT_VECTOR_NEAR(actual_values[i], expected_values[i], Defaults::absolute_tolerance)
    }
}

KRATOS_TEST_CASE_IN_SUITE(TestExtrapolationProcess_ExtrapolatesCorrectlyWhenNodesAreSharedBetweenModelParts,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    //   This test uses the following two-element system.
    //   4------3------6
    //   |  El1 |  El2 |
    //   1------2------5
    //
    Model model;

    auto& r_left_model_part = model.CreateModelPart("Left");
    r_left_model_part.AddNodalSolutionStepVariable(HYDRAULIC_HEAD);
    const auto nodes_of_element_1 = ModelSetupUtilities::CreateNodes(
        r_left_model_part,
        {{1, {0.0, 0.0, 0.0}}, {2, {1.0, 0.0, 0.0}}, {3, {1.0, 1.0, 0.0}}, {4, {0.0, 1.0, 0.0}}});
    r_left_model_part.AddElement(make_intrusive<StubElementForNodalExtrapolationTest>(
        1, std::make_shared<Quadrilateral2D4<Node>>(nodes_of_element_1)));

    auto& r_right_model_part = model.CreateModelPart("Right");
    r_right_model_part.AddNodalSolutionStepVariable(HYDRAULIC_HEAD);
    r_right_model_part.AddNode(r_left_model_part.pGetNode(2));
    r_right_model_part.AddNode(r_left_model_part.pGetNode(3));
    ModelSetupUtilities::CreateNodes(r_right_model_part, {{5, {2.0, 0.0, 0.0}}, {6, {2.0, 1.0, 0.0}}});
    // Unfortunately, it appeared that `PointerVector<Node>::insert` does not work properly, so
    // we're limited to using `PointerVector<Node>::push_back`
    auto nodes_of_element_2 = PointerVector<Node>{};
    nodes_of_element_2.push_back(r_right_model_part.pGetNode(2));
    nodes_of_element_2.push_back(r_right_model_part.pGetNode(5));
    nodes_of_element_2.push_back(r_right_model_part.pGetNode(6));
    nodes_of_element_2.push_back(r_right_model_part.pGetNode(3));
    r_right_model_part.AddElement(make_intrusive<StubElementForNodalExtrapolationTest>(
        2, std::make_shared<Quadrilateral2D4<Node>>(nodes_of_element_2)));

    dynamic_cast<StubElementForNodalExtrapolationTest&>(r_left_model_part.Elements()[1]).mIntegrationDoubleValues = {
        1.0, 1.0, 1.0, 1.0};
    dynamic_cast<StubElementForNodalExtrapolationTest&>(r_right_model_part.Elements()[2]).mIntegrationDoubleValues = {
        1.0, 1.0, 1.0, 1.0};

    const auto parameters = Parameters(R"(
    {
        "model_part_name_list"       : ["Left", "Right"],
        "echo_level"                 : 0,
        "list_of_variables"          : ["HYDRAULIC_HEAD"]
    })");

    GeoExtrapolateIntegrationPointValuesToNodesProcess process(model, parameters);
    process.ExecuteBeforeSolutionLoop();
    process.ExecuteFinalizeSolutionStep();

    for (auto& r_node : r_left_model_part.Nodes()) {
        KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(HYDRAULIC_HEAD), 1.0, Defaults::absolute_tolerance);
    }

    for (auto& r_node : r_right_model_part.Nodes()) {
        KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(HYDRAULIC_HEAD), 1.0, Defaults::absolute_tolerance);
    }
}

KRATOS_TEST_CASE_IN_SUITE(TestExtrapolationProcess_ExtrapolatesCorrectlyWhenModelPartsWithSharedNodesAreUsedByDifferentProcesses,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    //   This test uses the following three-element system.
    //          7------8
    //          |  El3 |
    //   4------3------6
    //   |  El1 |  El2 |
    //   1------2------5
    //

    Model model;
    auto& r_bottom_model_part = model.CreateModelPart("Bottom");
    r_bottom_model_part.AddNodalSolutionStepVariable(HYDRAULIC_HEAD);
    r_bottom_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto node_1 = r_bottom_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto node_2 = r_bottom_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto node_3 = r_bottom_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    auto node_4 = r_bottom_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
    auto node_5 = r_bottom_model_part.CreateNewNode(5, 2.0, 0.0, 0.0);
    auto node_6 = r_bottom_model_part.CreateNewNode(6, 2.0, 1.0, 0.0);

    auto geometry_1 = std::make_shared<Quadrilateral2D4<Node>>(node_1, node_2, node_3, node_4);
    auto element_1  = Kratos::make_intrusive<StubElementForNodalExtrapolationTest>(1, geometry_1);
    r_bottom_model_part.AddElement(element_1);

    auto geometry_2 = std::make_shared<Quadrilateral2D4<Node>>(node_2, node_5, node_6, node_3);
    auto element_2  = Kratos::make_intrusive<StubElementForNodalExtrapolationTest>(2, geometry_2);
    r_bottom_model_part.AddElement(element_2);

    auto& r_top_model_part = model.CreateModelPart("Top");
    r_top_model_part.AddNodalSolutionStepVariable(HYDRAULIC_HEAD);
    r_top_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    r_top_model_part.AddNode(node_3);
    r_top_model_part.AddNode(node_6);
    auto node_7     = r_top_model_part.CreateNewNode(7, 1.0, 2.0, 0.0);
    auto node_8     = r_top_model_part.CreateNewNode(8, 2.0, 2.0, 0.0);
    auto geometry_3 = std::make_shared<Quadrilateral2D4<Node>>(node_3, node_6, node_8, node_7);
    auto element_3  = Kratos::make_intrusive<StubElementForNodalExtrapolationTest>(3, geometry_3);
    r_top_model_part.AddElement(element_3);

    dynamic_cast<StubElementForNodalExtrapolationTest&>(r_bottom_model_part.Elements()[1]).mIntegrationDoubleValues = {
        1.0, 1.0, 1.0, 1.0};
    dynamic_cast<StubElementForNodalExtrapolationTest&>(r_bottom_model_part.Elements()[2]).mIntegrationDoubleValues = {
        1.0, 1.0, 1.0, 1.0};
    dynamic_cast<StubElementForNodalExtrapolationTest&>(r_top_model_part.Elements()[3]).mIntegrationDoubleValues = {
        1.0, 1.0, 1.0, 1.0};

    const auto                                         parameters_1 = Parameters(R"(
    {
        "model_part_name"       : "Bottom",
        "echo_level"                 : 0,
        "list_of_variables"          : ["HYDRAULIC_HEAD"]
    })");
    GeoExtrapolateIntegrationPointValuesToNodesProcess process_1(model, parameters_1);
    process_1.ExecuteBeforeSolutionLoop();

    const auto                                         parameters_2 = Parameters(R"(
    {
        "model_part_name"       : "Top",
        "echo_level"                 : 0,
        "list_of_variables"          : ["HYDRAULIC_DISCHARGE"]
    })");
    GeoExtrapolateIntegrationPointValuesToNodesProcess process_2(model, parameters_2);
    process_2.ExecuteBeforeSolutionLoop();

    process_1.ExecuteFinalizeSolutionStep();

    for (auto& node : r_bottom_model_part.Nodes()) {
        EXPECT_NEAR(node.FastGetSolutionStepValue(HYDRAULIC_HEAD), 1.0, Defaults::absolute_tolerance)
            << "Hydraulic head at node " << node.Id();
    }

    process_2.ExecuteFinalizeSolutionStep();

    for (auto& node : r_top_model_part.Nodes()) {
        EXPECT_NEAR(node.FastGetSolutionStepValue(HYDRAULIC_DISCHARGE), 1.0, Defaults::absolute_tolerance)
            << "Hydraulic discharge at node " << node.Id();
    }
}

KRATOS_TEST_CASE_IN_SUITE(CheckInfoGeoExtrapolateIntegrationPointValuesToNodesProcess, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    model.CreateModelPart("foo");

    const auto parameters = Parameters(R"(
     {
         "model_part_name"            : "foo",
         "list_of_variables"          : ["FLUID_FLUX_VECTOR"]
     })");

    const GeoExtrapolateIntegrationPointValuesToNodesProcess process(model, parameters);

    // Act & assert
    KRATOS_EXPECT_EQ(process.Info(), "GeoExtrapolateIntegrationPointValuesToNodesProcess");
}

} // namespace Kratos::Testing
