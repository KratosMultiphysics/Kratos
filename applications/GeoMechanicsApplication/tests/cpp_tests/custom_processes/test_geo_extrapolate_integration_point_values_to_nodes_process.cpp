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

namespace
{

using namespace Kratos;

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

ModelPart& CreateModelPartWithTwoStubElements(Model& model, const VariableData& rVariable)
{
    auto& r_result = model.CreateModelPart("MainModelPart");
    r_result.AddNodalSolutionStepVariable(rVariable);

    auto node_1 = r_result.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto node_2 = r_result.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto node_3 = r_result.CreateNewNode(3, 1.0, 1.0, 0.0);
    auto node_4 = r_result.CreateNewNode(4, 0.0, 1.0, 0.0);
    auto node_5 = r_result.CreateNewNode(5, 2.0, 0.0, 0.0);
    auto node_6 = r_result.CreateNewNode(6, 2.0, 1.0, 0.0);

    auto geometry_1 = std::make_shared<Quadrilateral2D4<Node>>(node_1, node_2, node_3, node_4);
    auto element_1  = Kratos::make_intrusive<StubElementForNodalExtrapolationTest>(1, geometry_1);
    r_result.AddElement(element_1);

    auto geometry_2 = std::make_shared<Quadrilateral2D4<Node>>(node_2, node_5, node_6, node_3);
    auto element_2  = Kratos::make_intrusive<StubElementForNodalExtrapolationTest>(2, geometry_2);
    r_result.AddElement(element_2);

    return r_result;
}

using ConstModelPartReferenceVector = std::vector<std::reference_wrapper<const ModelPart>>;

std::string MakeModelPartNameFrom(const ModelPart& rModelPart)
{
    return R"(        "model_part_name": ")" + rModelPart.Name() + '"';
}

std::string MakeModelPartNameListFrom(const ConstModelPartReferenceVector& rModelPartRefs)
{
    auto ss = std::stringstream{};
    ss << R"(        "model_part_name_list": [)";
    if (!rModelPartRefs.empty()) {
        // The names need to be separated by commas, but we must avoid a trailing comma. Therefore,
        // using `std::ranges::copy` with an `std::ostream_iterator` with ", " as the delimiter is
        // not good enough.
        auto quoted_name = [](const auto& rModelPartRef) {
            return '"' + rModelPartRef.get().Name() + '"';
        };
        ss << quoted_name(rModelPartRefs.front());
        for (auto it = rModelPartRefs.begin() + 1; it != rModelPartRefs.end(); ++it) {
            ss << ", " << quoted_name(*it);
        }
    }
    ss << "]";

    return ss.str();
}

std::string MakeVariableListFrom(const VariableData& rVariableData)
{
    return R"(        "list_of_variables": [")" + rVariableData.Name() + R"("])";
}

Parameters CreateExtrapolationProcessSettings(const std::string& rModelPartsSetting, const VariableData& rVariableData)
{
    // clang-format off
    return Parameters{
        "    {\n"s +
            rModelPartsSetting + ",\n"s +
            R"(        "echo_level": 0)" + ",\n" +
            MakeVariableListFrom(rVariableData) + '\n' +
        "    }\n"s
    };
    // clang-format on
}

Parameters CreateExtrapolationProcessSettings(const ModelPart& rModelPart, const VariableData& rVariableData)
{
    return CreateExtrapolationProcessSettings(MakeModelPartNameFrom(rModelPart), rVariableData);
}

Parameters CreateExtrapolationProcessSettings(const ConstModelPartReferenceVector& rModelPartRefs,
                                              const VariableData&                  rVariableData)
{
    return CreateExtrapolationProcessSettings(MakeModelPartNameListFrom(rModelPartRefs), rVariableData);
}

void SetIntegrationPointValues(Element& rElement, const std::vector<double>& rValues)
{
    dynamic_cast<StubElementForNodalExtrapolationTest&>(rElement).mIntegrationDoubleValues = rValues;
}

void SetIntegrationPointValues(Element& rElement, const std::vector<Vector>& rValues)
{
    dynamic_cast<StubElementForNodalExtrapolationTest&>(rElement).mIntegrationVectorValues = rValues;
}

void SetIntegrationPointValues(Element& rElement, const std::vector<Matrix>& rValues)
{
    dynamic_cast<StubElementForNodalExtrapolationTest&>(rElement).mIntegrationMatrixValues = rValues;
}

void BuildAndRunExtrapolationProcess(Model& rModel, const Parameters& rProcessSettings)
{
    auto extrapolation_process = GeoExtrapolateIntegrationPointValuesToNodesProcess{rModel, rProcessSettings};
    extrapolation_process.ExecuteBeforeSolutionLoop();
    extrapolation_process.ExecuteFinalizeSolutionStep();
}

template <typename NodeContainerType>
void AssertNodalValues(const NodeContainerType&   rNodes,
                       const Variable<double>&    rVariable,
                       const std::vector<double>& rExpectedNodalValues)
{
    auto actual_nodal_values = std::vector<double>{};
    actual_nodal_values.reserve(rNodes.size());
    auto get_nodal_value = [&rVariable](const auto& rNode) {
        return rNode.FastGetSolutionStepValue(rVariable);
    };
    std::ranges::transform(rNodes, std::back_inserter(actual_nodal_values), get_nodal_value);

    KRATOS_EXPECT_VECTOR_NEAR(actual_nodal_values, rExpectedNodalValues, Testing::Defaults::absolute_tolerance);
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TestExtrapolationProcess_ExtrapolatesCorrectlyForConstantField,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    //   This test uses the following two-element system.
    //   4------3------6
    //   |  El1 |  El2 |
    //   1------2------5

    auto        model           = Model{};
    const auto& r_test_variable = HYDRAULIC_HEAD;
    auto&       r_model_part    = CreateModelPartWithTwoStubElements(model, r_test_variable);

    SetIntegrationPointValues(r_model_part.Elements()[1], std::vector(4, 1.0));
    SetIntegrationPointValues(r_model_part.Elements()[2], std::vector(4, 1.0));

    BuildAndRunExtrapolationProcess(model, CreateExtrapolationProcessSettings(r_model_part, r_test_variable));

    AssertNodalValues(r_model_part.Nodes(), r_test_variable, std::vector(6, 1.0));
}

KRATOS_TEST_CASE_IN_SUITE(TestExtrapolationProcess_ExtrapolatesCorrectlyForConstantFieldWithInactiveElement,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    //   This test uses the following two-element system, where the second element is inactive.
    //   4------3------6
    //   |  El1 |  El2 |
    //   1------2------5

    auto        model           = Model{};
    const auto& r_test_variable = HYDRAULIC_HEAD;
    auto&       r_model_part    = CreateModelPartWithTwoStubElements(model, r_test_variable);

    SetIntegrationPointValues(r_model_part.Elements()[1], std::vector(4, 1.0));
    SetIntegrationPointValues(r_model_part.Elements()[2], std::vector(4, 1.0));

    r_model_part.Elements()[2].Set(ACTIVE, false);

    BuildAndRunExtrapolationProcess(model, CreateExtrapolationProcessSettings(r_model_part, r_test_variable));

    AssertNodalValues(r_model_part.Elements()[1].GetGeometry(), r_test_variable, std::vector(4, 1.0));
}

KRATOS_TEST_CASE_IN_SUITE(TestExtrapolationProcess_ExtrapolatesCorrectlyForTwoConstantButDifferentFields,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    //   This test uses the following two-element system.
    //   4------3------6
    //   |  El1 |  El2 |
    //   1------2------5

    auto        model           = Model{};
    const auto& r_test_variable = HYDRAULIC_HEAD;
    auto&       r_model_part    = CreateModelPartWithTwoStubElements(model, r_test_variable);

    SetIntegrationPointValues(r_model_part.Elements()[1], std::vector(4, 1.0));
    SetIntegrationPointValues(r_model_part.Elements()[2], std::vector(4, 2.0));

    BuildAndRunExtrapolationProcess(model, CreateExtrapolationProcessSettings(r_model_part, r_test_variable));

    AssertNodalValues(r_model_part.Nodes(), r_test_variable, {1.0, 1.5, 1.5, 1.0, 2.0, 2.0});
}

KRATOS_TEST_CASE_IN_SUITE(TestExtrapolationProcess_ExtrapolatesCorrectlyForLinearFields,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    //   This test uses the following two-element system.
    //   4------3------6
    //   |  El1 |  El2 |
    //   1------2------5

    auto        model           = Model{};
    const auto& r_test_variable = HYDRAULIC_HEAD;
    auto&       r_model_part    = CreateModelPartWithTwoStubElements(model, r_test_variable);

    // Linear field in x between -1 and 1
    SetIntegrationPointValues(r_model_part.Elements()[1], {-inv_sqrt3, inv_sqrt3, inv_sqrt3, -inv_sqrt3});

    // Linear field in y between -1 and 1
    SetIntegrationPointValues(r_model_part.Elements()[2], {-inv_sqrt3, -inv_sqrt3, inv_sqrt3, inv_sqrt3});

    BuildAndRunExtrapolationProcess(model, CreateExtrapolationProcessSettings(r_model_part, r_test_variable));

    AssertNodalValues(r_model_part.Nodes(), r_test_variable, {-1.0, 0.0, 1.0, -1.0, -1.0, 1.0});
}

KRATOS_TEST_CASE_IN_SUITE(TestExtrapolationProcess_ExtrapolatesCorrectlyForLinearFields_EvenIfUnrelatedEmptyModelPartsAreSupplied,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    //   This test uses the following two-element system.
    //   4------3------6
    //   |  El1 |  El2 |
    //   1------2------5

    auto        model             = Model{};
    const auto& r_test_variable   = HYDRAULIC_HEAD;
    auto&       r_main_model_part = CreateModelPartWithTwoStubElements(model, r_test_variable);
    const auto& r_foo_model_part  = model.CreateModelPart("foo"s);
    const auto& r_bar_model_part  = model.CreateModelPart("bar"s);

    // Linear field in x between -1 and 1
    SetIntegrationPointValues(r_main_model_part.Elements()[1], {-inv_sqrt3, inv_sqrt3, inv_sqrt3, -inv_sqrt3});

    // Linear field in y between -1 and 1
    SetIntegrationPointValues(r_main_model_part.Elements()[2], {-inv_sqrt3, -inv_sqrt3, inv_sqrt3, inv_sqrt3});

    BuildAndRunExtrapolationProcess(
        model, CreateExtrapolationProcessSettings(
                   {std::cref(r_main_model_part), std::cref(r_foo_model_part), std::cref(r_bar_model_part)},
                   r_test_variable));

    AssertNodalValues(r_main_model_part.Nodes(), r_test_variable, {-1.0, 0.0, 1.0, -1.0, -1.0, 1.0});
}

KRATOS_TEST_CASE_IN_SUITE(TestExtrapolationProcess_ExtrapolatesMatrixCorrectlyForLinearFields,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    //   This test uses the following two-element system.
    //   4------3------6
    //   |  El1 |  El2 |
    //   1------2------5

    auto        model           = Model{};
    const auto& r_test_variable = CAUCHY_STRESS_TENSOR;
    auto&       r_model_part    = CreateModelPartWithTwoStubElements(model, r_test_variable);

    // Linear field in x between -1 and 1
    SetIntegrationPointValues(r_model_part.Elements()[1],
                              {ScalarMatrix(3, 3, -inv_sqrt3), ScalarMatrix(3, 3, inv_sqrt3),
                               ScalarMatrix(3, 3, inv_sqrt3), ScalarMatrix(3, 3, -inv_sqrt3)});

    // Linear field in y between -1 and 1
    SetIntegrationPointValues(r_model_part.Elements()[2],
                              {ScalarMatrix(3, 3, -inv_sqrt3), ScalarMatrix(3, 3, -inv_sqrt3),
                               ScalarMatrix(3, 3, inv_sqrt3), ScalarMatrix(3, 3, inv_sqrt3)});

    BuildAndRunExtrapolationProcess(model, CreateExtrapolationProcessSettings(r_model_part, r_test_variable));

    std::vector<Matrix> expected_values = {ScalarMatrix(3, 3, -1), ScalarMatrix(3, 3, 0),
                                           ScalarMatrix(3, 3, 1),  ScalarMatrix(3, 3, -1),
                                           ScalarMatrix(3, 3, -1), ScalarMatrix(3, 3, 1)};
    std::vector<Matrix> actual_values;
    actual_values.reserve(r_model_part.Nodes().size());
    std::transform(r_model_part.Nodes().begin(), r_model_part.Nodes().end(),
                   std::back_inserter(actual_values), [&r_test_variable](const auto& node) {
        return node.FastGetSolutionStepValue(r_test_variable);
    });

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

    auto        model           = Model{};
    const auto& r_test_variable = CAUCHY_STRESS_VECTOR;
    auto&       r_model_part    = CreateModelPartWithTwoStubElements(model, r_test_variable);

    // Linear field in x between -1 and 1
    SetIntegrationPointValues(r_model_part.Elements()[1],
                              {ScalarVector(6, -inv_sqrt3), ScalarVector(6, inv_sqrt3),
                               ScalarVector(6, inv_sqrt3), ScalarVector(6, -inv_sqrt3)});

    // Linear field in y between -1 and 1
    SetIntegrationPointValues(r_model_part.Elements()[2],
                              {ScalarVector(6, -inv_sqrt3), ScalarVector(6, -inv_sqrt3),
                               ScalarVector(6, inv_sqrt3), ScalarVector(6, inv_sqrt3)});

    BuildAndRunExtrapolationProcess(model, CreateExtrapolationProcessSettings(r_model_part, r_test_variable));

    std::vector<Vector> expected_values = {ScalarVector(6, -1), ScalarVector(6, 0),
                                           ScalarVector(6, 1),  ScalarVector(6, -1),
                                           ScalarVector(6, -1), ScalarVector(6, 1)};
    std::vector<Vector> actual_values;
    actual_values.reserve(r_model_part.Nodes().size());
    std::transform(r_model_part.Nodes().begin(), r_model_part.Nodes().end(),
                   std::back_inserter(actual_values), [&r_test_variable](const auto& node) {
        return node.FastGetSolutionStepValue(r_test_variable);
    });

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

    auto        model           = Model{};
    const auto& r_test_variable = FLUID_FLUX_VECTOR;
    auto&       r_model_part    = CreateModelPartWithTwoStubElements(model, r_test_variable);

    // Linear field in x between -1 and 1
    dynamic_cast<StubElementForNodalExtrapolationTest&>(r_model_part.Elements()[1]).mIntegrationArrayValues = {
        ScalarVector(3, -inv_sqrt3), ScalarVector(3, inv_sqrt3), ScalarVector(3, inv_sqrt3),
        ScalarVector(3, -inv_sqrt3)};

    // Linear field in y between -1 and 1
    dynamic_cast<StubElementForNodalExtrapolationTest&>(r_model_part.Elements()[2]).mIntegrationArrayValues = {
        ScalarVector(3, -inv_sqrt3), ScalarVector(3, -inv_sqrt3), ScalarVector(3, inv_sqrt3),
        ScalarVector(3, inv_sqrt3)};

    BuildAndRunExtrapolationProcess(model, CreateExtrapolationProcessSettings(r_model_part, r_test_variable));

    std::vector<Vector>              expected_values = {ScalarVector(3, -1), ScalarVector(3, 0),
                                                        ScalarVector(3, 1),  ScalarVector(3, -1),
                                                        ScalarVector(3, -1), ScalarVector(3, 1)};
    std::vector<array_1d<double, 3>> actual_values;
    actual_values.reserve(r_model_part.Nodes().size());
    std::transform(r_model_part.Nodes().begin(), r_model_part.Nodes().end(),
                   std::back_inserter(actual_values), [&r_test_variable](const auto& node) {
        return node.FastGetSolutionStepValue(r_test_variable);
    });

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

    auto model = Model{};

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

    SetIntegrationPointValues(r_left_model_part.Elements()[1], std::vector(4, 1.0));
    SetIntegrationPointValues(r_right_model_part.Elements()[2], std::vector(4, 1.0));

    BuildAndRunExtrapolationProcess(
        model, CreateExtrapolationProcessSettings(
                   {std::cref(r_left_model_part), std::cref(r_right_model_part)}, HYDRAULIC_HEAD));

    AssertNodalValues(r_left_model_part.Nodes(), HYDRAULIC_HEAD, std::vector(4, 1.0));
    AssertNodalValues(r_right_model_part.Nodes(), HYDRAULIC_HEAD, std::vector(4, 1.0));
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

    auto model = Model{};

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

    SetIntegrationPointValues(r_bottom_model_part.Elements()[1], std::vector(4, 1.0));
    SetIntegrationPointValues(r_bottom_model_part.Elements()[2], std::vector(4, 1.0));
    SetIntegrationPointValues(r_top_model_part.Elements()[3], std::vector(4, 1.0));

    GeoExtrapolateIntegrationPointValuesToNodesProcess process_1(
        model, CreateExtrapolationProcessSettings(r_bottom_model_part, HYDRAULIC_HEAD));
    process_1.ExecuteBeforeSolutionLoop();

    GeoExtrapolateIntegrationPointValuesToNodesProcess process_2(
        model, CreateExtrapolationProcessSettings(r_top_model_part, HYDRAULIC_DISCHARGE));
    process_2.ExecuteBeforeSolutionLoop();

    process_1.ExecuteFinalizeSolutionStep();

    AssertNodalValues(r_bottom_model_part.Nodes(), HYDRAULIC_HEAD, std::vector(6, 1.0));

    process_2.ExecuteFinalizeSolutionStep();

    AssertNodalValues(r_top_model_part.Nodes(), HYDRAULIC_DISCHARGE, std::vector(4, 1.0));
}

KRATOS_TEST_CASE_IN_SUITE(CheckInfoGeoExtrapolateIntegrationPointValuesToNodesProcess, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto        model        = Model{};
    const auto& r_model_part = model.CreateModelPart("foo");

    const auto process = GeoExtrapolateIntegrationPointValuesToNodesProcess{
        model, CreateExtrapolationProcessSettings(r_model_part, FLUID_FLUX_VECTOR)};

    // Act & assert
    KRATOS_EXPECT_EQ(process.Info(), "GeoExtrapolateIntegrationPointValuesToNodesProcess");
}

} // namespace Kratos::Testing
