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

    Testing::ModelSetupUtilities::CreateNodes(r_result, {{1, {0.0, 0.0, 0.0}},
                                                         {2, {1.0, 0.0, 0.0}},
                                                         {3, {1.0, 1.0, 0.0}},
                                                         {4, {0.0, 1.0, 0.0}},
                                                         {5, {2.0, 0.0, 0.0}},
                                                         {6, {2.0, 1.0, 0.0}}});

    const auto nodes_of_element_1 = Testing::ModelSetupUtilities::GetNodesFromIds(r_result, {1, 2, 3, 4});
    r_result.AddElement(make_intrusive<StubElementForNodalExtrapolationTest>(
        1, std::make_shared<Quadrilateral2D4<Node>>(nodes_of_element_1)));

    const auto nodes_of_element_2 = Testing::ModelSetupUtilities::GetNodesFromIds(r_result, {2, 5, 6, 3});
    r_result.AddElement(make_intrusive<StubElementForNodalExtrapolationTest>(
        2, std::make_shared<Quadrilateral2D4<Node>>(nodes_of_element_2)));

    return r_result;
}

std::string MakeModelPartNameFrom(const ModelPart& rModelPart)
{
    return R"(        "model_part_name": ")" + rModelPart.Name() + '"';
}

using ConstModelPartReferenceVector = std::vector<std::reference_wrapper<const ModelPart>>;

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

void SetIntegrationPointValues(Element& rElement, const std::vector<array_1d<double, 3>>& rValues)
{
    dynamic_cast<StubElementForNodalExtrapolationTest&>(rElement).mIntegrationArrayValues = rValues;
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

template <typename NodeContainerType, typename DataType>
std::vector<DataType> GetNodalValues(const NodeContainerType& rNodes, const Variable<DataType>& rVariable)
{
    auto result = std::vector<DataType>{};
    result.reserve(rNodes.size());
    auto get_nodal_value = [&rVariable](const auto& rNode) {
        return rNode.FastGetSolutionStepValue(rVariable);
    };
    std::ranges::transform(rNodes, std::back_inserter(result), get_nodal_value);

    return result;
}

template <typename NodeContainerType>
void AssertNodalValues(const NodeContainerType&   rNodes,
                       const Variable<double>&    rVariable,
                       const std::vector<double>& rExpectedNodalValues)
{
    KRATOS_EXPECT_VECTOR_NEAR(GetNodalValues(rNodes, rVariable), rExpectedNodalValues,
                              Testing::Defaults::absolute_tolerance);
}

template <typename NodeContainerType>
void AssertNodalValues(const NodeContainerType&   rNodes,
                       const Variable<Vector>&    rVariable,
                       const std::vector<Vector>& rExpectedNodalValues)
{
    const auto actual_nodal_values = GetNodalValues(rNodes, rVariable);

    ASSERT_EQ(actual_nodal_values.size(), rExpectedNodalValues.size());
    for (auto i = std::size_t{0}; i < actual_nodal_values.size(); ++i) {
        KRATOS_EXPECT_VECTOR_NEAR(actual_nodal_values[i], rExpectedNodalValues[i],
                                  Testing::Defaults::absolute_tolerance);
    }
}

template <typename NodeContainerType>
void AssertNodalValues(const NodeContainerType&                rNodes,
                       const Variable<array_1d<double, 3>>&    rVariable,
                       const std::vector<array_1d<double, 3>>& rExpectedNodalValues)
{
    const auto actual_nodal_values = GetNodalValues(rNodes, rVariable);

    ASSERT_EQ(actual_nodal_values.size(), rExpectedNodalValues.size());
    for (auto i = std::size_t{0}; i < actual_nodal_values.size(); ++i) {
        KRATOS_EXPECT_VECTOR_NEAR(actual_nodal_values[i], rExpectedNodalValues[i],
                                  Testing::Defaults::absolute_tolerance);
    }
}

template <typename NodeContainerType>
void AssertNodalValues(const NodeContainerType&   rNodes,
                       const Variable<Matrix>&    rVariable,
                       const std::vector<Matrix>& rExpectedNodalValues)
{
    const auto actual_nodal_values = GetNodalValues(rNodes, rVariable);

    ASSERT_EQ(actual_nodal_values.size(), rExpectedNodalValues.size());
    for (auto i = std::size_t{0}; i < actual_nodal_values.size(); ++i) {
        KRATOS_EXPECT_MATRIX_NEAR(actual_nodal_values[i], rExpectedNodalValues[i],
                                  Testing::Defaults::absolute_tolerance);
    }
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
                              {Matrix(3, 3, -inv_sqrt3), Matrix(3, 3, inv_sqrt3),
                               Matrix(3, 3, inv_sqrt3), Matrix(3, 3, -inv_sqrt3)});

    // Linear field in y between -1 and 1
    SetIntegrationPointValues(r_model_part.Elements()[2],
                              {Matrix(3, 3, -inv_sqrt3), Matrix(3, 3, -inv_sqrt3),
                               Matrix(3, 3, inv_sqrt3), Matrix(3, 3, inv_sqrt3)});

    BuildAndRunExtrapolationProcess(model, CreateExtrapolationProcessSettings(r_model_part, r_test_variable));

    const auto expected_values =
        std::vector{Matrix(3, 3, -1.0), Matrix(3, 3, 0.0),  Matrix(3, 3, 1.0),
                    Matrix(3, 3, -1.0), Matrix(3, 3, -1.0), Matrix(3, 3, 1.0)};
    AssertNodalValues(r_model_part.Nodes(), r_test_variable, expected_values);
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
                              std::vector{Vector(6, -inv_sqrt3), Vector(6, inv_sqrt3),
                                          Vector(6, inv_sqrt3), Vector(6, -inv_sqrt3)});

    // Linear field in y between -1 and 1
    SetIntegrationPointValues(r_model_part.Elements()[2],
                              std::vector{Vector(6, -inv_sqrt3), Vector(6, -inv_sqrt3),
                                          Vector(6, inv_sqrt3), Vector(6, inv_sqrt3)});

    BuildAndRunExtrapolationProcess(model, CreateExtrapolationProcessSettings(r_model_part, r_test_variable));

    const auto expected_values = std::vector{Vector(6, -1.0), Vector(6, 0.0),  Vector(6, 1.0),
                                             Vector(6, -1.0), Vector(6, -1.0), Vector(6, 1.0)};
    AssertNodalValues(r_model_part.Nodes(), r_test_variable, expected_values);
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
    SetIntegrationPointValues(
        r_model_part.Elements()[1],
        std::vector{array_1d<double, 3>(3, -inv_sqrt3), array_1d<double, 3>(3, inv_sqrt3),
                    array_1d<double, 3>(3, inv_sqrt3), array_1d<double, 3>(3, -inv_sqrt3)});

    // Linear field in y between -1 and 1
    SetIntegrationPointValues(
        r_model_part.Elements()[2],
        std::vector{array_1d<double, 3>(3, -inv_sqrt3), array_1d<double, 3>(3, -inv_sqrt3),
                    array_1d<double, 3>(3, inv_sqrt3), array_1d<double, 3>(3, inv_sqrt3)});

    BuildAndRunExtrapolationProcess(model, CreateExtrapolationProcessSettings(r_model_part, r_test_variable));

    const auto expected_values = std::vector{
        array_1d<double, 3>(3, -1.0), array_1d<double, 3>(3, 0.0),  array_1d<double, 3>(3, 1.0),
        array_1d<double, 3>(3, -1.0), array_1d<double, 3>(3, -1.0), array_1d<double, 3>(3, 1.0)};
    AssertNodalValues(r_model_part.Nodes(), r_test_variable, expected_values);
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
    const auto nodes_of_element_2 = ModelSetupUtilities::GetNodesFromIds(r_right_model_part, {2, 5, 6, 3});
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

    auto        model             = Model{};
    const auto& r_test_variable_1 = HYDRAULIC_HEAD;
    const auto& r_test_variable_2 = HYDRAULIC_DISCHARGE;

    auto& r_bottom_model_part = model.CreateModelPart("Bottom");
    r_bottom_model_part.AddNodalSolutionStepVariable(r_test_variable_1);
    r_bottom_model_part.AddNodalSolutionStepVariable(r_test_variable_2);
    ModelSetupUtilities::CreateNodes(r_bottom_model_part, {{1, {0.0, 0.0, 0.0}},
                                                           {2, {1.0, 0.0, 0.0}},
                                                           {3, {1.0, 1.0, 0.0}},
                                                           {4, {0.0, 1.0, 0.0}},
                                                           {5, {2.0, 0.0, 0.0}},
                                                           {6, {2.0, 1.0, 0.0}}});
    const auto nodes_of_element_1 = ModelSetupUtilities::GetNodesFromIds(r_bottom_model_part, {1, 2, 3, 4});
    r_bottom_model_part.AddElement(make_intrusive<StubElementForNodalExtrapolationTest>(
        1, std::make_shared<Quadrilateral2D4<Node>>(nodes_of_element_1)));
    const auto nodes_of_element_2 = ModelSetupUtilities::GetNodesFromIds(r_bottom_model_part, {2, 5, 6, 3});
    r_bottom_model_part.AddElement(make_intrusive<StubElementForNodalExtrapolationTest>(
        2, std::make_shared<Quadrilateral2D4<Node>>(nodes_of_element_2)));

    auto& r_top_model_part = model.CreateModelPart("Top");
    r_top_model_part.AddNodalSolutionStepVariable(r_test_variable_1);
    r_top_model_part.AddNodalSolutionStepVariable(r_test_variable_2);
    r_top_model_part.AddNode(r_bottom_model_part.pGetNode(3));
    r_top_model_part.AddNode(r_bottom_model_part.pGetNode(6));
    ModelSetupUtilities::CreateNodes(r_top_model_part, {{7, {1.0, 2.0, 0.0}}, {8, {2.0, 2.0, 0.0}}});
    const auto nodes_of_element_3 = ModelSetupUtilities::GetNodesFromIds(r_top_model_part, {3, 6, 8, 7});
    r_top_model_part.AddElement(make_intrusive<StubElementForNodalExtrapolationTest>(
        3, std::make_shared<Quadrilateral2D4<Node>>(nodes_of_element_3)));

    SetIntegrationPointValues(r_bottom_model_part.Elements()[1], std::vector(4, 1.0));
    SetIntegrationPointValues(r_bottom_model_part.Elements()[2], std::vector(4, 1.0));
    SetIntegrationPointValues(r_top_model_part.Elements()[3], std::vector(4, 1.0));

    auto process_1 = GeoExtrapolateIntegrationPointValuesToNodesProcess{
        model, CreateExtrapolationProcessSettings(r_bottom_model_part, r_test_variable_1)};
    process_1.ExecuteBeforeSolutionLoop();
    auto process_2 = GeoExtrapolateIntegrationPointValuesToNodesProcess{
        model, CreateExtrapolationProcessSettings(r_top_model_part, r_test_variable_2)};
    process_2.ExecuteBeforeSolutionLoop();
    process_1.ExecuteFinalizeSolutionStep();

    AssertNodalValues(r_bottom_model_part.Nodes(), r_test_variable_1, std::vector(6, 1.0));

    process_2.ExecuteFinalizeSolutionStep();

    AssertNodalValues(r_top_model_part.Nodes(), r_test_variable_2, std::vector(4, 1.0));
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
