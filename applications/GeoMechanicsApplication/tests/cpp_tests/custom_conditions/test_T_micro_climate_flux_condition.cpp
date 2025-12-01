// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Anne van de Graaf

#include "containers/model.h"
#include "custom_conditions/T_microclimate_flux_condition.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_2d_4.h"
#include "geometries/line_2d_5.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "includes/expect.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"

#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;
using namespace boost::numeric::ublas;

namespace
{

std::shared_ptr<Properties> CreateDummyConditionProperties(ModelPart& rModelPart)
{
    constexpr auto properties_id = ModelPart::IndexType{1};
    auto           p_result      = rModelPart.CreateNewProperties(properties_id);
    p_result->SetValue(ALPHA_COEFFICIENT, 0.0);
    p_result->SetValue(A1_COEFFICIENT, 0.0);
    p_result->SetValue(A2_COEFFICIENT, 0.0);
    p_result->SetValue(A3_COEFFICIENT, 0.0);
    p_result->SetValue(QF_COEFFICIENT, 0.0);
    p_result->SetValue(SMIN_COEFFICIENT, 0.0);
    p_result->SetValue(SMAX_COEFFICIENT, 0.0);
    p_result->SetValue(DENSITY_WATER, 1.0e3);
    return p_result;
}

void CreateNodesForLineCondition(ModelPart& rModelPart, std::size_t NumberOfNodes)
{
    for (auto node_index = std::size_t{0}; node_index < NumberOfNodes; ++node_index) {
        const auto x = static_cast<double>(node_index);
        const auto y = 0.5 * x;
        const auto z = 0.0;
        rModelPart.CreateNewNode(node_index + 1, x, y, z);
    }
}

void CreateNodesForTriangle(ModelPart& rModelPart, bool WantMidSideNodes)
{
    rModelPart.CreateNewNode(ModelPart::IndexType{1}, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(ModelPart::IndexType{2}, 2.0, 0.0, 0.0);
    rModelPart.CreateNewNode(ModelPart::IndexType{3}, 0.0, 2.0, 0.0);
    if (WantMidSideNodes) {
        rModelPart.CreateNewNode(ModelPart::IndexType{4}, 1.0, 0.0, 0.0);
        rModelPart.CreateNewNode(ModelPart::IndexType{5}, 1.0, 1.0, 0.0);
        rModelPart.CreateNewNode(ModelPart::IndexType{6}, 0.0, 1.0, 0.0);
    }
}

enum class HigherOrderElementConfiguration { None, Serendipity, LaGrange };

void CreateNodesForQuadrilateral(ModelPart& rModelPart, HigherOrderElementConfiguration configuration)
{
    rModelPart.CreateNewNode(ModelPart::IndexType{1}, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(ModelPart::IndexType{2}, 2.0, 0.0, 0.0);
    rModelPart.CreateNewNode(ModelPart::IndexType{3}, 2.0, 2.0, 0.0);
    rModelPart.CreateNewNode(ModelPart::IndexType{4}, 0.0, 2.0, 0.0);
    if (configuration == HigherOrderElementConfiguration::None) return;

    rModelPart.CreateNewNode(ModelPart::IndexType{5}, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(ModelPart::IndexType{6}, 2.0, 1.0, 0.0);
    rModelPart.CreateNewNode(ModelPart::IndexType{7}, 1.0, 2.0, 0.0);
    rModelPart.CreateNewNode(ModelPart::IndexType{8}, 0.0, 1.0, 0.0);
    if (configuration == HigherOrderElementConfiguration::LaGrange) {
        rModelPart.CreateNewNode(ModelPart::IndexType{9}, 1.0, 1.0, 0.0);
    }
}

void AddSolutionStepVariablesToModelPart(ModelPart&                                  rModelPart,
                                         const std::vector<const Variable<double>*>& rVariables)
{
    for (auto p_variable : rVariables) {
        rModelPart.AddNodalSolutionStepVariable(*p_variable);
    }
}

void AddSolutionStepValuesToNodes(ModelPart::NodesContainerType& rNodes, const Variable<double>& rVariable)
{
    for (auto& node : rNodes) {
        node.FastGetSolutionStepValue(rVariable, 0) = 1.0;
        node.FastGetSolutionStepValue(rVariable, 1) = 0.0;
    }
}

void AddSolutionStepValuesToNodes(ModelPart::NodesContainerType&              rNodes,
                                  const std::vector<const Variable<double>*>& rVariables)
{
    for (auto p_variable : rVariables) {
        AddSolutionStepValuesToNodes(rNodes, *p_variable);
    }
}

ModelPart& CreateDummyModelPartWithNodes(Model& rModel, const std::function<void(ModelPart&)>& rCreateNodesFunc)
{
    constexpr auto buffer_size = Model::IndexType{2};
    auto&          r_result    = rModel.CreateModelPart("dummy", buffer_size);

    const auto variables = std::vector<const Variable<double>*>{
        &AIR_TEMPERATURE, &AIR_HUMIDITY, &SOLAR_RADIATION, &WIND_SPEED,
        &PRECIPITATION,   &TEMPERATURE,  &DT_TEMPERATURE};
    AddSolutionStepVariablesToModelPart(r_result, variables);

    rCreateNodesFunc(r_result);
    AddSolutionStepValuesToNodes(r_result.Nodes(), variables);

    r_result.GetProcessInfo()[DELTA_TIME] = 0.5;

    return r_result;
}

intrusive_ptr<Condition> CreateMicroClimateCondition(ModelPart&             rModelPart,
                                                     shared_ptr<Properties> pProperties,
                                                     std::size_t            DimensionSize)
{
    auto r_nodes  = rModelPart.Nodes();
    auto node_ids = std::vector<ModelPart::IndexType>{};
    node_ids.reserve(r_nodes.size());
    std::transform(r_nodes.begin(), r_nodes.end(), std::back_inserter(node_ids),
                   [](const auto& node) { return node.Id(); });

    PointerVector<Node> node_pointers(node_ids.size());
    std::transform(node_ids.begin(), node_ids.end(), node_pointers.ptr_begin(),
                   [&rModelPart](auto node_id) { return rModelPart.pGetNode(node_id); });

    Condition::Pointer condition = nullptr;
    if (DimensionSize == 2 && node_ids.size() == 2) {
        condition = make_intrusive<GeoTMicroClimateFluxCondition<2, 2>>(
            1, Kratos::make_shared<Line2D2<Node>>(node_pointers), pProperties);
    } else if (DimensionSize == 2 && node_ids.size() == 3) {
        condition = make_intrusive<GeoTMicroClimateFluxCondition<2, 3>>(
            1, Kratos::make_shared<Line2D3<Node>>(node_pointers), pProperties);
    } else if (DimensionSize == 2 && node_ids.size() == 4) {
        condition = make_intrusive<GeoTMicroClimateFluxCondition<2, 4>>(
            1, Kratos::make_shared<Line2D4<Node>>(node_pointers), pProperties);
    } else if (DimensionSize == 2 && node_ids.size() == 5) {
        condition = make_intrusive<GeoTMicroClimateFluxCondition<2, 5>>(
            1, Kratos::make_shared<Line2D5<Node>>(node_pointers), pProperties);
    } else if (DimensionSize == 3 && node_ids.size() == 3) {
        condition = make_intrusive<GeoTMicroClimateFluxCondition<3, 3>>(
            1, Kratos::make_shared<Triangle3D3<Node>>(node_pointers), pProperties);
    } else if (DimensionSize == 3 && node_ids.size() == 6) {
        condition = make_intrusive<GeoTMicroClimateFluxCondition<3, 6>>(
            1, Kratos::make_shared<Triangle3D6<Node>>(node_pointers), pProperties);
    } else if (DimensionSize == 3 && node_ids.size() == 4) {
        condition = make_intrusive<GeoTMicroClimateFluxCondition<3, 4>>(
            1, Kratos::make_shared<Quadrilateral3D4<Node>>(node_pointers), pProperties);
    } else if (DimensionSize == 3 && node_ids.size() == 8) {
        condition = make_intrusive<GeoTMicroClimateFluxCondition<3, 8>>(
            1, Kratos::make_shared<Quadrilateral3D8<Node>>(node_pointers), pProperties);
    } else if (DimensionSize == 3 && node_ids.size() == 9) {
        condition = make_intrusive<GeoTMicroClimateFluxCondition<3, 9>>(
            1, Kratos::make_shared<Quadrilateral3D9<Node>>(node_pointers), pProperties);
    }

    KRATOS_ERROR_IF_NOT(condition) << "Condition could not be created" << std::endl;

    rModelPart.AddCondition(condition);
    return condition;
}

std::string ExecuteInitializeSolutionStep(intrusive_ptr<Condition> pCondition, const ProcessInfo& rProcessInfo)
{
    try {
        pCondition->InitializeSolutionStep(rProcessInfo);
    } catch (const Exception& e) {
        return e.what();
    }

    return {}; // no error message
}

constexpr auto relative_tolerance = 1.0e-3;
constexpr auto absolute_tolerance = 1.0e-3;

} // namespace

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, NoThrowWhenInitializingThermalMicroClimateCondition)
{
    Model test_model;
    auto  create_nodes_func = [](ModelPart& rModelPart) {
        constexpr auto number_of_nodes = std::size_t{2};
        CreateNodesForLineCondition(rModelPart, number_of_nodes);
    };
    auto&          r_model_part   = CreateDummyModelPartWithNodes(test_model, create_nodes_func);
    auto           p_properties   = CreateDummyConditionProperties(r_model_part);
    constexpr auto dimension_size = std::size_t{2};
    auto p_condition = CreateMicroClimateCondition(r_model_part, p_properties, dimension_size);

    auto has_thrown = false;
    try {
        p_condition->Initialize(r_model_part.GetProcessInfo());
    } catch (const Exception&) {
        has_thrown = true;
    }

    KRATOS_EXPECT_FALSE(has_thrown)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, NoErrorWhenInitializingSolutionStepOnThermalMicroClimateCondition2D2N)
{
    Model test_model;
    auto  create_nodes_func = [](ModelPart& rModelPart) {
        constexpr auto number_of_nodes = std::size_t{2};
        CreateNodesForLineCondition(rModelPart, number_of_nodes);
    };
    auto&          r_model_part   = CreateDummyModelPartWithNodes(test_model, create_nodes_func);
    auto           p_properties   = CreateDummyConditionProperties(r_model_part);
    constexpr auto dimension_size = std::size_t{2};
    auto p_condition = CreateMicroClimateCondition(r_model_part, p_properties, dimension_size);

    const auto error_text = ExecuteInitializeSolutionStep(p_condition, r_model_part.GetProcessInfo());

    KRATOS_EXPECT_STREQ(error_text.data(), "")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, NoErrorWhenInitializingSolutionStepOnThermalMicroClimateCondition2D3N)
{
    Model test_model;
    auto  create_nodes_func = [](ModelPart& rModelPart) {
        constexpr auto number_of_nodes = std::size_t{3};
        CreateNodesForLineCondition(rModelPart, number_of_nodes);
    };
    auto&          r_model_part   = CreateDummyModelPartWithNodes(test_model, create_nodes_func);
    auto           p_properties   = CreateDummyConditionProperties(r_model_part);
    constexpr auto dimension_size = std::size_t{2};
    auto p_condition = CreateMicroClimateCondition(r_model_part, p_properties, dimension_size);

    const auto error_text = ExecuteInitializeSolutionStep(p_condition, r_model_part.GetProcessInfo());

    KRATOS_EXPECT_STREQ(error_text.data(), "")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, NoErrorWhenInitializingSolutionStepOnThermalMicroClimateCondition2D4N)
{
    Model test_model;
    auto  create_nodes_func = [](ModelPart& rModelPart) {
        constexpr auto number_of_nodes = std::size_t{4};
        CreateNodesForLineCondition(rModelPart, number_of_nodes);
    };
    auto&          r_model_part   = CreateDummyModelPartWithNodes(test_model, create_nodes_func);
    auto           p_properties   = CreateDummyConditionProperties(r_model_part);
    constexpr auto dimension_size = std::size_t{2};
    auto p_condition = CreateMicroClimateCondition(r_model_part, p_properties, dimension_size);

    const auto error_text = ExecuteInitializeSolutionStep(p_condition, r_model_part.GetProcessInfo());

    KRATOS_EXPECT_STREQ(error_text.data(), "")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, NoErrorWhenInitializingSolutionStepOnThermalMicroClimateCondition2D5N)
{
    Model test_model;
    auto  create_nodes_func = [](ModelPart& rModelPart) {
        constexpr auto number_of_nodes = std::size_t{5};
        CreateNodesForLineCondition(rModelPart, number_of_nodes);
    };
    auto&          r_model_part   = CreateDummyModelPartWithNodes(test_model, create_nodes_func);
    auto           p_properties   = CreateDummyConditionProperties(r_model_part);
    constexpr auto dimension_size = std::size_t{2};
    auto p_condition = CreateMicroClimateCondition(r_model_part, p_properties, dimension_size);

    const auto error_text = ExecuteInitializeSolutionStep(p_condition, r_model_part.GetProcessInfo());

    KRATOS_EXPECT_STREQ(error_text.data(), "")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, NoErrorWhenInitializingSolutionStepOnThermalMicroClimateCondition3D3N)
{
    Model test_model;
    auto  create_nodes_func = [](ModelPart& rModelPart) {
        constexpr auto want_mid_side_nodes = false;
        CreateNodesForTriangle(rModelPart, want_mid_side_nodes);
    };
    auto&          r_model_part   = CreateDummyModelPartWithNodes(test_model, create_nodes_func);
    auto           p_properties   = CreateDummyConditionProperties(r_model_part);
    constexpr auto dimension_size = std::size_t{3};
    auto p_condition = CreateMicroClimateCondition(r_model_part, p_properties, dimension_size);

    const auto error_text = ExecuteInitializeSolutionStep(p_condition, r_model_part.GetProcessInfo());

    KRATOS_EXPECT_STREQ(error_text.data(), "")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, NoErrorWhenInitializingSolutionStepOnThermalMicroClimateCondition3D6N)
{
    Model test_model;
    auto  create_nodes_func = [](ModelPart& rModelPart) {
        constexpr auto want_mid_side_nodes = true;
        CreateNodesForTriangle(rModelPart, want_mid_side_nodes);
    };
    auto&          r_model_part   = CreateDummyModelPartWithNodes(test_model, create_nodes_func);
    auto           p_properties   = CreateDummyConditionProperties(r_model_part);
    constexpr auto dimension_size = std::size_t{3};
    auto p_condition = CreateMicroClimateCondition(r_model_part, p_properties, dimension_size);

    const auto error_text = ExecuteInitializeSolutionStep(p_condition, r_model_part.GetProcessInfo());

    KRATOS_EXPECT_STREQ(error_text.data(), "")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, NoErrorWhenInitializingSolutionStepOnThermalMicroClimateCondition3D4N)
{
    Model test_model;
    auto  create_nodes_func = [](ModelPart& rModelPart) {
        CreateNodesForQuadrilateral(rModelPart, HigherOrderElementConfiguration::None);
    };
    auto&          r_model_part   = CreateDummyModelPartWithNodes(test_model, create_nodes_func);
    auto           p_properties   = CreateDummyConditionProperties(r_model_part);
    constexpr auto dimension_size = std::size_t{3};
    auto p_condition = CreateMicroClimateCondition(r_model_part, p_properties, dimension_size);

    const auto error_text = ExecuteInitializeSolutionStep(p_condition, r_model_part.GetProcessInfo());

    KRATOS_EXPECT_STREQ(error_text.data(), "")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, NoErrorWhenInitializingSolutionStepOnThermalMicroClimateCondition3D8N)
{
    Model test_model;
    auto  create_nodes_func = [](ModelPart& rModelPart) {
        CreateNodesForQuadrilateral(rModelPart, HigherOrderElementConfiguration::Serendipity);
    };
    auto&          r_model_part   = CreateDummyModelPartWithNodes(test_model, create_nodes_func);
    auto           p_properties   = CreateDummyConditionProperties(r_model_part);
    constexpr auto dimension_size = std::size_t{3};
    auto p_condition = CreateMicroClimateCondition(r_model_part, p_properties, dimension_size);

    const auto error_text = ExecuteInitializeSolutionStep(p_condition, r_model_part.GetProcessInfo());

    KRATOS_EXPECT_STREQ(error_text.data(), "")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, NoErrorWhenInitializingSolutionStepOnThermalMicroClimateCondition3D9N)
{
    Model test_model;
    auto  create_nodes_func = [](ModelPart& rModelPart) {
        CreateNodesForQuadrilateral(rModelPart, HigherOrderElementConfiguration::LaGrange);
    };
    auto&          r_model_part   = CreateDummyModelPartWithNodes(test_model, create_nodes_func);
    auto           p_properties   = CreateDummyConditionProperties(r_model_part);
    constexpr auto dimension_size = std::size_t{3};
    auto p_condition = CreateMicroClimateCondition(r_model_part, p_properties, dimension_size);

    const auto error_text = ExecuteInitializeSolutionStep(p_condition, r_model_part.GetProcessInfo());

    KRATOS_EXPECT_STREQ(error_text.data(), "")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CalculateLocalSystemForThermalMicroClimateCondition2D3N)
{
    Model test_model;
    auto  create_nodes_func = [](ModelPart& rModelPart) {
        constexpr auto number_of_nodes = std::size_t{3};
        CreateNodesForLineCondition(rModelPart, number_of_nodes);
    };
    auto&          r_model_part   = CreateDummyModelPartWithNodes(test_model, create_nodes_func);
    auto           p_properties   = CreateDummyConditionProperties(r_model_part);
    constexpr auto dimension_size = std::size_t{2};
    auto p_condition = CreateMicroClimateCondition(r_model_part, p_properties, dimension_size);
    p_condition->Initialize(r_model_part.GetProcessInfo());
    p_condition->InitializeSolutionStep(r_model_part.GetProcessInfo());

    Matrix lhs_matrix;
    Vector rhs_vector;
    p_condition->CalculateLocalSystem(lhs_matrix, rhs_vector, r_model_part.GetProcessInfo());

    auto expected_lhs_matrix = Matrix{3, 3, 0.0};
    // clang-format off
    expected_lhs_matrix <<= 21.2568,  -8.50271, 25.5081,
                            -8.50271, 12.7541,   8.50271,
                            25.5081,   8.50271, 68.0217;
    // clang-format on
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(expected_lhs_matrix, lhs_matrix, relative_tolerance)

    auto expected_rhs_vector = Vector{3, 0.0};
    expected_rhs_vector <<= -76.6356, -25.5452, -204.362;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_rhs_vector, rhs_vector, relative_tolerance)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CalculateLocalSystemForThermalMicroClimateCondition3D6N)
{
    Model test_model;
    auto  create_nodes_func = [](ModelPart& rModelPart) {
        constexpr auto want_mid_side_nodes = true;
        CreateNodesForTriangle(rModelPart, want_mid_side_nodes);
    };
    auto&          r_model_part   = CreateDummyModelPartWithNodes(test_model, create_nodes_func);
    auto           p_properties   = CreateDummyConditionProperties(r_model_part);
    constexpr auto dimension_size = std::size_t{3};
    auto p_condition = CreateMicroClimateCondition(r_model_part, p_properties, dimension_size);
    p_condition->Initialize(r_model_part.GetProcessInfo());
    p_condition->InitializeSolutionStep(r_model_part.GetProcessInfo());

    Matrix lhs_matrix;
    Vector rhs_vector;
    p_condition->CalculateLocalSystem(lhs_matrix, rhs_vector, r_model_part.GetProcessInfo());

    auto expected_lhs_matrix = Matrix{6, 6, 0.0};
    // clang-format off
    expected_lhs_matrix <<= 1.95146,  -0.975729, -0.975729,  0.975729, -1.95146,   0.975729,
                           -0.975729,  1.95146,  -0.975729,  0.975729,  0.975729, -1.95146,
                           -0.975729, -0.975729,  1.95146,  -1.95146,   0.975729,  0.975729,
                            0.975729,  0.975729, -1.95146,  10.733,     7.80583,   7.80583,
                           -1.95146,   0.975729,  0.975729,  7.80583,  10.733,     7.80583,
                            0.975729, -1.95146,   0.975729,  7.80583,   7.80583,  10.733;
    // clang-format on
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(expected_lhs_matrix, lhs_matrix, relative_tolerance)

    auto expected_rhs_vector = Vector{6, 0.0};
    expected_rhs_vector <<= 0.0, 0.0, 0.0, -52.7659, -52.7659, -52.7659;
    // To compare computed zeros (the first three elements of 'rhs_vector') use an absolute_tolerance
    KRATOS_EXPECT_VECTOR_NEAR(expected_rhs_vector, rhs_vector, absolute_tolerance)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CalculateLocalSystemForThermalMicroClimateCondition3D8N)
{
    Model test_model;
    auto  create_nodes_func = [](ModelPart& rModelPart) {
        CreateNodesForQuadrilateral(rModelPart, HigherOrderElementConfiguration::Serendipity);
    };
    auto&          r_model_part   = CreateDummyModelPartWithNodes(test_model, create_nodes_func);
    auto           p_properties   = CreateDummyConditionProperties(r_model_part);
    constexpr auto dimension_size = std::size_t{3};
    auto p_condition = CreateMicroClimateCondition(r_model_part, p_properties, dimension_size);
    p_condition->Initialize(r_model_part.GetProcessInfo());
    p_condition->InitializeSolutionStep(r_model_part.GetProcessInfo());

    Matrix lhs_matrix;
    Vector rhs_vector;
    p_condition->CalculateLocalSystem(lhs_matrix, rhs_vector, r_model_part.GetProcessInfo());

    auto expected_lhs_matrix = Matrix{8, 8, 0.0};
    // clang-format off
    expected_lhs_matrix <<= 5.26894,  1.75631,  2.63447,  1.75631, -5.26894, -7.02525, -7.02525, -5.26894,
                            1.75631,  5.26894,  1.75631,  2.63447, -5.26894, -5.26894, -7.02525, -7.02525,
                            2.63447,  1.75631,  5.26894,  1.75631, -7.02525, -5.26894, -5.26894, -7.02525,
                            1.75631,  2.63447,  1.75631,  5.26894, -7.02525, -7.02525, -5.26894, -5.26894,
                           -5.26894, -5.26894, -7.02525, -7.02525, 28.101,   17.5631,  14.0505,  17.5631,
                           -7.02525, -5.26894, -5.26894, -7.02525, 17.5631,  28.101,   17.5631,  14.0505,
                           -7.02525, -7.02525, -5.26894, -5.26894, 14.0505,  17.5631,  28.101,   17.5631,
                           -5.26894, -7.02525, -7.02525, -5.26894, 17.5631,  14.0505,  17.5631,  28.101;
    // clang-format on

    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(expected_lhs_matrix, lhs_matrix, relative_tolerance)

    auto expected_rhs_vector = Vector{8, 0.0};
    expected_rhs_vector <<= 26.383, 26.383, 26.383, 26.383, -105.532, -105.532, -105.532, -105.532;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_rhs_vector, rhs_vector, relative_tolerance)
}

} // namespace Kratos::Testing
