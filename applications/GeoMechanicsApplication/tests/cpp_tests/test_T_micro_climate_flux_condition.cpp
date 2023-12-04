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
#include "custom_conditions/T_microclimate_flux_condition.hpp"
#include "testing/testing.h"

using namespace Kratos;

namespace {

std::shared_ptr<Properties> CreateDummyConditionProperties(ModelPart& rModelPart)
{
    constexpr auto properties_id = ModelPart::IndexType{1};
    auto p_result = rModelPart.CreateNewProperties(properties_id);
    p_result->SetValue(ALPHA_COEFFICIENT, 0.0);
    p_result->SetValue(A1_COEFFICIENT, 0.0);
    p_result->SetValue(A2_COEFFICIENT, 0.0);
    p_result->SetValue(A3_COEFFICIENT, 0.0);
    p_result->SetValue(QF_COEFFICIENT, 0.0);
    p_result->SetValue(SMIN_COEFFICIENT, 0.0);
    p_result->SetValue(SMAX_COEFFICIENT, 0.0);
    return p_result;
}

void CreateNodesInModelPart(ModelPart& rModelPart,
                            std::size_t NumberOfNodes)
{
    for (auto node_index = 0; node_index < NumberOfNodes; ++node_index) {
        const auto x = static_cast<double>(node_index);
        const auto y = 0.5 * x;
        const auto z = 0.0;
        rModelPart.CreateNewNode(node_index+1, x, y, z);
    }
}

void AddSolutionStepVariablesToModelPart(ModelPart& rModelPart,
                                         const std::vector<const Variable<double>*>& rVariables)
{
    for (auto p_variable : rVariables) {
        rModelPart.AddNodalSolutionStepVariable(*p_variable);
    }
}

void AddSolutionStepValuesToNodes(ModelPart::NodesContainerType& rNodes,
                                  const Variable<double>& rVariable)
{
    for (auto& node : rNodes) {
        node.FastGetSolutionStepValue(rVariable, 0) = 1.0;
        node.FastGetSolutionStepValue(rVariable, 1) = 0.0;
    }
}

void AddSolutionStepValuesToNodes(ModelPart::NodesContainerType& rNodes,
                                  const std::vector<const Variable<double>*>& rVariables)
{
    for (auto p_variable : rVariables) {
        AddSolutionStepValuesToNodes(rNodes, *p_variable);
    }
}

ModelPart& CreateDummyModelPartWithNodes(Model& rModel,
                                         std::size_t NumberOfNodes)
{
    constexpr auto buffer_size = Model::IndexType{2};
    auto& r_result = rModel.CreateModelPart("dummy", buffer_size);

    const auto variables = std::vector<const Variable<double>*>{&AIR_TEMPERATURE,
                                                                &SOLAR_RADIATION,
                                                                &WIND_SPEED,
                                                                &TEMPERATURE};
    AddSolutionStepVariablesToModelPart(r_result, variables);

    CreateNodesInModelPart(r_result, NumberOfNodes);
    AddSolutionStepValuesToNodes(r_result.Nodes(), variables);

    r_result.GetProcessInfo()[DELTA_TIME] = 0.5;

    return r_result;
}

intrusive_ptr<Condition> CreateMicroClimateCondition(ModelPart& rModelPart,
                                                     shared_ptr<Properties> pProperties)
{
    constexpr auto condition_id = std::size_t{1};
    auto node_ids = std::vector<ModelPart::IndexType>{};
    for (const auto& node : rModelPart.Nodes()) {
        node_ids.emplace_back(node.Id());
    }
    const auto condition_name = std::string{"GeoTMicroClimateFluxCondition2D"} + std::to_string(node_ids.size()) + "N";
    return rModelPart.CreateNewCondition(condition_name, condition_id, node_ids, pProperties);
}

std::string ExecuteInitializeSolutionStep(intrusive_ptr<Condition> pCondition,
                                          const ProcessInfo& rProcessInfo)
{
    try {
        pCondition->InitializeSolutionStep(rProcessInfo);
    }
    catch (const Exception& e) {
        return e.what();
    }

    return {}; // no error message
}

}


namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(NoThrowWhenInitializingThermalMicroClimateCondition, KratosGeoMechanicsFastSuite)
{
    Model test_model;
    constexpr auto buffer_size = Model::IndexType{2};
    auto& r_model_part = test_model.CreateModelPart("dummy", buffer_size);

    constexpr auto number_of_dimensions = 2u;
    constexpr auto number_of_nodes = 2u;
    auto condition = TMicroClimateFluxCondition<number_of_dimensions, number_of_nodes>{};

    auto has_thrown = false;
    try {
        condition.Initialize(r_model_part.GetProcessInfo());
    }
    catch (const Exception&) {
        has_thrown = true;
    }

    KRATOS_EXPECT_FALSE(has_thrown)
}

KRATOS_TEST_CASE_IN_SUITE(NoErrorWhenInitializingSolutionStepOnThermalMicroClimateCondition2D2N, KratosGeoMechanicsFastSuite)
{
    Model test_model;
    constexpr auto number_of_nodes = std::size_t{2};
    auto& r_model_part = CreateDummyModelPartWithNodes(test_model, number_of_nodes);
    auto  p_properties = CreateDummyConditionProperties(r_model_part);
    auto  p_condition  = CreateMicroClimateCondition(r_model_part, p_properties);

    const auto error_text = ExecuteInitializeSolutionStep(p_condition, r_model_part.GetProcessInfo());

    KRATOS_EXPECT_STREQ(error_text.data(), "")
}

KRATOS_TEST_CASE_IN_SUITE(NoErrorWhenInitializingSolutionStepOnThermalMicroClimateCondition2D3N, KratosGeoMechanicsFastSuite)
{
    Model test_model;
    constexpr auto number_of_nodes = std::size_t{3};
    auto& r_model_part = CreateDummyModelPartWithNodes(test_model, number_of_nodes);
    auto  p_properties = CreateDummyConditionProperties(r_model_part);
    auto  p_condition  = CreateMicroClimateCondition(r_model_part, p_properties);

    const auto error_text = ExecuteInitializeSolutionStep(p_condition, r_model_part.GetProcessInfo());

    KRATOS_EXPECT_STREQ(error_text.data(), "")
}

KRATOS_TEST_CASE_IN_SUITE(NoErrorWhenInitializingSolutionStepOnThermalMicroClimateCondition2D4N, KratosGeoMechanicsFastSuite)
{
    Model test_model;
    constexpr auto number_of_nodes = std::size_t{4};
    auto& r_model_part = CreateDummyModelPartWithNodes(test_model, number_of_nodes);
    auto  p_properties = CreateDummyConditionProperties(r_model_part);
    auto  p_condition  = CreateMicroClimateCondition(r_model_part, p_properties);

    const auto error_text = ExecuteInitializeSolutionStep(p_condition, r_model_part.GetProcessInfo());

    KRATOS_EXPECT_STREQ(error_text.data(), "")
}

}
