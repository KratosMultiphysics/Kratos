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

KRATOS_TEST_CASE_IN_SUITE(NoThrowWhenInitializingSolutionStepOnThermalMicroClimateCondition, KratosGeoMechanicsFastSuite)
{
    Model test_model;
    constexpr auto buffer_size = Model::IndexType{2};
    auto& r_model_part = test_model.CreateModelPart("dummy", buffer_size);
    r_model_part.AddNodalSolutionStepVariable(AIR_TEMPERATURE);
    r_model_part.AddNodalSolutionStepVariable(SOLAR_RADIATION);
    r_model_part.AddNodalSolutionStepVariable(WIND_SPEED);
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_geometry = r_model_part.CreateNewGeometry("Line2D2", std::vector<ModelPart::IndexType>{1, 2});

    for (auto& node : r_model_part.Nodes()) {
        node.FastGetSolutionStepValue(AIR_TEMPERATURE, 0) = 1.0;
        node.FastGetSolutionStepValue(AIR_TEMPERATURE, 1) = 0.0;
        node.FastGetSolutionStepValue(SOLAR_RADIATION, 0) = 1.0;
        node.FastGetSolutionStepValue(SOLAR_RADIATION, 1) = 0.0;
        node.FastGetSolutionStepValue(WIND_SPEED, 0) = 1.0;
        node.FastGetSolutionStepValue(WIND_SPEED, 1) = 0.0;
        node.FastGetSolutionStepValue(TEMPERATURE, 0) = 1.0;
        node.FastGetSolutionStepValue(TEMPERATURE, 1) = 0.0;
    }

    auto p_properties = CreateDummyConditionProperties(r_model_part);

    constexpr auto number_of_dimensions = 2u;
    constexpr auto number_of_nodes = 2u;
    constexpr auto condition_id = size_t{1};
    auto condition = TMicroClimateFluxCondition<number_of_dimensions, number_of_nodes>{condition_id, p_geometry, p_properties};

    r_model_part.GetProcessInfo()[DELTA_TIME] = 0.5;

    auto error_text = std::string{};
    auto has_thrown = false;
    try {
        condition.InitializeSolutionStep(r_model_part.GetProcessInfo());
    }
    catch (const Exception& e) {
        error_text = e.what();
        has_thrown = true;
    }

    KRATOS_EXPECT_STREQ(error_text.data(), std::string{}.data())
    KRATOS_EXPECT_FALSE(has_thrown)
}

}
