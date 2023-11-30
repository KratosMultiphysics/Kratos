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

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(NoThrowWhenInitializingThermalMicroClimateCondition, KratosGeoMechanicsFastSuite)
{
    Model test_model;
    constexpr auto buffer_size = Model::IndexType{2};
    auto& r_model_part = test_model.CreateModelPart("dummy", buffer_size);

    constexpr auto number_of_dimensions = 2u;
    constexpr auto number_of_nodes = 3u;
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

}
