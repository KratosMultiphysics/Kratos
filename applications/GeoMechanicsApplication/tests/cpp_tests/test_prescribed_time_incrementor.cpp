// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//                   Anne van de Graaf
//

#include "testing/testing.h"
#include "custom_workflows/prescribed_time_incrementor.h"
#include "custom_workflows/time_step_end_state.hpp"

using namespace Kratos;


namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(NoNextTimeStepWhenPrescribedTimeIncrementorHasEmptyList, KratosGeoMechanicsFastSuite)
{
    std::vector<double> increments;
    PrescribedTimeIncrementor incrementor{increments};
    TimeStepEndState previous_state;

    KRATOS_EXPECT_FALSE(incrementor.WantNextStep(previous_state));
}

KRATOS_TEST_CASE_IN_SUITE(WantFirstTimeStepWhenPrescribedTimeIncrementorHasNonEmptyList, KratosGeoMechanicsFastSuite)
{
    std::vector<double> increments{0.4, 0.6};
    PrescribedTimeIncrementor incrementor{increments};
    TimeStepEndState previous_state;

    KRATOS_EXPECT_TRUE(incrementor.WantNextStep(previous_state));
}

}
