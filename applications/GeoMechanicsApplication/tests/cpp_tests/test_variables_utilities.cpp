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
#include "custom_utilities/variables_utilities.hpp"
#include "geo_mechanics_fast_suite.h"
#include "includes/variables.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TestVariablesUtilitiesGetsCorrectComponents, KratosGeoMechanicsFastSuite)
{
    const auto& component = VariablesUtilities::GetComponentFromVectorVariable(ACCELERATION.Name(), "X");

    KRATOS_EXPECT_EQ(component, ACCELERATION_X);
}

KRATOS_TEST_CASE_IN_SUITE(TestVariablesUtilitiesThrowsWhenComponentDoesNotExist, KratosGeoMechanicsFastSuite)
{
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        VariablesUtilities::GetComponentFromVectorVariable(ACCELERATION.Name(), "?"),
        "Error: The component \"ACCELERATION_?\" is not registered!")
}

} // namespace Kratos::Testing