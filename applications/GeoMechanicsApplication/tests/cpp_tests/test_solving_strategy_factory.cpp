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

#include "testing/testing.h"
#include "custom_utilities/solving_strategy_factory.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(Create_ReturnsNullptr_WhenNoCreatorWasAddedForRequestedStrategy, WorkInProgress)
{
    SolvingStrategyFactory factory;

    const auto created_strategy = factory.Create({R"({"strategy_type":"Unknown"})"});

    KRATOS_EXPECT_EQ(created_strategy, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(Create_Throws_WhenCallbackFunctionThrowsAndRequestIsInvalid, WorkInProgress)
{

}

KRATOS_TEST_CASE_IN_SUITE(Create_ReturnsSolvingStrategy_ForLinearStrategy, WorkInProgress)
{
    SolvingStrategyFactory factory;

    const auto created_strategy = factory.Create({R"({"strategy_type":"linear_strategy"})"});

    KRATOS_EXPECT_NE(created_strategy, nullptr);
}

}