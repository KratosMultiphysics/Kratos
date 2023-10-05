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

namespace
{

TimeStepEndState MakeConvergedEndState()
{
    auto result = TimeStepEndState{};
    result.convergence_state = TimeStepEndState::ConvergenceState::converged;
    return result;
}

}


namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(NoNextTimeStepWhenPrescribedTimeIncrementorHasEmptyList, KratosGeoMechanicsFastSuite)
{
    std::vector<double> increments;
    PrescribedTimeIncrementor incrementor{increments};
    const auto previous_state = MakeConvergedEndState();

    KRATOS_EXPECT_FALSE(incrementor.WantNextStep(previous_state))
}

KRATOS_TEST_CASE_IN_SUITE(WantFirstTimeStepWhenPrescribedTimeIncrementorHasNonEmptyList, KratosGeoMechanicsFastSuite)
{
    std::vector<double> increments{0.4, 0.6};
    PrescribedTimeIncrementor incrementor{increments};
    const auto previous_state = MakeConvergedEndState();

    KRATOS_EXPECT_TRUE(incrementor.WantNextStep(previous_state))
}

KRATOS_TEST_CASE_IN_SUITE(WantFirstTwoTimesNewStepWhenPrescribedTimeIncrementorHasTwoItems, KratosGeoMechanicsFastSuite)
{
    std::vector<double> increments{0.4, 0.6};
    PrescribedTimeIncrementor incrementor{increments};
    const auto previous_state = MakeConvergedEndState();

    KRATOS_EXPECT_TRUE(incrementor.WantNextStep(previous_state))
    incrementor.PostTimeStepExecution(previous_state);
    KRATOS_EXPECT_TRUE(incrementor.WantNextStep(previous_state))
    incrementor.PostTimeStepExecution(previous_state);
    KRATOS_EXPECT_FALSE(incrementor.WantNextStep(previous_state))
    KRATOS_EXPECT_FALSE(incrementor.WantNextStep(previous_state))
}

KRATOS_TEST_CASE_IN_SUITE(PrescribedTimeIncrementorThrowsIfAnyIncrementIsNegative, KratosGeoMechanicsFastSuite)
{
    std::vector<double> increments{0.4, -0.6};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(PrescribedTimeIncrementor{increments},
                                      "All prescribed increments must not be negative")
}

KRATOS_TEST_CASE_IN_SUITE(WantRetryStepAlwaysReturnsTrueOnFirstCycle, KratosGeoMechanicsFastSuite)
{
    std::vector<double> increments{0.4, 0.6};
    PrescribedTimeIncrementor incrementor{increments};
    const auto previous_state = MakeConvergedEndState();
    const auto cycle_number = std::size_t{0};

    KRATOS_EXPECT_TRUE(incrementor.WantRetryStep(cycle_number, previous_state))
}

KRATOS_TEST_CASE_IN_SUITE(WantRetryStepAlwaysReturnsFalseOnAnySubsequentCycle, KratosGeoMechanicsFastSuite)
{
    std::vector<double> increments{0.4, 0.6};
    PrescribedTimeIncrementor incrementor{increments};
    const auto previous_state = MakeConvergedEndState();
    const auto cycle_number = std::size_t{1};

    KRATOS_EXPECT_FALSE(incrementor.WantRetryStep(cycle_number, previous_state))
}

KRATOS_TEST_CASE_IN_SUITE(PrescribedTimeIncrementsMustMatchInput, KratosGeoMechanicsFastSuite)
{
    std::vector<double> increments{0.4, 0.6};
    PrescribedTimeIncrementor incrementor{increments};
    const auto previous_state = MakeConvergedEndState();

    KRATOS_EXPECT_DOUBLE_EQ(increments.front(), incrementor.GetIncrement());
    incrementor.PostTimeStepExecution(previous_state);
    KRATOS_EXPECT_DOUBLE_EQ(increments.back(), incrementor.GetIncrement());
}

KRATOS_TEST_CASE_IN_SUITE(PrescribedTimeIncrementorThrowsWhenAskingForIncrementBeyondEnd, KratosGeoMechanicsFastSuite)
{
    std::vector<double> increments{0.4, 0.6};
    PrescribedTimeIncrementor incrementor{increments};
    const auto previous_state = MakeConvergedEndState();
    incrementor.PostTimeStepExecution(previous_state);
    incrementor.PostTimeStepExecution(previous_state);
    // Now we are beyond the end of the increment vector

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(incrementor.GetIncrement(), "Out of increment range");
}

KRATOS_TEST_CASE_IN_SUITE(WithoutPostTimeStepExecutionAlwaysGetSameIncrement, KratosGeoMechanicsFastSuite)
{
    std::vector<double> increments{0.4, 0.6};
    PrescribedTimeIncrementor incrementor{increments};

    KRATOS_EXPECT_DOUBLE_EQ(increments.front(), incrementor.GetIncrement());
    KRATOS_EXPECT_DOUBLE_EQ(increments.front(), incrementor.GetIncrement());

    const auto previous_state = MakeConvergedEndState();
    incrementor.PostTimeStepExecution(previous_state);

    KRATOS_EXPECT_DOUBLE_EQ(increments.back(), incrementor.GetIncrement());
    KRATOS_EXPECT_DOUBLE_EQ(increments.back(), incrementor.GetIncrement());
}

KRATOS_TEST_CASE_IN_SUITE(NoNextTimeStepWhenPreviousEndStateDidNotConverge, KratosGeoMechanicsFastSuite)
{
    std::vector<double> increments{0.4, 0.6};
    PrescribedTimeIncrementor incrementor{increments};
    TimeStepEndState previous_state;  // contains non-converged state

    KRATOS_EXPECT_FALSE(incrementor.WantNextStep(previous_state));
}

}
