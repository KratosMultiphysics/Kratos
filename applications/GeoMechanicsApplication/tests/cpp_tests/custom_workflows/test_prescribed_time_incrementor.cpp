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

#include "custom_workflows/prescribed_time_incrementor.h"
#include "custom_workflows/time_step_end_state.hpp"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace
{

TimeStepEndState MakeConvergedEndState()
{
    auto result              = TimeStepEndState{};
    result.convergence_state = TimeStepEndState::ConvergenceState::converged;
    return result;
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(NoNextTimeStepWhenPrescribedTimeIncrementorHasEmptyList, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto increments     = std::vector<double>{};
    const auto incrementor    = PrescribedTimeIncrementor{increments};
    const auto previous_state = MakeConvergedEndState();

    KRATOS_EXPECT_FALSE(incrementor.WantNextStep(previous_state))
}

KRATOS_TEST_CASE_IN_SUITE(WantFirstTimeStepWhenPrescribedTimeIncrementorHasNonEmptyList,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto increments     = std::vector<double>{0.4, 0.6};
    const auto incrementor    = PrescribedTimeIncrementor{increments};
    const auto previous_state = MakeConvergedEndState();

    KRATOS_EXPECT_TRUE(incrementor.WantNextStep(previous_state))
}

KRATOS_TEST_CASE_IN_SUITE(WantFirstTwoTimesNewStepWhenPrescribedTimeIncrementorHasTwoItems,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto increments     = std::vector<double>{0.4, 0.6};
    auto       incrementor    = PrescribedTimeIncrementor{increments};
    const auto previous_state = MakeConvergedEndState();

    KRATOS_EXPECT_TRUE(incrementor.WantNextStep(previous_state))
    incrementor.PostTimeStepExecution(previous_state);
    KRATOS_EXPECT_TRUE(incrementor.WantNextStep(previous_state))
    incrementor.PostTimeStepExecution(previous_state);
    KRATOS_EXPECT_FALSE(incrementor.WantNextStep(previous_state))
    KRATOS_EXPECT_FALSE(incrementor.WantNextStep(previous_state))
}

KRATOS_TEST_CASE_IN_SUITE(PrescribedTimeIncrementorThrowsIfAnyIncrementIsNegative, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto increments = std::vector<double>{0.4, -0.6};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(PrescribedTimeIncrementor{increments},
                                      "All prescribed increments must not be negative")
}

KRATOS_TEST_CASE_IN_SUITE(WantRetryStepAlwaysReturnsTrueOnFirstCycle, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto increments     = std::vector<double>{0.4, 0.6};
    const auto incrementor    = PrescribedTimeIncrementor{increments};
    const auto previous_state = MakeConvergedEndState();
    const auto cycle_number   = std::size_t{0};

    KRATOS_EXPECT_TRUE(incrementor.WantRetryStep(cycle_number, previous_state))
}

KRATOS_TEST_CASE_IN_SUITE(WantRetryStepAlwaysReturnsFalseOnAnySubsequentCycle, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto increments     = std::vector<double>{0.4, 0.6};
    const auto incrementor    = PrescribedTimeIncrementor{increments};
    const auto previous_state = MakeConvergedEndState();
    const auto cycle_number   = std::size_t{1};

    KRATOS_EXPECT_FALSE(incrementor.WantRetryStep(cycle_number, previous_state))
}

KRATOS_TEST_CASE_IN_SUITE(PrescribedTimeIncrementsMustMatchInput, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto increments     = std::vector<double>{0.4, 0.6};
    auto       incrementor    = PrescribedTimeIncrementor{increments};
    const auto previous_state = MakeConvergedEndState();

    KRATOS_EXPECT_DOUBLE_EQ(increments.front(), incrementor.GetIncrement());
    incrementor.PostTimeStepExecution(previous_state);
    KRATOS_EXPECT_DOUBLE_EQ(increments.back(), incrementor.GetIncrement());
}

KRATOS_TEST_CASE_IN_SUITE(PrescribedTimeIncrementorThrowsWhenAskingForIncrementBeyondEnd,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto increments     = std::vector<double>{0.4, 0.6};
    auto       incrementor    = PrescribedTimeIncrementor{increments};
    const auto previous_state = MakeConvergedEndState();
    incrementor.PostTimeStepExecution(previous_state);
    incrementor.PostTimeStepExecution(previous_state);
    // Now we are beyond the end of the increment vector

    // Note: avoid a warning triggered by the `[[nodiscard]]` attribute of the `GetIncrement()`
    // member function by assigning the return value to a dummy variable. In turn, the dummy
    // variable needs to be marked `[[maybe_unused]]` to avoid a warning about an unused variable.
    KRATOS_EXPECT_EXCEPTION_IS_THROWN([[maybe_unused]] const auto increment = incrementor.GetIncrement(),
                                      "Out of increment range")
}

KRATOS_TEST_CASE_IN_SUITE(WithoutPostTimeStepExecutionAlwaysGetSameIncrement, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto     increments    = std::vector<double>{0.4, 0.6};
    auto           incrementor   = PrescribedTimeIncrementor{increments};

    KRATOS_EXPECT_DOUBLE_EQ(increments.front(), incrementor.GetIncrement());
    KRATOS_EXPECT_DOUBLE_EQ(increments.front(), incrementor.GetIncrement());

    const auto previous_state = MakeConvergedEndState();
    incrementor.PostTimeStepExecution(previous_state);

    KRATOS_EXPECT_DOUBLE_EQ(increments.back(), incrementor.GetIncrement());
    KRATOS_EXPECT_DOUBLE_EQ(increments.back(), incrementor.GetIncrement());
}

KRATOS_TEST_CASE_IN_SUITE(NoNextTimeStepWhenPreviousEndStateDidNotConverge, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto increments     = std::vector<double>{0.4, 0.6};
    const auto incrementor    = PrescribedTimeIncrementor{increments};
    const auto previous_state = TimeStepEndState{}; // contains non-converged state

    KRATOS_EXPECT_FALSE(incrementor.WantNextStep(previous_state))
}

} // namespace Kratos::Testing
