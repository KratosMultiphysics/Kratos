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
#include "custom_workflows/adaptive_time_incrementor.h"
#include "custom_workflows/time_step_end_state.hpp"

using namespace Kratos;


namespace
{

struct AdaptiveTimeIncrementorSettings
{
    double      StartTime{0.0};
    double      EndTime{8.0};
    double      StartIncrement{0.5};
    std::size_t MaxNumOfCycles{8};
    double      ReductionFactor{0.5};
    double      IncreaseFactor{2.0};
    double      MaxDeltaTimeFactor{1000.0};
    std::size_t MinNumOfIterations{3};
    std::size_t MaxNumOfIterations{15};
};

AdaptiveTimeIncrementor MakeAdaptiveTimeIncrementor(const AdaptiveTimeIncrementorSettings& rSettings)
{
    return AdaptiveTimeIncrementor{rSettings.StartTime,
                                   rSettings.EndTime,
                                   rSettings.StartIncrement,
                                   rSettings.MaxNumOfCycles,
                                   rSettings.ReductionFactor,
                                   rSettings.IncreaseFactor,
                                   rSettings.MaxDeltaTimeFactor,
                                   rSettings.MinNumOfIterations,
                                   rSettings.MaxNumOfIterations};
}

}


namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(AdaptiveTimeIncrementorThrowsWhenStartTimeExceedsEndTime, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.EndTime = settings.StartTime - 2.0;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Start time (0) must be smaller than end time (-2)")
}

KRATOS_TEST_CASE_IN_SUITE(AdaptiveTimeIncrementorThrowsWhenStartTimeEqualsEndTime, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.EndTime = settings.StartTime;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Start time (0) must be smaller than end time (0)")
}

KRATOS_TEST_CASE_IN_SUITE(AdaptiveTimeIncrementorThrowsWhenStartIncrementIsNegative, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.StartIncrement = -0.5;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Start increment must be positive, but got -0.5")
}

KRATOS_TEST_CASE_IN_SUITE(AdaptiveTimeIncrementorThrowsWhenStartIncrementEqualsZero, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.StartIncrement = 0.0;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Start increment must be positive, but got 0")
}

KRATOS_TEST_CASE_IN_SUITE(AdaptiveTimeIncrementorThrowsWhenMaxNumberOfCyclesEqualsZero, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.MaxNumOfCycles = std::size_t{0};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Maximum number of cycles must be positive")
}

KRATOS_TEST_CASE_IN_SUITE(AdaptiveTimeIncrementorThrowsWhenReductionFactorIsGreaterThanOne, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.ReductionFactor = 2.0;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Reduction factor must not be greater than 1, but got 2")
}

KRATOS_TEST_CASE_IN_SUITE(AdaptiveTimeIncrementorThrowsWhenReductionFactorIsNegative, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.ReductionFactor = -0.5;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Reduction factor must be positive, but got -0.5")
}

KRATOS_TEST_CASE_IN_SUITE(AdaptiveTimeIncrementorThrowsWhenReductionFactorEqualsZero, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.ReductionFactor = 0.0;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Reduction factor must be positive, but got 0")
}

KRATOS_TEST_CASE_IN_SUITE(AdaptiveTimeIncrementorDoesNotThrowWhenReductionFactorIsInRange, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.ReductionFactor = 0.5;
    auto has_thrown = false;

    try
    {
        MakeAdaptiveTimeIncrementor(settings);
    }
    catch (const Exception&)
    {
        has_thrown = true;
    }

    KRATOS_EXPECT_FALSE(has_thrown) // No other way to check that the constructor does not throw
}

KRATOS_TEST_CASE_IN_SUITE(AdaptiveTimeIncrementorDoesNotThrowWhenReductionFactorEqualsOne, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.ReductionFactor = 1.0;
    auto has_thrown = false;

    try
    {
        MakeAdaptiveTimeIncrementor(settings);
    }
    catch (const Exception&)
    {
        has_thrown = true;
    }

    KRATOS_EXPECT_FALSE(has_thrown) // No other way to check that the constructor does not throw
}

KRATOS_TEST_CASE_IN_SUITE(AdaptiveTimeIncrementorThrowsWhenIncreaseFactorIsSmallerThanOne, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.IncreaseFactor = 0.5;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Increase factor must be greater than or equal to 1, but got 0.5")
}

KRATOS_TEST_CASE_IN_SUITE(AdaptiveTimeIncrementorThrowsWhenMaxDeltaTimeFactorIsSmallerThanOne, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.MaxDeltaTimeFactor = 0.9;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Max_delta_time_factor must be greater than or equal to 1, but got 0.9")
}

KRATOS_TEST_CASE_IN_SUITE(AdaptiveTimeIncrementorDoesNotThrowWhenIncreaseFactorEqualsOne, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.IncreaseFactor = 1.0;
    auto has_thrown = false;

    try
    {
        MakeAdaptiveTimeIncrementor(settings);
    }
    catch (const Exception&)
    {
        has_thrown = true;
    }

    KRATOS_EXPECT_FALSE(has_thrown) // No other way to check that the constructor does not throw
}

KRATOS_TEST_CASE_IN_SUITE(AdaptiveTimeIncrementorThrowsWhenMaxNumberOfIterationsEqualsZero, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.MinNumOfIterations = std::size_t{0};
    settings.MaxNumOfIterations = std::size_t{0};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Minimum number of iterations (0) is not less than maximum number of iterations (0)")
}

KRATOS_TEST_CASE_IN_SUITE(AdaptiveTimeIncrementorThrowsWhenMinNumberOfIterationsIsGreaterThanMaxNumberOfIterations, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.MaxNumOfIterations = std::size_t{10};
    settings.MinNumOfIterations = settings.MaxNumOfIterations + 1;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Minimum number of iterations (11) is not less than maximum number of iterations (10)")
}

KRATOS_TEST_CASE_IN_SUITE(WantNextTimeStepWhenMoreThanStartIncrementLeftUntilEndTime, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    const auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto       previous_state   = TimeStepEndState{};
    previous_state.time         = settings.EndTime - 2 * settings.StartIncrement;

    KRATOS_EXPECT_TRUE(time_incrementor.WantNextStep(previous_state))
}

KRATOS_TEST_CASE_IN_SUITE(WantNextTimeStepWhenLessThanStartIncrementLeftUntilEndTime, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    const auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto       previous_state   = TimeStepEndState{};
    previous_state.time         = settings.EndTime - 0.5 * settings.StartIncrement;

    KRATOS_EXPECT_TRUE(time_incrementor.WantNextStep(previous_state))
}

KRATOS_TEST_CASE_IN_SUITE(NoNextTimeStepWhenEndTimeIsExceeded, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    const auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto       previous_state   = TimeStepEndState{};
    previous_state.time         = settings.EndTime + 1.0;

    KRATOS_EXPECT_FALSE(time_incrementor.WantNextStep(previous_state))
}

KRATOS_TEST_CASE_IN_SUITE(NoNextTimeStepWhenAtEndTime, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    const auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto       previous_state   = TimeStepEndState{};
    previous_state.time         = settings.EndTime;

    KRATOS_EXPECT_FALSE(time_incrementor.WantNextStep(previous_state))
}

KRATOS_TEST_CASE_IN_SUITE(RetryWhenNotAttemptedYet, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.MaxNumOfCycles     = std::size_t{4};
    const auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    const auto cycle_number     = std::size_t{0};
    const auto previous_state   = TimeStepEndState{};

    KRATOS_EXPECT_TRUE(time_incrementor.WantRetryStep(cycle_number, previous_state))
}

KRATOS_TEST_CASE_IN_SUITE(RetryWhenPreviousAttemptDidNotConvergeButAtLeastOneMoreLeft, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.MaxNumOfCycles     = std::size_t{4};
    const auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    const auto cycle_number     = settings.MaxNumOfCycles - 1;
    auto previous_state         = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::non_converged;

    KRATOS_EXPECT_TRUE(time_incrementor.WantRetryStep(cycle_number, previous_state))
}

KRATOS_TEST_CASE_IN_SUITE(DontRetryWhenPreviousAttemptDidNotConvergeAndNoAttemptsLeft, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.MaxNumOfCycles     = std::size_t{4};
    const auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto previous_state         = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::non_converged;

    KRATOS_EXPECT_FALSE(time_incrementor.WantRetryStep(settings.MaxNumOfCycles, previous_state))
}

KRATOS_TEST_CASE_IN_SUITE(DontRetryWhenPreviousAttemptConverged, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.MaxNumOfCycles     = std::size_t{4};
    const auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    const auto cycle_number     = std::size_t{1};
    auto previous_state         = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::converged;

    KRATOS_EXPECT_FALSE(time_incrementor.WantRetryStep(cycle_number, previous_state))
}

KRATOS_TEST_CASE_IN_SUITE(GetStartIncrementWhenItWouldNotResultInExceedingTheEndTime, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.StartIncrement = 0.6;
    const auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);

    KRATOS_EXPECT_DOUBLE_EQ(settings.StartIncrement, time_incrementor.GetIncrement());
}

KRATOS_TEST_CASE_IN_SUITE(ReduceStartIncrementWhenItWouldResultInExceedingTheEndTime, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.StartTime      = 0.0;
    settings.StartIncrement = 1.0;
    settings.EndTime        = settings.StartTime + 0.5 * settings.StartIncrement; // EndTime would be exceeded if StartIncrement is applied
    const auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);

    KRATOS_EXPECT_DOUBLE_EQ(settings.EndTime - settings.StartTime, time_incrementor.GetIncrement());
}

KRATOS_TEST_CASE_IN_SUITE(ReduceIncrementWhenPreviousAttemptDidNotConverge, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto previous_state   = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::non_converged;

    time_incrementor.PostTimeStepExecution(previous_state); // process previous non-converged state
    KRATOS_EXPECT_DOUBLE_EQ(settings.StartIncrement * settings.ReductionFactor, time_incrementor.GetIncrement());
}

KRATOS_TEST_CASE_IN_SUITE(ReduceIncrementEvenMoreWhenPreviousTwoAttemptsDidNotConverge, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto previous_state   = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::non_converged;

    time_incrementor.PostTimeStepExecution(previous_state); // process first non-converged state
    time_incrementor.PostTimeStepExecution(previous_state); // process second non-converged state
    KRATOS_EXPECT_DOUBLE_EQ(settings.StartIncrement * settings.ReductionFactor * settings.ReductionFactor,
                            time_incrementor.GetIncrement());
}

KRATOS_TEST_CASE_IN_SUITE(ReduceIncrementWhenStepConvergedAndMaxNumberOfIterationsWasAttained, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto previous_state   = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    previous_state.num_of_iterations = settings.MaxNumOfIterations;

    time_incrementor.PostTimeStepExecution(previous_state); // previous attempt converged, but required maximum number of iterations
    KRATOS_EXPECT_DOUBLE_EQ(settings.StartIncrement * settings.ReductionFactor, time_incrementor.GetIncrement());
}

KRATOS_TEST_CASE_IN_SUITE(IncreaseIncrementWhenStepRequiredLessThanMinNumberOfIterations, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto previous_state   = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    previous_state.num_of_iterations = settings.MinNumOfIterations - 1;

    time_incrementor.PostTimeStepExecution(previous_state); // previous attempt converged and required less than minimum number of iterations
    KRATOS_EXPECT_DOUBLE_EQ(settings.StartIncrement * settings.IncreaseFactor, time_incrementor.GetIncrement());
}

KRATOS_TEST_CASE_IN_SUITE(ReduceIncrementToAvoidExceedingEndTime, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto previous_state   = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    previous_state.num_of_iterations = settings.MinNumOfIterations + 1; // don't attempt to increase the increment
    previous_state.time              = settings.EndTime - 0.5 * settings.StartIncrement; // only half of StartIncrement left before reaching EndTime

    time_incrementor.PostTimeStepExecution(previous_state);
    KRATOS_EXPECT_DOUBLE_EQ(0.5 * settings.StartIncrement, time_incrementor.GetIncrement());
}

KRATOS_TEST_CASE_IN_SUITE(ReduceUpscaledIncrementToAvoidExceedingEndTime, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto previous_state   = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    previous_state.num_of_iterations = settings.MinNumOfIterations - 1; // this should normally increase the increment
    previous_state.time              = settings.EndTime - 0.5 * settings.StartIncrement; // only half of StartIncrement left before reaching EndTime

    time_incrementor.PostTimeStepExecution(previous_state);
    KRATOS_EXPECT_DOUBLE_EQ(0.5 * settings.StartIncrement, time_incrementor.GetIncrement());
}

KRATOS_TEST_CASE_IN_SUITE(ReduceUpscaledIncrementToAvoidExceedingMaxDeltaTimeFactor, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.MaxDeltaTimeFactor = 1.5; // lower than the increase factor, this should lead to truncation
    auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto previous_state   = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    previous_state.num_of_iterations = settings.MinNumOfIterations - 1; // this should normally increase the increment

    time_incrementor.PostTimeStepExecution(previous_state);
    KRATOS_EXPECT_DOUBLE_EQ(settings.MaxDeltaTimeFactor * settings.StartIncrement, time_incrementor.GetIncrement());
}

KRATOS_TEST_CASE_IN_SUITE(ScaleIncrementToAvoidExtraSmallTimeStep, KratosGeoMechanicsFastSuite)
{
    AdaptiveTimeIncrementorSettings settings; // with EndTime = 8.0
    settings.StartIncrement = 7.9999;
    auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto previous_state   = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::converged;

    time_incrementor.PostTimeStepExecution(previous_state);
    KRATOS_EXPECT_DOUBLE_EQ(8.0, time_incrementor.GetIncrement());
}

}
