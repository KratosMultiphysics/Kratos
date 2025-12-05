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
//                   Richard Faasse
//                   Gennady Markelov
//

#include "custom_workflows/adaptive_time_incrementor.h"
#include "custom_workflows/time_step_end_state.hpp"
#include "includes/expect.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"

using namespace Kratos;

namespace
{

struct AdaptiveTimeIncrementorSettings {
    double                StartTime{0.0};
    double                EndTime{8.0};
    double                StartIncrement{0.5};
    std::optional<double> UserMinDeltaTime{std::make_optional(1e-06)};
    std::size_t           MaxNumOfCycles{8};
    double                ReductionFactor{0.5};
    double                IncreaseFactor{2.0};
    double                MaxDeltaTimeFactor{1000.0};
    std::size_t           MinNumOfIterations{3};
    std::size_t           MaxNumOfIterations{15};
};

AdaptiveTimeIncrementor MakeAdaptiveTimeIncrementor(const AdaptiveTimeIncrementorSettings& rSettings)
{
    return AdaptiveTimeIncrementor{rSettings.StartTime,          rSettings.EndTime,
                                   rSettings.StartIncrement,     rSettings.MaxNumOfCycles,
                                   rSettings.ReductionFactor,    rSettings.IncreaseFactor,
                                   rSettings.UserMinDeltaTime,   rSettings.MaxDeltaTimeFactor,
                                   rSettings.MinNumOfIterations, rSettings.MaxNumOfIterations};
}

} // namespace

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, AdaptiveTimeIncrementorThrowsWhenStartTimeExceedsEndTime)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.EndTime = settings.StartTime - 2.0;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Start time (0) must be smaller than end time (-2)")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, AdaptiveTimeIncrementorThrowsWhenStartTimeEqualsEndTime)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.EndTime = settings.StartTime;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Start time (0) must be smaller than end time (0)")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, AdaptiveTimeIncrementorThrowsWhenStartIncrementIsNegative)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.StartIncrement = -0.5;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Start increment must be positive, but got -0.5")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, AdaptiveTimeIncrementorThrowsWhenStartIncrementEqualsZero)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.StartIncrement = 0.0;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Start increment must be positive, but got 0")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, AdaptiveTimeIncrementorThrowsWhenMaxNumberOfCyclesEqualsZero)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.MaxNumOfCycles = std::size_t{0};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Maximum number of cycles must be positive")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, AdaptiveTimeIncrementorThrowsWhenReductionFactorIsGreaterThanOne)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.ReductionFactor = 2.0;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Reduction factor must not be greater than 1, but got 2")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, AdaptiveTimeIncrementorThrowsWhenReductionFactorIsNegative)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.ReductionFactor = -0.5;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Reduction factor must be positive, but got -0.5")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, AdaptiveTimeIncrementorThrowsWhenReductionFactorEqualsZero)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.ReductionFactor = 0.0;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MakeAdaptiveTimeIncrementor(settings),
                                      "Reduction factor must be positive, but got 0")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, AdaptiveTimeIncrementorDoesNotThrowWhenReductionFactorIsInRange)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.ReductionFactor = 0.5;
    auto has_thrown          = false;

    try {
        MakeAdaptiveTimeIncrementor(settings);
    } catch (const Exception&) {
        has_thrown = true;
    }

    KRATOS_EXPECT_FALSE(has_thrown) // No other way to check that the constructor does not throw
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, AdaptiveTimeIncrementorDoesNotThrowWhenReductionFactorEqualsOne)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.ReductionFactor = 1.0;
    auto has_thrown          = false;

    try {
        MakeAdaptiveTimeIncrementor(settings);
    } catch (const Exception&) {
        has_thrown = true;
    }

    KRATOS_EXPECT_FALSE(has_thrown) // No other way to check that the constructor does not throw
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, AdaptiveTimeIncrementorThrowsWhenIncreaseFactorIsSmallerThanOne)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.IncreaseFactor = 0.5;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MakeAdaptiveTimeIncrementor(settings),
        "Increase factor must be greater than or equal to 1, but got 0.5")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, AdaptiveTimeIncrementorThrowsWhenMaxDeltaTimeFactorIsSmallerThanOne)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.MaxDeltaTimeFactor = 0.9;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MakeAdaptiveTimeIncrementor(settings),
        "Max_delta_time_factor must be greater than or equal to 1, but got 0.9")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, AdaptiveTimeIncrementorDoesNotThrowWhenIncreaseFactorEqualsOne)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.IncreaseFactor = 1.0;
    auto has_thrown         = false;

    try {
        MakeAdaptiveTimeIncrementor(settings);
    } catch (const Exception&) {
        has_thrown = true;
    }

    KRATOS_EXPECT_FALSE(has_thrown) // No other way to check that the constructor does not throw
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, AdaptiveTimeIncrementorThrowsWhenMaxNumberOfIterationsEqualsZero)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.MinNumOfIterations = std::size_t{0};
    settings.MaxNumOfIterations = std::size_t{0};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MakeAdaptiveTimeIncrementor(settings),
        "Minimum number of iterations (0) is not less than maximum number of iterations (0)")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       AdaptiveTimeIncrementorThrowsWhenMinNumberOfIterationsIsGreaterThanMaxNumberOfIterations)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.MaxNumOfIterations = std::size_t{10};
    settings.MinNumOfIterations = settings.MaxNumOfIterations + 1;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        MakeAdaptiveTimeIncrementor(settings),
        "Minimum number of iterations (11) is not less than maximum number of iterations (10)")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, WantNextTimeStepWhenMoreThanStartIncrementLeftUntilEndTime)
{
    AdaptiveTimeIncrementorSettings settings;
    const auto                      time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto                            previous_state   = TimeStepEndState{};
    previous_state.time = settings.EndTime - 2 * settings.StartIncrement;

    KRATOS_EXPECT_TRUE(time_incrementor.WantNextStep(previous_state))
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, WantNextTimeStepWhenLessThanStartIncrementLeftUntilEndTime)
{
    AdaptiveTimeIncrementorSettings settings;
    const auto                      time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto                            previous_state   = TimeStepEndState{};
    previous_state.time = settings.EndTime - 0.5 * settings.StartIncrement;

    KRATOS_EXPECT_TRUE(time_incrementor.WantNextStep(previous_state))
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, NoNextTimeStepWhenEndTimeIsExceeded)
{
    AdaptiveTimeIncrementorSettings settings;
    const auto                      time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto                            previous_state   = TimeStepEndState{};
    previous_state.time                              = settings.EndTime + 1.0;

    KRATOS_EXPECT_FALSE(time_incrementor.WantNextStep(previous_state))
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, NoNextTimeStepWhenAtEndTime)
{
    AdaptiveTimeIncrementorSettings settings;
    const auto                      time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto                            previous_state   = TimeStepEndState{};
    previous_state.time                              = settings.EndTime;

    KRATOS_EXPECT_FALSE(time_incrementor.WantNextStep(previous_state))
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, RetryWhenNotAttemptedYet)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.MaxNumOfCycles     = std::size_t{4};
    const auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    const auto cycle_number     = std::size_t{0};
    const auto previous_state   = TimeStepEndState{};

    KRATOS_EXPECT_TRUE(time_incrementor.WantRetryStep(cycle_number, previous_state))
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, RetryWhenPreviousAttemptDidNotConvergeButAtLeastOneMoreLeft)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.MaxNumOfCycles          = std::size_t{4};
    const auto time_incrementor      = MakeAdaptiveTimeIncrementor(settings);
    const auto cycle_number          = settings.MaxNumOfCycles - 1;
    auto       previous_state        = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::non_converged;

    KRATOS_EXPECT_TRUE(time_incrementor.WantRetryStep(cycle_number, previous_state))
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, DontRetryWhenPreviousAttemptDidNotConvergeAndNoAttemptsLeft)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.MaxNumOfCycles          = std::size_t{4};
    const auto time_incrementor      = MakeAdaptiveTimeIncrementor(settings);
    auto       previous_state        = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::non_converged;

    KRATOS_EXPECT_FALSE(time_incrementor.WantRetryStep(settings.MaxNumOfCycles, previous_state))
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, DontRetryWhenPreviousAttemptConverged)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.MaxNumOfCycles          = std::size_t{4};
    const auto time_incrementor      = MakeAdaptiveTimeIncrementor(settings);
    const auto cycle_number          = std::size_t{1};
    auto       previous_state        = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::converged;

    KRATOS_EXPECT_FALSE(time_incrementor.WantRetryStep(cycle_number, previous_state))
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GetStartIncrementWhenItWouldNotResultInExceedingTheEndTime)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.StartIncrement     = 0.6;
    const auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);

    KRATOS_EXPECT_DOUBLE_EQ(settings.StartIncrement, time_incrementor.GetIncrement());
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, ReduceStartIncrementWhenItWouldResultInExceedingTheEndTime)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.StartTime      = 0.0;
    settings.StartIncrement = 1.0;
    settings.EndTime = settings.StartTime + 0.5 * settings.StartIncrement; // EndTime would be exceeded if StartIncrement is applied
    const auto time_incrementor = MakeAdaptiveTimeIncrementor(settings);

    KRATOS_EXPECT_DOUBLE_EQ(settings.EndTime - settings.StartTime, time_incrementor.GetIncrement());
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, ReduceIncrementWhenPreviousAttemptDidNotConverge)
{
    AdaptiveTimeIncrementorSettings settings;
    auto                            time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto                            previous_state   = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::non_converged;

    time_incrementor.PostTimeStepExecution(previous_state); // process previous non-converged state
    KRATOS_EXPECT_DOUBLE_EQ(settings.StartIncrement * settings.ReductionFactor, time_incrementor.GetIncrement());
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, ReduceIncrementEvenMoreWhenPreviousTwoAttemptsDidNotConverge)
{
    AdaptiveTimeIncrementorSettings settings;
    auto                            time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto                            previous_state   = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::non_converged;

    time_incrementor.PostTimeStepExecution(previous_state); // process first non-converged state
    time_incrementor.PostTimeStepExecution(previous_state); // process second non-converged state
    KRATOS_EXPECT_DOUBLE_EQ(settings.StartIncrement * settings.ReductionFactor * settings.ReductionFactor,
                            time_incrementor.GetIncrement());
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, ReduceIncrementWhenStepConvergedAndMaxNumberOfIterationsWasAttained)
{
    AdaptiveTimeIncrementorSettings settings;
    auto                            time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto                            previous_state   = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    previous_state.num_of_iterations = settings.MaxNumOfIterations;

    time_incrementor.PostTimeStepExecution(previous_state); // previous attempt converged, but required maximum number of iterations
    KRATOS_EXPECT_DOUBLE_EQ(settings.StartIncrement * settings.ReductionFactor, time_incrementor.GetIncrement());
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, IncreaseIncrementWhenStepRequiredLessThanMinNumberOfIterations)
{
    AdaptiveTimeIncrementorSettings settings;
    auto                            time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto                            previous_state   = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    previous_state.num_of_iterations = settings.MinNumOfIterations - 1;

    time_incrementor.PostTimeStepExecution(previous_state); // previous attempt converged and required less than minimum number of iterations
    KRATOS_EXPECT_DOUBLE_EQ(settings.StartIncrement * settings.IncreaseFactor, time_incrementor.GetIncrement());
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, ReduceIncrementToAvoidExceedingEndTime)
{
    AdaptiveTimeIncrementorSettings settings;
    auto                            time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto                            previous_state   = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    previous_state.num_of_iterations = settings.MinNumOfIterations + 1; // don't attempt to increase the increment
    previous_state.time =
        settings.EndTime - 0.5 * settings.StartIncrement; // only half of StartIncrement left before reaching EndTime

    time_incrementor.PostTimeStepExecution(previous_state);
    KRATOS_EXPECT_DOUBLE_EQ(0.5 * settings.StartIncrement, time_incrementor.GetIncrement());
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, ReduceUpscaledIncrementToAvoidExceedingEndTime)
{
    AdaptiveTimeIncrementorSettings settings;
    auto                            time_incrementor = MakeAdaptiveTimeIncrementor(settings);
    auto                            previous_state   = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    previous_state.num_of_iterations =
        settings.MinNumOfIterations - 1; // this should normally increase the increment
    previous_state.time =
        settings.EndTime - 0.5 * settings.StartIncrement; // only half of StartIncrement left before reaching EndTime

    time_incrementor.PostTimeStepExecution(previous_state);
    KRATOS_EXPECT_DOUBLE_EQ(0.5 * settings.StartIncrement, time_incrementor.GetIncrement());
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, ReduceUpscaledIncrementToAvoidExceedingMaxDeltaTimeFactor)
{
    AdaptiveTimeIncrementorSettings settings;
    settings.MaxDeltaTimeFactor = 1.5; // lower than the increase factor, this should lead to truncation
    auto time_incrementor            = MakeAdaptiveTimeIncrementor(settings);
    auto previous_state              = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    previous_state.num_of_iterations =
        settings.MinNumOfIterations - 1; // this should normally increase the increment

    time_incrementor.PostTimeStepExecution(previous_state);
    KRATOS_EXPECT_DOUBLE_EQ(settings.MaxDeltaTimeFactor * settings.StartIncrement,
                            time_incrementor.GetIncrement());
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, ScaleIncrementToAvoidExtraSmallTimeStep)
{
    AdaptiveTimeIncrementorSettings settings; // with EndTime = 8.0
    settings.StartIncrement          = 7.9999;
    auto time_incrementor            = MakeAdaptiveTimeIncrementor(settings);
    auto previous_state              = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::converged;

    time_incrementor.PostTimeStepExecution(previous_state);
    KRATOS_EXPECT_DOUBLE_EQ(8.0, time_incrementor.GetIncrement());
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, ThrowExceptionWhenDeltaTimeSmallerThanTheLimit)
{
    AdaptiveTimeIncrementorSettings settings; // with EndTime = 8.0
    settings.StartIncrement          = 8.0;   // We jump to the end time right away
    auto time_incrementor            = MakeAdaptiveTimeIncrementor(settings);
    auto previous_state              = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    previous_state.time              = 7.9999999; // to have a zero time step

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        time_incrementor.PostTimeStepExecution(previous_state),
        "Delta time (1e-07) is smaller than given minimum allowable value 1e-06");

    previous_state.convergence_state = TimeStepEndState::ConvergenceState::non_converged;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        time_incrementor.PostTimeStepExecution(previous_state),
        "Delta time (5e-08) is smaller than given minimum allowable value 1e-06");
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, HalfTimeStepAtNonConverged)
{
    AdaptiveTimeIncrementorSettings settings; // with EndTime = 8.0
    settings.StartIncrement          = 8.0;   // We jump to the end time right away
    auto time_incrementor            = MakeAdaptiveTimeIncrementor(settings);
    auto previous_state              = TimeStepEndState{};
    previous_state.convergence_state = TimeStepEndState::ConvergenceState::non_converged;
    previous_state.time              = 8.0;

    time_incrementor.PostTimeStepExecution(previous_state);
    // The increment should be halved, since the step didn't converge
    KRATOS_EXPECT_DOUBLE_EQ(4.0, time_incrementor.GetIncrement());
}

} // namespace Kratos::Testing
