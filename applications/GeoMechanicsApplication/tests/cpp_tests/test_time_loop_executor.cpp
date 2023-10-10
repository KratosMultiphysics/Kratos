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
#include "custom_workflows/time_loop_executor.hpp"
#include "custom_workflows/time_step_end_state.hpp"
#include "custom_workflows/prescribed_time_incrementor.h"
#include "solving_strategies/strategies/solving_strategy.h"

#include <numeric>

using namespace Kratos;

namespace
{

class DummySolverStrategy : public StrategyWrapper
{
public:
    explicit DummySolverStrategy(TimeStepEndState::ConvergenceState AConvergenceState) : mConvergenceState{AConvergenceState} {}
    [[nodiscard]] TimeStepEndState::ConvergenceState GetConvergenceState()         override
    {
        return mConvergenceState;
    }
    [[nodiscard]] std::size_t                        GetNumberOfIterations() const override
    {
        return 1;
    }
    [[nodiscard]] double                             GetEndTime() const            override
    {
        return 0.0;
    }
    void Initialize()             override {}
    void InitializeSolutionStep() override {}
    void Predict()                override {}
    bool SolveSolutionStep()      override { return false; }
    void FinalizeSolutionStep()   override {}

private:
    TimeStepEndState::ConvergenceState mConvergenceState{TimeStepEndState::ConvergenceState::converged};
};

class FixedCyclesTimeIncrementor : public TimeIncrementor
{
public:
    explicit FixedCyclesTimeIncrementor(std::size_t DesiredNumOfCyclesPerStep)
    : mNumCyclesPerStep{DesiredNumOfCyclesPerStep}
    {}

    [[nodiscard]] bool WantNextStep(const TimeStepEndState& rPreviousState) const override
    {
        return rPreviousState.time < mEndTime;
    }
    [[nodiscard]] bool WantRetryStep(std::size_t             CycleNumber,
                                     const TimeStepEndState& rPreviousState) const override
    {
        return CycleNumber < mNumCyclesPerStep;
    }
    [[nodiscard]] double GetIncrement() const override
    {
        return 0.5;
    }
    void PostTimeStepExecution(const TimeStepEndState& rResultantState) override
    {}

private:
    std::size_t mNumCyclesPerStep;
    double      mEndTime{1.0};
};

}

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TimeLoopReturnsPerformedStatesAfterRunningAnAlwaysConvergingSolverStrategy, KratosGeoMechanicsFastSuite)
{
    TimeLoopExecutor executor;
    const auto increments = std::vector<double>{1.1, 1.2, 1.3, 1.4};
    executor.SetTimeIncrementor(std::make_unique<PrescribedTimeIncrementor>(increments));
    executor.SetSolverStrategyTimeIncrementor(std::make_unique<DummySolverStrategy>(TimeStepEndState::ConvergenceState::converged));
    TimeStepEndState start_state;
    start_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    const auto step_states = executor.Run(start_state);

    KRATOS_EXPECT_EQ(increments.size(), step_states.size());
    KRATOS_EXPECT_TRUE(std::all_of(step_states.begin(), step_states.end(), [](const auto& step_state) { return step_state.Converged(); }))
}

KRATOS_TEST_CASE_IN_SUITE(TimeLoopReturnsOneNonConvergedPerformedStatesAfterRunningAnNeverConvergingSolverStrategy, KratosGeoMechanicsFastSuite)
{
    TimeLoopExecutor executor;
    const auto increments = std::vector<double>{1.1, 1.2, 1.3, 1.4};
    executor.SetTimeIncrementor(std::make_unique<PrescribedTimeIncrementor>(increments));
    executor.SetSolverStrategyTimeIncrementor(std::make_unique<DummySolverStrategy>(TimeStepEndState::ConvergenceState::non_converged));
    TimeStepEndState start_state;
    start_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    const auto step_states = executor.Run(start_state);

    KRATOS_EXPECT_EQ(1, step_states.size());
    KRATOS_EXPECT_TRUE(step_states[0].NonConverged());
}

KRATOS_TEST_CASE_IN_SUITE(TimeLoopReturnsEndTimesAfterRunningAnAlwaysConvergingSolverStrategy, KratosGeoMechanicsFastSuite)
{
    TimeLoopExecutor executor;
    const auto increments = std::vector<double>{1.1, 1.2, 1.3, 1.4};
    executor.SetTimeIncrementor(std::make_unique<PrescribedTimeIncrementor>(increments));
    executor.SetSolverStrategyTimeIncrementor(std::make_unique<DummySolverStrategy>(TimeStepEndState::ConvergenceState::converged));
    TimeStepEndState start_state;
    start_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    const auto step_states = executor.Run(start_state);

    KRATOS_EXPECT_DOUBLE_EQ(std::accumulate(increments.begin(), increments.end(), 0.), step_states.back().time);
}

KRATOS_TEST_CASE_IN_SUITE(GetThreeCyclesWhenThreeCyclesAreNeeded, KratosGeoMechanicsFastSuite)
{
    TimeLoopExecutor executor;
    const auto wanted_num_of_cycles_per_step = 3;
    executor.SetTimeIncrementor(std::make_unique<FixedCyclesTimeIncrementor>(wanted_num_of_cycles_per_step));
    executor.SetSolverStrategyTimeIncrementor(std::make_unique<DummySolverStrategy>(TimeStepEndState::ConvergenceState::converged));
    const auto step_states = executor.Run(TimeStepEndState{});

    // The test time incrementor assumes fixed time increments of 0.5 and an end time of 1.0
    KRATOS_EXPECT_EQ(2, step_states.size());
    KRATOS_EXPECT_DOUBLE_EQ(0.5, step_states[0].time);
    KRATOS_EXPECT_EQ(wanted_num_of_cycles_per_step, step_states[0].num_of_cycles);
    KRATOS_EXPECT_DOUBLE_EQ(1.0, step_states[1].time);
    KRATOS_EXPECT_EQ(wanted_num_of_cycles_per_step, step_states[1].num_of_cycles);
}

KRATOS_TEST_CASE_IN_SUITE(GetFiveCyclesWhenFiveCyclesAreNeeded, KratosGeoMechanicsFastSuite)
{
    TimeLoopExecutor executor;
    const auto wanted_num_of_cycles_per_step = 5;
    executor.SetTimeIncrementor(std::make_unique<FixedCyclesTimeIncrementor>(wanted_num_of_cycles_per_step));
    executor.SetSolverStrategyTimeIncrementor(std::make_unique<DummySolverStrategy>(TimeStepEndState::ConvergenceState::converged));
    const auto step_states = executor.Run(TimeStepEndState{});

    // The test time incrementor assumes fixed time increments of 0.5 and an end time of 1.0
    KRATOS_EXPECT_EQ(2, step_states.size());
    KRATOS_EXPECT_DOUBLE_EQ(0.5, step_states[0].time);
    KRATOS_EXPECT_EQ(wanted_num_of_cycles_per_step, step_states[0].num_of_cycles);
    KRATOS_EXPECT_DOUBLE_EQ(1.0, step_states[1].time);
    KRATOS_EXPECT_EQ(wanted_num_of_cycles_per_step, step_states[1].num_of_cycles);
}

}