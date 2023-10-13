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
#include "containers/model.h"

#include <numeric>

using namespace Kratos;

namespace
{

class DummySolverStrategy : public StrategyWrapper
{
public:
    explicit DummySolverStrategy(TimeStepEndState::ConvergenceState AConvergenceState)
    : mConvergenceState{AConvergenceState}
    {
        mModel.Reset();
        mpModelPart = &mModel.CreateModelPart("TestModelPart");
        mpModelPart->Reset();
    }

    [[nodiscard]] TimeStepEndState::ConvergenceState GetConvergenceState() override
    {
        return mConvergenceState;
    }

    [[nodiscard]] std::size_t GetNumberOfIterations() const override
    {
        return 1;
    }

    [[nodiscard]] double GetEndTime() const override
    {
        return mpModelPart->GetProcessInfo()[TIME];
    }

    void SetEndTime(double EndTime) override
    {
        mpModelPart->GetProcessInfo()[TIME] = EndTime;
    }

    [[nodiscard]] double GetTimeIncrement() const override
    {
        return mpModelPart->GetProcessInfo()[DELTA_TIME];
    }

    void SetTimeIncrement(double TimeIncrement) override
    {
        mpModelPart->GetProcessInfo()[DELTA_TIME] = TimeIncrement;
    }

    [[nodiscard]] std::size_t GetStepNumber() const override
    {
        return static_cast<std::size_t>(mpModelPart->GetProcessInfo()[STEP]);
    }

    void IncrementStepNumber() override
    {
        mpModelPart->GetProcessInfo()[STEP] += 1;
    }

    void CloneTimeStep() override
    {
        mIsCloned = true;
    }

    [[nodiscard]] bool IsCloned() const
    {
        return mIsCloned;
    }

    void RestorePositionsAndDOFVectorToStartOfStep() override
    {
        mIsRestoreCalled = true;
    }

    [[nodiscard]] bool IsRestoreCalled() const
    {
        return mIsRestoreCalled;
    }

    void SaveTotalDisplacementFieldAtStartOfStage() override
    {
        mIsSaveTotalDisplacementFieldCalled = true;
    }

    [[nodiscard]] bool IsSaveFieldCalled() const
    {
        return mIsSaveTotalDisplacementFieldCalled;
    }

    void AccumulateTotalDisplacementField() override
    {
        mCountAccumulateTotalDisplacementFieldCalled += 1;
    }

    [[nodiscard]] std::size_t GetCountAccumulateTotalDisplacementFieldCalled() const
    {
        return mCountAccumulateTotalDisplacementFieldCalled;
    }

    void OutputProcess() override
    {
        mCountOutputProcessCalled += 1;
    }

    [[nodiscard]] std::size_t GetCountOutputProcessCalled() const
    {
        return mCountOutputProcessCalled;
    }

    [[nodiscard]] std::size_t GetCountFinalizeSolutionStepCalled() const
    {
        return mCountFinalizeSolutionStepCalled;
    }

    void Initialize()             override {
        // intentionally empty
    }
    void InitializeSolutionStep() override {
        // intentionally empty
    }
    void Predict()                override {
        // intentionally empty
    }
    bool SolveSolutionStep()      override { return false; }
    void FinalizeSolutionStep()   override {mCountFinalizeSolutionStepCalled += 1;}

private:
    TimeStepEndState::ConvergenceState mConvergenceState{TimeStepEndState::ConvergenceState::converged};
    Model mModel;
    ModelPart * mpModelPart = nullptr;
    bool mIsCloned = false;
    bool mIsRestoreCalled = false;
    bool mIsSaveTotalDisplacementFieldCalled = false;
    std::size_t mCountAccumulateTotalDisplacementFieldCalled = 0;
    std::size_t mCountOutputProcessCalled = 0;
    std::size_t mCountFinalizeSolutionStepCalled = 0;
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
    {
        // intentionally empty
    }

private:
    std::size_t mNumCyclesPerStep;
    double      mEndTime{1.0};
};

TimeStepEndState MakeConvergedStepState()
{
    auto result = TimeStepEndState{};
    result.convergence_state = TimeStepEndState::ConvergenceState::converged;
    return result;
}

std::shared_ptr<DummySolverStrategy> RunFixedCycleTimeLoop(std::size_t WantedNumOfCyclesPerStep)
{
    TimeLoopExecutor executor;
    executor.SetTimeIncrementor(std::make_unique<FixedCyclesTimeIncrementor>(WantedNumOfCyclesPerStep));
    auto solver_strategy = std::make_shared<DummySolverStrategy>(TimeStepEndState::ConvergenceState::converged);
    executor.SetSolverStrategyTimeStepExecutor(solver_strategy);
    const auto step_states = executor.Run(TimeStepEndState{});

    return solver_strategy;
}

}

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TimeLoopReturnsPerformedStatesAfterRunningAnAlwaysConvergingSolverStrategy, KratosGeoMechanicsFastSuite)
{
    TimeLoopExecutor executor;
    const auto increments = std::vector<double>{1.1, 1.2, 1.3, 1.4};
    executor.SetTimeIncrementor(std::make_unique<PrescribedTimeIncrementor>(increments));
    executor.SetSolverStrategyTimeStepExecutor(
            std::make_unique<DummySolverStrategy>(TimeStepEndState::ConvergenceState::converged));
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
    executor.SetSolverStrategyTimeStepExecutor(
            std::make_unique<DummySolverStrategy>(TimeStepEndState::ConvergenceState::non_converged));
    TimeStepEndState start_state;
    start_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    const auto step_states = executor.Run(start_state);

    KRATOS_EXPECT_EQ(1, step_states.size());
    KRATOS_EXPECT_TRUE(step_states[0].NonConverged())
}

KRATOS_TEST_CASE_IN_SUITE(TimeLoopReturnsEndTimesAfterRunningAnAlwaysConvergingSolverStrategy, KratosGeoMechanicsFastSuite)
{
    TimeLoopExecutor executor;
    const auto increments = std::vector<double>{1.1, 1.2, 1.3, 1.4};
    executor.SetTimeIncrementor(std::make_unique<PrescribedTimeIncrementor>(increments));
    executor.SetSolverStrategyTimeStepExecutor(
            std::make_unique<DummySolverStrategy>(TimeStepEndState::ConvergenceState::converged));
    TimeStepEndState start_state;
    start_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
    const auto step_states = executor.Run(start_state);

    KRATOS_EXPECT_DOUBLE_EQ(std::accumulate(increments.begin(), increments.end(), 0.), step_states.back().time);
}

KRATOS_TEST_CASE_IN_SUITE(GetFiveCyclesWhenFiveCyclesAreNeeded, KratosGeoMechanicsFastSuite)
{
    TimeLoopExecutor executor;
    const auto wanted_num_of_cycles_per_step = 5;
    executor.SetTimeIncrementor(std::make_unique<FixedCyclesTimeIncrementor>(wanted_num_of_cycles_per_step));
    executor.SetSolverStrategyTimeStepExecutor(
            std::make_unique<DummySolverStrategy>(TimeStepEndState::ConvergenceState::converged));
    const auto step_states = executor.Run(TimeStepEndState{});

    // The test time incrementor assumes fixed time increments of 0.5 and an end time of 1.0
    KRATOS_EXPECT_EQ(2, step_states.size());
    KRATOS_EXPECT_DOUBLE_EQ(0.5, step_states[0].time);
    KRATOS_EXPECT_EQ(wanted_num_of_cycles_per_step, step_states[0].num_of_cycles);
    KRATOS_EXPECT_DOUBLE_EQ(1.0, step_states[1].time);
    KRATOS_EXPECT_EQ(wanted_num_of_cycles_per_step, step_states[1].num_of_cycles);
}

KRATOS_TEST_CASE_IN_SUITE(ExpectOneCyclePerStepWhenUsingAPrescribedTimeIncrementor, KratosGeoMechanicsFastSuite)
{
    TimeLoopExecutor executor;
    const auto increments = std::vector<double>{0.1, 0.2, 0.3, 0.4};
    executor.SetTimeIncrementor(std::make_unique<PrescribedTimeIncrementor>(increments));
    executor.SetSolverStrategyTimeStepExecutor(
            std::make_unique<DummySolverStrategy>(TimeStepEndState::ConvergenceState::converged));
    const auto step_states = executor.Run(MakeConvergedStepState());

    KRATOS_EXPECT_EQ(increments.size(), step_states.size());
    KRATOS_EXPECT_TRUE(std::all_of(step_states.begin(), step_states.end(),
                       [](const auto& state) {return state.num_of_cycles == 1;}))
    auto time = 0.0;
    for (auto index = std::size_t{0}; index < increments.size(); ++index) {
        time += increments[index];
        KRATOS_EXPECT_DOUBLE_EQ(time, step_states[index].time);
    }
}

KRATOS_TEST_CASE_IN_SUITE(ExpectEndTimeToBeSetOnSolverStrategyAfterRunningAStepLoop, KratosGeoMechanicsFastSuite)
{
    auto solver_strategy = RunFixedCycleTimeLoop(1);

    KRATOS_EXPECT_DOUBLE_EQ(1.0, solver_strategy->GetEndTime());
    KRATOS_EXPECT_DOUBLE_EQ(0.5, solver_strategy->GetTimeIncrement());
    KRATOS_EXPECT_EQ(2, solver_strategy->GetStepNumber());
}

KRATOS_TEST_CASE_IN_SUITE(ExpectATimeStepCloneAtStartOfStep, KratosGeoMechanicsFastSuite)
{
    auto solver_strategy = RunFixedCycleTimeLoop(1);

    KRATOS_EXPECT_TRUE(solver_strategy->IsCloned())
}

KRATOS_TEST_CASE_IN_SUITE(ExpectNoRestoreCalledForOneCycle, KratosGeoMechanicsFastSuite)
{
    auto solver_strategy = RunFixedCycleTimeLoop(1);

    KRATOS_EXPECT_FALSE(solver_strategy->IsRestoreCalled())
}

KRATOS_TEST_CASE_IN_SUITE(ExpectRestoreCalledForTwoCycles, KratosGeoMechanicsFastSuite)
{
    auto solver_strategy = RunFixedCycleTimeLoop(2);

    KRATOS_EXPECT_TRUE(solver_strategy->IsRestoreCalled())
}

KRATOS_TEST_CASE_IN_SUITE(ExpectDisplacementFieldStoredForResetDisplacements, KratosGeoMechanicsFastSuite)
{
    auto solver_strategy = RunFixedCycleTimeLoop(1);

    KRATOS_EXPECT_TRUE(solver_strategy->IsSaveFieldCalled())
}

KRATOS_TEST_CASE_IN_SUITE(ExpectDisplacementFieldUpdateForEveryStep, KratosGeoMechanicsFastSuite)
{
    auto solver_strategy = RunFixedCycleTimeLoop(1);

    KRATOS_EXPECT_EQ(2, solver_strategy->GetCountAccumulateTotalDisplacementFieldCalled());
}

KRATOS_TEST_CASE_IN_SUITE(ExpectOutputProcessCalledForEveryStep, KratosGeoMechanicsFastSuite)
{
    TimeLoopExecutor executor;
    const auto wanted_num_of_cycles_per_step = 1;
    executor.SetTimeIncrementor(std::make_unique<FixedCyclesTimeIncrementor>(wanted_num_of_cycles_per_step));
    auto solver_strategy = std::make_shared<DummySolverStrategy>(TimeStepEndState::ConvergenceState::converged);
    executor.SetSolverStrategyTimeStepExecutor(solver_strategy);
    const auto step_states = executor.Run(TimeStepEndState{});
    KRATOS_EXPECT_EQ(step_states.size(), solver_strategy->GetCountOutputProcessCalled());
}

KRATOS_TEST_CASE_IN_SUITE(ExpectFinalizeSolutionStepCalledOnceForEveryStep, KratosGeoMechanicsFastSuite)
{
    TimeLoopExecutor executor;
    const auto wanted_num_of_cycles_per_step = 3;
    executor.SetTimeIncrementor(std::make_unique<FixedCyclesTimeIncrementor>(wanted_num_of_cycles_per_step));
    auto solver_strategy = std::make_shared<DummySolverStrategy>(TimeStepEndState::ConvergenceState::converged);
    executor.SetSolverStrategyTimeStepExecutor(solver_strategy);
    const auto step_states = executor.Run(TimeStepEndState{});
    KRATOS_EXPECT_EQ(step_states.size(), solver_strategy->GetCountFinalizeSolutionStepCalled());
}

}