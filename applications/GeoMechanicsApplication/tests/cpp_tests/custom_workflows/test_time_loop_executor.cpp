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

#include "containers/model.h"
#include "custom_workflows/adaptive_time_incrementor.h"
#include "custom_workflows/prescribed_time_incrementor.h"
#include "custom_workflows/time_loop_executor.hpp"
#include "custom_workflows/time_step_end_state.hpp"
#include "solving_strategies/strategies/solving_strategy.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <numeric>

using namespace Kratos;

namespace
{

class DummySolverStrategy : public StrategyWrapper
{
public:
    DummySolverStrategy()
    {
        mModel.Reset();
        mpModelPart = &mModel.CreateModelPart("TestModelPart");
        mpModelPart->Reset();
        ON_CALL(*this, SolveSolutionStep()).WillByDefault(testing::Return(TimeStepEndState::ConvergenceState::converged));
    }

    [[nodiscard]] double GetEndTime() const override { return mpModelPart->GetProcessInfo()[TIME]; }

    void SetEndTime(double EndTime) override { mpModelPart->GetProcessInfo()[TIME] = EndTime; }

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

    void IncrementStepNumber() override { ++mpModelPart->GetProcessInfo()[STEP]; }

    void SetOutputProcessCallback(std::function<void()> Callback)
    {
        mOutputProcessCallback = std::move(Callback);
    }

    void OutputProcess() override
    {
        ++mCountOutputProcessCalled;
        if (mOutputProcessCallback) mOutputProcessCallback();
    }

    [[nodiscard]] std::size_t GetCountOutputProcessCalled() const
    {
        return mCountOutputProcessCalled;
    }

    MOCK_METHOD(void, Initialize, (), (override));
    MOCK_METHOD(void, InitializeSolutionStep, (), (override));
    MOCK_METHOD(TimeStepEndState::ConvergenceState, SolveSolutionStep, (), (override));
    MOCK_METHOD(void, FinalizeSolutionStep, (), (override));

    MOCK_METHOD(void, Predict, (), (override));

    MOCK_METHOD(void, ComputeIncrementalDisplacementField, (), (override));
    MOCK_METHOD(void, AccumulateTotalDisplacementField,(), (override));

    MOCK_METHOD(void, CloneTimeStep, (), (override));
    MOCK_METHOD(void, RestorePositionsAndDOFVectorToStartOfStep, (), (override));
    MOCK_METHOD(std::size_t, GetNumberOfIterations, (), (const, override));

    MOCK_METHOD(void, InitializeOutput,(), (override));
    MOCK_METHOD(void, FinalizeOutput,(), (override));

private:
    Model                 mModel;
    ModelPart*            mpModelPart                                  = nullptr;
    std::size_t           mCountOutputProcessCalled                    = 0;
    std::function<void()> mOutputProcessCallback;
};

class MockProcess : public Process
{
public:
    MOCK_METHOD(void, ExecuteBeforeSolutionLoop, (), (override));
};

class FixedCyclesTimeIncrementor : public TimeIncrementor
{
public:
    explicit FixedCyclesTimeIncrementor(std::size_t DesiredNumOfCyclesPerStep)
        : mNumCyclesPerStep{DesiredNumOfCyclesPerStep}
    {
    }

    [[nodiscard]] bool WantNextStep(const TimeStepEndState& rPreviousState) const override
    {
        return rPreviousState.time < mEndTime;
    }

    [[nodiscard]] bool WantRetryStep(std::size_t CycleNumber, const TimeStepEndState& rPreviousState) const override
    {
        return CycleNumber < mNumCyclesPerStep;
    }

    [[nodiscard]] double GetIncrement() const override { return 0.5; }

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
    auto result              = TimeStepEndState{};
    result.convergence_state = TimeStepEndState::ConvergenceState::converged;
    return result;
}

void RunFixedCycleTimeLoop(std::size_t WantedNumOfCyclesPerStep, const std::shared_ptr<DummySolverStrategy>& rpSolver)
{
    TimeLoopExecutor executor;
    executor.SetTimeIncrementor(std::make_unique<FixedCyclesTimeIncrementor>(WantedNumOfCyclesPerStep));
    executor.SetSolverStrategyWrapper(rpSolver);

    const auto step_states = executor.Run(TimeStepEndState{});
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TimeLoopReturnsPerformedStatesAfterRunningAnAlwaysConvergingSolverStrategy,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TimeLoopExecutor executor;
    const auto       increments = std::vector<double>{1.1, 1.2, 1.3, 1.4};
    executor.SetTimeIncrementor(std::make_unique<PrescribedTimeIncrementor>(increments));
    auto converging_strategy = std::make_shared<DummySolverStrategy>();
    EXPECT_CALL(*converging_strategy, SolveSolutionStep())
        .WillRepeatedly(testing::Return(TimeStepEndState::ConvergenceState::converged));
    executor.SetSolverStrategyWrapper(converging_strategy);

    const auto step_states = executor.Run(MakeConvergedStepState());

    KRATOS_EXPECT_EQ(increments.size(), step_states.size());
    KRATOS_EXPECT_TRUE(std::all_of(step_states.begin(), step_states.end(),
                                   [](const auto& step_state) { return step_state.Converged(); }))
}

KRATOS_TEST_CASE_IN_SUITE(TimeLoopExecutorThrowsAfterRunningANeverConvergingSolverStrategy,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TimeLoopExecutor executor;
    const auto       increments = std::vector<double>{1.1, 1.2, 1.3, 1.4};
    executor.SetTimeIncrementor(std::make_unique<PrescribedTimeIncrementor>(increments));
    auto p_solving_strategy = std::make_shared<DummySolverStrategy>();
    EXPECT_CALL(*p_solving_strategy, SolveSolutionStep()).WillOnce(testing::Return(TimeStepEndState::ConvergenceState::non_converged));
    executor.SetSolverStrategyWrapper(p_solving_strategy);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(executor.Run(MakeConvergedStepState()),
                                      "The calculation exited without converging.");
}

KRATOS_TEST_CASE_IN_SUITE(TimeLoopExecutorThrowsAfterRunningAdaptiveTimeIncrementorWithNeverConvergingSolverStrategy,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TimeLoopExecutor executor;
    executor.SetTimeIncrementor(std::make_unique<AdaptiveTimeIncrementor>(0.0, 1.0, 1.0, 2, 0.5));
    auto p_solving_strategy = std::make_shared<DummySolverStrategy>();
    EXPECT_CALL(*p_solving_strategy, SolveSolutionStep())
        .WillRepeatedly(testing::Return(TimeStepEndState::ConvergenceState::non_converged));
    executor.SetSolverStrategyWrapper(p_solving_strategy);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(executor.Run(MakeConvergedStepState()),
                                      "The calculation exited without converging.");
}

KRATOS_TEST_CASE_IN_SUITE(TimeLoopReturnsEndTimesAfterRunningAnAlwaysConvergingSolverStrategy,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TimeLoopExecutor executor;
    const auto       increments = std::vector<double>{1.1, 1.2, 1.3, 1.4};
    executor.SetTimeIncrementor(std::make_unique<PrescribedTimeIncrementor>(increments));
    auto converging_strategy = std::make_shared<DummySolverStrategy>();
    EXPECT_CALL(*converging_strategy, SolveSolutionStep())
        .WillRepeatedly(testing::Return(TimeStepEndState::ConvergenceState::converged));
    executor.SetSolverStrategyWrapper(converging_strategy);

    const auto step_states = executor.Run(MakeConvergedStepState());

    KRATOS_EXPECT_DOUBLE_EQ(std::accumulate(increments.begin(), increments.end(), 0.),
                            step_states.back().time);
}

KRATOS_TEST_CASE_IN_SUITE(GetFiveCyclesWhenFiveCyclesAreNeeded, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TimeLoopExecutor executor;
    const auto       wanted_num_of_cycles_per_step = 5;
    executor.SetTimeIncrementor(std::make_unique<FixedCyclesTimeIncrementor>(wanted_num_of_cycles_per_step));
    auto converging_strategy = std::make_shared<DummySolverStrategy>();
    executor.SetSolverStrategyWrapper(converging_strategy);

    const auto step_states = executor.Run(TimeStepEndState{});

    // The test time incrementor assumes fixed time increments of 0.5 and an end time of 1.0
    KRATOS_EXPECT_EQ(2, step_states.size());
    KRATOS_EXPECT_DOUBLE_EQ(0.5, step_states[0].time);
    KRATOS_EXPECT_EQ(wanted_num_of_cycles_per_step, step_states[0].num_of_cycles);
    KRATOS_EXPECT_DOUBLE_EQ(1.0, step_states[1].time);
    KRATOS_EXPECT_EQ(wanted_num_of_cycles_per_step, step_states[1].num_of_cycles);
}

KRATOS_TEST_CASE_IN_SUITE(ExpectOneCyclePerStepWhenUsingAPrescribedTimeIncrementor, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TimeLoopExecutor executor;
    const auto       increments = std::vector<double>{0.1, 0.2, 0.3, 0.4};
    executor.SetTimeIncrementor(std::make_unique<PrescribedTimeIncrementor>(increments));
    auto converging_strategy = std::make_shared<DummySolverStrategy>();
    executor.SetSolverStrategyWrapper(converging_strategy);

    const auto step_states = executor.Run(MakeConvergedStepState());

    KRATOS_EXPECT_EQ(increments.size(), step_states.size());
    KRATOS_EXPECT_TRUE(std::all_of(step_states.begin(), step_states.end(),
                                   [](const auto& state) { return state.num_of_cycles == 1; }))
    auto time = 0.0;
    for (auto index = std::size_t{0}; index < increments.size(); ++index) {
        time += increments[index];
        KRATOS_EXPECT_DOUBLE_EQ(time, step_states[index].time);
    }
}

KRATOS_TEST_CASE_IN_SUITE(ExpectEndTimeToBeSetOnSolverStrategyAfterRunningAStepLoop, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto solver_strategy = std::make_shared<DummySolverStrategy>();
    RunFixedCycleTimeLoop(1, solver_strategy);

    KRATOS_EXPECT_DOUBLE_EQ(1.0, solver_strategy->GetEndTime());
    KRATOS_EXPECT_DOUBLE_EQ(0.5, solver_strategy->GetTimeIncrement());
    KRATOS_EXPECT_EQ(2, solver_strategy->GetStepNumber());
}

KRATOS_TEST_CASE_IN_SUITE(ExpectATimeStepCloneAtStartOfStep, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto solver_strategy = std::make_shared<DummySolverStrategy>();

    EXPECT_CALL(*solver_strategy, CloneTimeStep).Times(testing::AtLeast(1));
    RunFixedCycleTimeLoop(1, solver_strategy);
}

KRATOS_TEST_CASE_IN_SUITE(ExpectNoRestoreCalledForOneCycle, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto solver_strategy = std::make_shared<DummySolverStrategy>();

    EXPECT_CALL(*solver_strategy, RestorePositionsAndDOFVectorToStartOfStep).Times(0);
    RunFixedCycleTimeLoop(1, solver_strategy);
}

KRATOS_TEST_CASE_IN_SUITE(ExpectRestoreCalledForTwoCycles, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto solver_strategy = std::make_shared<DummySolverStrategy>();

    EXPECT_CALL(*solver_strategy, RestorePositionsAndDOFVectorToStartOfStep).Times(testing::AtLeast(1));
    RunFixedCycleTimeLoop(2, solver_strategy);
}

KRATOS_TEST_CASE_IN_SUITE(ExpectDisplacementFieldUpdateForEveryStep, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto solver_strategy = std::make_shared<DummySolverStrategy>();

    EXPECT_CALL(*solver_strategy, AccumulateTotalDisplacementField).Times(2);
    RunFixedCycleTimeLoop(1, solver_strategy);
}

KRATOS_TEST_CASE_IN_SUITE(ExpectOutputProcessCalledForEveryStep, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TimeLoopExecutor executor;
    const auto       wanted_num_of_cycles_per_step = 1;
    executor.SetTimeIncrementor(std::make_unique<FixedCyclesTimeIncrementor>(wanted_num_of_cycles_per_step));
    auto solver_strategy = std::make_shared<DummySolverStrategy>();
    executor.SetSolverStrategyWrapper(solver_strategy);

    const auto step_states = executor.Run(TimeStepEndState{});
    KRATOS_EXPECT_EQ(step_states.size(), solver_strategy->GetCountOutputProcessCalled());
    KRATOS_EXPECT_EQ(2, step_states.size());
}

KRATOS_TEST_CASE_IN_SUITE(ExpectFinalizeSolutionStepCalledOnceForEveryStep, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TimeLoopExecutor executor;
    const auto       wanted_num_of_cycles_per_step = 3;
    executor.SetTimeIncrementor(std::make_unique<FixedCyclesTimeIncrementor>(wanted_num_of_cycles_per_step));
    auto solver_strategy = std::make_shared<DummySolverStrategy>();
    executor.SetSolverStrategyWrapper(solver_strategy);

    EXPECT_CALL(*solver_strategy, FinalizeSolutionStep).Times(2);
    const auto step_states = executor.Run(TimeStepEndState{});
    KRATOS_EXPECT_EQ(2, step_states.size());
}

KRATOS_TEST_CASE_IN_SUITE(ExpectOutputIsInitializedAndFinalizedWhenRunCompletesOk, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TimeLoopExecutor executor;
    const auto       wanted_num_of_cycles_per_step = std::size_t{1};
    executor.SetTimeIncrementor(std::make_unique<FixedCyclesTimeIncrementor>(wanted_num_of_cycles_per_step));
    auto solver_strategy = std::make_shared<DummySolverStrategy>();
    executor.SetSolverStrategyWrapper(solver_strategy);

    EXPECT_CALL(*solver_strategy, InitializeOutput).Times(1);
    EXPECT_CALL(*solver_strategy, FinalizeOutput).Times(1);
    executor.Run(TimeStepEndState{});
}

KRATOS_TEST_CASE_IN_SUITE(ExpectOutputIsInitializedAndFinalizedWhenRunThrows, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TimeLoopExecutor executor;
    const auto       wanted_num_of_cycles_per_step = std::size_t{1};
    executor.SetTimeIncrementor(std::make_unique<FixedCyclesTimeIncrementor>(wanted_num_of_cycles_per_step));
    auto solver_strategy = std::make_shared<DummySolverStrategy>();
    solver_strategy->SetOutputProcessCallback([]() { throw Exception{"Test exception"}; });
    executor.SetSolverStrategyWrapper(solver_strategy);

    EXPECT_CALL(*solver_strategy, InitializeOutput).Times(1);
    EXPECT_CALL(*solver_strategy, FinalizeOutput).Times(1);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(executor.Run(MakeConvergedStepState()), "Test exception");
}

KRATOS_TEST_CASE_IN_SUITE(TimeLoopExecutor_CallsProcessExecuteBeforeSolutionLoop_AfterInitialize,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    TimeLoopExecutor executor;
    const auto       increments = std::vector<double>{1.1, 1.2, 1.3, 1.4};
    executor.SetTimeIncrementor(std::make_unique<PrescribedTimeIncrementor>(increments));

    auto solving_strategy = std::make_shared<DummySolverStrategy>();
    executor.SetSolverStrategyWrapper(solving_strategy);

    auto                                p_process = std::make_shared<MockProcess>();
    std::vector<std::weak_ptr<Process>> process_observables{p_process};
    executor.SetProcessObservables(process_observables);

    EXPECT_CALL(*p_process, ExecuteBeforeSolutionLoop).After(EXPECT_CALL(*solving_strategy, Initialize));
    executor.Run(MakeConvergedStepState());
}

} // namespace Kratos::Testing
