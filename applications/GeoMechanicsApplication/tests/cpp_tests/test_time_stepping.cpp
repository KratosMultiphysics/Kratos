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
#include "custom_workflows/time_stepping.hpp"
#include "custom_workflows/strategy_wrapper.hpp"
#include "solving_strategies/strategies/solving_strategy.h"

using namespace Kratos;


namespace
{

class DummySparseType
{
public:
    struct DataType{};
    struct MatrixType{};
    struct VectorType{};
    struct MatrixPointerType{};
    struct VectorPointerType{};
};

class DummyDenseType
{
public:
    struct MatrixType{};
    struct VectorType{};
};

class NeverConvergingSolverStrategy : public SolvingStrategy<DummySparseType, DummyDenseType>
{
public:
    bool IsConverged()            override {return false;}
};

class AlwaysConvergingSolverStrategy : public SolvingStrategy<DummySparseType, DummyDenseType>
{
public:
    bool IsConverged()            override {return true;}
    void Initialize()             override {mSolverStrategyInitializeCalls += 1;}
    void InitializeSolutionStep() override {mSolverStrategyInitializeSolutionStepCalls += 1;}
    void Predict()                override {mSolverStrategyPredictCalls += 1;}
    bool SolveSolutionStep()      override {mSolverStrategySolveSolutionsStepCalls += 1;
                                            return true;}
    void FinalizeSolutionStep()   override {mSolverStrategyFinalizeSolutionStepCalls += 1;}

private:
    unsigned int mSolverStrategyInitializeCalls             = 0;
    unsigned int mSolverStrategyInitializeSolutionStepCalls = 0;
    unsigned int mSolverStrategyPredictCalls                = 0;
    unsigned int mSolverStrategySolveSolutionsStepCalls     = 0;
    unsigned int mSolverStrategyFinalizeSolutionStepCalls   = 0;
};

class ProcessSpy : public Process
{
public:
    void ExecuteInitializeSolutionStep() override {mSolutionStepInitializedCalls += 1;}
    void ExecuteFinalizeSolutionStep()   override {mSolutionStepFinalizedCalls   += 1;}

    unsigned int NumberOfExecuteInitializeSolutionStepCalls() const {return mSolutionStepInitializedCalls;}
    unsigned int NumberOfExecuteFinalizeSolutionStepCalls()   const {return mSolutionStepFinalizedCalls;}

private:
    unsigned int mSolutionStepInitializedCalls = 0;
    unsigned int mSolutionStepFinalizedCalls   = 0;
};

class DummyAlwaysConvergingStrategyWrapper : public StrategyWrapper
{
public:
    [[nodiscard]] bool        GetConvergenceState(const TimeStepEndState& rEndState)   override { return IsConverged();};
    [[nodiscard]] std::size_t GetNumberOfIterations()                                  const override { return 4;};
    [[nodiscard]] double      GetEndTime()                                             const override { return 10.;};
    bool IsConverged() override {return true; }
    void Initialize()             override {mSolverStrategyInitializeCalls += 1;}
    void InitializeSolutionStep() override {mSolverStrategyInitializeSolutionStepCalls += 1;}
    void Predict()                override {mSolverStrategyPredictCalls += 1;}
    bool SolveSolutionStep()      override {mSolverStrategySolveSolutionsStepCalls += 1;
                                            return true;}
    void FinalizeSolutionStep()   override {mSolverStrategyFinalizeSolutionStepCalls += 1;}
    unsigned int NumberOfSolverStrategyInitializeCalls()             const {return mSolverStrategyInitializeCalls;}
    unsigned int NumberOfSolverStrategyInitializeSolutionStepCalls() const {return mSolverStrategyInitializeSolutionStepCalls;}
    unsigned int NumberOfSolverStrategyPredictCalls()                const {return mSolverStrategyPredictCalls;}
    unsigned int NumberOfSolverStrategySolveSolutionStepCalls()      const {return mSolverStrategySolveSolutionsStepCalls;}
    unsigned int NumberOfSolverStrategyFinalizeSolutionStepCalls()   const {return mSolverStrategyFinalizeSolutionStepCalls;}

private:
    unsigned int mSolverStrategyInitializeCalls             = 0;
    unsigned int mSolverStrategyInitializeSolutionStepCalls = 0;
    unsigned int mSolverStrategyPredictCalls                = 0;
    unsigned int mSolverStrategySolveSolutionsStepCalls     = 0;
    unsigned int mSolverStrategyFinalizeSolutionStepCalls   = 0;

};


class DummyNeverConvergingStrategyWrapper : public StrategyWrapper
{
public:
    [[nodiscard]] bool        GetConvergenceState(const TimeStepEndState& rEndState)   override { return IsConverged();};
    [[nodiscard]] std::size_t GetNumberOfIterations()                                  const override { return 2;};
    [[nodiscard]] double      GetEndTime()                                             const override { return 10.;};
    bool IsConverged() override {return false; }
    void Initialize() override {
        // Deliberately empty method for dummy
    }
    void InitializeSolutionStep() override {
        // Deliberately empty method for dummy
    }
    void Predict() override {
        // Deliberately empty method for dummy
    }
    bool SolveSolutionStep() override {return true;}
    void FinalizeSolutionStep() override {
        // Deliberately empty method for dummy
    }
};

}

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(RunReturnsNonConvergedWhenStrategyDoesNotConverge, KratosGeoMechanicsFastSuite)
{
    TimeStepExecutor executor;
    auto nonconverging_strategy = std::make_shared<DummyNeverConvergingStrategyWrapper>();
    executor.SetSolverStrategy(nonconverging_strategy);
    const auto time = 0.0;
    KRATOS_EXPECT_EQ(TimeStepEndState::ConvergenceState::non_converged, executor.Run(time).convergence_state);
}

KRATOS_TEST_CASE_IN_SUITE(RunReturnsConvergedWhenStrategyConverged, KratosGeoMechanicsFastSuite)
{
    TimeStepExecutor executor;
    auto converging_strategy = std::make_shared<DummyAlwaysConvergingStrategyWrapper>();
    executor.SetSolverStrategy(converging_strategy);
    const auto time = 0.0;
    KRATOS_EXPECT_EQ(TimeStepEndState::ConvergenceState::converged, executor.Run(time).convergence_state);
}

KRATOS_TEST_CASE_IN_SUITE(ProcessMemberFunctionsAllCalledOnce, KratosGeoMechanicsFastSuite)
{
    TimeStepExecutor executor;
    auto converging_strategy = std::make_shared<DummyAlwaysConvergingStrategyWrapper>();
    executor.SetSolverStrategy(converging_strategy);
    ProcessSpy spy;
    TimeStepExecutor::ProcessRefVec process_refs{spy};
    executor.SetProcessReferences(process_refs);
    const auto time = 0.0;

    executor.Run(time);

    KRATOS_EXPECT_EQ(1,spy.NumberOfExecuteInitializeSolutionStepCalls());
    KRATOS_EXPECT_EQ(1,spy.NumberOfExecuteFinalizeSolutionStepCalls());
}

KRATOS_TEST_CASE_IN_SUITE(SolverStrategyMemberFunctionsAllCalledOnce, KratosGeoMechanicsFastSuite)
{
    TimeStepExecutor executor;
    auto converging_strategy = std::make_shared<DummyAlwaysConvergingStrategyWrapper>();
    executor.SetSolverStrategy(converging_strategy);
    ProcessSpy spy;
    TimeStepExecutor::ProcessRefVec process_refs{spy};
    executor.SetProcessReferences(process_refs);
    const auto time = 0.0;

    executor.Run(time);

    KRATOS_EXPECT_EQ(1,converging_strategy->NumberOfSolverStrategyInitializeCalls());
    KRATOS_EXPECT_EQ(1,converging_strategy->NumberOfSolverStrategyInitializeSolutionStepCalls());
    KRATOS_EXPECT_EQ(1,converging_strategy->NumberOfSolverStrategyPredictCalls());
    KRATOS_EXPECT_EQ(1,converging_strategy->NumberOfSolverStrategySolveSolutionStepCalls());
    KRATOS_EXPECT_EQ(1,converging_strategy->NumberOfSolverStrategyFinalizeSolutionStepCalls());
}

KRATOS_TEST_CASE_IN_SUITE(ConvergingTimeStepExecutionReturnsGivenTime, KratosGeoMechanicsFastSuite)
{
    TimeStepExecutor executor;
    auto converging_strategy = std::make_shared<DummyAlwaysConvergingStrategyWrapper>();
    executor.SetSolverStrategy(converging_strategy);
    const auto time = 2.0;
    executor.Run(time);
    KRATOS_EXPECT_DOUBLE_EQ(time, executor.Run(time).time);
}

KRATOS_TEST_CASE_IN_SUITE(NonConvergingTimeStepExecutionReturnsGivenTime, KratosGeoMechanicsFastSuite)
{
    TimeStepExecutor executor;
    auto non_converging_strategy = std::make_shared<DummyNeverConvergingStrategyWrapper>();
    executor.SetSolverStrategy(non_converging_strategy);
    const auto time = 2.0;
    KRATOS_EXPECT_DOUBLE_EQ(time, executor.Run(time).time);
}

KRATOS_TEST_CASE_IN_SUITE(TimeStepExecutionReturnsNumberOfIterations, KratosGeoMechanicsFastSuite)
{
    TimeStepExecutor executor;
    auto converging_strategy = std::make_shared<DummyAlwaysConvergingStrategyWrapper>();
    executor.SetSolverStrategy(converging_strategy);
    const auto time = 2.0;
    KRATOS_EXPECT_EQ(4, executor.Run(time).number_of_nonlinear_iterations);
}

}
