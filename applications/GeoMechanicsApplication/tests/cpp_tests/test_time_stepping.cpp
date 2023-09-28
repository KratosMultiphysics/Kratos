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
    bool IsConverged() override {return false;}
};

class AlwaysConvergingSolverStrategy : public SolvingStrategy<DummySparseType, DummyDenseType>
{
public:
    bool IsConverged() override {return true;}
};

class ProcessSpy : public Process
{
public:
    void ExecuteInitializeSolutionStep() override {mWasSolutionStepInitialized = true;}
    void ExecuteFinalizeSolutionStep()   override {mWasSolutionStepFinalized   = true;}

    bool WasSolutionStepInitialized() const {return mWasSolutionStepInitialized;}
    bool WasSolutionStepFinalized()   const {return mWasSolutionStepFinalized;}

private:
    bool mWasSolutionStepInitialized = false;
    bool mWasSolutionStepFinalized   = false;
};

}


namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(RunReturnsNonConvergedWhenStrategyDoesNotConverge, KratosGeoMechanicsWip)
{
    TimeStepExecuter<NeverConvergingSolverStrategy> executer;
    auto nonconverging_strategy = std::make_shared<NeverConvergingSolverStrategy>();
    executer.SetSolverStrategy(nonconverging_strategy);
    const auto time = 0.0;
    KRATOS_EXPECT_EQ(TimeStepEndState::ConvergenceState::non_converged, executer.Run(time).convergence_state);
}

KRATOS_TEST_CASE_IN_SUITE(RunReturnsConvergedWhenStrategyConverged, KratosGeoMechanicsWip)
{
    TimeStepExecuter<AlwaysConvergingSolverStrategy> executer;
    auto converging_strategy = std::make_shared<AlwaysConvergingSolverStrategy>();
    executer.SetSolverStrategy(converging_strategy);
    const auto time = 0.0;
    KRATOS_EXPECT_EQ(TimeStepEndState::ConvergenceState::converged, executer.Run(time).convergence_state);
}

KRATOS_TEST_CASE_IN_SUITE(ProcessSolutionStepWasInitializedAndFinalized, KratosGeoMechanicsWip)
{
    TimeStepExecuter<AlwaysConvergingSolverStrategy> executer;
    auto converging_strategy = std::make_shared<AlwaysConvergingSolverStrategy>();
    executer.SetSolverStrategy(converging_strategy);
    ProcessSpy spy;
    TimeStepExecuter<AlwaysConvergingSolverStrategy>::ProcessRefVec process_refs{spy};
    executer.SetProcessReferences(process_refs);
    const auto time = 0.0;

    executer.Run(time);

    KRATOS_EXPECT_TRUE(spy.WasSolutionStepInitialized());
    KRATOS_EXPECT_TRUE(spy.WasSolutionStepFinalized());
}

KRATOS_TEST_CASE_IN_SUITE(TimeStepExecutionReturnsGivenTime, KratosGeoMechanicsWip)
{
    TimeStepExecuter<AlwaysConvergingSolverStrategy> executer;
    auto converging_strategy = std::make_shared<AlwaysConvergingSolverStrategy>();
    executer.SetSolverStrategy(converging_strategy);
    const auto time = 2.0;
    KRATOS_EXPECT_DOUBLE_EQ(time, executer.Run(time).time);
}

}
