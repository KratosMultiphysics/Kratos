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

}


namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(RunReturnsNonConvergedWhenStrategyDoesNotConverge, KratosGeoMechanicsFastSuite)
{
    TimeStepExecuter<NeverConvergingSolverStrategy> executer;
    auto nonconverging_strategy = std::make_shared<NeverConvergingSolverStrategy>();
    executer.SetSolverStrategy(nonconverging_strategy);
    KRATOS_EXPECT_EQ(TimeStepExecuter<NeverConvergingSolverStrategy>::ConvergenceState::non_converged, executer.Run());
}

KRATOS_TEST_CASE_IN_SUITE(RunReturnsConvergedWhenStrategyConverged, KratosGeoMechanicsFastSuite)
{
    TimeStepExecuter<AlwaysConvergingSolverStrategy> executer;
    auto converging_strategy = std::make_shared<AlwaysConvergingSolverStrategy>();
    executer.SetSolverStrategy(converging_strategy);
    KRATOS_EXPECT_EQ(TimeStepExecuter<AlwaysConvergingSolverStrategy>::ConvergenceState::converged, executer.Run());
}

}
