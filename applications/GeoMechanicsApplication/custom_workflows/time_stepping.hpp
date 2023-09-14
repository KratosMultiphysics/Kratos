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

#pragma once

#include <memory>


namespace Kratos
{

template <typename StrategyType>
class TimeStepExecuter
{
public:
    enum class ConvergenceState {converged, non_converged};

    void SetSolverStrategy(std::shared_ptr<StrategyType> SolverStrategy)
    {
        mStrategy = std::move(SolverStrategy);
    }

    ConvergenceState Run()
    {
        mStrategy->Initialize();
        mStrategy->InitializeSolutionStep();

        mStrategy->Predict();
        mStrategy->SolveSolutionStep();

        mStrategy->FinalizeSolutionStep();

        return mStrategy->IsConverged() ? ConvergenceState::converged : ConvergenceState::non_converged;
    }

private:
    std::shared_ptr<StrategyType> mStrategy;
};

}
