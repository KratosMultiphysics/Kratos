// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//                   Anne van de Graaf
//                   Wijtze Pieter Kikstra
//

#pragma once

#include <vector>
#include <memory>

#include "time_step_end_state.hpp"
#include "time_incrementor.hpp"
#include "processes/process.h"

namespace Kratos
{

class Process;

template <typename StrategyType>
class TimeLoopExecutor{
public :
    virtual ~TimeLoopExecutor() = default;
    virtual void SetProcessReferences(const std::vector<std::weak_ptr<Process>>& rProcessRefs) = 0;

    void SetTimeIncrementor(std::unique_ptr<TimeIncrementor> ATimeIncrementor)
    {
        mTimeIncrementor = std::move(ATimeIncrementor);
    }

    std::vector<TimeStepEndState> Run()
    {
        std::vector<TimeStepEndState> result;
        TimeStepEndState              step_state;
        if (no_previous_step_state_available) {
            step_state.time = 0.0; // mrModelPart->GetProcessInfo()[TIME];
            step_state.convergence_state = TimeStepEndState::ConvergenceState::converged;
            step_state.num_of_iterations = 0;
            result.push_back(step_state);
        }
        while (mTimeIncrementor->WantNextStep(step_state)) {
            //mrModelPart->CloneTimeStep();
            // nodig voor TimeStepExecutor: solving strategy
            mTimeStepExecutor.run(step_state.time);
            result.push_back(step_state);
        }

        return result;
    }

private:
    std::unique_ptr<TimeIncrementor> mTimeIncrementor;
    TimeStepExecutor<StrategyType>   mTimeStepExecutor;
    ModelPart*                       mrModelPart = nullptr;
};

}