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
#include <functional>

#include "time_step_end_state.hpp"
#include "time_incrementor.h"
#include "strategy_wrapper.hpp"

namespace Kratos
{

class Process;

class TimeLoopExecutorInterface{
public :
    virtual ~TimeLoopExecutorInterface() = default;
    virtual void SetCancelDelegate(const std::function<bool()>& rCancelDelegate) = 0;
    virtual void SetProgressDelegate(const std::function<void(double)>& rProgressDelegate) = 0;
    virtual void SetProcessObservables(const std::vector<std::weak_ptr<Process>>& rProcessObservables) = 0;
    virtual void SetTimeIncrementor(std::unique_ptr<TimeIncrementor> pTimeIncrementor) = 0;
    virtual void SetSolverStrategyWrapper(std::shared_ptr<StrategyWrapper> pStrategyWrapper) = 0;
    virtual std::vector<TimeStepEndState> Run(const TimeStepEndState& EndState) = 0;
};

}