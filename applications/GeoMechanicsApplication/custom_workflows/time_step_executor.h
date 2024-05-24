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

#include "processes/process.h"
#include "strategy_wrapper.hpp"
#include "time_step_end_state.hpp"

#include <functional>
#include <memory>

namespace Kratos
{

class TimeStepExecutor
{
public:
    using ProcessRef    = std::reference_wrapper<Process>;
    using ProcessRefVec = std::vector<ProcessRef>;

    void SetSolverStrategy(std::shared_ptr<StrategyWrapper> SolverStrategy);
    void SetProcessObservables(const std::vector<std::weak_ptr<Process>>& rProcessObservables);
    TimeStepEndState Run(double Time);

private:
    std::shared_ptr<StrategyWrapper>    mStrategyWrapper;
    std::vector<std::weak_ptr<Process>> mProcessObservables;
};

} // namespace Kratos
