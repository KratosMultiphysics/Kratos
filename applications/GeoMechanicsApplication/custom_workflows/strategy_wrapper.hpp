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

#include <cstddef>

#include "time_step_end_state.hpp"

namespace Kratos
{
class StrategyWrapper
{
public:
    virtual ~StrategyWrapper() = default;

    [[nodiscard]] virtual std::size_t GetNumberOfIterations() const                 = 0;
    [[nodiscard]] virtual double      GetEndTime() const                            = 0;
    virtual void                      SetEndTime(double EndTime)                    = 0;
    [[nodiscard]] virtual double      GetTimeIncrement() const                      = 0;
    virtual void                      SetTimeIncrement(double TimeIncrement)        = 0;
    [[nodiscard]] virtual std::size_t GetStepNumber() const                         = 0;
    virtual void                      IncrementStepNumber()                         = 0;
    virtual void                      CloneTimeStep()                               = 0;
    virtual void                      RestorePositionsAndDOFVectorToStartOfStep()   = 0;
    virtual void                      SaveTotalDisplacementFieldAtStartOfTimeLoop() = 0;
    virtual void                      AccumulateTotalDisplacementField()            = 0;
    virtual void                      OutputProcess()                               = 0;

    virtual void                               Initialize()             = 0;
    virtual void                               InitializeOutput()       = 0;
    virtual void                               InitializeSolutionStep() = 0;
    virtual void                               Predict()                = 0;
    virtual TimeStepEndState::ConvergenceState SolveSolutionStep()      = 0;
    virtual void                               FinalizeSolutionStep()   = 0;
    virtual void                               FinalizeOutput()         = 0;
};
} // namespace Kratos