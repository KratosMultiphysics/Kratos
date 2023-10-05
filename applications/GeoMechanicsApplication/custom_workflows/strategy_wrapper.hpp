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

        [[nodiscard]] virtual TimeStepEndState::ConvergenceState GetConvergenceState()       = 0;
        [[nodiscard]] virtual std::size_t      GetNumberOfIterations()                 const = 0;
        [[nodiscard]] virtual double           GetEndTime()                            const = 0;
        virtual void Initialize()             = 0;
        virtual void InitializeSolutionStep() = 0;
        virtual void Predict()                = 0;
        virtual bool SolveSolutionStep()      = 0;
        virtual void FinalizeSolutionStep()   = 0;
    };
}