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

#include "time_step_end_state.hpp"

#include <optional>

namespace Kratos
{

class TimeIncrementor
{
public:
    virtual ~TimeIncrementor() = default;

    [[nodiscard]] virtual bool   WantNextStep(const TimeStepEndState& rPreviousState) const     = 0;
    [[nodiscard]] virtual bool   WantRetryStep(std::size_t             CycleNumber,
                                               const TimeStepEndState& rPreviousState) const    = 0;
    [[nodiscard]] virtual double GetIncrement() const                                           = 0;
    virtual void                 PostTimeStepExecution(const TimeStepEndState& rResultantState) = 0;
};

} // namespace Kratos
