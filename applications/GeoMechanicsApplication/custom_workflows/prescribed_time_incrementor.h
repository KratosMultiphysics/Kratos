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

#include "includes/kratos_export_api.h"
#include "time_incrementor.h"

#include <vector>

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) PrescribedTimeIncrementor : public TimeIncrementor
{
public:
    explicit PrescribedTimeIncrementor(const std::vector<double>& increments);

    [[nodiscard]] bool WantNextStep(const TimeStepEndState& rPreviousState) const override;
    [[nodiscard]] bool WantRetryStep(std::size_t CycleNumber, const TimeStepEndState& rPreviousState) const override;
    [[nodiscard]] double GetIncrement() const override;
    void                 PostTimeStepExecution(const TimeStepEndState& rResultantState) override;

private:
    std::vector<double>                 mIncrements;
    std::vector<double>::const_iterator mPos;
};

} // namespace Kratos
