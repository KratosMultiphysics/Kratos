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

#include "prescribed_time_incrementor.h"

#include "includes/exception.h"

#include <algorithm>

namespace Kratos
{

PrescribedTimeIncrementor::PrescribedTimeIncrementor(const std::vector<double>& increments)
    : mIncrements(increments), mPos(mIncrements.begin())
{
    KRATOS_ERROR_IF(std::any_of(mIncrements.begin(), mIncrements.end(), [](auto value) {
        return value < 0.0;
    })) << "All prescribed increments must not be negative";
}

bool PrescribedTimeIncrementor::WantNextStep(const Kratos::TimeStepEndState& rPreviousState) const
{
    if (rPreviousState.NonConverged()) return false;

    return mPos != mIncrements.end();
}

bool PrescribedTimeIncrementor::WantRetryStep(std::size_t CycleNumber, const TimeStepEndState& rPreviousState) const
{
    return CycleNumber == 0;
}

double PrescribedTimeIncrementor::GetIncrement() const
{
    KRATOS_ERROR_IF(mPos == mIncrements.end()) << "Out of increment range";

    return *mPos;
}

void PrescribedTimeIncrementor::PostTimeStepExecution(const Kratos::TimeStepEndState& rResultantState)
{
    ++mPos;
}

} // namespace Kratos
