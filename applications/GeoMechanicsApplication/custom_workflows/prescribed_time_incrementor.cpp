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


namespace Kratos
{

PrescribedTimeIncrementor::PrescribedTimeIncrementor(const std::vector<double>& increments)
    : mIncrements(increments)
{}

bool PrescribedTimeIncrementor::WantNextStep(const Kratos::TimeStepEndState& rPreviousState) const
{
    return !mIncrements.empty();
}

bool PrescribedTimeIncrementor::WantRetryStep(std::size_t             CycleNumber,
                                              const TimeStepEndState& rPreviousState) const
{
    return false;
}

double PrescribedTimeIncrementor::GetIncrement() const
{
    return 0.0;
}

void PrescribedTimeIncrementor::PostTimeStepExecution(const Kratos::TimeStepEndState& rResultantState)
{
}

}
