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
//                   Marjan Fathian
//                   Richard Faasse
//                   Anne van de Graaf
//

#include "time_step_end_state.h"

namespace Kratos
{
bool TimeStepEndState::Converged() const
{
    return convergence_state == ConvergenceState::converged;
}

bool TimeStepEndState::NonConverged() const
{
    return convergence_state == ConvergenceState::non_converged;
}

} // namespace Kratos
