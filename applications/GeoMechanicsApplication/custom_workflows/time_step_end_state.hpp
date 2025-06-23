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

#pragma once

#include <cstddef>

namespace Kratos
{

struct TimeStepEndState {
    enum class ConvergenceState { converged, non_converged };

    double           time              = 0.0;
    ConvergenceState convergence_state = ConvergenceState::non_converged;
    std::size_t      num_of_cycles     = 0;
    std::size_t      num_of_iterations = 0;

    [[nodiscard]] bool Converged() const
    {
        return convergence_state == ConvergenceState::converged;
    }

    [[nodiscard]] bool NonConverged() const
    {
        return convergence_state == ConvergenceState::non_converged;
    }
};

} // namespace Kratos
