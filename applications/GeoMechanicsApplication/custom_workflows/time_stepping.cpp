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

#include "time_stepping.h"


namespace Kratos
{

TimeStepExecuter::ConvergenceState TimeStepExecuter::Run()
{
    return ConvergenceState::non_converged;
}

}
