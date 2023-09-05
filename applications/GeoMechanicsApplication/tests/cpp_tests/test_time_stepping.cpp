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

#include "testing/testing.h"
#include "custom_workflows/time_stepping.h"

using namespace Kratos;


namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TimeStepNeverConverges, KratosGeoMechanicsFastSuite)
{
    TimeStepExecuter executer;
    KRATOS_EXPECT_EQ(TimeStepExecuter::ConvergenceState::non_converged, executer.Run());
}

}
