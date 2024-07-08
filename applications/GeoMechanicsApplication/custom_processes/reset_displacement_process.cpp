// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//
#include "reset_displacement_process.h"
#include "includes/mat_variables.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "processes/set_initial_state_process.h"

namespace Kratos
{

ResetDisplacementProcess::ResetDisplacementProcess(ModelPart& rModelPart, const Parameters&)
    : mrModelPart(rModelPart)
{
}

void ResetDisplacementProcess::ExecuteInitializeSolutionStep()
{
    SetInitialStateProcess<3> set_initial_state_process(mrModelPart);
    set_initial_state_process.ExecuteInitializeSolutionStep();
}

void ResetDisplacementProcess::ExecuteFinalize()
{
    block_for_each(mrModelPart.Elements(), [this](Element& rElement) {
        std::vector<Vector> test;
        rElement.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, test, mrModelPart.GetProcessInfo());
        rElement.GetGeometry().SetValue(INITIAL_STRESS_VECTOR, test[0]);
    });
}

} // namespace Kratos
