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

#include <boost/format/free_funcs.hpp>

namespace Kratos
{

ResetDisplacementProcess::ResetDisplacementProcess(ModelPart& rModelPart, const Parameters&)
    : mrModelPart(rModelPart)
{
}

void ResetDisplacementProcess::ExecuteInitialize()
{
    block_for_each(mrModelPart.Elements(), [this](Element& rElement) {
        std::vector<Vector> stresses_on_integration_points;
        rElement.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stresses_on_integration_points,
                                              mrModelPart.GetProcessInfo());
        KRATOS_ERROR_IF(stresses_on_integration_points.size() != 1);
        rElement.GetGeometry().SetValue(INITIAL_STRESS_VECTOR, stresses_on_integration_points[0]);
    });

    SetInitialStateProcess<3> set_initial_state_process(mrModelPart);
    set_initial_state_process.ExecuteInitializeSolutionStep();
}

} // namespace Kratos