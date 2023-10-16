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
#include "apply_normal_load_table_process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{

ApplyNormalLoadTableProcess::ApplyNormalLoadTableProcess(ModelPart&        rModelPart,
                                                         const Parameters& rProcessSettings)
    : mrModelPart{rModelPart}
{
    if (rProcessSettings["active"][0].GetBool()) {
        auto normal_parameters = Parameters{"{}"};
        normal_parameters.AddValue("model_part_name", rProcessSettings["model_part_name"]);
        normal_parameters.AddValue("variable_name", rProcessSettings["variable_name"]);

        if (rProcessSettings["fluid_pressure_type"].GetString() == "Uniform") {

        }
    }
}

void ApplyNormalLoadTableProcess::ExecuteInitialize()
{
    for (const auto& process : mProcesses) {
        process->ExecuteInitialize();
    }
}

void ApplyNormalLoadTableProcess::ExecuteInitializeSolutionStep()
{
    for (const auto& process : mProcesses) {
        process->ExecuteInitializeSolutionStep();
    }
}

}
