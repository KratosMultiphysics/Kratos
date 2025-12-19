// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

#include "set_absorbing_boundary_parameters_process.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{
using namespace std::string_literals;

SetAbsorbingBoundaryParametersProcess::SetAbsorbingBoundaryParametersProcess(ModelPart& model_part, Parameters rParameters)
    : Process(Flags()), mrModelPart(model_part)
{
    KRATOS_TRY

    // only include validation with c++11 since raw_literals do not exist in c++03
    Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "absorbing_factors": [1.0,1.0],
                "virtual_thickness": 1e10
            }  )");

    // Some values need to be mandatory prescribed since no meaningful default value exist. For
    // this reason try accessing to them So that an error is thrown if they don't exist
    rParameters["model_part_name"];

    // Now validate against defaults -- this also ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);

    // get absorbing factors
    mAbsorbingFactors.resize(2, false);
    mAbsorbingFactors(0) = rParameters["absorbing_factors"][0].GetDouble();
    mAbsorbingFactors(1) = rParameters["absorbing_factors"][1].GetDouble();

    // get virtual thickness
    mVirtualThickness = rParameters["virtual_thickness"].GetDouble();

    KRATOS_CATCH("")
}

void SetAbsorbingBoundaryParametersProcess::ExecuteInitialize()
{
    KRATOS_TRY

    block_for_each(mrModelPart.Conditions(), [&](Condition& rCondition) {
        rCondition.SetValue(ABSORBING_FACTORS, mAbsorbingFactors);
        rCondition.SetValue(VIRTUAL_THICKNESS, mVirtualThickness);
    });

    KRATOS_CATCH("")
}

std::string SetAbsorbingBoundaryParametersProcess::Info() const
{
    return "SetAbsorbingBoundaryParametersProcess"s;
}

} // namespace Kratos