// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Lorenzo Gracia,
//                   Aron Noordam,
//                   Vahid Galavi
//                   Marjan Fathian
//
#include "apply_excavation_process.hpp"
#include "includes/model_part.h"
#include "includes/element.h"
#include "utilities/variable_utils.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{

ApplyExcavationProcess::ApplyExcavationProcess(ModelPart&        rModelPart,
                                               const Parameters& rProcessSettings)
           : Process(Flags()),
             mrModelPart{rModelPart}
{
    Parameters params = Parameters("{}");
    params.AddValue("model_part_name", rProcessSettings["model_part_name"]);
    params.AddValue("variable_name", rProcessSettings["variable_name"]);

    if (rProcessSettings.Has("deactivate_soil_part")) {
        params.AddValue("deactivate_soil_part", rProcessSettings["deactivate_soil_part"]);
        KRATOS_TRY
        mDeactivateSoilPart = rProcessSettings["deactivate_soil_part"].GetBool();
        KRATOS_CATCH("")
    }
}

void ApplyExcavationProcess::ExecuteInitialize()
{
    KRATOS_TRY
    VariableUtils{}.SetFlag(ACTIVE, !mDeactivateSoilPart, mrModelPart.Elements());

    if (mDeactivateSoilPart) {
        block_for_each(mrModelPart.Elements(), [](Element &rElement) {
            rElement.ResetConstitutiveLaw();
        });
    } else {
        VariableUtils{}.SetFlag(ACTIVE, true, mrModelPart.Nodes());
    }

    VariableUtils{}.SetFlag(ACTIVE, !mDeactivateSoilPart, mrModelPart.Conditions());
    KRATOS_CATCH("")
}

ApplyExcavationProcess::~ApplyExcavationProcess() = default;

}
