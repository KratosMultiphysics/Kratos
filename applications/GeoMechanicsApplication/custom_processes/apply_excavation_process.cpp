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
//                   Vahid Galavi,
//                   Marjan Fathian
//
#include "apply_excavation_process.h"
#include "containers/model.h"
#include "custom_utilities/process_utilities.h"
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"

#include <string>

namespace Kratos
{
using namespace std::string_literals;

ApplyExcavationProcess::ApplyExcavationProcess(Model& rModel, const Parameters& rProcessSettings)
    : mDeactivateSoilPart{rProcessSettings["deactivate_soil_part"].GetBool()}
{
    mrModelParts = ProcessUtilities::GetModelPartsFromSettings(rModel, rProcessSettings,
                                                               ApplyExcavationProcess::Info());
}

void ApplyExcavationProcess::ExecuteInitialize()
{
    KRATOS_TRY
    for (const auto& r_model_part : mrModelParts) {
        VariableUtils{}.SetFlag(ACTIVE, !mDeactivateSoilPart, r_model_part.get().Elements());
    }

    if (mDeactivateSoilPart) {
        for (const auto& r_model_part : mrModelParts) {
            block_for_each(r_model_part.get().Elements(),
                           [](Element& rElement) { rElement.ResetConstitutiveLaw(); });
        }
    } else {
        for (const auto& r_model_part : mrModelParts) {
            VariableUtils{}.SetFlag(ACTIVE, true, r_model_part.get().Nodes());
        }
    }

    for (const auto& r_model_part : mrModelParts) {
        VariableUtils{}.SetFlag(ACTIVE, !mDeactivateSoilPart, r_model_part.get().Conditions());
    }
    KRATOS_CATCH("")
}

std::string ApplyExcavationProcess::Info() const { return "ApplyExcavationProcess"s; }

} // namespace Kratos
