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
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

ApplyExcavationProcess::ApplyExcavationProcess(Model& rModel, const Parameters& rProcessSettings)
    : mDeactivateSoilPart{rProcessSettings["deactivate_soil_part"].GetBool()}
{
    KRATOS_ERROR_IF_NOT(rProcessSettings.Has("model_part_name") ||
                        rProcessSettings.Has("model_part_name_list"))
        << "Please specify 'model_part_name' or 'model_part_name_list' for "
        << ApplyExcavationProcess::Info();

    KRATOS_ERROR_IF(rProcessSettings.Has("model_part_name") &&
                    rProcessSettings.Has("model_part_name_list"))
        << "The parameters 'model_part_name' and 'model_part_name_list' are mutually exclusive for "
        << ApplyExcavationProcess::Info();

    const auto model_part_names = rProcessSettings.Has("model_part_name_list")
                                      ? rProcessSettings["model_part_name_list"].GetStringArray()
                                      : std::vector{rProcessSettings["model_part_name"].GetString()};
    KRATOS_ERROR_IF(model_part_names.empty()) << "The parameters 'model_part_name_list' needs "
                                                 "to contain at least one model part name for "
                                              << ApplyExcavationProcess::Info();
    std::ranges::transform(
        model_part_names, std::back_inserter(mrModelParts),
        [&rModel](const auto& rName) -> ModelPart& { return rModel.GetModelPart(rName); });
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

ApplyExcavationProcess::~ApplyExcavationProcess() = default;

} // namespace Kratos
