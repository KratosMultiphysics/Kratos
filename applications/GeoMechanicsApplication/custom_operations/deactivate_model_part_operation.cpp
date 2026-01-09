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
//                   Ruben Zorrilla
//

// Project includes
#include "containers/model.h"
#include "custom_utilities/process_utilities.h"
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "utilities/variable_utils.h"

// Application includes
#include "deactivate_model_part_operation.h"

namespace Kratos
{
DeactivateModelPartOperation::DeactivateModelPartOperation(Model& rModel, const Parameters& rSettings)
{
    mrModelParts = ProcessUtilities::GetModelPartsFromSettings(rModel, rSettings, "DeactivateModelPartOperation");
}

Operation::Pointer DeactivateModelPartOperation::Create(Model& rModel, Parameters Settings) const
{
    return Kratos::make_shared<DeactivateModelPartOperation>(rModel, Settings);
}

void DeactivateModelPartOperation::Execute()
{
    KRATOS_TRY

    for (const auto& r_model_part : mrModelParts) {
        // Deactivate the elements of the model part
        VariableUtils().SetFlag(ACTIVE, false, r_model_part.get().Elements());

        // Reset the elements' constitutive law (e.g., excavation)
        block_for_each(r_model_part.get().Elements(),
                       [](Element& rElement) { rElement.ResetConstitutiveLaw(); });

        // Deactivate the nodes of the model part
        VariableUtils().SetFlag(ACTIVE, false, r_model_part.get().Nodes());

        // Deactivate the conditions of the model part
        VariableUtils().SetFlag(ACTIVE, false, r_model_part.get().Conditions());
    }

    KRATOS_CATCH("")
}

} // namespace Kratos
