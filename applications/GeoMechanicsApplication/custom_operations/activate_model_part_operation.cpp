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
#include "activate_model_part_operation.h"

namespace Kratos
{
ActivateModelPartOperation::ActivateModelPartOperation() : Operation() {}

ActivateModelPartOperation::ActivateModelPartOperation(Model& rModel, const Parameters& rSettings)
{
    mrModelParts = ProcessUtilities::GetModelPartsFromSettings(rModel, rSettings, "ActivateModelPartOperation");
}

Operation::Pointer ActivateModelPartOperation::Create(Model& rModel, Parameters Settings) const
{
    return Kratos::make_shared<ActivateModelPartOperation>(rModel, Settings);
}

void ActivateModelPartOperation::Execute()
{
    KRATOS_TRY

    for (const auto& r_model_part : mrModelParts) {
        // Activate the elements of the model part
        VariableUtils().SetFlag(ACTIVE, true, r_model_part.get().Elements());

        // Activate the nodes of the model part
        VariableUtils().SetFlag(ACTIVE, true, r_model_part.get().Nodes());

        // Activate the conditions of the model part
        VariableUtils().SetFlag(ACTIVE, true, r_model_part.get().Conditions());
    }

    KRATOS_CATCH("")
}
} // namespace Kratos
