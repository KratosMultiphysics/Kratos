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
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "utilities/variable_utils.h"

// Application includes
#include "deactivate_model_part_operation.h"

namespace Kratos
{

DeactivateModelPartOperation::DeactivateModelPartOperation(
    Model& rModel,
    const Parameters rSettings)
    : Operation()
    , mpModelPart(&rModel.GetModelPart(rSettings["model_part_name"].GetString()))
{}

Operation::Pointer DeactivateModelPartOperation::Create(
    Model &rModel,
    Parameters Parameters) const
{
    return Kratos::make_shared<DeactivateModelPartOperation>(rModel, Parameters);
}

void DeactivateModelPartOperation::Execute()
{
    KRATOS_TRY

    // Deactivate the elements of the model part
    VariableUtils().SetFlag(ACTIVE, false, mpModelPart->Elements());

    // Reset the elements' constitutive law (e.g., excavation)
    block_for_each(mpModelPart->Elements(), [](Element &rElement) {
        rElement.ResetConstitutiveLaw();
    });

    // Deactivate the nodes of the model part
    VariableUtils().SetFlag(ACTIVE, false, mpModelPart->Nodes());

    // Deactivate the conditions of the model part
    VariableUtils().SetFlag(ACTIVE, false, mpModelPart->Conditions());

    KRATOS_CATCH("")
}

}
