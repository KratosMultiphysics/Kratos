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

// System includes


// External includes


// Project includes
#include "containers/model.h"
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "utilities/variable_utils.h"

// Application includes
#include "apply_excavation_operation.h"

namespace Kratos
{

ApplyExcavationOperation::ApplyExcavationOperation(
    Model& rModel,
    const Parameters rSettings)
    : Operation()
    , mDeactivateSoilPart(rSettings["deactivate_soil_part"].GetBool())
    , mpModelPart(&rModel.GetModelPart(rSettings["model_part_name"].GetString()))
{}

Operation::Pointer ApplyExcavationOperation::Create(
    Model &rModel,
    Parameters Parameters) const
{
    return Kratos::make_shared<ApplyExcavationOperation>(rModel, Parameters);
}

void ApplyExcavationOperation::Execute()
{
    KRATOS_TRY

    // Activate/deactivate the elements of the model part according to the excavation status
    VariableUtils().SetFlag(ACTIVE, !mDeactivateSoilPart, mpModelPart->Elements());

    // Reset the elements' constitutive law in case of excavation
    if (mDeactivateSoilPart) {
        block_for_each(mpModelPart->Elements(), [](Element &rElement) {
            rElement.ResetConstitutiveLaw();
        });
    } else {
        //TODO: Ask Deltares team why this is needed
        VariableUtils().SetFlag(ACTIVE, true, mpModelPart->Nodes());
    }

    // Activate/deactivate the conditions of the model part according to the excavation status
    VariableUtils().SetFlag(ACTIVE, !mDeactivateSoilPart, mpModelPart->Conditions());

    KRATOS_CATCH("")
}

}
