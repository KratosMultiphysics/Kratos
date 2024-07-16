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
#include "activate_model_part_operation.h"

namespace Kratos
{

ActivateModelPartOperation::ActivateModelPartOperation(
    Model& rModel,
    const Parameters rSettings)
    : Operation()
    , mpModelPart(&rModel.GetModelPart(rSettings["model_part_name"].GetString()))
{}

Operation::Pointer ActivateModelPartOperation::Create(
    Model &rModel,
    Parameters Parameters) const
{
    return Kratos::make_shared<ActivateModelPartOperation>(rModel, Parameters);
}

void ActivateModelPartOperation::Execute()
{
    KRATOS_TRY

    // Activate the elements of the model part
    VariableUtils().SetFlag(ACTIVE, true, mpModelPart->Elements());

    // Activate the nodes of the model part
    VariableUtils().SetFlag(ACTIVE, true, mpModelPart->Nodes());

    // Activate the conditions of the model part
    VariableUtils().SetFlag(ACTIVE, true, mpModelPart->Conditions());

    KRATOS_CATCH("")
}

}
