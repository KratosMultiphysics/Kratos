//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/kratos_flags.h"
#include "processes/entity_erase_process.h"
#include "utilities/variable_utils.h"
#include "utilities/auxiliar_model_part_utilities.h"

namespace Kratos
{
template<class TEntity>
EntitiesEraseProcess<TEntity>::EntitiesEraseProcess(
    ModelPart& rModelPart,
    const bool RemoveFromAllLevels,
    const bool AssignFlag,
    const bool RemoveBelongings
    )
    : mrModelPart(rModelPart),
     mRemoveFromAllLevels(RemoveFromAllLevels),
     mAssignFlag(AssignFlag),
     mRemoveBelongings(RemoveBelongings)
{
    KRATOS_TRY

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
EntitiesEraseProcess<TEntity>::EntitiesEraseProcess(
    ModelPart& rModelPart,
    Parameters ThisParameters
    )
    : mrModelPart(rModelPart)
{
    KRATOS_TRY

    // The default parameters
    Parameters default_parameters = GetDefaultParameters();
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Assign parameters
    mRemoveFromAllLevels = ThisParameters["remove_from_all_levels"].GetBool();
    mAssignFlag = ThisParameters["assign_flag"].GetBool();
    mRemoveBelongings = ThisParameters["remove_belongings"].GetBool();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void EntitiesEraseProcess<Node<3>>::Execute()
{
    KRATOS_TRY;

    // If we assign flags
    if (mAssignFlag) VariableUtils().SetFlag(TO_ERASE, true, mrModelPart.Nodes());

    // Remove nodes
    if (mRemoveFromAllLevels) {
        mrModelPart.RemoveNodesFromAllLevels(TO_ERASE);
    } else {
        mrModelPart.RemoveNodes(TO_ERASE);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void EntitiesEraseProcess<Element>::Execute()
{
    KRATOS_TRY;

    // If we assign flags
    if (mAssignFlag) VariableUtils().SetFlag(TO_ERASE, true, mrModelPart.Elements());

    // Remove elements
    if (mRemoveBelongings) {
        if (mRemoveFromAllLevels) {
            AuxiliarModelPartUtilities(mrModelPart).RemoveElementsAndBelongings(TO_ERASE);
        } else {
            AuxiliarModelPartUtilities(mrModelPart).RemoveElementsAndBelongingsFromAllLevels(TO_ERASE);
        }
    } else {
        if (mRemoveFromAllLevels) {
            mrModelPart.RemoveElementsFromAllLevels(TO_ERASE);
        } else {
            mrModelPart.RemoveElements(TO_ERASE);
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void EntitiesEraseProcess<Condition>::Execute()
{
    KRATOS_TRY;

    // If we assign flags
    if (mAssignFlag) VariableUtils().SetFlag(TO_ERASE, true, mrModelPart.Conditions());

    // Remove conditions
    if (mRemoveBelongings) {
        if (mRemoveFromAllLevels) {
            AuxiliarModelPartUtilities(mrModelPart).RemoveConditionsAndBelongings(TO_ERASE);
        } else {
            AuxiliarModelPartUtilities(mrModelPart).RemoveConditionsAndBelongingsFromAllLevels(TO_ERASE);
        }
    } else {
        if (mRemoveFromAllLevels) {
            mrModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
        } else {
            mrModelPart.RemoveConditions(TO_ERASE);
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
Parameters EntitiesEraseProcess<TEntity>::GetDefaultParameters()
{
    Parameters default_parameters = Parameters(R"(
    {
        "remove_from_all_levels" : false,
        "assign_flag"            : false,
        "remove_belongings"      : false
    })" );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class EntitiesEraseProcess<Node<3>>;
template class EntitiesEraseProcess<Element>;
template class EntitiesEraseProcess<Condition>;

}




