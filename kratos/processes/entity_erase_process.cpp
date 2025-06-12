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

KRATOS_CREATE_LOCAL_FLAG(EntitiesEraseProcessFlags, REMOVE_FROM_ALL_LEVELS, 0);
KRATOS_CREATE_LOCAL_FLAG(EntitiesEraseProcessFlags, ERASE_ALL_ENTITIES,     1);

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
EntitiesEraseProcess<TEntity>::EntitiesEraseProcess(
    ModelPart& rModelPart,
    const Flags Options
    )
    : mrModelPart(rModelPart),
     mrOptions(Options)
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
    mrOptions.Set(EntitiesEraseProcessFlags::REMOVE_FROM_ALL_LEVELS, ThisParameters["remove_from_all_levels"].GetBool());
    mrOptions.Set(EntitiesEraseProcessFlags::ERASE_ALL_ENTITIES, ThisParameters["remove_all_entities"].GetBool());

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void EntitiesEraseProcess<Node>::Execute()
{
    KRATOS_TRY;

    // If we assign flags
    if (mrOptions.Is(EntitiesEraseProcessFlags::ERASE_ALL_ENTITIES)) VariableUtils().SetFlag(TO_ERASE, true, mrModelPart.Nodes());

    // Remove nodes
    if (mrOptions.Is(EntitiesEraseProcessFlags::REMOVE_FROM_ALL_LEVELS)) {
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
    if (mrOptions.Is(EntitiesEraseProcessFlags::ERASE_ALL_ENTITIES)) VariableUtils().SetFlag(TO_ERASE, true, mrModelPart.Elements());

    // Remove elements
    if (mrOptions.Is(EntitiesEraseProcessFlags::REMOVE_FROM_ALL_LEVELS)) {
        mrModelPart.RemoveElementsFromAllLevels(TO_ERASE);
    } else {
        mrModelPart.RemoveElements(TO_ERASE);
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
    if (mrOptions.Is(EntitiesEraseProcessFlags::ERASE_ALL_ENTITIES)) VariableUtils().SetFlag(TO_ERASE, true, mrModelPart.Conditions());

    // Remove conditions
    if (mrOptions.Is(EntitiesEraseProcessFlags::REMOVE_FROM_ALL_LEVELS)) {
        mrModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
    } else {
        mrModelPart.RemoveConditions(TO_ERASE);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void EntitiesEraseProcess<MasterSlaveConstraint>::Execute()
{
    KRATOS_TRY;

    // If we assign flags
    if (mrOptions.Is(EntitiesEraseProcessFlags::ERASE_ALL_ENTITIES)) VariableUtils().SetFlag(TO_ERASE, true, mrModelPart.MasterSlaveConstraints());

    // Remove conditions
    if (mrOptions.Is(EntitiesEraseProcessFlags::REMOVE_FROM_ALL_LEVELS)) {
        mrModelPart.RemoveMasterSlaveConstraintsFromAllLevels(TO_ERASE);
    } else {
        mrModelPart.RemoveMasterSlaveConstraints(TO_ERASE);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
const Parameters EntitiesEraseProcess<TEntity>::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "remove_from_all_levels" : false,
        "remove_all_entities"    : false
    })" );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class EntitiesEraseProcess<Node>;
template class EntitiesEraseProcess<Element>;
template class EntitiesEraseProcess<Condition>;
template class EntitiesEraseProcess<MasterSlaveConstraint>;
}
