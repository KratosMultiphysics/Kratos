//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//                   Philipp Bucher
//

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "processes/replace_elements_and_condition_process.h"
#include "utilities/parallel_utilities.h"
#include "utilities/string_utilities.h"
#include "utilities/geometry_utilities.h"

namespace Kratos {
namespace {

/**
 * @brief Replace entities in a given container if the entity id is present in a list of ids.
 * @param pEntitityIdentifier Entity identifier
 * @param rEntityContainer Container of elements susceptible to be replaces
 * @param rSetOfIds Set of entities ids we want to replace
 */
template <class TContainer>
void ReplaceEntities(
    typename EntitiesUtilities::EntitityIdentifier<typename TContainer::value_type>& rEntitityIdentifier,
    TContainer& rEntityContainer,
    const std::unordered_set<std::size_t>& rSetOfIds
    )
{
    KRATOS_TRY;

    IndexPartition<std::size_t>(rEntityContainer.size()).for_each([&rEntitityIdentifier, &rEntityContainer, &rSetOfIds](std::size_t Index){
        auto it_entity = rEntityContainer.begin() + Index;
        const std::size_t id = it_entity->Id();
        if (rSetOfIds.find(id) != rSetOfIds.end()) {
            auto p_geometry = it_entity->pGetGeometry();
            const auto& r_reference_entity = rEntitityIdentifier.GetPrototypeEntity(p_geometry);
            auto p_new_entity = r_reference_entity.Create(id, p_geometry, it_entity->pGetProperties());
            // Deep copy data and flags
            p_new_entity->Set(Flags(*it_entity));

            // Reassign the pointer
            (*it_entity.base()) = p_new_entity;
        }
    });

    KRATOS_CATCH("");
}

/**
 * @brief Replace elements in a given submodelpart using the elements from the root model part if the element id is present in a given set of ids
 * @param rModelPart Model part whose elements we want to replace
 * @param rRootModelPart Root model part with the replaced elements
 * @param rSetOfElementsIds Set of elements ids we want to replace
 */
void UpdateElementsInSubModelPart(
    ModelPart& rModelPart,
    ModelPart& rRootModelPart,
    std::unordered_set<std::size_t>& rSetOfElementsIds
    )
{
    KRATOS_TRY;

    if(!rRootModelPart.Elements().IsSorted()) {
        rRootModelPart.Elements().Sort();
    }
    IndexPartition<std::size_t>(rModelPart.Elements().size()).for_each([&](std::size_t Index){
        auto it_elem = rModelPart.ElementsBegin() + Index;
        if (rSetOfElementsIds.find(it_elem->Id()) != rSetOfElementsIds.end()) {
            (*it_elem.base()) = rRootModelPart.Elements()(it_elem->Id());
        }
    });

    // Change the submodelparts
    for (auto& r_sub_model_part : rModelPart.SubModelParts()) {
        UpdateElementsInSubModelPart(r_sub_model_part, rRootModelPart, rSetOfElementsIds);
    }

    KRATOS_CATCH("");
}

/**
 * @brief Replace conditions in a given submodelpart using the conditions from the root model part if the condition id is present in a given set of ids
 * @param rModelPart Model part whose conditions we want to replace
 * @param rRootModelPart Root model part with the replaced conditions
 * @param rSetOfConditions Set of conditions ids we want to replace
 */
void UpdateConditionsInSubModelPart(
    ModelPart& rModelPart,
    ModelPart& rRootModelPart,
    std::unordered_set<std::size_t>& rSetOfConditions
    )
{
    KRATOS_TRY;

    if(!rRootModelPart.Conditions().IsSorted()) {
        rRootModelPart.Conditions().Sort();
    }
    IndexPartition<std::size_t>(rModelPart.Conditions().size()).for_each([&](std::size_t Index){
        auto it_cond = rModelPart.ConditionsBegin() + Index;
        if (rSetOfConditions.find(it_cond->Id()) != rSetOfConditions.end()) {
            (*it_cond.base()) = rRootModelPart.Conditions()(it_cond->Id());
        }
    });

    // Change the submodelparts
    for (auto& r_sub_model_part : rModelPart.SubModelParts()) {
        UpdateConditionsInSubModelPart(r_sub_model_part, rRootModelPart, rSetOfConditions);
    }

    KRATOS_CATCH("");
}

}

/***********************************************************************************/
/***********************************************************************************/

void ReplaceElementsAndConditionsProcess::InitializeMemberVariables()
{
    KRATOS_TRY;

    // We can provide a list of replacements instead of a single one

    /* Elements */
    if (mSettings.Has("element_name")) {
        const std::string& r_element_name = mSettings["element_name"].GetString();
        mElementIdentifier = EntitiesUtilities::EntitityIdentifier<Element>(r_element_name);
    }

    /* Conditions */
    if (mSettings.Has("condition_name")) {
        const std::string& r_condition_name = mSettings["condition_name"].GetString();
        mConditionIdentifier = EntitiesUtilities::EntitityIdentifier<Condition>(r_condition_name);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void ReplaceElementsAndConditionsProcess::Execute()
{
    KRATOS_TRY;

    ModelPart& r_root_model_part = mrModelPart.GetRootModelPart();

    /* Elements */
    if (mElementIdentifier.IsInitialized()) {
        std::unordered_set<std::size_t> set_element_ids (mrModelPart.NumberOfElements());
        for(const auto& r_elem : mrModelPart.Elements()) {
            set_element_ids.insert(r_elem.Id());
        }
        ReplaceEntities(mElementIdentifier, r_root_model_part.Elements(), set_element_ids);
        for (auto& r_sub_model_part : r_root_model_part.SubModelParts()) {
            UpdateElementsInSubModelPart(r_sub_model_part, r_root_model_part, set_element_ids);
        }
    }

    /* Conditions */
    if (mConditionIdentifier.IsInitialized()) {
        std::unordered_set<std::size_t> set_conditions_ids (mrModelPart.NumberOfConditions());
        for(const auto& r_cond : mrModelPart.Conditions()) {
            set_conditions_ids.insert(r_cond.Id());
        }
        ReplaceEntities(mConditionIdentifier, r_root_model_part.Conditions(), set_conditions_ids);
        for (auto& r_sub_model_part : r_root_model_part.SubModelParts()) {
            UpdateConditionsInSubModelPart(r_sub_model_part, r_root_model_part, set_conditions_ids);
        }
    }

    KRATOS_CATCH("");
}

}  // namespace Kratos.