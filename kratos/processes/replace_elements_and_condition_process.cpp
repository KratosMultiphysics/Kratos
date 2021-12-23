//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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

namespace Kratos {
namespace {

/// Replace entities in a given container if the entity id is present in a list of ids.
    /*  @param rReferenceEntity New type of entity that will replace old one
     *  @param rEntityContainer Container of elements susceptible to be replaces
     *  @param rSetOfIds Set of entities ids we want to replace
     */
template <class TEntity, class TEntityContainer>
void ReplaceEntities(
    const TEntity& rReferenceEntity,
    TEntityContainer& rEntityContainer,
    std::unordered_set<std::size_t>& rSetOfIds)
{
    IndexPartition<std::size_t>(rEntityContainer.size()).for_each([&](std::size_t Index){
        auto it_entity = rEntityContainer.begin() + Index;
        if (rSetOfIds.find(it_entity->Id()) != rSetOfIds.end()) {
            auto p_new_entity = rReferenceEntity.Create(it_entity->Id(), it_entity->pGetGeometry(), it_entity->pGetProperties());
            // Deep copy data and flags
            p_new_entity->Set(Flags(*it_entity));

            (*it_entity.base()) = p_new_entity;
        }
    });
}

/// Replace elements in a given submodelpart using the elements from the root model part
/// if the element id is present in a given set of ids
    /*  @param rModelPart Model part whose elements we want to replace
     *  @param rRootModelPart Root model part with the replaced elements
     *  @param rSetOfElementsIds Set of elements ids we want to replace
     */
void UpdateElementsInSubModelPart(
    ModelPart& rModelPart,
    ModelPart& rRootModelPart,
    std::unordered_set<std::size_t>& rSetOfElementsIds)
{
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
}


/// Replace conditions in a given submodelpart using the conditions from the root model part
/// if the condition id is present in a given set of ids
    /*  @param rModelPart Model part whose conditions we want to replace
     *  @param rRootModelPart Root model part with the replaced conditions
     *  @param rSetOfConditions Set of conditions ids we want to replace
     */
void UpdateConditionsInSubModelPart(
    ModelPart& rModelPart,
    ModelPart& rRootModelPart,
    std::unordered_set<std::size_t>& rSetOfConditions)
{
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
}

}

/***********************************************************************************/
/***********************************************************************************/

void ReplaceElementsAndConditionsProcess::Execute()
{
    const std::string element_name = mSettings["element_name"].GetString();
    const std::string condition_name = mSettings["condition_name"].GetString();
    ModelPart& root_model_part = mrModelPart.GetRootModelPart();

    if (element_name != "") {
        std::unordered_set<std::size_t> set_element_ids (mrModelPart.NumberOfElements());
        for(auto& ielem : mrModelPart.Elements()) {
            set_element_ids.insert(ielem.Id());
        }
        ReplaceEntities(KratosComponents<Element>::Get(element_name), root_model_part.Elements(), set_element_ids);
        for (auto& r_sub_model_part : root_model_part.SubModelParts()) {
            UpdateElementsInSubModelPart(r_sub_model_part, root_model_part, set_element_ids);
        }
    }

    if (condition_name != "") {
        std::unordered_set<std::size_t> set_conditions_ids (mrModelPart.NumberOfConditions());
        for(auto& icond : mrModelPart.Conditions()) {
            set_conditions_ids.insert(icond.Id());
        }
        ReplaceEntities(KratosComponents<Condition>::Get(condition_name), root_model_part.Conditions(), set_conditions_ids);
        for (auto& r_sub_model_part : root_model_part.SubModelParts()) {
            UpdateConditionsInSubModelPart(r_sub_model_part, root_model_part, set_conditions_ids);
        }
    }

}

}  // namespace Kratos.



