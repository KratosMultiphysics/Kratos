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
    std::unordered_set<std::size_t>& rSetOfIds
    )
{
    const auto& r_reference_geometry = rReferenceEntity.GetGeometry();
    const auto& r_reference_geometry_type = r_reference_geometry.GetGeometryType();
    IndexPartition<std::size_t>(rEntityContainer.size()).for_each([&](std::size_t Index){
        auto it_entity = rEntityContainer.begin() + Index;
        if (rSetOfIds.find(it_entity->Id()) != rSetOfIds.end()) {
            auto p_geometry = it_entity->pGetGeometry();
            KRATOS_DEBUG_ERROR_IF_NOT(p_geometry->GetGeometryType() == r_reference_geometry_type) << "Trying to replace an element with a different geometry type. Reference entity " << r_reference_geometry.Info() << " vs  " << p_geometry->Info() << "\n Entity info: " << rReferenceEntity.Info() << std::endl;
            auto p_new_entity = rReferenceEntity.Create(it_entity->Id(), p_geometry, it_entity->pGetProperties());
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
    std::unordered_set<std::size_t>& rSetOfElementsIds
    )
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
    std::unordered_set<std::size_t>& rSetOfConditions
    )
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
    const std::string& r_element_name = mSettings["element_name"].GetString();
    const std::string& r_condition_name = mSettings["condition_name"].GetString();
    ModelPart& r_root_model_part = mrModelPart.GetRootModelPart();

    if (r_element_name != "") {
        std::unordered_set<std::size_t> set_element_ids (mrModelPart.NumberOfElements());
        for(const auto& r_elem : mrModelPart.Elements()) {
            set_element_ids.insert(r_elem.Id());
        }
        ReplaceEntities(KratosComponents<Element>::Get(r_element_name), r_root_model_part.Elements(), set_element_ids);
        for (auto& r_sub_model_part : r_root_model_part.SubModelParts()) {
            UpdateElementsInSubModelPart(r_sub_model_part, r_root_model_part, set_element_ids);
        }
    }

    if (r_condition_name != "") {
        std::unordered_set<std::size_t> set_conditions_ids (mrModelPart.NumberOfConditions());
        for(const auto& r_cond : mrModelPart.Conditions()) {
            set_conditions_ids.insert(r_cond.Id());
        }
        ReplaceEntities(KratosComponents<Condition>::Get(r_condition_name), r_root_model_part.Conditions(), set_conditions_ids);
        for (auto& r_sub_model_part : r_root_model_part.SubModelParts()) {
            UpdateConditionsInSubModelPart(r_sub_model_part, r_root_model_part, set_conditions_ids);
        }
    }
}

}  // namespace Kratos.