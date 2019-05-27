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

// External includes

// Project includes
#include "processes/replace_elements_and_condition_process.h"

namespace Kratos {
namespace {

template <class TEntity, class TEntityContainer>
void ReplaceEntities(
    const TEntity& rReferenceEntity,
    TEntityContainer& rEntityContainer,
    std::unordered_set<std::size_t>& rSetOfIds)
{
    #pragma omp parallel for
    for (int i=0; i<static_cast<int>(rEntityContainer.size()); ++i) {
        auto it_entity = rEntityContainer.begin() + i;
        if (rSetOfIds.find(it_entity->Id()) != rSetOfIds.end()) {
            auto p_new_entity = rReferenceEntity.Create(it_entity->Id(), it_entity->pGetGeometry(), it_entity->pGetProperties());
            // Deep copy data and flags
            p_new_entity->Data() = it_entity->Data();
            p_new_entity->Set(Flags(*it_entity));

            (*it_entity.base()) = p_new_entity;
        }        
    }
}

void UpdateElementsInSubModelPart(
    ModelPart& rModelPart,
    ModelPart& rRootModelPart,
    std::unordered_set<std::size_t>& rSetOfElementsIds)
{
    #pragma omp parallel for
    for(int i=0; i< static_cast<int>(rModelPart.Elements().size()); i++) {
        auto it_elem = rModelPart.ElementsBegin() + i;
        if (rSetOfElementsIds.find(it_elem->Id()) != rSetOfElementsIds.end()) {
            (*it_elem.base()) = rRootModelPart.Elements()(it_elem->Id());
        }
    }
    // Change the submodelparts
    for (auto& r_sub_model_part : rModelPart.SubModelParts()) {
        UpdateElementsInSubModelPart(r_sub_model_part, rRootModelPart, rSetOfElementsIds);
    }
}

void UpdateConditionsInSubModelPart(
    ModelPart& rModelPart,
    ModelPart& rRootModelPart,
    std::unordered_set<std::size_t>& rSetOfConditions)
{
    #pragma omp parallel for
    for(int i=0; i< static_cast<int>(rModelPart.Conditions().size()); i++) {
        auto it_cond = rModelPart.ConditionsBegin() + i;
        if (rSetOfConditions.find(it_cond->Id()) != rSetOfConditions.end()) {
            (*it_cond.base()) = rRootModelPart.Conditions()(it_cond->Id());
        }
    }
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



