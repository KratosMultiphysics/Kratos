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
void ReplaceEntities(const TEntity& rReferenceEntity,
                     TEntityContainer& rEntityContainer)
{
    #pragma omp parallel for
    for (int i=0; i<static_cast<int>(rEntityContainer.size()); ++i) {
        auto it_entity = rEntityContainer.begin() + i;

        auto p_new_entity = rReferenceEntity.Create(it_entity->Id(), it_entity->pGetGeometry(), it_entity->pGetProperties());

        // Deep copy data and flags
        p_new_entity->Data() = it_entity->Data();
        p_new_entity->Set(Flags(*it_entity));

        (*it_entity.base()) = p_new_entity;
    }
}

void UpdateSubModelPart(
    ModelPart& rModelPart,
    ModelPart& rRootModelPart,
    const bool UpdateElements,
    const bool UpdateConditions
    )
{
    // Change the model part itself
    if (UpdateElements) {
        #pragma omp parallel for
        for(int i=0; i< static_cast<int>(rModelPart.Elements().size()); i++) {
            auto it_elem = rModelPart.ElementsBegin() + i;

            (*it_elem.base()) = rRootModelPart.Elements()(it_elem->Id());
        }
    }

    if (UpdateConditions) {
        #pragma omp parallel for
        for(int i=0; i< static_cast<int>(rModelPart.Conditions().size()); i++) {
            auto it_cond = rModelPart.ConditionsBegin() + i;

            (*it_cond.base()) = rRootModelPart.Conditions()(it_cond->Id());
        }
    }

    // Change the submodelparts
    for (auto& r_sub_model_part : rModelPart.SubModelParts()) {
        UpdateSubModelPart(r_sub_model_part, rRootModelPart, UpdateElements, UpdateConditions);
    }
}

}

/***********************************************************************************/
/***********************************************************************************/

void ReplaceElementsAndConditionsProcess::Execute()
{
    ModelPart& r_root_model_part = mrModelPart.GetRootModelPart();

    const std::string element_name = mSettings["element_name"].GetString();
    const std::string condition_name = mSettings["condition_name"].GetString();

    if (element_name != "") {
        ReplaceEntities(KratosComponents<Element>::Get(element_name), r_root_model_part.Elements());
    }

    if (condition_name != "") {
        ReplaceEntities(KratosComponents<Condition>::Get(condition_name), r_root_model_part.Conditions());
    }

    // Change the submodelparts
    for (auto& r_sub_model_part : r_root_model_part.SubModelParts()) {
        UpdateSubModelPart(r_sub_model_part, r_root_model_part,
                           element_name != "", condition_name != "");
    }
}

}  // namespace Kratos.



