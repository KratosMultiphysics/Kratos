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

/**
 * @brief Replace entities in a given container if the entity id is present in a list of ids.
 * @param rReferenceEntity New type of entity that will replace old one
 * @param rEntityContainer Container of elements susceptible to be replaces
 * @param rSetOfIds Set of entities ids we want to replace
 */
template <class TEntity>
void ReplaceEntities(
    const TEntity& rReferenceEntity,
    PointerVectorSet<TEntity, IndexedObject, std::less<typename IndexedObject::result_type>, std::equal_to<typename IndexedObject::result_type>, typename TEntity::Pointer, std::vector< typename TEntity::Pointer>>& rEntityContainer,
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

/**
 * @brief Replace entities in a given container if the entity id is present in a list of ids.
 * @param rReferenceEntity New type of entity that will replace old one
 * @param rEntityContainer Container of elements susceptible to be replaces
 * @param rSetOfIds Set of entities ids we want to replace
 */
template <class TEntity>
void ReplaceEntities(
    const Parameters ListReferenceEntity,
    const std::unordered_map<GeometryData::KratosGeometryType, std::string>& rGeometryTypesToStrings,
    PointerVectorSet<TEntity, IndexedObject, std::less<typename IndexedObject::result_type>, std::equal_to<typename IndexedObject::result_type>, typename TEntity::Pointer, std::vector< typename TEntity::Pointer>>& rEntityContainer,
    std::unordered_set<std::size_t>& rSetOfIds
    )
{
    typename TEntity::Pointer p_reference_entity = nullptr;
    GeometryData::KratosGeometryType current_geometry_type = GeometryData::KratosGeometryType::Kratos_generic_type;
    GeometryData::KratosGeometryType reference_geometry_type = GeometryData::KratosGeometryType::Kratos_generic_type;
    IndexPartition<std::size_t>(rEntityContainer.size()).for_each([&](std::size_t Index){
        auto it_entity = rEntityContainer.begin() + Index;
        if (rSetOfIds.find(it_entity->Id()) != rSetOfIds.end()) {
            auto p_geometry = it_entity->pGetGeometry();
            const auto& r_geometry_type = p_geometry->GetGeometryType();
            // Checking if geometry type is the same
            if (r_geometry_type != current_geometry_type) {
                auto it_find = rGeometryTypesToStrings.find(r_geometry_type);
                KRATOS_ERROR_IF(it_find == rGeometryTypesToStrings.end()) << "Trying to replace an element with a different geometry type. No compatible geometry type: " << static_cast<int>(r_geometry_type) << std::endl;
                const std::string& r_type = it_find->second;
                KRATOS_ERROR_IF_NOT(ListReferenceEntity.Has(r_type)) << "Trying to replace an element with a different geometry type. No reference entity found for geometry type: " << r_type << std::endl;
                const auto& r_reference_entity = KratosComponents<TEntity>::Get(ListReferenceEntity[r_type].GetString());
                p_reference_entity = r_reference_entity.Create(it_entity->Id(), p_geometry, it_entity->pGetProperties());;
                current_geometry_type = r_geometry_type;
                reference_geometry_type = r_reference_entity.GetGeometry().GetGeometryType();
            }
            KRATOS_DEBUG_ERROR_IF_NOT(r_geometry_type == reference_geometry_type) << "Trying to replace an element with a different geometry type. Reference entity " << p_reference_entity->GetGeometry().Info() << " vs  " << p_geometry->Info() << "\n Entity info: " << p_reference_entity->Info() << std::endl;
            auto p_new_entity = p_reference_entity->Create(it_entity->Id(), p_geometry, it_entity->pGetProperties());
            // Deep copy data and flags
            p_new_entity->Set(Flags(*it_entity));

            (*it_entity.base()) = p_new_entity;
        }
    });
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

void ReplaceElementsAndConditionsProcess::InitializeMemberVariables()
{
    // We can provide a list of replacements instead of a single one

    /* Elements */
    if (mSettings["element_name"].IsString()) {
        const std::string& r_element_name = mSettings["element_name"].GetString();
        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them so that an error is thrown if they don't exist
        KRATOS_ERROR_IF(r_element_name != "" && !KratosComponents<Element>::Has(r_element_name)) << "Element name not found in KratosComponents< Element > -- name is " << r_element_name << std::endl;
        mOnlyOneElementOrCondition[0] = true;
    } else if (mSettings["element_name"].IsSubParameter()) {
        for (auto& r_sub_parameter : mSettings["element_name"]) {
            const std::string& r_element_name = r_sub_parameter.GetString();
            KRATOS_ERROR_IF(r_element_name != "" && !KratosComponents<Element>::Has(r_element_name)) << "Element name not found in KratosComponents< Element > -- name is " << r_element_name << std::endl;
        }
        mOnlyOneElementOrCondition[0] = false;
    } else {
        KRATOS_ERROR << "The element_name should be either a string or a list of strings" << std::endl;
    }

    /* Conditions */
    if (mSettings["condition_name"].IsString()) {
        const std::string& r_condition_name = mSettings["condition_name"].GetString();
        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them so that an error is thrown if they don't exist
        KRATOS_ERROR_IF(r_condition_name != "" && !KratosComponents<Condition>::Has(r_condition_name)) << "Element name not found in KratosComponents< Condition > -- name is " << r_condition_name << std::endl;
        mOnlyOneElementOrCondition[1] = true;
    } else if (mSettings["condition_name"].IsSubParameter()) {
        for (auto& r_sub_parameter : mSettings["condition_name"]) {
            const std::string& r_condition_name = r_sub_parameter.GetString();
            KRATOS_ERROR_IF(r_condition_name != "" && !KratosComponents<Condition>::Has(r_condition_name)) << "Condition name not found in KratosComponents< Condition > -- name is " << r_condition_name << std::endl;
        }
        mOnlyOneElementOrCondition[1] = false;
    } else {
        KRATOS_ERROR << "The condition_name should be either a string or a list of strings" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ReplaceElementsAndConditionsProcess::Execute()
{
    ModelPart& r_root_model_part = mrModelPart.GetRootModelPart();
    /* Elements */
    bool replace_elements = false;
    if (mOnlyOneElementOrCondition[0]) {
        if (mSettings["element_name"].GetString() != "") replace_elements = true;
    } else {
        std::size_t counter = 0;
        for (auto& r_sub_parameter : mSettings["element_name"]) {
            if (r_sub_parameter.GetString() != "") counter += 1;
        }
        if (counter > 0) replace_elements = true;
    }
    if (replace_elements) {
        std::unordered_set<std::size_t> set_element_ids (mrModelPart.NumberOfElements());
        for(const auto& r_elem : mrModelPart.Elements()) {
            set_element_ids.insert(r_elem.Id());
        }
        if (mOnlyOneElementOrCondition[0]) {
            const std::string& r_element_name = mSettings["element_name"].GetString();
            ReplaceEntities(KratosComponents<Element>::Get(r_element_name), r_root_model_part.Elements(), set_element_ids);
        } else {
            ReplaceEntities<Element>(mSettings["element_name"], mGeometryTypesToStrings, r_root_model_part.Elements(), set_element_ids);
        }
        for (auto& r_sub_model_part : r_root_model_part.SubModelParts()) {
            UpdateElementsInSubModelPart(r_sub_model_part, r_root_model_part, set_element_ids);
        }
    }

    /* Conditions */
    bool replace_conditions = false;
    if (mOnlyOneElementOrCondition[1]) {
        if (mSettings["condition_name"].GetString() != "") replace_conditions = true;
    } else {
        std::size_t counter = 0;
        for (auto& r_sub_parameter : mSettings["condition_name"]) {
            if (r_sub_parameter.GetString() != "") counter += 1;
        }
        if (counter > 0) replace_conditions = true;
    }
    if (replace_conditions) {
        std::unordered_set<std::size_t> set_conditions_ids (mrModelPart.NumberOfConditions());
        for(const auto& r_cond : mrModelPart.Conditions()) {
            set_conditions_ids.insert(r_cond.Id());
        }
        if (mOnlyOneElementOrCondition[1]) {
            const std::string& r_condition_name = mSettings["condition_name"].GetString();
            ReplaceEntities(KratosComponents<Condition>::Get(r_condition_name), r_root_model_part.Conditions(), set_conditions_ids);
        } else {
            ReplaceEntities<Condition>(mSettings["condition_name"], mGeometryTypesToStrings, r_root_model_part.Conditions(), set_conditions_ids);
        }
        for (auto& r_sub_model_part : r_root_model_part.SubModelParts()) {
            UpdateConditionsInSubModelPart(r_sub_model_part, r_root_model_part, set_conditions_ids);
        }
    }
}

}  // namespace Kratos.