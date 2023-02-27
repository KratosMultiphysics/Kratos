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
 * @param rReferenceEntity New type of entity that will replace old one
 * @param rEntityContainer Container of elements susceptible to be replaces
 * @param rSetOfIds Set of entities ids we want to replace
 */
template <class TContainer>
void ReplaceEntities(
    const typename TContainer::value_type& rReferenceEntity,
    TContainer& rEntityContainer,
    const std::unordered_set<std::size_t>& rSetOfIds
    )
{
    KRATOS_TRY;

    const auto& r_reference_geometry = rReferenceEntity.GetGeometry();
    const auto& r_reference_geometry_type = r_reference_geometry.GetGeometryType();
    IndexPartition<std::size_t>(rEntityContainer.size()).for_each([&rReferenceEntity, &rEntityContainer, &rSetOfIds, &r_reference_geometry, &r_reference_geometry_type](std::size_t Index){
        auto it_entity = rEntityContainer.begin() + Index;
        const std::size_t id = it_entity->Id();
        if (rSetOfIds.find(id) != rSetOfIds.end()) {
            auto p_geometry = it_entity->pGetGeometry();
            KRATOS_DEBUG_ERROR_IF_NOT(p_geometry->GetGeometryType() == r_reference_geometry_type) << "Trying to replace an entity with a different geometry type. Reference entity " << r_reference_geometry.Info() << " vs  " << p_geometry->Info() << "\n Entity info: " << rReferenceEntity.Info() << std::endl;
            auto p_new_entity = rReferenceEntity.Create(id, p_geometry, it_entity->pGetProperties());
            // Deep copy data and flags
            p_new_entity->Set(Flags(*it_entity));

            (*it_entity.base()) = p_new_entity;
        }
    });

    KRATOS_CATCH("");
}

/**
 * @brief Replace entities in a given container if the entity id is present in a list of ids.
 * @param rReferenceEntity New type of entity that will replace old one
 * @param rEntityContainer Container of elements susceptible to be replaces
 * @param rSetOfIds Set of entities ids we want to replace
 */
template <class TContainer>
void ReplaceEntities(
    const Parameters ListReferenceEntity,
    TContainer& rEntityContainer,
    const std::unordered_set<std::size_t>& rSetOfIds
    )
{
    KRATOS_TRY;

    struct TLS {
        typename TContainer::value_type::Pointer p_reference_entity = nullptr;
        GeometryData::KratosGeometryType current_geometry_type = GeometryData::KratosGeometryType::Kratos_generic_type;
        GeometryData::KratosGeometryType reference_geometry_type = GeometryData::KratosGeometryType::Kratos_generic_type;
    };
    IndexPartition<std::size_t>(rEntityContainer.size()).for_each(TLS(), [&ListReferenceEntity, &rEntityContainer, &rSetOfIds](std::size_t Index, TLS& rTLS){
        auto it_entity = rEntityContainer.begin() + Index;
        const std::size_t id = it_entity->Id();
        if (rSetOfIds.find(id) != rSetOfIds.end()) {
            auto p_geometry = it_entity->pGetGeometry();
            const auto& r_geometry_type = p_geometry->GetGeometryType();
            // Checking if geometry type is the same
            if (r_geometry_type != rTLS.current_geometry_type) {
                const std::string& r_type = GeometryUtils::GetGeometryName(r_geometry_type);
                KRATOS_ERROR_IF_NOT(ListReferenceEntity.Has(r_type)) << "Trying to replace an entity with a different geometry type. No reference entity found for geometry type: " << r_type << "\nReference list: " << ListReferenceEntity << std::endl;
                const auto& r_reference_entity = KratosComponents<typename TContainer::value_type>::Get(ListReferenceEntity[r_type].GetString());
                rTLS.p_reference_entity = r_reference_entity.Create(id, p_geometry, it_entity->pGetProperties());;
                rTLS.current_geometry_type = r_geometry_type;
                rTLS.reference_geometry_type = r_reference_entity.GetGeometry().GetGeometryType();
            }
            KRATOS_DEBUG_ERROR_IF_NOT(r_geometry_type == rTLS.reference_geometry_type) << "Trying to replace an entity with a different geometry type. Reference entity " << rTLS.p_reference_entity->GetGeometry().Info() << " vs  " << p_geometry->Info() << "\n Entity info: " << rTLS.p_reference_entity->Info() << std::endl;
            auto p_new_entity = rTLS.p_reference_entity->Create(id, p_geometry, it_entity->pGetProperties());
            // Deep copy data and flags
            p_new_entity->Set(Flags(*it_entity));

            (*it_entity.base()) = p_new_entity;
        }
    });

    KRATOS_CATCH("");
}

/**
 * @brief Replace entities in a given container if the entity id is present in a list of ids.
 * @param rName The name of the entity to be replaced
 * @param rEntityContainer Container of elements susceptible to be replaces
 * @param rSetOfIds Set of entities ids we want to replace
 */
template <class TContainer>
void ReplaceEntities(
    const std::string& rName,
    TContainer& rEntityContainer,
    const std::unordered_set<std::size_t>& rSetOfIds
    )
{
    KRATOS_TRY;

    struct TLS {
        typename TContainer::value_type::Pointer p_reference_entity = nullptr;
        GeometryData::KratosGeometryType current_geometry_type = GeometryData::KratosGeometryType::Kratos_generic_type;
        GeometryData::KratosGeometryType reference_geometry_type = GeometryData::KratosGeometryType::Kratos_generic_type;
    };
    IndexPartition<std::size_t>(rEntityContainer.size()).for_each(TLS(), [&rName, &rEntityContainer, &rSetOfIds](std::size_t Index, TLS& rTLS){
        auto it_entity = rEntityContainer.begin() + Index;
        const std::size_t id = it_entity->Id();
        if (rSetOfIds.find(id) != rSetOfIds.end()) {
            auto p_geometry = it_entity->pGetGeometry();
            const auto& r_geometry_type = p_geometry->GetGeometryType();
            // Checking if geometry type is the same
            if (r_geometry_type != rTLS.current_geometry_type) {
                const std::size_t dimension = p_geometry->WorkingSpaceDimension();
                const std::size_t number_of_nodes = p_geometry->size();
                const std::string replace_dimension = StringUtilities::ReplaceAllSubstrings(rName, "#D", std::to_string(dimension) + "D");
                const std::string replace_number_of_nodes = StringUtilities::ReplaceAllSubstrings(replace_dimension, "#N", std::to_string(number_of_nodes) + "N");
                KRATOS_ERROR_IF_NOT(KratosComponents<typename TContainer::value_type>::Has(replace_number_of_nodes)) << "Entity not registered: " << replace_number_of_nodes << std::endl;
                const auto& r_reference_entity = KratosComponents<typename TContainer::value_type>::Get(replace_number_of_nodes);
                rTLS.p_reference_entity = r_reference_entity.Create(id, p_geometry, it_entity->pGetProperties());;
                rTLS.current_geometry_type = r_geometry_type;
                rTLS.reference_geometry_type = r_reference_entity.GetGeometry().GetGeometryType();
            }
            KRATOS_DEBUG_ERROR_IF_NOT(r_geometry_type == rTLS.reference_geometry_type) << "Trying to replace an entity with a different geometry type. Reference entity " << rTLS.p_reference_entity->GetGeometry().Info() << " vs  " << p_geometry->Info() << "\n Entity info: " << rTLS.p_reference_entity->Info() << std::endl;
            auto p_new_entity = rTLS.p_reference_entity->Create(it_entity->Id(), p_geometry, it_entity->pGetProperties());
            // Deep copy data and flags
            p_new_entity->Set(Flags(*it_entity));

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
        if (r_element_name.find(";") == std::string::npos) {
            if (r_element_name.find("#") == std::string::npos) {
                // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them so that an error is thrown if they don't exist
                KRATOS_ERROR_IF(r_element_name != "" && !KratosComponents<Element>::Has(r_element_name)) << "Element name not found in KratosComponents< Element > -- name is " << r_element_name << std::endl;
                mDefinitionElementCondition[0] = DefinitionType::Single;
            } else {
                mDefinitionElementCondition[0] = DefinitionType::Templated;
            }
        } else {
            const std::vector<std::string> splitted_names = StringUtilities::SplitStringByDelimiter(r_element_name, ';');
            mSettings.RemoveValue("element_name");
            mSettings.AddEmptyValue("element_name");
            for (auto& r_element_name : splitted_names) {
                KRATOS_ERROR_IF(r_element_name != "" && !KratosComponents<Element>::Has(r_element_name)) << "Element name not found in KratosComponents< Element > -- name is " << r_element_name << std::endl;
                const auto& r_ref_element = KratosComponents<Element>::Get(r_element_name);
                const auto& r_reference_geometry = r_ref_element.GetGeometry();
                const auto& r_reference_geometry_type = r_reference_geometry.GetGeometryType();
                mSettings["element_name"].AddString(GeometryUtils::GetGeometryName(r_reference_geometry_type), r_element_name);
            }
            mDefinitionElementCondition[0] = DefinitionType::Multiple;
        }
    }

    /* Conditions */
    if (mSettings.Has("condition_name")) {
        const std::string& r_condition_name = mSettings["condition_name"].GetString();
        if (r_condition_name.find(";") == std::string::npos) {
            if (r_condition_name.find("#") == std::string::npos) {
                // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them so that an error is thrown if they don't exist
                KRATOS_ERROR_IF(r_condition_name != "" && !KratosComponents<Condition>::Has(r_condition_name)) << "Element name not found in KratosComponents< Condition > -- name is " << r_condition_name << std::endl;
                mDefinitionElementCondition[1] = DefinitionType::Single;
            } else {
                mDefinitionElementCondition[1] = DefinitionType::Templated;
            }
        } else {
            const std::vector<std::string> splitted_names = StringUtilities::SplitStringByDelimiter(r_condition_name, ';');
            mSettings.RemoveValue("condition_name");
            mSettings.AddEmptyValue("condition_name");
            for (auto& r_condition_name : splitted_names) {
                KRATOS_ERROR_IF(r_condition_name != "" && !KratosComponents<Condition>::Has(r_condition_name)) << "Condition name not found in KratosComponents< Condition > -- name is " << r_condition_name << std::endl;
                const auto& r_ref_element = KratosComponents<Condition>::Get(r_condition_name);
                const auto& r_reference_geometry = r_ref_element.GetGeometry();
                const auto& r_reference_geometry_type = r_reference_geometry.GetGeometryType();
                mSettings["condition_name"].AddString(GeometryUtils::GetGeometryName(r_reference_geometry_type), r_condition_name);
            }
            mDefinitionElementCondition[1] = DefinitionType::Multiple;
        }
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
    bool replace_elements = false;
    if (mDefinitionElementCondition[0] == DefinitionType::Multiple) {
        std::size_t counter = 0;
        for (auto& r_sub_parameter : mSettings["element_name"]) {
            if (r_sub_parameter.GetString() != "") counter += 1;
        }
        if (counter > 0) replace_elements = true;
    } else {
        if (mSettings["element_name"].GetString() != "") replace_elements = true;
    }
    if (replace_elements) {
        std::unordered_set<std::size_t> set_element_ids (mrModelPart.NumberOfElements());
        for(const auto& r_elem : mrModelPart.Elements()) {
            set_element_ids.insert(r_elem.Id());
        }
        if (mDefinitionElementCondition[0] == DefinitionType::Single) {
            const std::string& r_element_name = mSettings["element_name"].GetString();
            ReplaceEntities(KratosComponents<Element>::Get(r_element_name), r_root_model_part.Elements(), set_element_ids);
        } else if (mDefinitionElementCondition[0] == DefinitionType::Multiple) {
            ReplaceEntities(mSettings["element_name"], r_root_model_part.Elements(), set_element_ids);
        } else {
            const std::string& r_element_name = mSettings["element_name"].GetString();
            ReplaceEntities(r_element_name, r_root_model_part.Elements(), set_element_ids);
        }
        for (auto& r_sub_model_part : r_root_model_part.SubModelParts()) {
            UpdateElementsInSubModelPart(r_sub_model_part, r_root_model_part, set_element_ids);
        }
    }

    /* Conditions */
    bool replace_conditions = false;
    if (mDefinitionElementCondition[1] == DefinitionType::Multiple) {
        std::size_t counter = 0;
        for (auto& r_sub_parameter : mSettings["condition_name"]) {
            if (r_sub_parameter.GetString() != "") counter += 1;
        }
        if (counter > 0) replace_conditions = true;
    } else {
        if (mSettings["condition_name"].GetString() != "") replace_conditions = true;
    }
    if (replace_conditions) {
        std::unordered_set<std::size_t> set_conditions_ids (mrModelPart.NumberOfConditions());
        for(const auto& r_cond : mrModelPart.Conditions()) {
            set_conditions_ids.insert(r_cond.Id());
        }
        if (mDefinitionElementCondition[1] == DefinitionType::Single) {
            const std::string& r_condition_name = mSettings["condition_name"].GetString();
            ReplaceEntities(KratosComponents<Condition>::Get(r_condition_name), r_root_model_part.Conditions(), set_conditions_ids);
        } else if (mDefinitionElementCondition[1] == DefinitionType::Multiple) {
            ReplaceEntities(mSettings["condition_name"], r_root_model_part.Conditions(), set_conditions_ids);
        } else {
            const std::string& r_condition_name = mSettings["condition_name"].GetString();
            ReplaceEntities(r_condition_name, r_root_model_part.Conditions(), set_conditions_ids);
        }
        for (auto& r_sub_model_part : r_root_model_part.SubModelParts()) {
            UpdateConditionsInSubModelPart(r_sub_model_part, r_root_model_part, set_conditions_ids);
        }
    }

    KRATOS_CATCH("");
}

}  // namespace Kratos.