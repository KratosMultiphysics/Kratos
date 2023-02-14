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

namespace Kratos {
namespace {

std::string GetElementName(const GeometryData::KratosGeometryType ElementType)
{
    KRATOS_TRY;

    // Using switch over map as the compiler warns if some enum values are not handled in the switch
    switch(ElementType) {
        case GeometryData::KratosGeometryType::Kratos_Hexahedra3D20:
            return "Hexahedra3D20";
        case GeometryData::KratosGeometryType::Kratos_Hexahedra3D27:
            return "Hexahedra3D27";
        case GeometryData::KratosGeometryType::Kratos_Hexahedra3D8:
            return "Hexahedra3D8";
        case GeometryData::KratosGeometryType::Kratos_Prism3D15:
            return "Prism3D15";
        case GeometryData::KratosGeometryType::Kratos_Prism3D6:
            return "Prism3D6";
        case GeometryData::KratosGeometryType::Kratos_Pyramid3D13:
            return "Pyramid3D13";
        case GeometryData::KratosGeometryType::Kratos_Pyramid3D5:
            return "Pyramid3D5";
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4:
            return "Quadrilateral2D4";
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8:
            return "Quadrilateral2D8";
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9:
            return "Quadrilateral2D9";
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4:
            return "Quadrilateral3D4";
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8:
            return "Quadrilateral3D8";
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9:
            return "Quadrilateral3D9";
        case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10:
            return "Tetrahedra3D10";
        case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
            return "Tetrahedra3D4";
        case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
            return "Triangle3D3";
        case GeometryData::KratosGeometryType::Kratos_Triangle2D6:
            return "Triangle2D6";
        case GeometryData::KratosGeometryType::Kratos_Triangle2D10:
            return "Triangle2D10";
        case GeometryData::KratosGeometryType::Kratos_Triangle2D15:
            return "Triangle2D15";
        case GeometryData::KratosGeometryType::Kratos_Triangle3D3:
            return "Triangle3D3";
        case GeometryData::KratosGeometryType::Kratos_Triangle3D6:
            return "Triangle3D6";
        case GeometryData::KratosGeometryType::Kratos_Line2D2:
            return "Line2D2";
        case GeometryData::KratosGeometryType::Kratos_Line2D3:
            return "Line2D3";
        case GeometryData::KratosGeometryType::Kratos_Line2D4:
            return "Line2D4";
        case GeometryData::KratosGeometryType::Kratos_Line2D5:
            return "Line2D5";
        case GeometryData::KratosGeometryType::Kratos_Line3D2:
            return "Line3D2";
        case GeometryData::KratosGeometryType::Kratos_Line3D3:
            return "Line3D3";
        case GeometryData::KratosGeometryType::Kratos_Point2D:
            return "Point2D";
        case GeometryData::KratosGeometryType::Kratos_Point3D:
            return "Point3D";
        case GeometryData::KratosGeometryType::Kratos_Sphere3D1:
            return "Sphere3D1";
        case GeometryData::KratosGeometryType::Kratos_Nurbs_Curve:
            return "Nurbs_Curve";
        case GeometryData::KratosGeometryType::Kratos_Nurbs_Surface:
            return "Nurbs_Surface";
        case GeometryData::KratosGeometryType::Kratos_Nurbs_Volume:
            return "Nurbs_Volume";
        case GeometryData::KratosGeometryType::Kratos_Nurbs_Curve_On_Surface:
            return "Nurbs_Curve_On_Surface";
        case GeometryData::KratosGeometryType::Kratos_Surface_In_Nurbs_Volume:
            return "Surface_In_Nurbs_Volume";
        case GeometryData::KratosGeometryType::Kratos_Brep_Curve:
            return "Brep_Curve";
        case GeometryData::KratosGeometryType::Kratos_Brep_Surface:
            return "Brep_Surface";
        case GeometryData::KratosGeometryType::Kratos_Brep_Curve_On_Surface:
            return "Brep_Curve_On_Surface";
        case GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Geometry:
            return "Quadrature_Point_Geometry";
        case GeometryData::KratosGeometryType::Kratos_Coupling_Geometry:
            return "Coupling_Geometry";
        case GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Curve_On_Surface_Geometry:
            return "Quadrature_Point_Curve_On_Surface_Geometry";
        case GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Surface_In_Volume_Geometry:
            return "Quadrature_Point_Surface_In_Volume_Geometry";
        case GeometryData::KratosGeometryType::NumberOfGeometryTypes:
            KRATOS_ERROR << "Geometry type not supported" << std::endl;
            return "NumberOfGeometryTypes";
        default:
            KRATOS_ERROR << "Geometry type not supported: " << static_cast<int>(ElementType) << std::endl;
    };

    KRATOS_CATCH("");
}

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
    KRATOS_TRY;

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

    KRATOS_CATCH("");
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
    PointerVectorSet<TEntity, IndexedObject, std::less<typename IndexedObject::result_type>, std::equal_to<typename IndexedObject::result_type>, typename TEntity::Pointer, std::vector< typename TEntity::Pointer>>& rEntityContainer,
    std::unordered_set<std::size_t>& rSetOfIds
    )
{
    KRATOS_TRY;

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
                const std::string& r_type = GetElementName(r_geometry_type);
                KRATOS_ERROR_IF_NOT(ListReferenceEntity.Has(r_type)) << "Trying to replace an element with a different geometry type. No reference entity found for geometry type: " << r_type << "\nReference list: " << ListReferenceEntity << std::endl;
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

    KRATOS_CATCH("");
}

/**
 * @brief Replace entities in a given container if the entity id is present in a list of ids.
 * @param rName The name of the entity to be replaced
 * @param rEntityContainer Container of elements susceptible to be replaces
 * @param rSetOfIds Set of entities ids we want to replace
 */
template <class TEntity>
void ReplaceEntities(
    const std::string& rName,
    PointerVectorSet<TEntity, IndexedObject, std::less<typename IndexedObject::result_type>, std::equal_to<typename IndexedObject::result_type>, typename TEntity::Pointer, std::vector< typename TEntity::Pointer>>& rEntityContainer,
    std::unordered_set<std::size_t>& rSetOfIds
    )
{
    KRATOS_TRY;

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
                const std::size_t dimension = p_geometry->WorkingSpaceDimension();
                const std::size_t number_of_nodes = p_geometry->size();
                const std::string replace_dimension = StringUtilities::ReplaceAllSubstrings(rName, "#D", std::to_string(dimension) + "D");
                const std::string replace_number_of_nodes = StringUtilities::ReplaceAllSubstrings(replace_dimension, "#N", std::to_string(number_of_nodes) + "N");
                KRATOS_ERROR_IF_NOT(KratosComponents<TEntity>::Has(replace_number_of_nodes)) << "Entity not registered: " << replace_number_of_nodes << std::endl;
                const auto& r_reference_entity = KratosComponents<TEntity>::Get(replace_number_of_nodes);
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
                mSettings["element_name"].AddString(GetElementName(r_reference_geometry_type), r_element_name);
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
                mSettings["condition_name"].AddString(GetElementName(r_reference_geometry_type), r_condition_name);
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
            ReplaceEntities<Element>(mSettings["element_name"], r_root_model_part.Elements(), set_element_ids);
        } else {
            const std::string& r_element_name = mSettings["element_name"].GetString();
            ReplaceEntities<Element>(r_element_name, r_root_model_part.Elements(), set_element_ids);
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
            ReplaceEntities<Condition>(mSettings["condition_name"], r_root_model_part.Conditions(), set_conditions_ids);
        } else {
            const std::string& r_condition_name = mSettings["condition_name"].GetString();
            ReplaceEntities<Condition>(r_condition_name, r_root_model_part.Conditions(), set_conditions_ids);
        }
        for (auto& r_sub_model_part : r_root_model_part.SubModelParts()) {
            UpdateConditionsInSubModelPart(r_sub_model_part, r_root_model_part, set_conditions_ids);
        }
    }

    KRATOS_CATCH("");
}

}  // namespace Kratos.