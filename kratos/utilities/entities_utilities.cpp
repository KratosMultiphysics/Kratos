//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes

// Project includes
#include "utilities/entities_utilities.h"
#include "utilities/string_utilities.h"
#include "utilities/geometry_utilities.h"

namespace Kratos::EntitiesUtilities
{

template<class TEntity>
EntitityIdentifier<TEntity>::EntitityIdentifier(const std::string& rName)
    : mName(rName)
{
    // Check if the name is empty
    if (rName != "") {
        // Check the type of entity considered
        std::string entity_name;
        if constexpr (std::is_same<TEntity, Element>::value) {
            entity_name = "Element";
        } else if constexpr (std::is_same<TEntity, Condition>::value) {
            entity_name = "Condition";
        } else {
            static_assert(std::is_same<TEntity, Element>::value || std::is_same<TEntity, Condition>::value, "Unsupported entity type. Only element and condition are supported.");
        }

        // Identify type of definition
        if (rName.find(";") == std::string::npos) { // Single or templated definition
            // Identified if simple or templated
            if (rName.find("#") == std::string::npos) {
                // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them so that an error is thrown if they don't exist
                KRATOS_ERROR_IF(rName != "" && !KratosComponents<TEntity>::Has(rName)) << entity_name << " name not found in KratosComponents<" << entity_name << "> -- name is " << rName << std::endl;
                mpPrototypeEntity = KratosComponents<TEntity>::Get(rName).Create(0, nullptr, nullptr);
            } else {
                mDefinitionType = DefinitionType::Templated;
                typename GeometryType::Pointer p_geometry_type = Kratos::make_shared<GeometryType>();
                mpPrototypeEntity = KratosComponents<TEntity>::Get(entity_name + "2D2N").Create(0, p_geometry_type, nullptr);
            }
        } else { // Multiple definition
            const auto splitted_names = StringUtilities::SplitStringByDelimiter(rName, ';');
            for (auto& r_entity_name : splitted_names) {
                // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them so that an error is thrown if they don't exist
                KRATOS_ERROR_IF(r_entity_name != "" && !KratosComponents<TEntity>::Has(r_entity_name)) << entity_name << " name not found in KratosComponents<" << entity_name << "> -- name is " << r_entity_name << std::endl;
                const auto& r_ref_entity = KratosComponents<TEntity>::Get(r_entity_name);
                const auto& r_reference_geometry = r_ref_entity.GetGeometry();
                const auto& r_reference_geometry_type = r_reference_geometry.GetGeometryType();
                mTypes.insert({r_reference_geometry_type, r_entity_name});
            }
            mDefinitionType = EntitiesUtilities::DefinitionType::Multiple;
            typename GeometryType::Pointer p_geometry_type = Kratos::make_shared<GeometryType>();
            mpPrototypeEntity = KratosComponents<TEntity>::Get(entity_name + "2D2N").Create(0, p_geometry_type, nullptr);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
bool EntitityIdentifier<TEntity>::IsInitialized()
{
    if (mName == "") {
        return false;
    } else {
        // Check if the class is initialized
        if (mDefinitionType == DefinitionType::Multiple) {
            // Count the number of non-empty types
            std::size_t counter = 0;
            for (auto& r_sub_replacement : mTypes) {
                if (r_sub_replacement.second != "") counter += 1;
            }
            // Check counter
            if (counter > 0) {
                return true;
            } else {
                return false;
            }
        } else { // Directly check the name
            return true;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
const TEntity& EntitityIdentifier<TEntity>::GetPrototypeEntity(typename GeometryType::Pointer pGeometry)
{
    switch (mDefinitionType) {
        case DefinitionType::Single:
            KRATOS_DEBUG_ERROR_IF_NOT(pGeometry->GetGeometryType() == KratosComponents<TEntity>::Get(mName).GetGeometry().GetGeometryType()) << "Trying to replace an entity with a different geometry type. Reference entity " << KratosComponents<TEntity>::Get(mName).GetGeometry().Info() << " vs  " << pGeometry->Info() << "\n Entity info: " << KratosComponents<TEntity>::Get(mName).Info() << std::endl;
            return *mpPrototypeEntity;
        case DefinitionType::Templated:
            // We check if the current type is correct or retrieving is required
            if (pGeometry->GetGeometryType() != mpPrototypeEntity->GetGeometry().GetGeometryType()) {
                const std::size_t dimension = pGeometry->WorkingSpaceDimension();
                const std::size_t number_of_nodes = pGeometry->size();
                const std::string replace_dimension = StringUtilities::ReplaceAllSubstrings(mName, "#D", std::to_string(dimension) + "D");
                const std::string replace_number_of_nodes = StringUtilities::ReplaceAllSubstrings(replace_dimension, "#N", std::to_string(number_of_nodes) + "N");
                KRATOS_ERROR_IF_NOT(KratosComponents<TEntity>::Has(replace_number_of_nodes)) << "Entity not registered: " << replace_number_of_nodes << std::endl;
                mpPrototypeEntity = KratosComponents<TEntity>::Get(replace_number_of_nodes).Create(0, pGeometry, nullptr);
            }
            return *mpPrototypeEntity;
        case DefinitionType::Multiple:
            // We check if the current type is correct or retrieving is required
            if (pGeometry->GetGeometryType() != mpPrototypeEntity->GetGeometry().GetGeometryType()) {
                mpPrototypeEntity = KratosComponents<TEntity>::Get(mTypes[pGeometry->GetGeometryType()]).Create(0, pGeometry, nullptr);
            }
            return *mpPrototypeEntity;
        default:
            KRATOS_ERROR << "Unsupported definition type" << std::endl;
            break;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void EntitityIdentifier<TEntity>::PrintData(std::ostream& rOStream) const
{
    rOStream << "Name: " << mName << "\n";
    std::string definition_type = mDefinitionType == DefinitionType::Single ? "Single" : mDefinitionType == DefinitionType::Multiple ? "Multiple" : "Templated";
    rOStream << "Definition type: " << definition_type << "\n";
    if (mDefinitionType == DefinitionType::Multiple) {
        rOStream << "Types: " << "\n";
        for (const auto& r_type : mTypes) {
            rOStream << "\t" << GeometryUtils::GetGeometryName(r_type.first) << " : " << r_type.second << "\n";
        }
    }
    rOStream << "Current prototype: " << *mpPrototypeEntity << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template struct EntitityIdentifier<Condition>;
template struct EntitityIdentifier<Element>;

/***********************************************************************************/
/***********************************************************************************/

template<> KRATOS_API(KRATOS_CORE) PointerVectorSet<Element, IndexedObject>& GetEntities<Element>(ModelPart& rModelPart);
template<> KRATOS_API(KRATOS_CORE) PointerVectorSet<Condition, IndexedObject>& GetEntities<Condition>(ModelPart& rModelPart);
template<> KRATOS_API(KRATOS_CORE) PointerVectorSet<MasterSlaveConstraint, IndexedObject>& GetEntities<MasterSlaveConstraint>(ModelPart& rModelPart);

void InitializeAllEntities(ModelPart& rModelPart)
{
    KRATOS_TRY

    InitializeEntities<Element>(rModelPart);
    InitializeEntities<Condition>(rModelPart);
    InitializeEntities<MasterSlaveConstraint>(rModelPart);

    KRATOS_CATCH("")
}

void InitializeSolutionStepAllEntities(ModelPart& rModelPart)
{
    KRATOS_TRY

    InitializeSolutionStepEntities<Element>(rModelPart);
    InitializeSolutionStepEntities<Condition>(rModelPart);
    InitializeSolutionStepEntities<MasterSlaveConstraint>(rModelPart);

    KRATOS_CATCH("")
}

void FinalizeSolutionStepAllEntities(ModelPart& rModelPart)
{
    KRATOS_TRY

    FinalizeSolutionStepEntities<Element>(rModelPart);
    FinalizeSolutionStepEntities<Condition>(rModelPart);
    FinalizeSolutionStepEntities<MasterSlaveConstraint>(rModelPart);

    KRATOS_CATCH("")
}

void InitializeNonLinearIterationAllEntities(ModelPart& rModelPart)
{
    KRATOS_TRY

    InitializeNonLinearIterationEntities<Element>(rModelPart);
    InitializeNonLinearIterationEntities<Condition>(rModelPart);
    InitializeNonLinearIterationEntities<MasterSlaveConstraint>(rModelPart);

    KRATOS_CATCH("")
}

void FinalizeNonLinearIterationAllEntities(ModelPart& rModelPart)
{
    KRATOS_TRY

    FinalizeNonLinearIterationEntities<Element>(rModelPart);
    FinalizeNonLinearIterationEntities<Condition>(rModelPart);
    FinalizeNonLinearIterationEntities<MasterSlaveConstraint>(rModelPart);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Element, IndexedObject>& GetEntities<Element>(ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>& GetEntities<Condition>(ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<MasterSlaveConstraint, IndexedObject>& GetEntities<MasterSlaveConstraint>(ModelPart& rModelPart)
{
    return rModelPart.MasterSlaveConstraints();
}

} // namespace Kratos::EntitiesUtilities
