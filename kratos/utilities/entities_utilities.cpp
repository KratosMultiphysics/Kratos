//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Philipp Bucher (https://github.com/philbucher)
//

// System includes
#include <regex>

// External includes

// Project includes
#include "utilities/entities_utilities.h"
#include "utilities/string_utilities.h"
#include "utilities/geometry_utilities.h"

namespace Kratos::EntitiesUtilities
{

template<class TEntity>
EntitityIdentifier<TEntity>::EntitityIdentifier(const std::string& rName)
{
    // Check if the name is empty (could be empty due to default values in settings)
    if (rName != "") {
        // Initialize the definition types
        mTypes.fill(nullptr);

        // Identify type of definition
        if (rName.find("#") != std::string::npos) { // Templated definition
            if (rName.find(";") == std::string::npos) { // Must be single definition
                GenerateTemplatedTypes(rName);
            } else {
                KRATOS_ERROR << "Unsupported definition type. Cannot use # and ; together" << std::endl;
            }
        } else { // Multiple definition or single definition
            if (rName.find(";") == std::string::npos) { // Must be single definition
                GenerateSingleType(rName);
            } else {
                GenerateMultipleTypes(rName);
            }
        }
        mIsInitialized = true;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
bool EntitityIdentifier<TEntity>::IsInitialized() const
{
    return mIsInitialized;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
bool EntitityIdentifier<TEntity>::HasPrototypeEntity(const GeometryType& rGeometry) const
{
    return mTypes[static_cast<std::size_t>(rGeometry.GetGeometryType())] != nullptr;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
const TEntity& EntitityIdentifier<TEntity>::GetPrototypeEntity(typename GeometryType::Pointer pGeometry) const
{
    KRATOS_DEBUG_ERROR_IF(mTypes[static_cast<std::size_t>(pGeometry->GetGeometryType())] == nullptr) << "Prototype not initialized for " << GeometryUtils::GetGeometryName(pGeometry->GetGeometryType()) << std::endl;
    return *mTypes[static_cast<std::size_t>(pGeometry->GetGeometryType())];
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
const TEntity& EntitityIdentifier<TEntity>::GetPrototypeEntity(const GeometryType& rGeometry) const
{
    KRATOS_DEBUG_ERROR_IF(mTypes[static_cast<std::size_t>(rGeometry.GetGeometryType())] == nullptr) << "Prototype not initialized for " << GeometryUtils::GetGeometryName(rGeometry.GetGeometryType()) << std::endl;
    return *mTypes[static_cast<std::size_t>(rGeometry.GetGeometryType())];
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void EntitityIdentifier<TEntity>::PrintData(std::ostream& rOStream) const
{
    rOStream << "Types: " << "\n";
    std::size_t counter = 1;
    for (const auto& r_type : mTypes) {
        if (r_type != nullptr) {
            rOStream << "\t" << GeometryUtils::GetGeometryName(static_cast<GeometryData::KratosGeometryType>(counter)) << " : " << r_type->Info() << "\n";
        }
        ++counter;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
std::string EntitityIdentifier<TEntity>::GetEntityTypeName() const
{
    std::string entity_type_name;
    if constexpr (std::is_same<TEntity, Element>::value) {
        entity_type_name = "Element";
    } else if constexpr (std::is_same<TEntity, Condition>::value) {
        entity_type_name = "Condition";
    } else {
        static_assert(std::is_same<TEntity, Element>::value || std::is_same<TEntity, Condition>::value, "Unsupported entity type. Only element and condition are supported.");
    }
    return entity_type_name;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void EntitityIdentifier<TEntity>::GenerateSingleType(const std::string& rName)
{
    // Check the type of entity considered
    const std::string entity_type_name = GetEntityTypeName();

    // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them so that an error is thrown if they don't exist
    KRATOS_ERROR_IF(rName != "" && !KratosComponents<TEntity>::Has(rName)) << entity_type_name << " name not found in KratosComponents<" << entity_type_name << "> -- name is " << rName << std::endl;
    const auto& r_ref_entity = KratosComponents<TEntity>::Get(rName);
    const auto& r_reference_geometry = r_ref_entity.GetGeometry();
    const auto& r_reference_geometry_type = r_reference_geometry.GetGeometryType();
    const std::size_t index = static_cast<std::size_t>(r_reference_geometry_type);
    mTypes[index] = &r_ref_entity;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void EntitityIdentifier<TEntity>::GenerateMultipleTypes(const std::string& rName)
{
    // Check the type of entity considered
    const std::string entity_type_name = GetEntityTypeName();

    // Split the name
    const auto splitted_names = StringUtilities::SplitStringByDelimiter(rName, ';');
    for (auto& r_entity_name : splitted_names) {
        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them so that an error is thrown if they don't exist
        KRATOS_ERROR_IF(r_entity_name != "" && !KratosComponents<TEntity>::Has(r_entity_name)) << entity_type_name << " name not found in KratosComponents<" << entity_type_name << "> -- name is " << r_entity_name << std::endl;
        const auto& r_ref_entity = KratosComponents<TEntity>::Get(r_entity_name);
        const auto& r_reference_geometry = r_ref_entity.GetGeometry();
        const auto& r_reference_geometry_type = r_reference_geometry.GetGeometryType();
        const std::size_t index = static_cast<std::size_t>(r_reference_geometry_type);
        mTypes[index] = &r_ref_entity;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void EntitityIdentifier<TEntity>::GenerateTemplatedTypes(const std::string& rName)
{
    // Prepare regex
    const std::string replace_dimension = StringUtilities::ReplaceAllSubstrings(rName, "#D", "(2D|3D)");
    const std::string replace_number_of_nodes = StringUtilities::ReplaceAllSubstrings(replace_dimension, "#N", "[0-9]+N");
    std::regex pattern(replace_number_of_nodes);

    // Geometry combinations
    for (const auto& r_pair : KratosComponents<TEntity>::GetComponents()) {
        const std::string& r_key = r_pair.first;
        if (std::regex_match(r_key, pattern)) {
            const auto& r_entity = KratosComponents<TEntity>::Get(r_key);
            const std::size_t index = static_cast<std::size_t>(r_entity.GetGeometry().GetGeometryType());
            mTypes[index] = &r_entity;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class EntitityIdentifier<Condition>;
template class EntitityIdentifier<Element>;

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
