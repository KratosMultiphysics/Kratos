//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//                   Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <type_traits>

// External includes

// Project includes
#include "includes/kratos_export_api.h"
#include "includes/model_part.h"
#include "includes/global_variables.h"

namespace Kratos
{

/// @brief Class collecting a set of free-floating utility functions for querying and mutating @ref ModelPart s.
class KRATOS_API(KRATOS_CORE) ModelPartUtils
{
public:
    /// @brief Templated interface for getting nodes, elements, conditions or @ref ProcessInfo from a @ref ModelPart.
    template <Globals::DataLocation TLocation>
    static const auto& GetContainer(const ModelPart& rModelPart)
    {
        if constexpr (TLocation == Globals::DataLocation::NodeHistorical || TLocation == Globals::DataLocation::NodeNonHistorical) {
            return rModelPart.Nodes();
        } else if constexpr (TLocation == Globals::DataLocation::Element) {
            return rModelPart.Elements();
        } else if constexpr (TLocation == Globals::DataLocation::Condition) {
            return rModelPart.Conditions();
        } else if constexpr (TLocation == Globals::DataLocation::ProcessInfo) {
            return rModelPart.GetProcessInfo();
        } else if constexpr (TLocation == Globals::DataLocation::ModelPart) {
            return rModelPart;
        }
    }

    /// @brief Templated interface for getting nodes, elements, conditions or @ref ProcessInfo from a @ref ModelPart.
    template <Globals::DataLocation TLocation>
    static auto& GetContainer(ModelPart& rModelPart)
    {
        if constexpr (TLocation == Globals::DataLocation::NodeHistorical || TLocation == Globals::DataLocation::NodeNonHistorical) {
            return rModelPart.Nodes();
        } else if constexpr (TLocation == Globals::DataLocation::Element) {
            return rModelPart.Elements();
        } else if constexpr (TLocation == Globals::DataLocation::Condition) {
            return rModelPart.Conditions();
        } else if constexpr (TLocation == Globals::DataLocation::ProcessInfo) {
            return rModelPart.GetProcessInfo();
        } else if constexpr (TLocation == Globals::DataLocation::ModelPart) {
            return rModelPart;
        }
    }

    /// @brief Templated interface to get nodes, elements and conditions from a @ref ModelPart
    template<class TContainerType>
    static const auto& GetContainer(const ModelPart& rModelPart)
    {
        if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
            return GetContainer<Globals::DataLocation::NodeNonHistorical>(rModelPart);
        } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
            return GetContainer<Globals::DataLocation::Condition>(rModelPart);
        } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
            return GetContainer<Globals::DataLocation::Element>(rModelPart);
        } else {
            static_assert(!std::is_same_v<TContainerType, TContainerType>, "Unsupported container type.");
            return 0;
        }
    }

    /// @brief Templated interface to get nodes, elements and conditions from a @ref ModelPart
    template<class TContainerType>
    static auto& GetContainer(ModelPart& rModelPart)
    {
        if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
            return GetContainer<Globals::DataLocation::NodeNonHistorical>(rModelPart);
        } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
            return GetContainer<Globals::DataLocation::Condition>(rModelPart);
        } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
            return GetContainer<Globals::DataLocation::Element>(rModelPart);
        } else {
            static_assert(!std::is_same_v<TContainerType, TContainerType>, "Unsupported container type.");
            return 0;
        }
    }

    /**
    * @brief Generates entities (e.g., Elements or Conditions) based on provided connectivities and adds them to a specified container, with an option to use default properties.
    * @details This function is designed to create entities of a specific type (denoted by TEntity) based on their node connectivities. It facilitates the construction of these entities by looking up nodes in the provided nodes container and either using a specified Properties object or defaulting to a provided one if not explicitly specified for each entity. The generated entities are then added to the specified container. The function ensures that the entity type is registered in Kratos and that all specified nodes are available. It allows specifying a common Properties object to be used for all entities, which can be particularly useful when all entities share the same physical properties.
    * @tparam TEntity The type of entities to be generated (typically Element or Condition). Must be a type registered with KratosComponents.
    * @param rEntityName Name of the entity type to generate, which must be registered in KratosComponents.
    * @param rEntitiesConnectivities A vector of vectors, where each sub-vector contains the node IDs defining the connectivity for one entity.
    * @param rThisNodes Reference to the container holding the nodes available in the current model part.
    * @param rThisEntities Reference to the container where the new entities will be stored.
    * @param pProperties Optional pointer to a Properties object to be used for all created entities. If not provided (nullptr), the function will use the first Properties object found in rThisProperties for all entities.
    * @throws std::runtime_error If the specified entity type is not registered in Kratos, or if any of the specified nodes do not exist.
    */
    template<class TEntity>
    static void GenerateEntitiesFromConnectivities(
        const std::string& rEntityName,
        const std::vector<std::vector<std::size_t>>& rEntitiesConnectivities,
        ModelPart::NodesContainerType& rThisNodes,
        PointerVectorSet<TEntity, IndexedObject, std::less<typename IndexedObject::result_type>, std::equal_to<typename IndexedObject::result_type>, typename TEntity::Pointer, std::vector<typename TEntity::Pointer>>& rThisEntities,
        const Properties::Pointer pProperties = nullptr
        );

    /**
     * @brief Generates entities (Elements or Conditions) based on the provided connectivities and adds them to a specified container.
     * @details This function creates entities of type TEntity (Element or Condition) based on their connectivity information. It looks up
     * the specified entities and properties by their IDs in the provided containers and constructs new entities accordingly. The new entities are then added to the specified container. It checks if the specified entity type is registered in Kratos and if the necessary nodes and properties are available. Throws an error if any requirement is not met.
     * @tparam TEntity The type of the entities to be generated (Element or Condition). Must be a type registered in Kratos.
     * @param rEntityName The name of the entity type to generate. Must be registered in KratosComponents.
     * @param rEntitiesIds A vector of IDs for the new entities to be created.
     * @param rPropertiesIds A vector of IDs for the properties to be associated with the new entities.
     * @param rEntitiesConnectivities A vector of vectors, each sub-vector provides the node IDs for the connectivity of one entity.
     * @param rThisNodes A reference to the container of nodes available in the current model part.
     * @param rThisProperties A reference to the container of properties available in the current model part.
     * @param rThisEntities A reference to the container where the new entities will be added.
     * @throws std::runtime_error If the specified entity type is not registered in Kratos, if any of the specified nodes do not exist,  or if any of the specified properties do not exist.
     */
    template<class TEntity>
    static void GenerateEntitiesFromConnectivities(
        const std::string& rEntityName,
        const std::vector<std::size_t>& rEntitiesIds,
        const std::vector<std::size_t>& rPropertiesIds,
        const std::vector<std::vector<std::size_t>>& rEntitiesConnectivities,
        ModelPart::NodesContainerType& rThisNodes,
        ModelPart::PropertiesContainerType& rThisProperties,
        PointerVectorSet<TEntity, IndexedObject, std::less<typename IndexedObject::result_type>, std::equal_to<typename IndexedObject::result_type>, typename TEntity::Pointer, std::vector<typename TEntity::Pointer>>& rThisEntities
        );

}; // class ModelPartUtils

} // namespace Kratos
