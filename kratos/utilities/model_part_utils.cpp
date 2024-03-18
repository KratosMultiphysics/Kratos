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

// System includes

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/model_part_utils.h"

namespace Kratos
{

template<class TEntity>
const TEntity& CheckEntity(const std::string& rEntityName)
{
    // Check if the entity is registered in Kratos
    if(!KratosComponents<TEntity>::Has(rEntityName)) {
        std::stringstream buffer;
        if constexpr(std::is_same_v<TEntity, Element>) {
            buffer << "Element " << rEntityName << " is not registered in Kratos.";
            buffer << " Please check the spelling of the element name and see if the application which containing it, is registered correctly.";
        } else if constexpr(std::is_same_v<TEntity, Condition>) {
            buffer << "Condition " << rEntityName << " is not registered in Kratos.";
            buffer << " Please check the spelling of the condition name and see if the application which containing it, is registered correctly.";
        } else {
            static_assert(!std::is_same_v<TEntity, TEntity>, "Unsupported entity type.");
        }
        KRATOS_ERROR << buffer.str() << std::endl;
    }
    return KratosComponents<TEntity>::Get(rEntityName);
}

template const Element& CheckEntity<Element>(const std::string& rEntityName);
template const Condition& CheckEntity<Condition>(const std::string& rEntityName);

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void ModelPartUtils::FromConnectivityGenerateEntities(
    const std::string& rEntityName,
    const std::vector<std::vector<std::size_t>>& rEntitiesConnectivities,
    ModelPart::NodesContainerType& rThisNodes,
    PointerVectorSet<TEntity, IndexedObject, std::less<typename IndexedObject::result_type>, std::equal_to<typename IndexedObject::result_type>, typename TEntity::Pointer, std::vector<typename TEntity::Pointer>>& rThisEntities,
    const Properties::Pointer pProperties
    )
{
    // Some definitions
    const std::size_t current_size = rThisEntities.size();
    const std::size_t number_of_entities = rEntitiesConnectivities.size();
    rThisEntities.reserve(current_size + number_of_entities);

    // We assume that the first entity in the list is the one we want to clone
    const TEntity& r_clone_entity = CheckEntity<TEntity>(rEntityName);
    const std::size_t number_of_nodes_per_entity = r_clone_entity.GetGeometry().size();

    // From connectivities generate entities
    const auto entities_vector = IndexPartition<IndexType>(number_of_entities).for_each<AccumReduction<typename TEntity::Pointer>>([&](const IndexType i) {
        typename TEntity::NodesArrayType temp_entity_nodes;
        temp_entity_nodes.reserve(number_of_nodes_per_entity);
        for(std::size_t j = 0; j < number_of_nodes_per_entity; ++j){
            const std::size_t node_id = rEntitiesConnectivities[i][j];
            auto it_node = rThisNodes.find(node_id);
            KRATOS_ERROR_IF(it_node == rThisNodes.end()) << "Node with Id " << node_id << " does not exist in the nodes array." << std::endl;
            temp_entity_nodes.push_back(*(it_node.base()));
        }

        // Add the entity to the list
        return r_clone_entity.Create(i + 1, temp_entity_nodes, pProperties);
    });

    // Add the entities to the list
    rThisEntities.insert(entities_vector.begin(), entities_vector.end());
}

template void ModelPartUtils::FromConnectivityGenerateEntities<Element>(
    const std::string& rEntityName,
    const std::vector<std::vector<std::size_t>>& rEntitiesConnectivities,
    ModelPart::NodesContainerType& rThisNodes,
    ModelPart::ElementsContainerType& rThisElements,
    const Properties::Pointer pProperties
    );

template void ModelPartUtils::FromConnectivityGenerateEntities<Condition>(
    const std::string& rEntityName,
    const std::vector<std::vector<std::size_t>>& rEntitiesConnectivities,
    ModelPart::NodesContainerType& rThisNodes,
    ModelPart::ConditionsContainerType& rThisConditions,
    const Properties::Pointer pProperties
    );

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void ModelPartUtils::FromConnectivityGenerateEntities(
    const std::string& rEntityName,
    const std::vector<std::size_t>& rEntitiesIds,
    const std::vector<std::size_t>& rPropertiesIds,
    const std::vector<std::vector<std::size_t>>& rEntitiesConnectivities,
    ModelPart::NodesContainerType& rThisNodes,
    ModelPart::PropertiesContainerType& rThisProperties,
    PointerVectorSet<TEntity, IndexedObject, std::less<typename IndexedObject::result_type>, std::equal_to<typename IndexedObject::result_type>, typename TEntity::Pointer, std::vector<typename TEntity::Pointer>>& rThisEntities
    )
{
    // Some definitions
    const std::size_t current_size = rThisEntities.size();
    const std::size_t number_of_entities = rEntitiesConnectivities.size();
    rThisEntities.reserve(current_size + number_of_entities);

    // We assume that the first entity in the list is the one we want to clone
    const TEntity& r_clone_entity = CheckEntity<TEntity>(rEntityName);
    const std::size_t number_of_nodes_per_entity = r_clone_entity.GetGeometry().size();

    // From connectivities generate entities
    // for(std::size_t i = 0; i < number_of_entities; ++i) {
    //     typename TEntity::NodesArrayType temp_entity_nodes;
    //     temp_entity_nodes.resize(number_of_nodes_per_entity);
    //     for(std::size_t j = 0 ; j < number_of_nodes_per_entity ; ++j){
    //         const std::size_t node_id = rEntitiesConnectivities[i][j];
    //         auto it_node = rThisNodes.find(node_id);
    //         KRATOS_ERROR_IF(it_node == rThisNodes.end()) << "Node with Id " << node_id << " does not exist in the nodes array." << std::endl;
    //         temp_entity_nodes.push_back(*(it_node.base()));
    //     }

    //     // Properties iterator
    //     auto it_properties = rThisProperties.find(rPropertiesIds[i]);
    //     KRATOS_ERROR_IF(it_properties == rThisProperties.end()) << "Properties with Id " << rPropertiesIds[i] << " does not exist in the properties array." << std::endl;
    //     Properties::Pointer p_properties = *(it_properties.base());

    //     // Add the entity to the list
    //     rThisEntities.push_back(r_clone_entity.Create(rEntitiesIds[i], temp_entity_nodes, p_properties));
    // }

    // From connectivities generate entities
    const auto entities_vector = IndexPartition<IndexType>(number_of_entities).for_each<AccumReduction<typename TEntity::Pointer>>([&](const IndexType i) {
        typename TEntity::NodesArrayType temp_entity_nodes;
        temp_entity_nodes.reserve(number_of_nodes_per_entity);
        for(std::size_t j = 0 ; j < number_of_nodes_per_entity; ++j){
            const std::size_t node_id = rEntitiesConnectivities[i][j];
            auto it_node = rThisNodes.find(node_id);
            KRATOS_ERROR_IF(it_node == rThisNodes.end()) << "Node with Id " << node_id << " does not exist in the nodes array." << std::endl;
            temp_entity_nodes.push_back(*(it_node.base()));
        }

        // Properties iterator
        auto it_properties = rThisProperties.find(rPropertiesIds[i]);
        KRATOS_ERROR_IF(it_properties == rThisProperties.end()) << "Properties with Id " << rPropertiesIds[i] << " does not exist in the properties array." << std::endl;
        Properties::Pointer p_properties = *(it_properties.base());

        // Add the entity to the list
        return r_clone_entity.Create(rEntitiesIds[i], temp_entity_nodes, p_properties);
    });

    // Add the entities to the list
    rThisEntities.insert(entities_vector.begin(), entities_vector.end());
}

template void ModelPartUtils::FromConnectivityGenerateEntities<Element>(
    const std::string& rEntityName,
    const std::vector<std::size_t>& rEntitiesIds,
    const std::vector<std::size_t>& rPropertiesIds,
    const std::vector<std::vector<std::size_t>>& rEntitiesConnectivities,
    ModelPart::NodesContainerType& rThisNodes,
    ModelPart::PropertiesContainerType& rThisProperties,
    ModelPart::ElementsContainerType& rThisElements
    );

template void ModelPartUtils::FromConnectivityGenerateEntities<Condition>(
    const std::string& rEntityName,
    const std::vector<std::size_t>& rEntitiesIds,
    const std::vector<std::size_t>& rPropertiesIds,
    const std::vector<std::vector<std::size_t>>& rEntitiesConnectivities,
    ModelPart::NodesContainerType& rThisNodes,
    ModelPart::PropertiesContainerType& rThisProperties,
    ModelPart::ConditionsContainerType& rThisConditions
    );

} // namespace Kratos
