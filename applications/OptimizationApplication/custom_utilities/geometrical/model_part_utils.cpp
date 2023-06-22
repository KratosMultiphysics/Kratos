//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "optimization_application_variables.h"

// Include base h
#include "model_part_utils.h"

namespace Kratos {

template<class TDataType>
std::set<TDataType> ModelPartUtils::SetReduction<TDataType>::GetValue() const
{
    return mValue;
}

template<class TDataType>
void ModelPartUtils::SetReduction<TDataType>::LocalReduce(const value_type& rValue)
{
    mValue.emplace(rValue);
}

template<class TDataType>
void ModelPartUtils::SetReduction<TDataType>::ThreadSafeReduce(SetReduction<TDataType>& rOther)
{
    KRATOS_CRITICAL_SECTION
    mValue.merge(rOther.mValue);
}

template<class TEntityType, class TMapValueType>
std::map<IndexType, TMapValueType> ModelPartUtils::ContainerEntityMapReduction<TEntityType, TMapValueType>::GetValue() const
{
    return mValue;
}

template<class TEntityType, class TMapValueType>
void ModelPartUtils::ContainerEntityMapReduction<TEntityType, TMapValueType>::LocalReduce(const value_type& rValue)
{
    if constexpr(std::is_same_v<TMapValueType, EntityPointerType<TEntityType>>) {
        for (const auto& r_item : rValue) {
            mValue.emplace(r_item);
        }
    } else if constexpr(std::is_same_v<TMapValueType, std::vector<EntityPointerType<TEntityType>>>) {
        for (const auto& r_item : rValue) {
            mValue[r_item.first].push_back(r_item.second);
        }
    } else {
        KRATOS_ERROR << "Unsupported type for TMapValueType";
    }
}

template<class TEntityType, class TMapValueType>
void ModelPartUtils::ContainerEntityMapReduction<TEntityType, TMapValueType>::ThreadSafeReduce(ContainerEntityMapReduction<TEntityType, TMapValueType>& rOther)
{
    KRATOS_CRITICAL_SECTION
    if constexpr(std::is_same_v<TMapValueType, EntityPointerType<TEntityType>>) {
        mValue.merge(rOther.mValue);
    } else if constexpr(std::is_same_v<TMapValueType, std::vector<EntityPointerType<TEntityType>>>) {
        for (const auto& it : rOther.mValue) {
            auto& r_current_vector = mValue[it.first];
            for (auto p_item : it.second) {
                r_current_vector.push_back(p_item);
            }
        }
    } else {
        KRATOS_ERROR << "Unsupported type for TMapValueType";
    }
}

void ModelPartUtils::AppendModelPartNames(
    std::stringstream& rOutputStream,
    const std::vector<ModelPart*>& rModelParts)
{
    std::vector<std::string> mp_names_list(rModelParts.size());
    for (IndexType i = 0; i < rModelParts.size(); ++i) {
        mp_names_list[i] = rModelParts[i]->FullName();
    }
    std::sort(mp_names_list.begin(), mp_names_list.end());

    for (const auto& r_mp_name : mp_names_list) {
        rOutputStream << r_mp_name << ";";
    }
}

template<class TContainerType>
void ModelPartUtils::UpdateEntityIdsSetFromContainer(
    std::set<IndexType>& rOutput,
    const TContainerType& rContainer)
{
    auto entity_ids_set = block_for_each<SetReduction<IndexType>>(rContainer, [&](const auto& rEntity) {
        return rEntity.Id();
    });

    rOutput.merge(entity_ids_set);
}

template<class TContainerType>
void ModelPartUtils::UpdateEntityGeometryNodeIdsSetFromContainer(
    std::set<NodeIdsType>& rOutput,
    const TContainerType& rContainer)
{
    auto entity_geometry_node_ids_set = block_for_each<SetReduction<NodeIdsType>>(rContainer, [&](const auto& rEntity) {
        const auto& r_geometry = rEntity.GetGeometry();
        NodeIdsType geometry_node_ids(r_geometry.size());
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            geometry_node_ids[i] = r_geometry[i].Id();
        }
        std::sort(geometry_node_ids.begin(), geometry_node_ids.end());
        return geometry_node_ids;
    });

    rOutput.merge(entity_geometry_node_ids_set);
}

template<class TContainerType>
void ModelPartUtils::UpdateEntityIdEntityPtrMapWithCommonEntitiesFromContainerAndEntityIdsSet(
    std::map<IndexType, ContainerEntityPointerType<TContainerType>>& rOutput,
    const std::set<IndexType>& rEntityIdsSet,
    TContainerType& rContainer)
{
    using map_type = std::map<IndexType, ContainerEntityPointerType<TContainerType>>;

    using reduction_type = MapReduction<map_type>;

    // need to use ptr_iterator here because, we need to increment the reference counter of the entity intrusive_ptr when push_back is used.
    auto entity_id_ptr_map = block_for_each<reduction_type>(rContainer, [&](auto& rEntity) {
        auto it = rEntityIdsSet.find(rEntity.Id());
        if (it != rEntityIdsSet.end()) {
            return std::make_pair(rEntity.Id(), &rEntity);
        } else {
            // return a dummy entry which is removed later.
            return std::make_pair(std::numeric_limits<IndexType>::max(), &rEntity);
        }
    });

    // remove the unwanted dummy entry of max.
    if (entity_id_ptr_map.find(std::numeric_limits<IndexType>::max()) != entity_id_ptr_map.end()) {
        entity_id_ptr_map.erase(std::numeric_limits<IndexType>::max());
    }

    rOutput.merge(entity_id_ptr_map);
}

template<class TContainerType>
void ModelPartUtils::UpdateEntityIdEntityPtrMapWithCommonEntitiesFromContainerAndEntityGeometryNodeIdsSet(
    std::map<IndexType, ContainerEntityPointerType<TContainerType>>& rOutput,
    const std::set<NodeIdsType>& rEntityGeometryNodeIdsSet,
    TContainerType& rContainer)
{
    using map_type = std::map<IndexType, ContainerEntityPointerType<TContainerType>>;

    using reduction_type = MapReduction<map_type>;

    // need to use ptr_iterator here because, we need to increment the reference counter of the entity intrusive_ptr when push_back is used.
    auto entity_id_ptr_map = block_for_each<reduction_type>(rContainer, [&](auto& rEntity) {
        auto& r_geometry = rEntity.GetGeometry();
        NodeIdsType geometry_node_ids(r_geometry.size());
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            geometry_node_ids[i] = r_geometry[i].Id();
        }
        std::sort(geometry_node_ids.begin(), geometry_node_ids.end());

        auto it = rEntityGeometryNodeIdsSet.find(geometry_node_ids);
        if (it != rEntityGeometryNodeIdsSet.end()) {
            return std::make_pair(rEntity.Id(), &rEntity);
        } else {
            // return a dummy entry which is removed later.
            return std::make_pair(std::numeric_limits<IndexType>::max(), &rEntity);
        }
    });

    // remove the unwanted dummy entry of max.
    if (entity_id_ptr_map.find(std::numeric_limits<IndexType>::max()) != entity_id_ptr_map.end()) {
        entity_id_ptr_map.erase(std::numeric_limits<IndexType>::max());
    }

    rOutput.merge(entity_id_ptr_map);
}

template<class TContainerType>
void ModelPartUtils::UpdateNeighbourMaps(
    std::map<IndexType, std::vector<ContainerEntityPointerType<TContainerType>>>& rOutput,
    const std::set<IndexType>& rNodeIdsSet,
    TContainerType& rContainer)
{
    using entity_pointer_type = ContainerEntityPointerType<TContainerType>;

    using reduction_type = ContainerEntityMapReduction<typename TContainerType::value_type, std::vector<entity_pointer_type>>;

    // need to use ptr_iterator here because, we need to increment the reference counter of the entity intrusive_ptr when push_back is used.
    auto entity_id_ptrs_map = block_for_each<reduction_type>(rContainer, [&](auto& rEntity) {
        std::vector<std::pair<IndexType, entity_pointer_type>> items;
        for (const auto& r_node : rEntity.GetGeometry()) {
            auto itr = rNodeIdsSet.find(r_node.Id());
            if (itr != rNodeIdsSet.end()) {
                items.push_back(std::make_pair(r_node.Id(), &rEntity));
            }
        }
        return items;
    });

    rOutput.merge(entity_id_ptrs_map);
}

template<class TEntityPointerType>
void ModelPartUtils::UpdateEntityIdEntityPtrMapFromNeighbourMap(
    std::map<IndexType, TEntityPointerType>& rOutput,
    const std::map<IndexType, std::vector<TEntityPointerType>>& rNodeIdNeighbourEntityPtrsMap)
{
    std::vector<IndexType> node_ids;
    for (const auto& it : rNodeIdNeighbourEntityPtrsMap) {
        node_ids.push_back(it.first);
    }

    auto entity_id_ptr_map = block_for_each<ContainerEntityMapReduction<typename TEntityPointerType::element_type, TEntityPointerType>>(node_ids, [&](const auto NodeId) {
        std::vector<std::pair<IndexType, TEntityPointerType>> items;
        const auto& r_neighbour_entities = rNodeIdNeighbourEntityPtrsMap.find(NodeId)->second;
        items.resize(r_neighbour_entities.size());
        for (IndexType i = 0; i < r_neighbour_entities.size(); ++i) {
            auto p_entity = r_neighbour_entities[i];
            items[i] = std::make_pair(p_entity->Id(), p_entity);
        }
        return items;
    });

    rOutput.merge(entity_id_ptr_map);
}

template<class TEntityPointerType>
void ModelPartUtils::UpdateNodeIdNodePtrMapFromEntityIdEntityPtrMap(
    std::map<IndexType, ModelPart::NodeType::Pointer>& rOutput,
    const std::map<IndexType, TEntityPointerType>& rInput)
{
    // here we use raw ptrs, because there is no need to increment the reference counts.
    std::vector<typename TEntityPointerType::element_type*> entities;
    for (const auto& it : rInput) {
        entities.push_back(&*(it.second));
    }

    auto node_id_ptr_map = block_for_each<ContainerEntityMapReduction<ModelPart::NodeType, ModelPart::NodeType::Pointer>>(entities, [&](auto& pEntity) {
        auto& r_geometry = pEntity->GetGeometry();
        std::vector<std::pair<IndexType, ModelPart::NodeType::Pointer>> items(r_geometry.size());
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            // this gives the intrusive_ptr of the node and increment the reference count.
            auto p_node = r_geometry(i);
            items[i] = std::make_pair(p_node->Id(), p_node);
        }
        return items;
    });

    rOutput.merge(node_id_ptr_map);
}

std::string ModelPartUtils::GetExaminedModelPartsInfo(
    const std::vector<ModelPart*>& rExaminedModelPartsList,
    const bool AreNodesConsidered,
    const bool AreConditionsConsidered,
    const bool AreElementsConsidered,
    const bool AreNeighboursConsidered)
{
    std::stringstream msg;
    msg << "sensitivity model parts [ ";

    for (const auto p_model_part : rExaminedModelPartsList) {
        msg << p_model_part->FullName() << ", ";
    }

    if (*msg.str().rbegin() == ' ') msg.seekp(-1, std::ios_base::end);
    if (*msg.str().rbegin() == ',') msg.seekp(-1, std::ios_base::end);

    msg << " ] with common [ ";
    msg << (AreNodesConsidered ? "nodes, " : "");
    msg << (AreConditionsConsidered ? "conditions, " : "");
    msg << (AreElementsConsidered ? "elements, "  : "");
    msg << (AreNeighboursConsidered ? "parents, " : "");

    if (*msg.str().rbegin() == ' ') msg.seekp(-1, std::ios_base::end);
    if (*msg.str().rbegin() == ',') msg.seekp(-1, std::ios_base::end);

    msg << " ]" << '\0';
    return msg.str();
}

void ModelPartUtils::GetModelParts(
    std::set<ModelPart*>& rOutput,
    ModelPart& rInput)
{
    // add the current model part
    rOutput.emplace(&rInput);

    for (auto& r_sub_model_part : rInput.SubModelParts()) {
        GetModelParts(rOutput, r_sub_model_part);
    }
}

void ModelPartUtils::ExamineModelParts(
    std::set<IndexType>& rExaminedNodeIds,
    std::set<NodeIdsType>& rExaminedConditionGeometryNodeIdsSet,
    std::set<NodeIdsType>& rExaminedElementGeometryNodeIdsSet,
    const std::vector<ModelPart*> rExaminedModelPartsList,
    const bool AreNodesConsidered,
    const bool AreConditionsConsidered,
    const bool AreElementsConsidered,
    const bool AreNeighboursConsidered)
{
    for (auto p_model_part : rExaminedModelPartsList) {
        // first generate examined model part node ids set. Same node can be shared
        // between rExamined model parts and rReferenceModelParts. If they are the
        // same node, then they should have the same id. Hence, only ids are stored.
        if (AreNodesConsidered || AreNeighboursConsidered) UpdateEntityIdsSetFromContainer(rExaminedNodeIds, p_model_part->Nodes());

        // now generate examined model part condition geometry sorted node ids set
        // to identify corresponding entities in reference model parts.
        if (AreConditionsConsidered) UpdateEntityGeometryNodeIdsSetFromContainer(rExaminedConditionGeometryNodeIdsSet, p_model_part->Conditions());

        // now generate examined model part element geometry sorted node ids set
        // to identify corresponding entities in reference model parts.
        if (AreElementsConsidered) UpdateEntityGeometryNodeIdsSetFromContainer(rExaminedElementGeometryNodeIdsSet, p_model_part->Elements());
    }
}

void ModelPartUtils::PopulateModelPart(
    ModelPart& rOutputModelPart,
    ModelPart& rReferenceModelPart,
    const bool AreNodesConsidered,
    const bool AreConditionsConsidered,
    const bool AreElementsConsidered,
    const bool AreNeighboursConsidered,
    const std::set<IndexType>& rExaminedNodeIds,
    const std::set<NodeIdsType>& rExaminedConditionGeometryNodeIdsSet,
    const std::set<NodeIdsType>& rExaminedElementGeometryNodeIdsSet)
{
    // check whether rOutputModelPart is empty.

    KRATOS_ERROR_IF(rOutputModelPart.NumberOfNodes() != 0)
        << rOutputModelPart.FullName() << " should be empty. It contains "
        << rOutputModelPart.NumberOfNodes() << " nodes.\n";

    KRATOS_ERROR_IF(rOutputModelPart.NumberOfConditions() != 0)
        << rOutputModelPart.FullName() << " should be empty. It contains "
        << rOutputModelPart.NumberOfConditions() << " conditions.\n";

    KRATOS_ERROR_IF(rOutputModelPart.NumberOfElements() != 0)
        << rOutputModelPart.FullName() << " should be empty. It contains "
        << rOutputModelPart.NumberOfElements() << " elements.\n";

    // create the necesary maps.
    std::map<IndexType, ModelPart::NodeType::Pointer> node_id_ptr_map;
    std::map<IndexType, ModelPart::ConditionType::Pointer> condition_id_ptr_map;
    std::map<IndexType, ModelPart::ElementType::Pointer> element_id_ptr_map;

    // now get the common nodes between rExaminedModelParts and rReferenceModelParts.
    if (AreNodesConsidered) UpdateEntityIdEntityPtrMapWithCommonEntitiesFromContainerAndEntityIdsSet(node_id_ptr_map, rExaminedNodeIds, rReferenceModelPart.Nodes());

    // now get the common conditions between rExaminedModelParts and rReferenceModelParts by checking geometry node ids
    if (AreConditionsConsidered) UpdateEntityIdEntityPtrMapWithCommonEntitiesFromContainerAndEntityGeometryNodeIdsSet(condition_id_ptr_map, rExaminedConditionGeometryNodeIdsSet, rReferenceModelPart.Conditions());

    // now get the common elements between rExaminedModelParts and rReferenceModelParts by checking geometry node ids
    if (AreElementsConsidered) UpdateEntityIdEntityPtrMapWithCommonEntitiesFromContainerAndEntityGeometryNodeIdsSet(element_id_ptr_map, rExaminedElementGeometryNodeIdsSet, rReferenceModelPart.Elements());

    // now get the neighbour entities for common nodes
    if (AreNeighboursConsidered) {
        // find and update neighbour conditions
        if (AreConditionsConsidered) {
            std::map<IndexType, std::vector<ModelPart::ConditionType::Pointer>> node_id_neighbour_condition_ptrs_map;
            UpdateNeighbourMaps(node_id_neighbour_condition_ptrs_map, rExaminedNodeIds, rReferenceModelPart.Conditions());
            UpdateEntityIdEntityPtrMapFromNeighbourMap(condition_id_ptr_map, node_id_neighbour_condition_ptrs_map);
        }

        // find and update neighbour elements
        if (AreElementsConsidered) {
            std::map<IndexType, std::vector<ModelPart::ElementType::Pointer>> node_id_neighbour_element_ptrs_map;
            UpdateNeighbourMaps(node_id_neighbour_element_ptrs_map, rExaminedNodeIds, rReferenceModelPart.Elements());
            UpdateEntityIdEntityPtrMapFromNeighbourMap(element_id_ptr_map, node_id_neighbour_element_ptrs_map);
        }
    }

    // get nodes from condition_id_ptr_map
    if (AreConditionsConsidered) UpdateNodeIdNodePtrMapFromEntityIdEntityPtrMap(node_id_ptr_map, condition_id_ptr_map);

    // get nodes from element_id_ptr_map
    if (AreElementsConsidered) UpdateNodeIdNodePtrMapFromEntityIdEntityPtrMap(node_id_ptr_map, element_id_ptr_map);

    // now we have to create the communicator for MPI communication.
    Communicator& r_reference_communicator = rReferenceModelPart.GetCommunicator();
    Communicator::Pointer p_output_communicator = r_reference_communicator.Create();
    p_output_communicator->SetNumberOfColors(r_reference_communicator.GetNumberOfColors());
    p_output_communicator->NeighbourIndices() = r_reference_communicator.NeighbourIndices();

    // finally we add all the nodes, conditions and elements from the maps, we don't need to call Unique in
    // here because, we are using a std::map which is sorted with entity ids. We also populate
    // communicator local meshes for conditions because conditions and elements are only present in the
    // local mesh, and not in interface or ghost meshes.
    for (auto& it : condition_id_ptr_map) {
        rOutputModelPart.Conditions().push_back(it.second);
        p_output_communicator->LocalMesh().Conditions().push_back(it.second);
    }

    for (auto& it : element_id_ptr_map) {
        rOutputModelPart.Elements().push_back(it.second);
        p_output_communicator->LocalMesh().Elements().push_back(it.second);
    }

    // now we have to add nodes and create corresponding nodal meshes in the communicator. This has to be done
    // irrespective of AreNodesConsidered true or false because, the nodes of the conditions and elements which was added
    // needs to be properly added and assigned to proper meshes in the communicator. Hence no checks for
    // AreNodesConsidered is done henceforth.

    // now we add all nodes
    for (auto& it : node_id_ptr_map) {
        rOutputModelPart.Nodes().push_back(it.second);
    }

    // get the local mesh nodes to set from rReferenceModelPart
    std::set<IndexType> local_node_ids;
    UpdateEntityIdsSetFromContainer(local_node_ids, r_reference_communicator.LocalMesh().Nodes());

    // get the interface mesh nodes to set from rReferenceModelPart
    std::set<IndexType> interface_node_ids;
    UpdateEntityIdsSetFromContainer(interface_node_ids, r_reference_communicator.InterfaceMesh().Nodes());

    // get the ghost mesh nodes to set from rReferenceModelPart
    std::set<IndexType> ghost_node_ids;
    UpdateEntityIdsSetFromContainer(ghost_node_ids, r_reference_communicator.GhostMesh().Nodes());

    for (auto p_node_itr = rOutputModelPart.Nodes().ptr_begin(); p_node_itr < rOutputModelPart.Nodes().ptr_end(); ++p_node_itr) {
        // populate the local mesh correctly in the output model part
        auto p_local_itr = local_node_ids.find((*p_node_itr)->Id());
        if (p_local_itr != local_node_ids.end()) {
            p_output_communicator->LocalMesh().Nodes().push_back(*p_node_itr);
        }

        // populate the interface mesh correctly in the output model part
        auto p_interface_itr = interface_node_ids.find((*p_node_itr)->Id());
        if (p_interface_itr != interface_node_ids.end()) {
            p_output_communicator->InterfaceMesh().Nodes().push_back(*p_node_itr);
        }

        // populate the ghost mesh correctly in the output model part
        auto p_ghost_itr = ghost_node_ids.find((*p_node_itr)->Id());
        if (p_ghost_itr != ghost_node_ids.end()) {
            p_output_communicator->GhostMesh().Nodes().push_back(*p_node_itr);
        }
    }

    // now set the communicator
    rOutputModelPart.SetCommunicator(p_output_communicator);

    // now set common information
    rOutputModelPart.SetProcessInfo(rReferenceModelPart.pGetProcessInfo());
    rOutputModelPart.PropertiesArray() = rReferenceModelPart.PropertiesArray();
    rOutputModelPart.Tables() = rReferenceModelPart.Tables();
}

std::vector<ModelPart*> ModelPartUtils::GetModelPartsWithCommonReferenceEntities(
    const std::vector<ModelPart*>& rExaminedModelPartsList,
    const std::vector<ModelPart*>& rReferenceModelParts,
    const bool AreNodesConsidered,
    const bool AreConditionsConsidered,
    const bool AreElementsConsidered,
    const bool AreNeighboursConsidered,
    const IndexType EchoLevel)
{
    std::stringstream mp_name_prefix;
    mp_name_prefix << "<OPTIMIZATION_APP_AUTO>"
                   << (AreNodesConsidered ? "_Nodes" : "_NoNodes")
                   << (AreConditionsConsidered ? "_Conditions" : "_NoConditions")
                   << (AreElementsConsidered ? "_Elements" : "_NoElements")
                   << (AreNeighboursConsidered ? "_Parents" : "_NoParents")
                   << "_ExaminedMPs_";

    AppendModelPartNames(mp_name_prefix, rExaminedModelPartsList);

    // generate the unique model part name
    std::string unique_mp_name = mp_name_prefix.str();
    std::replace(unique_mp_name.begin(), unique_mp_name.end(), '.', '>');

    IndexType total_number_of_entities{0};
    bool is_examined_model_parts_processed = false;

    // now create combined examined model part entity information sets.
    std::set<IndexType> examined_node_ids;
    std::set<NodeIdsType> examined_condition_geometry_node_ids;
    std::set<NodeIdsType> examined_element_geometry_node_ids;

    std::vector<ModelPart*> output_model_parts(rReferenceModelParts.size());
    for (IndexType i = 0; i < rReferenceModelParts.size(); ++i) {
        auto& r_reference_model_part = *rReferenceModelParts[i];

        if (!r_reference_model_part.HasSubModelPart(unique_mp_name)) {
            if (!is_examined_model_parts_processed) {
                // we only wants to examine the rExaminedModelPartsList once only for all the rReferenceModelParts.
                is_examined_model_parts_processed = true;

                ExamineModelParts(examined_node_ids, examined_condition_geometry_node_ids,
                                  examined_element_geometry_node_ids, rExaminedModelPartsList,
                                  AreNodesConsidered, AreConditionsConsidered,
                                  AreElementsConsidered, AreNeighboursConsidered);
            }

            // now we create the output sub model part with the unique model part name.
            auto& r_output_model_part = r_reference_model_part.CreateSubModelPart(unique_mp_name);
            PopulateModelPart(r_output_model_part, r_reference_model_part,
                              AreNodesConsidered, AreConditionsConsidered,
                              AreElementsConsidered, AreNeighboursConsidered,
                              examined_node_ids, examined_condition_geometry_node_ids,
                              examined_element_geometry_node_ids);

            KRATOS_INFO_IF("ModelPartUtils", EchoLevel > 0)
                << "Created common entity model part using "
                << r_reference_model_part.FullName() << " for "
                << GetExaminedModelPartsInfo(
                    rExaminedModelPartsList, AreNodesConsidered, AreConditionsConsidered,
                    AreElementsConsidered, AreNeighboursConsidered);
        }

        auto p_model_part = &r_reference_model_part.GetSubModelPart(unique_mp_name);
        output_model_parts[i] = p_model_part;

        KRATOS_INFO_IF("ModelPartUtils", EchoLevel > 1)
            << "Retrieved common entity model part using "
            << r_reference_model_part.FullName() << " for "
            << GetExaminedModelPartsInfo(
                   rExaminedModelPartsList, AreNodesConsidered, AreConditionsConsidered,
                   AreElementsConsidered, AreNeighboursConsidered);

        total_number_of_entities += (AreNodesConsidered ? p_model_part->GetCommunicator().GlobalNumberOfNodes() : 0);
        total_number_of_entities += (AreConditionsConsidered ? p_model_part->GetCommunicator().GlobalNumberOfConditions() : 0);
        total_number_of_entities += (AreElementsConsidered ? p_model_part->GetCommunicator().GlobalNumberOfElements() : 0);
    }

    if (total_number_of_entities == 0) {
        std::stringstream msg;
        msg << "No common entities found between the reference "
            "model parts and examined model parts.";

        msg << "\nFollowings are the reference model parts:";
        for (const auto p_model_part : rReferenceModelParts) {
            msg << "\n\t" << p_model_part->FullName();
        }

        msg << "\nFollowings are the examined model parts:";
        for (const auto p_model_part : rExaminedModelPartsList) {
            msg << "\n\t" << p_model_part->FullName();
        }

        KRATOS_ERROR << msg.str();
    }

    return output_model_parts;
}

void ModelPartUtils::RemoveModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(
    const std::vector<ModelPart*> rModelParts)
{
    std::set<ModelPart*> all_model_parts;
    for (auto p_model_part : rModelParts) {
        GetModelParts(all_model_parts, *p_model_part);
    }

    for (auto& it : all_model_parts) {
        if (it->Name().rfind("<OPTIMIZATION_APP_AUTO>", 0) == 0) {
            it->GetParentModelPart().RemoveSubModelPart(*it);
        }
    }
}

void ModelPartUtils::LogModelPartStatus(
    ModelPart& rModelPart,
    const std::string& rStatus)
{
    if (!rModelPart.Has(MODEL_PART_STATUS)) {
        rModelPart.SetValue(MODEL_PART_STATUS, {});
    }

    auto& r_statuses = rModelPart.GetValue(MODEL_PART_STATUS);
    const auto p_itr = std::find(r_statuses.begin(), r_statuses.end(), rStatus);
    if (p_itr == r_statuses.end()) {
        r_statuses.push_back(rStatus);
    }
}

std::vector<std::string> ModelPartUtils::GetModelPartStatusLog(ModelPart& rModelPart)
{
    if (!rModelPart.Has(MODEL_PART_STATUS)) {
        return {};
    } else {
        return rModelPart.GetValue(MODEL_PART_STATUS);
    }
}

bool ModelPartUtils::CheckModelPartStatus(
    const ModelPart& rModelPart,
    const std::string& rStatus)
{
    if (!rModelPart.Has(MODEL_PART_STATUS)) {
        return false;
    } else {
        const auto& r_statuses = rModelPart.GetValue(MODEL_PART_STATUS);
        const auto p_itr = std::find(r_statuses.begin(), r_statuses.end(), rStatus);
        return p_itr != r_statuses.end();
    }
}

} // namespace Kratos