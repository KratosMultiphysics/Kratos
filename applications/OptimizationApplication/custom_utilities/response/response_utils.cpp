//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl,
//                   Suneth Warnakulasuriya
//

// System includes
#include <map>
#include <sstream>
#include <limits>

// Project includes
#include "includes/communicator.h"
#include "includes/condition.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"

// Application includes

// Include base h
#include "response_utils.h"

namespace Kratos {

std::vector<ModelPart*> ResponseUtils::GetModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList(
    const std::vector<ModelPart*>& rExaminedModelPartsList,
    const std::vector<ModelPart*>& rReferenceModelParts,
    const bool AreNodesConsidered,
    const bool AreConditionsConsidered,
    const bool AreElementsConsidered,
    const bool AreParentsConsidered,
    const IndexType EchoLevel)
{
    KRATOS_TRY

    std::stringstream mp_name_prefix;
    mp_name_prefix << "<OPTIMIZATION_APP_AUTO>"
                   << (AreNodesConsidered
                            ? "_Nodes"
                            : "_NoNodes")
                   << (AreConditionsConsidered
                            ? "_Conditions"
                            : "_NoConditions")
                   << (AreElementsConsidered
                            ? "_Elements"
                            : "_NoElements")
                   << (AreParentsConsidered
                            ? "_Parents"
                            : "_NoParents")
                   << "_SensitivityMPs_";

    // now generate the unique name for model part
    const std::string& unique_mp_name = GetCombinedModelPartsName(mp_name_prefix.str(), rExaminedModelPartsList);

    IndexType total_number_of_entities{0};

    std::vector<ModelPart*> output_model_parts(rReferenceModelParts.size());

    for (IndexType i = 0; i < rReferenceModelParts.size(); ++i) {
        auto& r_reference_model_part = *(rReferenceModelParts[i]);

        // first check whether the required model part exists, if not create it.
        if (!r_reference_model_part.HasSubModelPart(unique_mp_name)) {
            CreateModelPartWithCommonReferenceEntitiesBetweenReferenceAndExamined(unique_mp_name, rExaminedModelPartsList,
                                       r_reference_model_part, AreNodesConsidered,
                                       AreConditionsConsidered, AreElementsConsidered,
                                       AreParentsConsidered, EchoLevel);
        }

        auto p_model_part = &r_reference_model_part.GetSubModelPart(unique_mp_name);
        output_model_parts[i] = p_model_part;

        KRATOS_INFO_IF("ResponseUtils", EchoLevel > 1)
            << "Retrieved sensitivity computation model part using "
            << r_reference_model_part.FullName() << " for "
            << GetSensitivityComputationModelPartsInfo(
                   rExaminedModelPartsList, AreNodesConsidered, AreConditionsConsidered,
                   AreElementsConsidered, AreParentsConsidered);

        total_number_of_entities += (AreNodesConsidered ? p_model_part->GetCommunicator().GlobalNumberOfNodes() : 0);
        total_number_of_entities += (AreConditionsConsidered ? p_model_part->GetCommunicator().GlobalNumberOfConditions() : 0);
        total_number_of_entities += (AreElementsConsidered ? p_model_part->GetCommunicator().GlobalNumberOfElements() : 0);
    }

    if (total_number_of_entities == 0) {
        std::stringstream msg;
        msg << "No common entities found between the reference "
            "model parts and sensitivity model parts for sensitivity computation.";

        msg << "\nFollowings are the reference model parts:";
        for (const auto p_model_part : rReferenceModelParts) {
            msg << "\n\t" << p_model_part->FullName();
        }

        msg << "\nFollowings are the sensitivity model parts:";
        for (const auto p_model_part : rExaminedModelPartsList) {
            msg << "\n\t" << p_model_part->FullName();
        }

        KRATOS_ERROR << msg.str();
    }

    return output_model_parts;

    KRATOS_CATCH("");
}

template<class TEntityType, class TMapValueType>
std::map<IndexType, TMapValueType> ResponseUtils::ContainerEntityMapReduction<TEntityType, TMapValueType>::GetValue() const
{
    return mValue;
}

template<class TEntityType, class TMapValueType>
void ResponseUtils::ContainerEntityMapReduction<TEntityType, TMapValueType>::LocalReduce(const value_type& rValue){
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
void ResponseUtils::ContainerEntityMapReduction<TEntityType, TMapValueType>::ThreadSafeReduce(ContainerEntityMapReduction<TEntityType, TMapValueType>& rOther)
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

std::string ResponseUtils::GetCombinedModelPartsName(
    const std::string& rPrefix,
    const std::vector<ModelPart*>& rModelParts)
{
    std::vector<std::string> mp_names_list(rModelParts.size());
    for (IndexType i = 0; i < rModelParts.size(); ++i) {
        mp_names_list[i] = rModelParts[i]->FullName();
    }
    std::sort(mp_names_list.begin(), mp_names_list.end());

    std::stringstream msg;
    msg << rPrefix;
    for (const auto& r_mp_name : mp_names_list) {
        msg << r_mp_name << ";";
    }

    std::string name = msg.str();
    std::replace(name.begin(), name.end(), '.', '>');

    return name;
}

std::string ResponseUtils::GetSensitivityComputationModelPartsInfo(
    const std::vector<ModelPart*>& rExaminedModelPartsList,
    const bool AreNodesConsidered,
    const bool AreConditionsConsidered,
    const bool AreElementsConsidered,
    const bool AreParentsConsidered)
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
    msg << (AreParentsConsidered ? "parents, " : "");

    if (*msg.str().rbegin() == ' ') msg.seekp(-1, std::ios_base::end);
    if (*msg.str().rbegin() == ',') msg.seekp(-1, std::ios_base::end);

    msg << " ]" << '\0';
    return msg.str();
}

template<class TContainerType>
void ResponseUtils::AddNeighbourEntitiesToFlaggedNodes(
    std::map<IndexType, std::vector<ContainerEntityPointerType<TContainerType>>>& rOutput,
    TContainerType& rContainer,
    const Flags& rFlag,
    const bool FlagValue)
{
    using reduction_type = ContainerEntityMapReduction<typename TContainerType::value_type, std::vector<ContainerEntityPointerType<TContainerType>>>;

    // need to use ptr_iterator here because, we need to increment the reference counter of the entity intrusive_ptr when push_back is used.
    auto node_id_neighbour_ptrs_map = BlockPartition<TContainerType, typename TContainerType::ptr_iterator>(rContainer.ptr_begin(), rContainer.ptr_end()).template for_each<reduction_type>([&](auto& pEntity) {
        std::vector<std::pair<IndexType, ContainerEntityPointerType<TContainerType>>> items;
        for (const auto& r_node : pEntity->GetGeometry()) {
            if (r_node.Is(rFlag) == FlagValue) {
                items.push_back(std::make_pair(r_node.Id(), pEntity));
            }
        }
        return items;
    });

    rOutput.merge(node_id_neighbour_ptrs_map);
}

template<class TEntityPointerType>
void ResponseUtils::UpdateEntityIdEntityPtrMapFromNodalNeighbourEntities(
    std::map<IndexType, TEntityPointerType>& rOutput,
    const std::map<IndexType, std::vector<TEntityPointerType>>& rNodeIdNeighbourEntityPtrsMap,
    const ModelPart::NodesContainerType& rNodes)
{
    auto entity_id_ptr_map = block_for_each<ContainerEntityMapReduction<typename TEntityPointerType::element_type, TEntityPointerType>>(rNodes, [&](auto& rNode) {
        std::vector<std::pair<IndexType, TEntityPointerType>> items;
        auto itr_item = rNodeIdNeighbourEntityPtrsMap.find(rNode.Id());
        if (itr_item != rNodeIdNeighbourEntityPtrsMap.end()) {
            const auto& r_neighbour_entities = itr_item->second;
            items.resize(r_neighbour_entities.size());
            for (IndexType i = 0; i < r_neighbour_entities.size(); ++i) {
                // the r_neighbour_entities(i) gives the intrusive_ptr which is incrementing the reference counter.
                auto p_entity = r_neighbour_entities[i];
                items[i] = std::make_pair(p_entity->Id(), p_entity);
            }
        }
        return items;
    });

    rOutput.merge(entity_id_ptr_map);
}

template<class TContainerType>
void ResponseUtils::UpdateEntityIdEntityPtrMapFromEntityContainer(
    std::map<IndexType, ContainerEntityPointerType<TContainerType>>& rOutput,
    TContainerType& rContainer)
{
    using map_type = std::map<IndexType, ContainerEntityPointerType<TContainerType>>;

    using reduction_type = MapReduction<map_type>;

    // need to use ptr_iterator here because, we need to increment the reference counter of the entity intrusive_ptr when push_back is used.
    auto entity_id_ptr_map = BlockPartition<TContainerType, typename TContainerType::ptr_iterator>(rContainer.ptr_begin(), rContainer.ptr_end()).template for_each<reduction_type>([&](auto pEntity) {
        return std::make_pair(pEntity->Id(), pEntity);
    });

    rOutput.merge(entity_id_ptr_map);
}

template<class TContainerType>
void ResponseUtils::UpdateNodeIdsEntityPtrMapFromEntityContainer(
    std::map<NodeIdsType, ContainerEntityPointerType<TContainerType>>& rOutput,
    TContainerType& rContainer)
{
    using map_type = std::map<NodeIdsType, ContainerEntityPointerType<TContainerType>>;

    using reduction_type = MapReduction<map_type>;

    // need to use ptr_iterator here because, we need to increment the reference counter of the entity intrusive_ptr when push_back is used.
    auto node_ids_entity_ptr_map = BlockPartition<TContainerType, typename TContainerType::ptr_iterator>(rContainer.ptr_begin(), rContainer.ptr_end()).template for_each<reduction_type>([&](auto pEntity) {
        const auto& r_geometry = pEntity->GetGeometry();
        NodeIdsType node_ids(r_geometry.size());
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            node_ids[i] = r_geometry[i].Id();
        }
        std::sort(node_ids.begin(), node_ids.end());
        return std::make_pair(node_ids, pEntity);
    });

    rOutput.merge(node_ids_entity_ptr_map);
}

template<class TContainerType>
void ResponseUtils::UpdateEntityIdEntityPtrMapFromNodeIdsEntityPtrMapAndEntityContainer(
    std::map<IndexType, ContainerEntityPointerType<TContainerType>>& rOutput,
    const std::map<NodeIdsType, ContainerEntityPointerType<TContainerType>>& rNodeIdsEntityPtrMap,
    const TContainerType& rContainer)
{
    using map_type = std::map<IndexType, ContainerEntityPointerType<TContainerType>>;

    using reduction_type = MapReduction<map_type>;

    auto entity_id_ptr_map = block_for_each<reduction_type>(rContainer, [&](auto& rEntity) {
        const auto& r_geometry = rEntity.GetGeometry();
        NodeIdsType node_ids(r_geometry.size());
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            node_ids[i] = r_geometry[i].Id();
        }
        std::sort(node_ids.begin(), node_ids.end());

        const auto p_it = rNodeIdsEntityPtrMap.find(node_ids);
        if (p_it != rNodeIdsEntityPtrMap.end()) {
            return std::make_pair(p_it->second->Id(), p_it->second);
        } else {
            KRATOS_ERROR << "Sensitivity entity with id " << rEntity.Id() << " not found in node ids entity map.";
        }

        // following line will not be executed in runtime.
        return std::make_pair(rEntity.Id(), p_it->second);
    });

    rOutput.merge(entity_id_ptr_map);
}

template<class TContainerType>
void ResponseUtils::UpdateEntityIdEntityPtrMapFromFlaggedEntityContainer(
    std::map<IndexType, ContainerEntityPointerType<TContainerType>>& rOutput,
    TContainerType& rContainer,
    const Flags& rFlag,
    const bool FlagValue)
{
    using map_type = std::map<IndexType, ContainerEntityPointerType<TContainerType>>;

    using reduction_type = MapReduction<map_type>;

    // need to use ptr_iterator here because, we need to increment the reference counter of the entity intrusive_ptr when push_back is used.
    auto entity_id_ptr_map = BlockPartition<TContainerType, typename TContainerType::ptr_iterator>(rContainer.ptr_begin(), rContainer.ptr_end()).template for_each<reduction_type>([&](auto pEntity) {
        if (pEntity->Is(rFlag) == FlagValue) {
            return std::make_pair(pEntity->Id(), pEntity);
        } else {
            // this is a dummy return which will removed eventually.
            return std::make_pair(std::numeric_limits<IndexType>::max(), pEntity);
        }
    });

    // remove the unwanted dummy entry of max
    if (entity_id_ptr_map.find(std::numeric_limits<IndexType>::max()) != entity_id_ptr_map.end()) {
        entity_id_ptr_map.erase(std::numeric_limits<IndexType>::max());
    }

    rOutput.merge(entity_id_ptr_map);
}

template<class TEntityPointerType>
void ResponseUtils::UpdateNodeIdNodePtrMapFromEntityIdEntityPtrMap(
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

void ResponseUtils::CreateModelPartWithCommonReferenceEntitiesBetweenReferenceAndExamined(
    const std::string& rOutputModelPartName,
    const std::vector<ModelPart*>& rExaminedModelPartsList,
    ModelPart& rReferenceModelPart,
    const bool AreNodesConsidered,
    const bool AreConditionsConsidered,
    const bool AreElementsConsidered,
    const bool AreParentsConsidered,
    const IndexType EchoLevel)
{
    KRATOS_TRY

    std::map<IndexType, ModelPart::NodeType::Pointer> node_id_ptr_map;
    std::map<IndexType, Condition::Pointer> condition_id_ptr_map;
    std::map<IndexType, Element::Pointer> element_id_ptr_map;

    if (AreNodesConsidered || AreParentsConsidered) {
        // clear flags in analysis model part
        VariableUtils().SetFlag(SELECTED, false, rReferenceModelPart.Nodes());

        // set flags for sensitivity model part nodes
        for (const auto p_model_part : rExaminedModelPartsList) {
            VariableUtils().SetFlag(SELECTED, true, p_model_part->Nodes());
        }
    }

    // update the common nodes from sensitivity model part and analysis model part
    if (AreNodesConsidered) UpdateEntityIdEntityPtrMapFromFlaggedEntityContainer(node_id_ptr_map, rReferenceModelPart.Nodes(), SELECTED);

    // first generate the analysis model part node_ids and entity ptr maps for later map sensitivity model part
    // entities with analysis model part
    std::map<NodeIdsType, Condition::Pointer> analysis_mp_node_ids_condition_ptr_map;
    std::map<NodeIdsType, Element::Pointer> analysis_mp_node_ids_element_ptr_map;

    // generate the analysis mp node ids and ptrs maps for conditions
    if (AreConditionsConsidered) UpdateNodeIdsEntityPtrMapFromEntityContainer(analysis_mp_node_ids_condition_ptr_map, rReferenceModelPart.Conditions());

    // generate the analysis mp node ids and ptrs maps for elements
    if (AreElementsConsidered) UpdateNodeIdsEntityPtrMapFromEntityContainer(analysis_mp_node_ids_element_ptr_map, rReferenceModelPart.Elements());

    for (const auto p_model_part : rExaminedModelPartsList) {
        // now we have to match node ids of each sensitivity model part conditions and add them to map
        if (AreConditionsConsidered) UpdateEntityIdEntityPtrMapFromNodeIdsEntityPtrMapAndEntityContainer(condition_id_ptr_map, analysis_mp_node_ids_condition_ptr_map, p_model_part->Conditions());

        // now we have to match node ids of each sensitivity model part elements and add them to map
        if (AreElementsConsidered) UpdateEntityIdEntityPtrMapFromNodeIdsEntityPtrMapAndEntityContainer(element_id_ptr_map, analysis_mp_node_ids_element_ptr_map, p_model_part->Elements());
    }

    if (AreParentsConsidered) {
        // create the maps to hold neighbours
        std::map<IndexType, std::vector<Condition::Pointer>> node_id_neighbour_condition_ptrs_map;
        std::map<IndexType, std::vector<Element::Pointer>> node_id_neighbour_element_ptrs_map;

        // now populate neighbour conditions for on sensitivity model parts' nodes from the analysis model part
        if (AreConditionsConsidered) AddNeighbourEntitiesToFlaggedNodes(node_id_neighbour_condition_ptrs_map, rReferenceModelPart.Conditions(), SELECTED);

        // now populate neighbour elements for on sensitivity model parts' nodes from the analysis model part
        if (AreElementsConsidered) AddNeighbourEntitiesToFlaggedNodes(node_id_neighbour_element_ptrs_map, rReferenceModelPart.Elements(), SELECTED);

        // now add the parent elements/conditions which are from the analysis model parts
        // we only iterate through nodal parent elements and nodal parent conditions
        // since this covers parent elements of conditions as well.
        for (const auto p_model_part : rExaminedModelPartsList) {
            // we update the map with condition ids and condition pointers in parallel
            if (AreConditionsConsidered) UpdateEntityIdEntityPtrMapFromNodalNeighbourEntities(condition_id_ptr_map, node_id_neighbour_condition_ptrs_map, p_model_part->Nodes());

            // we update the map with element ids and element pointers in parallel
            if (AreElementsConsidered) UpdateEntityIdEntityPtrMapFromNodalNeighbourEntities(element_id_ptr_map, node_id_neighbour_element_ptrs_map, p_model_part->Nodes());
        }
    }

    // now we create the submodel part in the reference model part because, this
    // model part contains nodes, conditions, elements from reference model part
    // only.
    auto& model_part = rReferenceModelPart.CreateSubModelPart(rOutputModelPartName);

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
        model_part.Conditions().push_back(it.second);
        p_output_communicator->LocalMesh().Conditions().push_back(it.second);
    }

    for (auto& it : element_id_ptr_map) {
        model_part.Elements().push_back(it.second);
        p_output_communicator->LocalMesh().Elements().push_back(it.second);
    }

    // now we have to add nodes and create corresponding nodal meshes in the communicator. This has to be done
    // irrespective of AreNodesConsidered true or false because, the nodes of the conditions and elements which was added
    // needs to be properly added and assigned to proper meshes in the communicator. Hence no checks for
    // AreNodesConsidered is done henceforth.

    // now we add all nodes
    for (auto& it : node_id_ptr_map) {
        model_part.Nodes().push_back(it.second);
    }

    // get the local mesh nodes to map from rReferenceModelPart
    std::map<IndexType, ModelPart::NodeType::Pointer> local_node_id_ptr_map;
    UpdateEntityIdEntityPtrMapFromEntityContainer(local_node_id_ptr_map, r_reference_communicator.LocalMesh().Nodes());

    // get the interface mesh nodes to map from rReferenceModelPart
    std::map<IndexType, ModelPart::NodeType::Pointer> interface_node_id_ptr_map;
    UpdateEntityIdEntityPtrMapFromEntityContainer(interface_node_id_ptr_map, r_reference_communicator.InterfaceMesh().Nodes());

    // get the ghost mesh nodes to map from rReferenceModelPart
    std::map<IndexType, ModelPart::NodeType::Pointer> ghost_node_id_ptr_map;
    UpdateEntityIdEntityPtrMapFromEntityContainer(ghost_node_id_ptr_map, r_reference_communicator.GhostMesh().Nodes());

    for (auto& r_node : model_part.Nodes()) {
        // populate the local mesh correctly in the output model part
        auto p_local_itr = local_node_id_ptr_map.find(r_node.Id());
        if (p_local_itr != local_node_id_ptr_map.end()) {
            p_output_communicator->LocalMesh().Nodes().push_back(p_local_itr->second);
        }

        // populate the interface mesh correctly in the output model part
        auto p_interface_itr = interface_node_id_ptr_map.find(r_node.Id());
        if (p_interface_itr != interface_node_id_ptr_map.end()) {
            p_output_communicator->InterfaceMesh().Nodes().push_back(p_interface_itr->second);
        }

        // populate the ghost mesh correctly in the output model part
        auto p_ghost_itr = ghost_node_id_ptr_map.find(r_node.Id());
        if (p_ghost_itr != ghost_node_id_ptr_map.end()) {
            p_output_communicator->GhostMesh().Nodes().push_back(p_ghost_itr->second);
        }
    }

    // now set the communicator
    model_part.SetCommunicator(p_output_communicator);

    // now set common information
    model_part.SetProcessInfo(rReferenceModelPart.pGetProcessInfo());
    model_part.PropertiesArray() = rReferenceModelPart.PropertiesArray();
    model_part.Tables() = rReferenceModelPart.Tables();

    KRATOS_INFO_IF("ResponseUtils", EchoLevel > 0)
        << "Created sensitivity computation model part using "
        << rReferenceModelPart.FullName() << " for "
        << GetSensitivityComputationModelPartsInfo(
               rExaminedModelPartsList, AreNodesConsidered, AreConditionsConsidered,
               AreElementsConsidered, AreParentsConsidered);

    KRATOS_CATCH("");
}

} // namespace Kratos