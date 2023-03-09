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
#include "containers/global_pointers_vector.h"
#include "containers/model.h"
#include "includes/condition.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"

// Application includes

// Include base h
#include "response_utils.h"

namespace Kratos {

template<class EntityType>
std::map<IndexType, EntityType*> ResponseUtils::ContainerEntityMapReduction<EntityType>::GetValue() const
{
    return mValue;
}

template<class EntityType>
void ResponseUtils::ContainerEntityMapReduction<EntityType>::LocalReduce(const value_type& rValue){
    for (const auto& r_item : rValue) {
        mValue.emplace(r_item);
    }
}

template<class EntityType>
void ResponseUtils::ContainerEntityMapReduction<EntityType>::ThreadSafeReduce(ContainerEntityMapReduction<EntityType>& rOther)
{
    KRATOS_CRITICAL_SECTION
    mValue.merge(rOther.mValue);
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

template<class TContainerType>
void ResponseUtils::AddNeighbourEntitiesToFlaggedNodes(
    TContainerType& rContainer,
    const Variable<GlobalPointersVector<typename TContainerType::value_type>>& rNeighbourEntitiesOutputVariable,
    const Flags& rFlag,
    const bool FlagValue)
{
    block_for_each(rContainer, [&](auto& rEntity) {
        for (auto& r_node : rEntity.GetGeometry()) {
            if (r_node.Is(rFlag) == FlagValue) {
                r_node.SetLock();
                auto& r_neighbour_entities_vector = r_node.GetValue(rNeighbourEntitiesOutputVariable);
                r_neighbour_entities_vector.push_back(&rEntity);
                r_node.UnSetLock();
            }
        }
    });
}

template<class TEntityType>
void ResponseUtils::UpdateEntityIdEntityPtrMapFromNodalNeighbourEntities(
    std::map<IndexType, TEntityType*>& rOutput,
    const ModelPart::NodesContainerType& rNodes,
    const Variable<GlobalPointersVector<TEntityType>>& rNeighbourEntitiesVariable)
{
    auto entity_id_ptr_map = block_for_each<ContainerEntityMapReduction<TEntityType>>(rNodes, [&](auto& rNode) {
        auto& r_neighbour_entities = rNode.GetValue(rNeighbourEntitiesVariable);
        std::vector<std::pair<IndexType, TEntityType*>> items(r_neighbour_entities.size());
        for (IndexType i = 0; i < r_neighbour_entities.size(); ++i) {
            auto& r_entity = r_neighbour_entities[i];
            items[i] = std::make_pair(r_entity.Id(), &r_entity);
        }
        return items;
    });

    rOutput.merge(entity_id_ptr_map);
}

template<class TContainerType>
void ResponseUtils::UpdateEntityIdEntityPtrMapFromEntityContainer(
    std::map<IndexType, typename TContainerType::value_type*>& rOutput,
    TContainerType& rContainer)
{
    auto entity_id_ptr_map = block_for_each<MapReduction<std::map<IndexType, typename TContainerType::value_type*>>>(rContainer, [](auto& rEntity) {
        return std::make_pair(rEntity.Id(), &rEntity);
    });

    rOutput.merge(entity_id_ptr_map);
}

template<class TContainerType>
void ResponseUtils::UpdateNodeIdsEntityPtrMapFromEntityContainer(
    std::map<NodeIdsType, typename TContainerType::value_type*>& rOutput,
    TContainerType& rContainer)
{
    auto node_ids_entity_ptr_map = block_for_each<MapReduction<std::map<NodeIdsType, typename TContainerType::value_type*>>>(rContainer, [](auto& rEntity) {
        const auto& r_geometry = rEntity.GetGeometry();
        NodeIdsType node_ids(r_geometry.size());
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            node_ids[i] = r_geometry[i].Id();
        }
        std::sort(node_ids.begin(), node_ids.end());
        return std::make_pair(node_ids, &rEntity);
    });

    rOutput.merge(node_ids_entity_ptr_map);
}

template<class TContainerType>
void ResponseUtils::UpdateEntityIdEntityPtrMapFromNodeIdsEntityPtrMapAndEntityContainer(
    std::map<IndexType, typename TContainerType::value_type*>& rOutput,
    const std::map<NodeIdsType, typename TContainerType::value_type*>& rNodeIdsEntityPtrMap,
    const TContainerType& rContainer)
{
    auto entity_id_ptr_map = block_for_each<MapReduction<std::map<IndexType, typename TContainerType::value_type*>>>(rContainer, [&](auto& rEntity) {
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

        return std::make_pair(rEntity.Id(), &rEntity);
    });

    rOutput.merge(entity_id_ptr_map);
}

template<class TContainerType>
void ResponseUtils::UpdateEntityIdEntityPtrMapFromFlaggedEntityContainer(
    std::map<IndexType, typename TContainerType::value_type*>& rOutput,
    TContainerType& rContainer,
    const Flags& rFlag,
    const bool FlagValue)
{
    auto entity_id_ptr_map = block_for_each<MapReduction<std::map<IndexType, typename TContainerType::value_type*>>>(rContainer, [&](auto& rEntity) {
        if (rEntity.Is(rFlag) == FlagValue) {
            return std::make_pair(rEntity.Id(), &rEntity);
        } else {
            return std::make_pair(std::numeric_limits<IndexType>::max(), &rEntity);
        }
    });

    // remove the unwanted entry of max
    if (entity_id_ptr_map.find(std::numeric_limits<IndexType>::max()) != entity_id_ptr_map.end()) {
        entity_id_ptr_map.erase(std::numeric_limits<IndexType>::max());
    }

    rOutput.merge(entity_id_ptr_map);
}

template<class TEntityType>
void ResponseUtils::UpdateNodeIdNodePtrMapFromEntityIdEntityPtrMap(
    std::map<IndexType, ModelPart::NodeType*>& rOutput,
    const std::map<IndexType, TEntityType*>& rInput)
{
    std::vector<TEntityType*> entities;
    for (const auto& it : rInput) {
        entities.push_back(it.second);
    }

    auto node_id_ptr_map = block_for_each<ContainerEntityMapReduction<ModelPart::NodeType>>(entities, [&](auto pEntity) {
        auto& r_geometry = pEntity->GetGeometry();
        std::vector<std::pair<IndexType, ModelPart::NodeType*>> items(r_geometry.size());
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            auto& r_node = r_geometry[i];
            items[i] = std::make_pair(r_node.Id(), &r_node);
        }
        return items;
    });

    rOutput.merge(node_id_ptr_map);
}

ModelPart& ResponseUtils::GetSensitivityModelPartForAdjointSensitivities(
    const std::vector<ModelPart*>& rSensitivityModelParts,
    ModelPart& rAnalysisModelPart,
    const bool AreSensitivityEntityParentsConsidered,
    const bool AreSensitivityEntitesConsidered,
    const bool ForceFindSensitivityEntitiesInAnalysisModelPart,
    const IndexType EchoLevel)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rSensitivityModelParts.size() == 0) << "No sensitivity model parts were provided.\n";

    std::stringstream mp_name_prefix;
    mp_name_prefix << "<OPTIMIZATION_APP_AUTO>_AnalysisMP_" << rAnalysisModelPart.FullName()
                   << (AreSensitivityEntityParentsConsidered
                            ? "_SensParentEntities_"
                            : "_NoSensParentEntities_")
                   << (AreSensitivityEntitesConsidered
                            ? "_SensEntities_"
                            : "_NoSensEntities_")
                   << (ForceFindSensitivityEntitiesInAnalysisModelPart
                            ? "_ForceFindInAnalysisMP_"
                            : "_NoForceFindInAnalysisMP_")
                   << "SensitivityMPs_";

    // now generate the unique name for model part
    const std::string& unique_mp_name = GetCombinedModelPartsName(mp_name_prefix.str(), rSensitivityModelParts);

    auto& r_model = rSensitivityModelParts[0]->GetModel();

    if (!r_model.HasModelPart(unique_mp_name)) {
        std::map<IndexType, Condition*> condition_id_ptr_map;
        std::map<IndexType, Element*> element_id_ptr_map;

        if (AreSensitivityEntityParentsConsidered) {
            // clear flags in analysis model part
            VariableUtils().SetFlag(SELECTED, false, rAnalysisModelPart.Nodes());

            // clear neighbours on nodes of sensitivity model parts
            for (const auto p_model_part : rSensitivityModelParts) {
                block_for_each(p_model_part->Nodes(), [](auto& rNode) {
                    rNode.SetValue(NEIGHBOUR_ELEMENTS, GlobalPointersVector<Element>());
                    rNode.SetValue(NEIGHBOUR_CONDITIONS, GlobalPointersVector<Condition>());
                    rNode.Set(SELECTED, true);
                });
            }

            // now populate neighbour elements for on sensitivity model parts' nodes from the analysis model part
            AddNeighbourEntitiesToFlaggedNodes(rAnalysisModelPart.Elements(), NEIGHBOUR_ELEMENTS, SELECTED);

            // now populate neighbour conditions for on sensitivity model parts' nodes from the analysis model part
            AddNeighbourEntitiesToFlaggedNodes(rAnalysisModelPart.Conditions(), NEIGHBOUR_CONDITIONS, SELECTED);

            // now add the parent elements/conditions which are from the analysis model parts
            // we only iterate through nodal parent elements and nodal parent conditions
            // since this covers parent elements of conditions as well.
            for (const auto p_model_part : rSensitivityModelParts) {
                // we update the map with condition ids and condition pointers in parallel
                UpdateEntityIdEntityPtrMapFromNodalNeighbourEntities(condition_id_ptr_map, p_model_part->Nodes(), NEIGHBOUR_CONDITIONS);

                // we update the map with element ids and element pointers in parallel
                UpdateEntityIdEntityPtrMapFromNodalNeighbourEntities(element_id_ptr_map, p_model_part->Nodes(), NEIGHBOUR_ELEMENTS);
            }
        }

        if (AreSensitivityEntitesConsidered) {
            bool is_analysis_mp_entity_ids_ptrs_maps_generted = false;

            std::map<NodeIdsType, Condition*> analysis_mp_node_ids_condition_ptr_map;
            std::map<NodeIdsType, Element*> analysis_mp_node_ids_element_ptr_map;

            for (const auto p_model_part : rSensitivityModelParts) {
                if (!ForceFindSensitivityEntitiesInAnalysisModelPart && (&(p_model_part->GetRootModelPart()) == &(rAnalysisModelPart.GetRootModelPart()))) {
                    // we update the map with condition ids and condition pointers in parallel
                    UpdateEntityIdEntityPtrMapFromEntityContainer(condition_id_ptr_map, p_model_part->Conditions());

                    // we update the map with element ids and element pointers in parallel
                    UpdateEntityIdEntityPtrMapFromEntityContainer(element_id_ptr_map, p_model_part->Elements());
                }  else {
                    // either ForceFindSensitivityEntitiesInAnalysisModelPart is true or the root model parts does not match.
                    // then we need to find for each sensitivity entity, the corresponding analysis model part entity.

                    if (!is_analysis_mp_entity_ids_ptrs_maps_generted) {
                        is_analysis_mp_entity_ids_ptrs_maps_generted = true;

                        // generate the analysis mp node ids and ptrs maps for conditions
                        UpdateNodeIdsEntityPtrMapFromEntityContainer(analysis_mp_node_ids_condition_ptr_map, rAnalysisModelPart.Conditions());

                        // generate the analysis mp node ids and ptrs maps for elements
                        UpdateNodeIdsEntityPtrMapFromEntityContainer(analysis_mp_node_ids_element_ptr_map, rAnalysisModelPart.Elements());
                    }

                    // now we have to match node ids of each sensitivity model part conditions and add them to map
                    UpdateEntityIdEntityPtrMapFromNodeIdsEntityPtrMapAndEntityContainer(condition_id_ptr_map, analysis_mp_node_ids_condition_ptr_map, p_model_part->Conditions());

                    // now we have to match node ids of each sensitivity model part elements and add them to map
                    UpdateEntityIdEntityPtrMapFromNodeIdsEntityPtrMapAndEntityContainer(element_id_ptr_map, analysis_mp_node_ids_element_ptr_map, p_model_part->Elements());
                }
            }
        }

        auto& model_part = r_model.CreateModelPart(unique_mp_name);

        std::map<IndexType, ModelPart::NodeType*> node_id_ptr_map;
        // get nodes from condition_id_ptr_map
        UpdateNodeIdNodePtrMapFromEntityIdEntityPtrMap(node_id_ptr_map, condition_id_ptr_map);

        // get nodes from element_id_ptr_map
        UpdateNodeIdNodePtrMapFromEntityIdEntityPtrMap(node_id_ptr_map, element_id_ptr_map);

        // finally we add all the nodes, conditions and elements from the maps, we don't need to call Unique in
        // here because, we are using a std::map which is sorted with entity ids.
        for (auto& it : node_id_ptr_map) {
            model_part.Nodes().push_back(it.second);
        }

        for (auto& it : condition_id_ptr_map) {
            model_part.Conditions().push_back(it.second);
        }

        for (auto& it : element_id_ptr_map) {
            model_part.Elements().push_back(it.second);
        }

        // check whether there is an overlap.
        // in adjoint based sensitivity analysis, we always require residuals
        // hence, there should be always at least one element to compute
        // residuals. Having conditions is optional.
        if (model_part.GetCommunicator().GlobalNumberOfElements() == 0) {

            std::stringstream msg;
            msg << "No common elements found between " << rAnalysisModelPart.FullName()
                << " and sensitivity model parts to be used with adjoint "
                   "sensitivity computation.";

            msg << "\nFollowings are the sensitivity model parts:";
            for (const auto p_model_part : rSensitivityModelParts) {
                msg << "\n\t" << p_model_part->FullName();
            }

            KRATOS_ERROR << msg.str();
        }

        if (EchoLevel > 0) {
            std::stringstream msg;
            msg << "Created sensitivity computation model part based on " << rAnalysisModelPart.FullName();
            msg << " for sensitivity model parts [ ";
            for (const auto p_model_part : rSensitivityModelParts) {
                msg << p_model_part->FullName() << " ";
            }
            msg << "]";
            msg << (AreSensitivityEntityParentsConsidered ? " with parents" : "");
            msg << (AreSensitivityEntitesConsidered ? " with sensitivity entities" : "");
            msg << ".";
            KRATOS_INFO("ResponseUtils") << msg.str() << std::endl;
        }
        return model_part;
    } else {
        if (EchoLevel > 1) {
            std::stringstream msg;
            msg << "Retrieved sensitivity computation model part based on " << rAnalysisModelPart.FullName();
            msg << " for sensitivity model parts [ ";
            for (const auto p_model_part : rSensitivityModelParts) {
                msg << p_model_part->FullName() << " ";
            }
            msg << "]";
            msg << (AreSensitivityEntityParentsConsidered ? " with parents" : "");
            msg << (AreSensitivityEntitesConsidered ? " with sensitivity entities" : "");
            msg << ".";
            KRATOS_INFO("ResponseUtils") << msg.str() << std::endl;
        }
        return r_model.GetModelPart(unique_mp_name);
    }

    KRATOS_CATCH("");
}

ModelPart& ResponseUtils::GetSensitivityModelPartForDirectSensitivities(
    const std::vector<ModelPart*>& rSensitivityModelParts,
    const std::vector<ModelPart*>& rEvaluatedModelParts,
    const bool AreNodesConsidered,
    const bool AreConditionsConsidered,
    const bool AreElementsConsidered,
    const IndexType EchoLevel)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rSensitivityModelParts.size() == 0) << "No sensitivity model parts were provided.\n";
    auto& r_model = rSensitivityModelParts[0]->GetModel();

    // now generate the unique name for model part
    std::stringstream mp_name_prefix;
    mp_name_prefix << GetCombinedModelPartsName("<OPTIMIZATION_APP_AUTO>_EvaluatedMPs_", rEvaluatedModelParts)
                   << (AreNodesConsidered
                            ? "_Nodes_"
                            : "_NoNodes_")
                   << (AreConditionsConsidered
                           ? "_Conditions_"
                           : "_NoConditions_")
                   << (AreElementsConsidered
                           ? "_Elements_"
                           : "_NoElements_")
                   << GetCombinedModelPartsName("SensitivityMPs_", rSensitivityModelParts);
    const std::string unique_mp_name = mp_name_prefix.str();

    if (!r_model.HasModelPart(unique_mp_name)) {

        // reset all entities
        for (const auto p_model_part : rSensitivityModelParts) {
            if (AreNodesConsidered) VariableUtils().SetFlag(SELECTED, false, p_model_part->Nodes());
            if (AreConditionsConsidered) VariableUtils().SetFlag(SELECTED, false, p_model_part->Conditions());
            if (AreElementsConsidered) VariableUtils().SetFlag(SELECTED, false, p_model_part->Elements());
        }

        // select evaluated model part quantities
        for (const auto p_model_part : rEvaluatedModelParts) {
            if (AreNodesConsidered) VariableUtils().SetFlag(SELECTED, true, p_model_part->Nodes());
            if (AreConditionsConsidered) VariableUtils().SetFlag(SELECTED, true, p_model_part->Conditions());
            if (AreElementsConsidered) VariableUtils().SetFlag(SELECTED, true, p_model_part->Elements());
        }

        std::map<IndexType, ModelPart::NodeType*> node_id_ptr_map;
        std::map<IndexType, ModelPart::ConditionType*> condition_id_ptr_map;
        std::map<IndexType, ModelPart::ElementType*> element_id_ptr_map;

        for (const auto p_model_part : rSensitivityModelParts) {
            if (AreNodesConsidered) UpdateEntityIdEntityPtrMapFromFlaggedEntityContainer(node_id_ptr_map, p_model_part->Nodes(), SELECTED);
            if (AreConditionsConsidered) UpdateEntityIdEntityPtrMapFromFlaggedEntityContainer(condition_id_ptr_map, p_model_part->Conditions(), SELECTED);
            if (AreElementsConsidered) UpdateEntityIdEntityPtrMapFromFlaggedEntityContainer(element_id_ptr_map, p_model_part->Elements(), SELECTED);
        }

        auto& model_part = r_model.CreateModelPart(unique_mp_name);

        // get nodes from condition_id_ptr_map
        UpdateNodeIdNodePtrMapFromEntityIdEntityPtrMap(node_id_ptr_map, condition_id_ptr_map);

        // get nodes from element_id_ptr_map
        UpdateNodeIdNodePtrMapFromEntityIdEntityPtrMap(node_id_ptr_map, element_id_ptr_map);

        // finally we add all the nodes, conditions and elements from the maps, we don't need to call Unique in
        // here because, we are using a std::map which is sorted with entity ids.
        for (auto& it : node_id_ptr_map) {
            model_part.Nodes().push_back(it.second);
        }

        for (auto& it : condition_id_ptr_map) {
            model_part.Conditions().push_back(it.second);
        }

        for (auto& it : element_id_ptr_map) {
            model_part.Elements().push_back(it.second);
        }

        // check whether there is an overlap.
        if (model_part.GetCommunicator().GlobalNumberOfNodes() == 0 &&
            model_part.GetCommunicator().GlobalNumberOfConditions() == 0 &&
            model_part.GetCommunicator().GlobalNumberOfElements() == 0) {

            std::stringstream msg;
            msg << "No common entities found between the evaluated "
                   "model parts and sensitivity model parts for direct "
                   "sensitivity computation.";

            msg << "\nFollowings are the evaluated model parts:";
            for (const auto p_model_part : rEvaluatedModelParts) {
                msg << "\n\t" << p_model_part->FullName();
            }

            msg << "\nFollowings are the sensitivity model parts:";
            for (const auto p_model_part : rSensitivityModelParts) {
                msg << "\n\t" << p_model_part->FullName();
            }

            KRATOS_ERROR << msg.str();
        }

        if (EchoLevel > 0) {
            std::stringstream msg;
            msg << "Created sensitivity computation model part based on evaluated model parts [ ";
            for (const auto p_model_part : rEvaluatedModelParts) {
                msg << p_model_part->FullName() << " ";
            }

            msg << "] for sensitivity model parts [ ";
            for (const auto p_model_part : rSensitivityModelParts) {
                msg << p_model_part->FullName() << " ";
            }
            msg << "]";
            msg << (AreNodesConsidered ? " with nodes" : "");
            msg << (AreConditionsConsidered ? " with conditions" : "");
            msg << (AreElementsConsidered ? " with elements" : "");
            msg << ".";
            KRATOS_INFO("ResponseUtils") << msg.str() << std::endl;
        }
        return model_part;
    } else {
        if (EchoLevel > 1) {
            std::stringstream msg;
            msg << "Retrieved sensitivity computation model part based on evaluated model parts [ ";
            for (const auto p_model_part : rEvaluatedModelParts) {
                msg << p_model_part->FullName() << " ";
            }

            msg << "] for sensitivity model parts [ ";
            for (const auto p_model_part : rSensitivityModelParts) {
                msg << p_model_part->FullName() << " ";
            }
            msg << "]";
            msg << (AreNodesConsidered ? " with nodes" : "");
            msg << (AreConditionsConsidered ? " with conditions" : "");
            msg << (AreElementsConsidered ? " with elements" : "");
            msg << ".";
            KRATOS_INFO("ResponseUtils") << msg.str() << std::endl;
        }
        return r_model.GetModelPart(unique_mp_name);
    }

    KRATOS_CATCH("");
}

} // namespace Kratos