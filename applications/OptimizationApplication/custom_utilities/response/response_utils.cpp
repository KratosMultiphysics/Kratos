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
#include <sstream>
#include <map>

// Project includes
#include "includes/condition.h"
#include "includes/element.h"
#include "containers/model.h"
#include "containers/global_pointers_vector.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
#include "utilities/reduction_utilities.h"

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
            return std::make_pair(0UL, &rEntity);
        }
    });

    // remove the nullptr if present
    if (entity_id_ptr_map.find(0) != entity_id_ptr_map.end()) {
        entity_id_ptr_map.erase(0);
    }

    rOutput.merge(entity_id_ptr_map);
}

ModelPart& ResponseUtils::GetSensitivityModelPartForAdjointSensitivities(
    const std::vector<ModelPart*>& rSensitivityModelParts,
    ModelPart& rAnalysisModelPart,
    const bool AreSensitivityEntityParentsConsidered,
    const bool AreSensitivityEntitesConsidered,
    const bool ForceFindSensitivityEntitiesInAnalysisModelPart)
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

        // finally we add all the conditions and elemenets from the maps, we don't need to call Unique in
        // here because, we are using a std::map which is sorted with entity ids.
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

        return model_part;
    } else {
        return r_model.GetModelPart(unique_mp_name);
    }

    KRATOS_CATCH("");
}

ModelPart& ResponseUtils::GetSensitivityModelPartForDirectSensitivities(
    const std::vector<ModelPart*>& rSensitivityModelParts,
    const std::vector<ModelPart*>& rEvaluatedModelParts,
    const bool AreNodesConsidered,
    const bool AreConditionsConsidered,
    const bool AreElementsConsidered)
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

        // finally we add all the conditions and elemenets from the maps, we don't need to call Unique in
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

        return model_part;
    } else {
        return r_model.GetModelPart(unique_mp_name);
    }

    KRATOS_CATCH("");
}

void ResponseUtils::CheckAndPrepareModelPartsForSensitivityComputation(
    const std::vector<ModelPart*>& rEvaluatedModelParts,
    const SensitivityModelPartVariablesListMap& rSensitivityModelPartVariableInfo,
    const Flags& rFlag,
    const std::vector<SensitivityFieldVariableTypes>& rUsedNodalSensitivityVariables)
{
    KRATOS_TRY

    // reset entity flags for sensitivity model parts
    for (auto& it : rSensitivityModelPartVariableInfo) {
        VariableUtils().SetFlag(rFlag, false, it.first->Elements());
    }

    // set entity flags for evaluated model parts
    for (auto& p_model_part : rEvaluatedModelParts) {
        VariableUtils().SetFlag(rFlag, true, p_model_part->Elements());
    }

    // get number of overlapping entities
    OptimizationUtils::SensitivityVariableModelPartsListMap reversed_map;
    OptimizationUtils::ReverseSensitivityModelPartVariablesListMap(reversed_map, rSensitivityModelPartVariableInfo);
    for (const auto& it : reversed_map) {

        IndexType number_of_common_entities = 0;
        for (const auto& p_model_part : it.second) {
            number_of_common_entities += OptimizationUtils::GetNumberOfContainerItemsWithFlag(p_model_part->Elements(), p_model_part->GetCommunicator().GetDataCommunicator(), rFlag);
        }

        std::visit([&](auto&& r_variable) {
            if (number_of_common_entities == 0) {
                std::stringstream msg;
                msg << "No common entities between evaluated and sensitivity "
                       "model parts found for sensitivity variable "
                    << r_variable->Name() << ".";

                msg << "\nFollowings are the evaluated model parts:";
                for (const auto& p_model_part : rEvaluatedModelParts) {
                    msg << "\n\t" << p_model_part->FullName();
                }

                msg << "\nFollowings are the sensitivity model parts:";
                for (const auto& p_model_part : it.second) {
                    msg << "\n\t" << p_model_part->FullName();
                }

                KRATOS_ERROR << msg.str();
            }
        }, it.first);
    }

    // clear all the sensitivity variables for nodes. Here we assume there are
    // no overlapping regions in Elements and/or Conditions between provided rSensitivityModelParts hence, SetValue is
    // used in Elements and/or Condtions. Nodal sensitivities are added so that common nodes between two model parts
    // will have correct sensitivities.
    for (const auto& it_var_1 : rUsedNodalSensitivityVariables) {
        std::visit([&](auto&& p_var_1) {
            for (const auto& it : reversed_map) {
                std::visit([&](auto&& p_var_2) {
                    if (*p_var_1 == *p_var_2) {
                        for (const auto p_model_part : it.second) {
                            VariableUtils().SetNonHistoricalVariableToZero(*p_var_1, p_model_part->Nodes());
                        }
                    }
                }, it.first);
            }
        }, it_var_1);
    }

    KRATOS_CATCH("");
}

} // namespace Kratos