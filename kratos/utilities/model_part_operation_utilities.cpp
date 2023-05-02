//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <set>
#include <limits>
#include <sstream>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"

// Include base h
#include "model_part_operation_utilities.h"

namespace Kratos {

namespace ModelPartOperationHelperUtilities {

using IndexType = std::size_t;

void FillNodesPointerSet(
    std::set<ModelPart::NodeType const*>& rOutput,
    const ModelPart::NodesContainerType& rNodes)
{
    for (auto& r_node : rNodes) {
        rOutput.insert(&r_node);
    }
}

template <class TContainerType>
void FillNodePointersForEntities(
    std::set<ModelPartOperationUtilities::CNodePointersType>& rOutput,
    const TContainerType& rEntityContainer)
{
    auto result = block_for_each<AccumReduction<ModelPartOperationUtilities::CNodePointersType, std::set<ModelPartOperationUtilities::CNodePointersType>>>(rEntityContainer, [](const auto& rEntity) {
        const auto& r_geometry = rEntity.GetGeometry();
        const IndexType number_of_nodes = r_geometry.size();

        ModelPartOperationUtilities::CNodePointersType nodes(number_of_nodes);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            nodes[i] = &r_geometry[i];
        }

        std::sort(nodes.begin(), nodes.end());
        return nodes;
    });

    rOutput.merge(result);
}

template<class TOutputContainerType>
IndexType FillCommonNodes(
    TOutputContainerType& rOutput,
    ModelPart::NodesContainerType& rMainNodes,
    const std::set<ModelPart::NodeType const*>& rNodesSet)
{
    const auto& result = block_for_each<AccumReduction<std::pair<IndexType, ModelPart::NodeType*>>>(rMainNodes, [&rNodesSet](auto& rNode) -> std::pair<IndexType, ModelPart::NodeType*> {
        auto p_itr = rNodesSet.find(&rNode);

        if (p_itr != rNodesSet.end()) {
            return std::make_pair(rNode.Id(), &rNode);
        } else {
            return std::make_pair(rNode.Id(), nullptr);
        }
    });

    IndexType found_nullptr = std::numeric_limits<IndexType>::max();
    std::for_each(result.begin(), result.end(), [&found_nullptr, &rOutput](auto& rPair) {
        if (rPair.second != nullptr) {
            rOutput.push_back(rPair.second);
        } else {
            found_nullptr = rPair.first;
        }
    });

    return found_nullptr;
}

template<class TEntityType>
void FillNodesFromEntities(
    std::vector<ModelPart::NodeType*>& rOutput,
    typename std::vector<TEntityType*>::iterator pEntityBegin,
    typename std::vector<TEntityType*>::iterator pEntityEnd)
{
    for (; pEntityBegin < pEntityEnd; ++pEntityBegin) {
        auto& r_geometry = (*pEntityBegin)->GetGeometry();
        for (auto& r_node : r_geometry) {
            rOutput.push_back(&r_node);
        }
    }
}

template<class TContainerType>
void FindNeighbourEntities(
    std::vector<typename TContainerType::value_type*>& rOutputEntities,
    const Flags& rNodeSelectionFlag,
    TContainerType& rMainEntities)
{
    using p_entity_type = typename TContainerType::value_type*;

    const auto& result = block_for_each<AccumReduction<p_entity_type>>(rMainEntities, [&rNodeSelectionFlag](auto& rEntity) -> p_entity_type {
        for (auto& r_node : rEntity.GetGeometry()) {
            if (r_node.Is(rNodeSelectionFlag)) {
                return &rEntity;
            }
        }
        return nullptr;
    });

    std::for_each(result.begin(), result.end(), [&rOutputEntities](auto pEntity) {
        if (pEntity != nullptr) {
            rOutputEntities.push_back(pEntity);
        }
    });
}

void AddNeighbours(
    std::vector<ModelPart::NodeType*>& rOutputNodes,
    std::vector<ModelPart::ConditionType*>& rOutputConditions,
    std::vector<ModelPart::ElementType*>& rOutputElements,
    ModelPart& rMainModelPart)
{
    // set the selection flags
    VariableUtils().SetFlag(VISITED, false, rMainModelPart.Nodes());
    block_for_each(rOutputNodes, [](auto p_node) {
        p_node->Set(VISITED, true);
    });

    // find the neighbour conditions
    FindNeighbourEntities(rOutputConditions, VISITED, rMainModelPart.Conditions());

    // find the neighbour elements
    FindNeighbourEntities(rOutputElements, VISITED, rMainModelPart.Elements());
}

void FillWithMainNodesFromSearchNodes(
    ModelPart::NodesContainerType& rOutput,
    ModelPart::NodesContainerType& rMainNodes,
    ModelPart::NodesContainerType& rSearchedNodes)
{
    std::set<ModelPart::NodeType const*> search_local_node_set;
    FillNodesPointerSet(search_local_node_set, rSearchedNodes);
    FillCommonNodes(rOutput, rMainNodes, search_local_node_set);
}

void SetCommunicator(
    ModelPart& rOutputModelPart,
    ModelPart& rMainModelPart)
{
    // create the proper communicators
    Communicator& r_reference_communicator = rMainModelPart.GetCommunicator();
    Communicator::Pointer p_output_communicator = r_reference_communicator.Create();
    p_output_communicator->SetNumberOfColors(r_reference_communicator.GetNumberOfColors());
    p_output_communicator->NeighbourIndices() = r_reference_communicator.NeighbourIndices();

    // there are no ghost or interface conditions and elements. Hence, we add all conditions
    // and elements to the local mesh
    std::for_each(rOutputModelPart.Conditions().begin(), rOutputModelPart.Conditions().end(), [&p_output_communicator](auto& rCondition) {
        p_output_communicator->LocalMesh().Conditions().push_back(Kratos::intrusive_ptr<ModelPart::ConditionType>(&rCondition));
    });

    std::for_each(rOutputModelPart.Elements().begin(), rOutputModelPart.Elements().end(), [&p_output_communicator](auto& rElement) {
        p_output_communicator->LocalMesh().Elements().push_back(Kratos::intrusive_ptr<ModelPart::ElementType>(&rElement));
    });

    // now set the local, interface and ghost meshes
    ModelPartOperationHelperUtilities::FillWithMainNodesFromSearchNodes(p_output_communicator->LocalMesh().Nodes(), r_reference_communicator.LocalMesh().Nodes(), rOutputModelPart.Nodes());
    ModelPartOperationHelperUtilities::FillWithMainNodesFromSearchNodes(p_output_communicator->InterfaceMesh().Nodes(), r_reference_communicator.InterfaceMesh().Nodes(), rOutputModelPart.Nodes());
    ModelPartOperationHelperUtilities::FillWithMainNodesFromSearchNodes(p_output_communicator->GhostMesh().Nodes(), r_reference_communicator.GhostMesh().Nodes(), rOutputModelPart.Nodes());

    // now set the communicator
    rOutputModelPart.SetCommunicator(p_output_communicator);
}

ModelPart& CreateOutputModelPart(
    const std::string& rOutputSubModelPartName,
    ModelPart& rMainModelPart,
    std::vector<ModelPart::NodeType*>& rOutputNodes,
    std::vector<ModelPart::ConditionType*>& rOutputConditions,
    std::vector<ModelPart::ElementType*>& rOutputElements)
{
    // create the output sub model part
    auto& r_output_model_part = rMainModelPart.CreateSubModelPart(rOutputSubModelPartName);

    // add unique conditions
    std::sort(rOutputConditions.begin(), rOutputConditions.end());
    const auto& condition_last = std::unique(rOutputConditions.begin(), rOutputConditions.end());
    std::for_each(rOutputConditions.begin(), condition_last, [&r_output_model_part](auto p_condition) {
        r_output_model_part.Conditions().push_back(Kratos::intrusive_ptr<ModelPart::ConditionType>(p_condition));
    });
    FillNodesFromEntities<ModelPart::ConditionType>(rOutputNodes, rOutputConditions.begin(), condition_last);

    // add uniqe elements
    std::sort(rOutputElements.begin(), rOutputElements.end());
    const auto& element_last = std::unique(rOutputElements.begin(), rOutputElements.end());
    std::for_each(rOutputElements.begin(), element_last, [&r_output_model_part](auto p_element) {
        r_output_model_part.Elements().push_back(Kratos::intrusive_ptr<ModelPart::ElementType>(p_element));
    });
    FillNodesFromEntities<ModelPart::ElementType>(rOutputNodes, rOutputElements.begin(), element_last);

    // populate the mesh with nodes.
    std::sort(rOutputNodes.begin(), rOutputNodes.end());
    const auto& node_last = std::unique(rOutputNodes.begin(), rOutputNodes.end());
    std::for_each(rOutputNodes.begin(), node_last, [&r_output_model_part](auto p_node) {
        r_output_model_part.Nodes().push_back(Kratos::intrusive_ptr<ModelPart::NodeType>(p_node));
    });

    // sets the communicator info.
    SetCommunicator(r_output_model_part, rMainModelPart);

    return r_output_model_part;
}

void CheckNodes(
    std::vector<IndexType>& rNodeIdsWithIssues,
    const ModelPart::NodesContainerType& rCheckNodes,
    const std::set<ModelPart::NodeType const*>& rMainNodes)
{
    const auto& result = block_for_each<AccumReduction<IndexType>>(rCheckNodes, [&rMainNodes](const auto& rCheckNode) -> IndexType {
        if (rMainNodes.find(&rCheckNode) == rMainNodes.end()) {
            return rCheckNode.Id();
        } else {
            return std::numeric_limits<IndexType>::max();
        }
    });

    std::for_each(result.begin(), result.end(), [&rNodeIdsWithIssues](auto NodeId) {
        if (NodeId != std::numeric_limits<IndexType>::max()) {
            rNodeIdsWithIssues.push_back(NodeId);
        }
    });
}

template<class TContainerType>
void CheckEntities(
    std::vector<IndexType>& rEntityIdsWithIssues,
    const TContainerType& rCheckEntities,
    const std::set<ModelPartOperationUtilities::CNodePointersType>& rMainEntities)
{
    const auto& result = block_for_each<AccumReduction<IndexType>>(rCheckEntities, ModelPartOperationUtilities::CNodePointersType(), [&rMainEntities](const auto& rCheckEntity, ModelPartOperationUtilities::CNodePointersType& rTLS) -> IndexType {
        const auto& r_geometry = rCheckEntity.GetGeometry();
        const IndexType number_of_nodes = r_geometry.size();

        if (rTLS.size() != number_of_nodes) {
            rTLS.resize(number_of_nodes);
        }

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rTLS[i] = &r_geometry[i];
        }

        std::sort(rTLS.begin(), rTLS.end());

        if (rMainEntities.find(rTLS) == rMainEntities.end()) {
            return rCheckEntity.Id();
        } else {
            return std::numeric_limits<IndexType>::max();
        }
    });

    std::for_each(result.begin(), result.end(), [&rEntityIdsWithIssues](auto NodeId) {
        if (NodeId != std::numeric_limits<IndexType>::max()) {
            rEntityIdsWithIssues.push_back(NodeId);
        }
    });
}

template<class TUnaryFunc>
void AddNodes(
    std::vector<ModelPart::NodeType*>& rNodesOutput,
    ModelPart::NodesContainerType& rNodesContainer,
    TUnaryFunc&& rIsValidEntity)
{
    using p_entity_type = ModelPart::NodeType*;

    const auto& result = block_for_each<AccumReduction<p_entity_type>>(rNodesContainer, [&rIsValidEntity](auto& rMainEntity) -> p_entity_type {
        ModelPart::NodeType const* p_entity = &rMainEntity;
        if (rIsValidEntity(p_entity)) {
            return &rMainEntity;
        } else {
            return nullptr;
        }
    });

    std::for_each(result.begin(), result.end(), [&rNodesOutput](auto pEntity) {
        if (pEntity != nullptr) {
            rNodesOutput.push_back(pEntity);
        }
    });
}

template<class TContainerType, class TUnaryFunc>
void AddEntities(
    std::vector<typename TContainerType::value_type*>& rEntitiesOutput,
    TContainerType& rMainEntityContainer,
    TUnaryFunc&& rIsValidEntity)
{
    using p_entity_type = typename TContainerType::value_type*;
    const auto& result = block_for_each<AccumReduction<p_entity_type>>(rMainEntityContainer, ModelPartOperationUtilities::CNodePointersType(), [&rIsValidEntity](auto& rMainEntity, ModelPartOperationUtilities::CNodePointersType& rTLS) -> p_entity_type {
        const auto& r_geometry = rMainEntity.GetGeometry();
        const IndexType number_of_nodes = r_geometry.size();

        if (rTLS.size() != number_of_nodes) {
            rTLS.resize(number_of_nodes);
        }

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rTLS[i] = &r_geometry[i];
        }

        std::sort(rTLS.begin(), rTLS.end());

        if (rIsValidEntity(rTLS)) {
            return &rMainEntity;
        } else {
            return nullptr;
        }
    });

    std::for_each(result.begin(), result.end(), [&rEntitiesOutput](auto pEntity) {
        if (pEntity != nullptr) {
            rEntitiesOutput.push_back(pEntity);
        }
    });
}

template<class TModelPartOperation>
ModelPart& ModelPartOperation(
    const std::string& rOutputSubModelPartName,
    ModelPart& rMainModelPart,
    const std::vector<ModelPart const*>& rModelPartOperationModelParts,
    const bool AddNeighbourEntities)
{
    KRATOS_ERROR_IF(rMainModelPart.HasSubModelPart(rOutputSubModelPartName))
        << "\"" << rOutputSubModelPartName << "\" already exists in the \""
        << rMainModelPart.FullName() << "\".\n";

    const IndexType number_of_operation_model_parts = rModelPartOperationModelParts.size();

    std::vector<std::set<ModelPart::NodeType const*>> set_operation_node_sets(number_of_operation_model_parts);
    std::vector<std::set<ModelPartOperationUtilities::CNodePointersType>> set_operation_condition_sets(number_of_operation_model_parts);
    std::vector<std::set<ModelPartOperationUtilities::CNodePointersType>> set_operation_element_sets(number_of_operation_model_parts);

    // now iterate through model parts' conditions/elements containers and get available main model part entities.
    for (IndexType i = 0; i < number_of_operation_model_parts; ++i) {
        auto p_model_part = rModelPartOperationModelParts[i];

        // fill the set_operation nodes
        FillNodesPointerSet(set_operation_node_sets[i], p_model_part->Nodes());

        // fill the set_operation conditions
        FillNodePointersForEntities(set_operation_condition_sets[i], p_model_part->Conditions());

        // fill the set_operation elements
        FillNodePointersForEntities(set_operation_element_sets[i], p_model_part->Elements());
    }

    std::vector<ModelPart::NodeType*> output_nodes;
    std::vector<ModelPart::ConditionType*> output_conditions;
    std::vector<ModelPart::ElementType*> output_elements;

    AddNodes(output_nodes, rMainModelPart.Nodes(), [&set_operation_node_sets](auto& rEntity) {
        return TModelPartOperation::IsValid(rEntity, set_operation_node_sets);
    });

    AddEntities(output_conditions, rMainModelPart.Conditions(), [&set_operation_condition_sets](auto& rEntity) {
        return TModelPartOperation::IsValid(rEntity, set_operation_condition_sets);
    });

    AddEntities(output_elements, rMainModelPart.Elements(), [&set_operation_element_sets](auto& rEntity) {
        return TModelPartOperation::IsValid(rEntity, set_operation_element_sets);
    });

    // now we have all the nodes to find and add neighbour entities.
    if (AddNeighbourEntities) {
        // we need to fill the boundary nodes for elements and conditions
        FillNodesFromEntities<ModelPart::ConditionType>(output_nodes, output_conditions.begin(), output_conditions.end());
        FillNodesFromEntities<ModelPart::ElementType>(output_nodes, output_elements.begin(), output_elements.end());

        // now add all the neighbours
        AddNeighbours(output_nodes, output_conditions, output_elements, rMainModelPart);
    }

    // now create the sub model part
    return CreateOutputModelPart(rOutputSubModelPartName, rMainModelPart, output_nodes, output_conditions, output_elements);
}

struct Merge
{
    template<class TCheckType>
    static bool IsValid(
        const TCheckType& rEntity,
        const std::vector<std::set<TCheckType>>& mUnionSets)
    {
        bool is_valid = false;

        for (const auto& r_set : mUnionSets) {
            if (r_set.find(rEntity) != r_set.end()) {
                is_valid = true;
                break;
            }
        }

        return is_valid;
    }
};

struct Substraction
{
    template<class TCheckType>
    static bool IsValid(
        const TCheckType& rEntity,
        const std::vector<std::set<TCheckType>>& mSubstractionSets)
    {
        bool is_valid = true;

        for (const auto& r_set : mSubstractionSets) {
            if (r_set.find(rEntity) != r_set.end()) {
                is_valid = false;
                break;
            }
        }

        return is_valid;
    }
};

struct Intersection
{
    template<class TCheckType>
    static bool IsValid(
        const TCheckType& rEntity,
        const std::vector<std::set<TCheckType>>& mIntersectionSets)
    {
        bool is_valid = true;

        for (const auto& r_set : mIntersectionSets) {
            if (r_set.find(rEntity) == r_set.end()) {
                is_valid = false;
                break;
            }
        }

        return is_valid;
    }
};

} // namespace ModelPartOperationHelperUtilities

bool ModelPartOperationUtilities::CheckValidityOfModelPartsForOperations(
    const ModelPart& rMainModelPart,
    const std::vector<ModelPart const*>& rCheckModelParts,
    const bool ThrowError)
{
    using namespace ModelPartOperationHelperUtilities;

    std::set<ModelPart::NodeType const*> main_node_sets;
    std::set<ModelPartOperationUtilities::CNodePointersType> main_condition_sets;
    std::set<ModelPartOperationUtilities::CNodePointersType> main_element_sets;

    // fill the check nodes
    FillNodesPointerSet(main_node_sets, rMainModelPart.Nodes());

    // fill the check conditions
    FillNodePointersForEntities(main_condition_sets, rMainModelPart.Conditions());

    // fill the check elements
    FillNodePointersForEntities(main_element_sets, rMainModelPart.Elements());

    std::vector<IndexType> issue_ids;

    for (auto p_model_part : rCheckModelParts) {
        CheckNodes(issue_ids, p_model_part->Nodes(), main_node_sets);

        if (issue_ids.size() > 0) {
            if (ThrowError) {
                std::stringstream msg;
                msg << "Following nodes with ids in \"" << p_model_part->FullName() << "\" not found in the main model part \"" << rMainModelPart.FullName() << "\":";
                for (const auto node_id : issue_ids) {
                    msg << "\n\t" << node_id;
                }

                KRATOS_ERROR << msg.str();
            }

            return false;
        }

        CheckEntities(issue_ids, p_model_part->Conditions(), main_condition_sets);

        if (issue_ids.size() > 0) {
            if (ThrowError) {
                std::stringstream msg;
                msg << "Following conditions with ids in \"" << p_model_part->FullName() << "\" not found in the main model part \"" << rMainModelPart.FullName() << "\":";
                for (const auto condition_id : issue_ids) {
                    msg << "\n\t" << condition_id;
                }

                KRATOS_ERROR << msg.str();
            }

            return false;
        }

        CheckEntities(issue_ids, p_model_part->Elements(), main_element_sets);

        if (issue_ids.size() > 0) {
            if (ThrowError) {
                std::stringstream msg;
                msg << "Following elements with ids in \"" << p_model_part->FullName() << "\" not found in the main model part \"" << rMainModelPart.FullName() << "\":";
                for (const auto element_id : issue_ids) {
                    msg << "\n\t" << element_id;
                }

                KRATOS_ERROR << msg.str();
            }

            return false;
        }
    }

    return true;
}

ModelPart& ModelPartOperationUtilities::Merge(
    const std::string& rOutputSubModelPartName,
    ModelPart& rMainModelPart,
    const std::vector<ModelPart const*>& rUnionModelParts,
    const bool AddNeighbourEntities)
{
    return ModelPartOperationHelperUtilities::ModelPartOperation<ModelPartOperationHelperUtilities::Merge>(
        rOutputSubModelPartName, rMainModelPart, rUnionModelParts,
        AddNeighbourEntities);
}

ModelPart& ModelPartOperationUtilities::Substract(
    const std::string& rOutputSubModelPartName,
    ModelPart& rMainModelPart,
    const std::vector<ModelPart const*>& rSubstractionModelParts,
    const bool AddNeighbourEntities)
{
    return ModelPartOperationHelperUtilities::ModelPartOperation<ModelPartOperationHelperUtilities::Substraction>(
        rOutputSubModelPartName, rMainModelPart, rSubstractionModelParts,
        AddNeighbourEntities);
}

ModelPart& ModelPartOperationUtilities::Intersect(
    const std::string& rOutputSubModelPartName,
    ModelPart& rMainModelPart,
    const std::vector<ModelPart const*>& rIntersectionModelParts,
    const bool AddNeighbourEntities)
{
    return ModelPartOperationHelperUtilities::ModelPartOperation<ModelPartOperationHelperUtilities::Intersection>(
        rOutputSubModelPartName, rMainModelPart, rIntersectionModelParts,
        AddNeighbourEntities);
}

} // namespace Kratos