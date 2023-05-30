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

void ModelPartOperationUtilities::FillNodesPointerSet(
    std::set<ModelPart::NodeType const*>& rOutput,
    const ModelPart::NodesContainerType& rNodes)
{
    for (auto& r_node : rNodes) {
        rOutput.insert(&r_node);
    }
}

template <class TContainerType>
void ModelPartOperationUtilities::FillNodePointersForEntities(
    std::set<CNodePointersType>& rOutput,
    const TContainerType& rEntityContainer)
{
    auto result = block_for_each<AccumReduction<CNodePointersType, std::set<CNodePointersType>>>(rEntityContainer, [](const auto& rEntity) {
        const auto& r_geometry = rEntity.GetGeometry();
        const IndexType number_of_nodes = r_geometry.size();

        CNodePointersType nodes(number_of_nodes);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            nodes[i] = &r_geometry[i];
        }

        std::sort(nodes.begin(), nodes.end());
        return nodes;
    });

    rOutput.merge(result);
}

void ModelPartOperationUtilities::FillCommonNodes(
    ModelPart::NodesContainerType& rOutput,
    ModelPart::NodesContainerType& rMainNodes,
    const std::set<ModelPart::NodeType const*>& rNodesSet)
{
    const auto& result = block_for_each<AccumReduction<ModelPart::NodeType*>>(rMainNodes, [&rNodesSet](auto& rNode) -> ModelPart::NodeType* {
        auto p_itr = rNodesSet.find(&rNode);

        if (p_itr != rNodesSet.end()) {
            return &rNode;
        } else {
            return nullptr;
        }
    });

    std::for_each(result.begin(), result.end(), [&rOutput](auto& p_node) {
        if (p_node != nullptr) {
            rOutput.push_back(p_node);
        }
    });
}

template<class TEntityType>
void ModelPartOperationUtilities::FillNodesFromEntities(
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
void ModelPartOperationUtilities::FindNeighbourEntities(
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

void ModelPartOperationUtilities::AddNeighbours(
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

void ModelPartOperationUtilities::FillWithMainNodesFromSearchNodes(
    ModelPart::NodesContainerType& rOutput,
    ModelPart::NodesContainerType& rMainNodes,
    ModelPart::NodesContainerType& rSearchedNodes)
{
    std::set<ModelPart::NodeType const*> search_local_node_set;
    FillNodesPointerSet(search_local_node_set, rSearchedNodes);
    FillCommonNodes(rOutput, rMainNodes, search_local_node_set);
}

void ModelPartOperationUtilities::SetCommunicator(
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
    FillWithMainNodesFromSearchNodes(p_output_communicator->LocalMesh().Nodes(), r_reference_communicator.LocalMesh().Nodes(), rOutputModelPart.Nodes());
    FillWithMainNodesFromSearchNodes(p_output_communicator->InterfaceMesh().Nodes(), r_reference_communicator.InterfaceMesh().Nodes(), rOutputModelPart.Nodes());
    FillWithMainNodesFromSearchNodes(p_output_communicator->GhostMesh().Nodes(), r_reference_communicator.GhostMesh().Nodes(), rOutputModelPart.Nodes());

    // now set the communicator
    rOutputModelPart.SetCommunicator(p_output_communicator);
}

ModelPart& ModelPartOperationUtilities::CreateOutputModelPart(
    const std::string& rOutputSubModelPartName,
    ModelPart& rMainModelPart,
    std::vector<ModelPart::NodeType*>& rOutputNodes,
    std::vector<ModelPart::ConditionType*>& rOutputConditions,
    std::vector<ModelPart::ElementType*>& rOutputElements)
{
    KRATOS_ERROR_IF(rMainModelPart.HasSubModelPart(rOutputSubModelPartName))
        << "\"" << rOutputSubModelPartName << "\" already exists in the \""
        << rMainModelPart.FullName() << "\".\n";

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

    // until now everything is sorted based on the memory location ptrs.
    // sorting them based on the ids.
    const auto& sort_mesh = [](ModelPart::MeshType& rMesh) {
        rMesh.Nodes().Sort();
        rMesh.Conditions().Sort();
        rMesh.Elements().Sort();
    };

    sort_mesh(r_output_model_part.GetMesh());
    sort_mesh(r_output_model_part.GetCommunicator().LocalMesh());
    sort_mesh(r_output_model_part.GetCommunicator().GhostMesh());
    sort_mesh(r_output_model_part.GetCommunicator().InterfaceMesh());

    return r_output_model_part;
}

void ModelPartOperationUtilities::CheckNodes(
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
void ModelPartOperationUtilities::CheckEntities(
    std::vector<IndexType>& rEntityIdsWithIssues,
    const TContainerType& rCheckEntities,
    const std::set<CNodePointersType>& rMainEntities)
{
    const auto& result = block_for_each<AccumReduction<IndexType>>(rCheckEntities, CNodePointersType(), [&rMainEntities](const auto& rCheckEntity, CNodePointersType& rTLS) -> IndexType {
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

bool ModelPartOperationUtilities::CheckValidityOfModelPartsForOperations(
    const ModelPart& rMainModelPart,
    const std::vector<ModelPart const*>& rCheckModelParts,
    const bool ThrowError)
{
    std::set<ModelPart::NodeType const*> main_node_sets;
    std::set<CNodePointersType> main_condition_sets;
    std::set<CNodePointersType> main_element_sets;

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

bool ModelPartOperationUtilities::HasIntersection(const std::vector<ModelPart*>& rIntersectionModelParts)
{
    KRATOS_ERROR_IF(rIntersectionModelParts.size() < 2) << "HasIntersection check requires atleast 2 model parts.\n";

    auto& r_main_model_part = rIntersectionModelParts[0]->GetRootModelPart();

    std::vector<ModelPart::NodeType*> output_nodes;
    std::vector<ModelPart::ConditionType*> output_conditions;
    std::vector<ModelPart::ElementType*> output_elements;

    // fill vectors
    std::vector<ModelPart const*> intersection_model_parts(rIntersectionModelParts.size());
    std::transform(rIntersectionModelParts.begin(), rIntersectionModelParts.end(), intersection_model_parts.begin(), [](ModelPart* p_model_part) -> ModelPart const* {return &*(p_model_part); });
    ModelPartOperation<ModelPartIntersectionOperator>(output_nodes, output_conditions, output_elements, r_main_model_part, intersection_model_parts, false);

    bool is_intersected = false;

    is_intersected = is_intersected || r_main_model_part.GetCommunicator().GetDataCommunicator().SumAll(static_cast<unsigned int>(output_nodes.size())) > 0;
    is_intersected = is_intersected || r_main_model_part.GetCommunicator().GetDataCommunicator().SumAll(static_cast<unsigned int>(output_conditions.size())) > 0;
    is_intersected = is_intersected || r_main_model_part.GetCommunicator().GetDataCommunicator().SumAll(static_cast<unsigned int>(output_elements.size())) > 0;

    return is_intersected;
}

// template instantiations
template KRATOS_API(KRATOS_CORE) void ModelPartOperationUtilities::FillNodePointersForEntities<ModelPart::ConditionsContainerType>(std::set<CNodePointersType>&, const ModelPart::ConditionsContainerType&);
template KRATOS_API(KRATOS_CORE) void ModelPartOperationUtilities::FillNodePointersForEntities<ModelPart::ElementsContainerType>(std::set<CNodePointersType>&, const ModelPart::ElementsContainerType&);

template KRATOS_API(KRATOS_CORE) void ModelPartOperationUtilities::FillNodesFromEntities<ModelPart::ConditionType>(std::vector<ModelPart::NodeType*>&, typename std::vector<ModelPart::ConditionType*>::iterator, typename std::vector<ModelPart::ConditionType*>::iterator);
template KRATOS_API(KRATOS_CORE) void ModelPartOperationUtilities::FillNodesFromEntities<ModelPart::ElementType>(std::vector<ModelPart::NodeType*>&, typename std::vector<ModelPart::ElementType*>::iterator, typename std::vector<ModelPart::ElementType*>::iterator);

template KRATOS_API(KRATOS_CORE) void ModelPartOperationUtilities::FindNeighbourEntities<ModelPart::ConditionsContainerType>(std::vector<ModelPart::ConditionType*>&, const Flags&, ModelPart::ConditionsContainerType&);
template KRATOS_API(KRATOS_CORE) void ModelPartOperationUtilities::FindNeighbourEntities<ModelPart::ElementsContainerType>(std::vector<ModelPart::ElementType*>&, const Flags&, ModelPart::ElementsContainerType&);

template KRATOS_API(KRATOS_CORE) void ModelPartOperationUtilities::CheckEntities<ModelPart::ConditionsContainerType>(std::vector<IndexType>&, const ModelPart::ConditionsContainerType&, const std::set<CNodePointersType>&);
template KRATOS_API(KRATOS_CORE) void ModelPartOperationUtilities::CheckEntities<ModelPart::ElementsContainerType>(std::vector<IndexType>&, const ModelPart::ElementsContainerType&, const std::set<CNodePointersType>&);

} // namespace Kratos