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
#include <algorithm>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"

// Include base h
#include "model_part_operation_utilities.h"

namespace Kratos {

namespace ModelPartOperationHelperUtilities
{
template <class TContainerType>
void FindNeighbourEntities(
    TContainerType& rOutputEntities,
    const Flags& rNodeSelectionFlag,
    const TContainerType& rMainEntities)
{
    TContainerType temp;
    temp.reserve(rMainEntities.size());

    for (auto itr_entity = rMainEntities.begin(); itr_entity != rMainEntities.end(); ++itr_entity) {
        for (const auto& r_node : itr_entity->GetGeometry()) {
            if (r_node.Is(rNodeSelectionFlag)) {
                // since the rMainEntities is a PVS, here we can use insert with the correct
                // hint without a problem which will use a push_back in the back end.
                temp.insert(temp.end(), *(itr_entity.base()));
            }
        }
    }

    rOutputEntities.insert(temp.begin(), temp.end());
}

void FillCommonNodes(
    ModelPart::NodesContainerType& rOutput,
    const ModelPart::NodesContainerType& rMainNodes,
    const ModelPart::NodesContainerType& rSearchNodes)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rOutput.empty()) << "Output container is not empty.";

    // checks whether rSearchNodes are within the rMainNodes, if so adds them to rOutput
    for (auto itr_search = rSearchNodes.begin(); itr_search != rSearchNodes.end(); ++itr_search) {
        if (rMainNodes.find(itr_search->Id()) != rMainNodes.end()) {
            rOutput.insert(rOutput.end(), *(itr_search.base()));
        }
    }

    KRATOS_CATCH("");
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
    p_output_communicator->LocalMesh().SetConditions(rOutputModelPart.pConditions());
    p_output_communicator->LocalMesh().SetElements(rOutputModelPart.pElements());

    // now set the local, interface and ghost meshes
    FillCommonNodes(p_output_communicator->LocalMesh().Nodes(), r_reference_communicator.LocalMesh().Nodes(), rOutputModelPart.Nodes());
    FillCommonNodes(p_output_communicator->InterfaceMesh().Nodes(), r_reference_communicator.InterfaceMesh().Nodes(), rOutputModelPart.Nodes());
    FillCommonNodes(p_output_communicator->GhostMesh().Nodes(), r_reference_communicator.GhostMesh().Nodes(), rOutputModelPart.Nodes());

    // now set the communicator
    rOutputModelPart.SetCommunicator(p_output_communicator);
}

template<class TContainerType>
void CheckEntities(
    std::vector<IndexType>& rEntityIdsWithIssues,
    const TContainerType& rCheckNodes,
    const TContainerType& rMainNodes)
{
    rEntityIdsWithIssues.clear();

    const auto& result = block_for_each<AccumReduction<IndexType>>(rCheckNodes, [&rMainNodes](const auto& rCheckNode) -> IndexType {
        if (rMainNodes.find(rCheckNode.Id()) == rMainNodes.end()) {
            return rCheckNode.Id();
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

} // namespace ModelPartOperationHelperUtilities

template<class TContainer>
void ModelPartOperationUtilities::AddNodesFromEntities(
    ModelPart::NodesContainerType& rOutputNodes,
    const TContainer& rInputEntities)
{
    // do not do anything if empty.
    if (rInputEntities.empty()) {
        return;
    }

    std::vector<ModelPart::NodeType::Pointer> temp;
    temp.reserve(rInputEntities.size() * rInputEntities.front().GetGeometry().size());

    for (const auto& r_entity : rInputEntities) {
        const auto& r_geometry = r_entity.GetGeometry();
        for (auto itr = r_geometry.begin(); itr != r_geometry.end(); ++itr) {
            temp.push_back(*(itr.base()));
        }
    }

    rOutputNodes.insert(std::move(temp));
}

void ModelPartOperationUtilities::AddNeighbours(
    ModelPart::NodesContainerType& rOutputNodes,
    ModelPart::ConditionsContainerType& rOutputConditions,
    ModelPart::ElementsContainerType& rOutputElements,
    ModelPart& rMainModelPart)
{
    // set the selection flags
    VariableUtils().SetFlag(VISITED, false, rMainModelPart.Nodes());
    block_for_each(rOutputNodes, [](auto& rNode) {
        rNode.Set(VISITED, true);
    });

    // find the neighbour conditions
    ModelPartOperationHelperUtilities::FindNeighbourEntities(rOutputConditions, VISITED, rMainModelPart.Conditions());

    // find the neighbour elements
    ModelPartOperationHelperUtilities::FindNeighbourEntities(rOutputElements, VISITED, rMainModelPart.Elements());
}

void ModelPartOperationUtilities::FillOutputSubModelPart(
    ModelPart& rOutputSubModelPart,
    ModelPart& rMainModelPart,
    ModelPart::NodesContainerType& rOutputNodes,
    ModelPart::ConditionsContainerType& rOutputConditions,
    ModelPart::ElementsContainerType& rOutputElements)
{
    // check if the sub model part is empty. Here we use the data communicator
    // of rMainModelPart assuming the communicators are not set yet in rOutputSubModelPart
    // since it is assumed to be empty in all contaienrs.
    const auto& r_data_communicator = rMainModelPart.GetCommunicator().GetDataCommunicator();
    const auto& entity_info = r_data_communicator.SumAll(std::vector<unsigned int>{
            static_cast<unsigned int>(rOutputSubModelPart.NumberOfNodes()),
            static_cast<unsigned int>(rOutputSubModelPart.NumberOfConditions()),
            static_cast<unsigned int>(rOutputSubModelPart.NumberOfElements())});

    KRATOS_ERROR_IF(entity_info[0] > 0 || entity_info[1] > 0 || entity_info[2] > 0)
        << rOutputSubModelPart.FullName() << " is not empty.";

    // add conditions
    rOutputSubModelPart.Conditions().insert(rOutputConditions.begin(), rOutputConditions.end());
    AddNodesFromEntities(rOutputNodes, rOutputConditions);

    // // add uniqe elements
    rOutputSubModelPart.Elements().insert(rOutputElements.begin(), rOutputElements.end());
    AddNodesFromEntities(rOutputNodes, rOutputElements);

    // populate the mesh with nodes.
    rOutputSubModelPart.Nodes().insert(rOutputNodes.begin(), rOutputNodes.end());

    // sets the communicator info.
    ModelPartOperationHelperUtilities::SetCommunicator(rOutputSubModelPart, rMainModelPart);
}

bool ModelPartOperationUtilities::CheckValidityOfModelPartsForOperations(
    const ModelPart& rMainModelPart,
    const std::vector<ModelPart const*>& rCheckModelParts,
    const bool ThrowError)
{
    std::vector<IndexType> issue_ids;

    for (auto p_model_part : rCheckModelParts) {
        ModelPartOperationHelperUtilities::CheckEntities(issue_ids, p_model_part->Nodes(), rMainModelPart.Nodes());

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

        ModelPartOperationHelperUtilities::CheckEntities(issue_ids, p_model_part->Conditions(), rMainModelPart.Conditions());

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

        ModelPartOperationHelperUtilities::CheckEntities(issue_ids, p_model_part->Elements(), rMainModelPart.Elements());

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

    ModelPart::NodesContainerType output_nodes;
    ModelPart::ConditionsContainerType output_conditions;
    ModelPart::ElementsContainerType output_elements;

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
template KRATOS_API(KRATOS_CORE) void ModelPartOperationUtilities::AddNodesFromEntities<ModelPart::ConditionsContainerType>(ModelPart::NodesContainerType&,const ModelPart::ConditionsContainerType&);
template KRATOS_API(KRATOS_CORE) void ModelPartOperationUtilities::AddNodesFromEntities<ModelPart::ElementsContainerType>(ModelPart::NodesContainerType&,const ModelPart::ElementsContainerType&);

} // namespace Kratos