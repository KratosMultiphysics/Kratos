//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Michael Andre, https://github.com/msandre
//                   Jordi Cotela Dalmau
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "mpi/utilities/gather_modelpart_utility.h"
#include "mpi/utilities/parallel_fill_communicator.h"
#include "includes/data_communicator.h"
#include "mpi/includes/mpi_communicator.h"

namespace Kratos
{

GatherModelPartUtility::GatherModelPartUtility(
    const int GatherRank,
    ModelPart& rOriginModelPart,
    const int MeshId,
    ModelPart& rDestinationModelPart
    ) : mrModelPart(rDestinationModelPart), mGatherRank(GatherRank)
{
    KRATOS_TRY;

    const DataCommunicator& r_comm = rOriginModelPart.GetCommunicator().GetDataCommunicator();
    const int mpi_rank = r_comm.Rank();
    const int mpi_size = r_comm.Size();

    rDestinationModelPart.GetNodalSolutionStepVariablesList() = rOriginModelPart.GetNodalSolutionStepVariablesList();
    if (r_comm.IsDistributed()) {
      VariablesList* pVariablesList = &rDestinationModelPart.GetNodalSolutionStepVariablesList();
      rDestinationModelPart.SetCommunicator(Communicator::Pointer(new MPICommunicator(pVariablesList, r_comm)));
    }
    rDestinationModelPart.SetBufferSize(rOriginModelPart.GetBufferSize());

    // Copy the mesh of interest to rDestinationModelPart
    // be careful to push back the pointer and not copy
    // construct the object
    auto& r_mesh = rOriginModelPart.GetMesh(MeshId);
    for (auto it = r_mesh.NodesBegin(); it != r_mesh.NodesEnd(); ++it) {
        rDestinationModelPart.Nodes().push_back(*it.base());
    }

    for (auto it = r_mesh.ElementsBegin();  it != r_mesh.ElementsEnd(); ++it) {
        rDestinationModelPart.Elements().push_back(*it.base());
    }

    for (auto it = r_mesh.ConditionsBegin(); it != r_mesh.ConditionsEnd(); ++it) {
        rDestinationModelPart.Conditions().push_back(*it.base());
    }

    // send everything to node with id "GatherRank"
    // transfer nodes
    std::vector<NodesContainerType> SendNodes(mpi_size);
    std::vector<NodesContainerType> RecvNodes(mpi_size);
    SendNodes[GatherRank].reserve(rDestinationModelPart.Nodes().size());
    if (r_comm.IsDistributed()) {
      for (auto it = rDestinationModelPart.NodesBegin(); it != rDestinationModelPart.NodesEnd(); ++it) {
          // only send the nodes owned by this partition
          if (it->FastGetSolutionStepValue(PARTITION_INDEX) == mpi_rank)
              SendNodes[GatherRank].push_back(*it.base());
      }
    } else {
      for (auto it = rDestinationModelPart.NodesBegin(); it != rDestinationModelPart.NodesEnd(); ++it) {
          SendNodes[GatherRank].push_back(*it.base());
      }
    }

    rDestinationModelPart.GetCommunicator().TransferObjects(SendNodes, RecvNodes);
    for (unsigned int i = 0; i < RecvNodes.size(); i++) {
        for (auto it = RecvNodes[i].begin(); it != RecvNodes[i].end(); ++it) {
            if (rDestinationModelPart.Nodes().find(it->Id()) == rDestinationModelPart.Nodes().end())
                rDestinationModelPart.Nodes().push_back(*it.base());
        }
    }
    int temp = rDestinationModelPart.Nodes().size();
    rDestinationModelPart.Nodes().Unique();
    KRATOS_ERROR_IF(temp != int(rDestinationModelPart.Nodes().size())) << "The rDestinationModelPart has repeated nodes" << std::endl;
    SendNodes.clear();
    RecvNodes.clear();

    // Transfer elements
    std::vector<ElementsContainerType> SendElements(mpi_size);
    std::vector<ElementsContainerType> RecvElements(mpi_size);
    SendElements[GatherRank].reserve(rDestinationModelPart.Elements().size());
    for (auto it = rDestinationModelPart.ElementsBegin(); it != rDestinationModelPart.ElementsEnd(); ++it) {
        SendElements[GatherRank].push_back(*it.base());
    }
    rDestinationModelPart.GetCommunicator().TransferObjects(SendElements, RecvElements);
    for (unsigned int i = 0; i < RecvElements.size(); i++) {
        for (auto it = RecvElements[i].begin(); it != RecvElements[i].end(); ++it) {
            // Replace the nodes copied with the element by nodes in the model part
            Element::GeometryType& rGeom = it->GetGeometry();
            unsigned int NumNodes = rGeom.PointsNumber();
            for (unsigned int iNode = 0; iNode < NumNodes; iNode++) {
                auto itNode = rDestinationModelPart.Nodes().find(rGeom(iNode)->Id());
                if (itNode != rDestinationModelPart.Nodes().end())
                    rGeom(iNode) = *itNode.base();
            }
            rDestinationModelPart.Elements().push_back(*it.base());
        }
    }
    SendElements.clear();
    RecvElements.clear();

    // Transfer conditions
    std::vector<ConditionsContainerType> SendConditions(mpi_size);
    std::vector<ConditionsContainerType> RecvConditions(mpi_size);
    SendConditions[GatherRank].reserve(rDestinationModelPart.Conditions().size());
    for (auto it = rDestinationModelPart.ConditionsBegin(); it != rDestinationModelPart.ConditionsEnd(); ++it) {
        SendConditions[GatherRank].push_back(*it.base());
    }
    rDestinationModelPart.GetCommunicator().TransferObjects(SendConditions, RecvConditions);
    for (unsigned int i = 0; i < RecvConditions.size(); i++) {
        for (auto it = RecvConditions[i].begin(); it != RecvConditions[i].end(); ++it) {
            // Replace the nodes copied with the condition by nodes in the model part
            Condition::GeometryType& rGeom = it->GetGeometry();
            unsigned int NumNodes = rGeom.PointsNumber();
            for (unsigned int iNode = 0; iNode < NumNodes; iNode++) {
                auto itNode = rDestinationModelPart.Nodes().find(rGeom(iNode)->Id());
                if (itNode != rDestinationModelPart.Nodes().end())
                    rGeom(iNode) = *itNode.base();
            }
            rDestinationModelPart.Conditions().push_back(*it.base());
        }
    }
    SendConditions.clear();
    RecvConditions.clear();

    if (r_comm.IsDistributed()) {
        ParallelFillCommunicator(rDestinationModelPart, r_comm).Execute();
    }

    KRATOS_CATCH("");
}

void GatherModelPartUtility::GatherOnMaster()
{
    KRATOS_TRY;
    mrModelPart.GetCommunicator().SynchronizeNodalSolutionStepsData();
    KRATOS_CATCH("");
}

template <class TDataType>
void GatherModelPartUtility::GatherOnMaster(const Variable<TDataType>& ThisVariable)
{
    KRATOS_TRY;
    mrModelPart.GetCommunicator().SynchronizeVariable(ThisVariable);
    KRATOS_CATCH("");
}

template <class TDataType>
void GatherModelPartUtility::ScatterFromMaster(const Variable<TDataType>& ThisVariable)
{
    KRATOS_TRY;

    Communicator& r_comm = mrModelPart.GetCommunicator();

    if (mGatherRank != r_comm.GetDataCommunicator().Rank()) {
        for (auto it = mrModelPart.NodesBegin(); it != mrModelPart.NodesEnd(); ++it)
            it->FastGetSolutionStepValue(ThisVariable) = ThisVariable.Zero();
    }
    r_comm.AssembleCurrentData(ThisVariable);

    KRATOS_CATCH("");
}

template void GatherModelPartUtility::GatherOnMaster(const Variable<double>&);
template void GatherModelPartUtility::GatherOnMaster(const Variable<array_1d<double, 3>>&);
template void GatherModelPartUtility::ScatterFromMaster(const Variable<double>&);
template void GatherModelPartUtility::ScatterFromMaster(const Variable<array_1d<double, 3>>&);

void GatherModelPartUtility::GatherEntitiesFromOtherPartitions(
    ModelPart& rModelPart,
    const std::map<int, std::vector<std::size_t>>& rNodesToBring,
    const std::map<int, std::vector<std::size_t>>& rElementsToBring,
    const std::map<int, std::vector<std::size_t>>& rConditionsToBring,
    const bool CallExecuteAfterBringingEntities,
    const int EchoLevel
    )
{
    KRATOS_TRY

    // Retrieving the model part and the communicator
    const auto& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();

    // Call auxiliary methods
    const std::size_t nodes_to_bring = r_data_communicator.SumAll(rNodesToBring.size());
    if (nodes_to_bring > 0) {
        GatherEntityFromOtherPartitions<Node>(rModelPart, rNodesToBring, EchoLevel);
    }
    const std::size_t elements_to_bring = r_data_communicator.SumAll(rElementsToBring.size());
    if (elements_to_bring > 0) {
        GatherEntityFromOtherPartitions<Element>(rModelPart, rElementsToBring, EchoLevel);
    }
    const std::size_t conditions_to_bring = r_data_communicator.SumAll(rConditionsToBring.size());
    if (conditions_to_bring > 0) {
        GatherEntityFromOtherPartitions<Condition>(rModelPart, rConditionsToBring, EchoLevel);
    }

    // Execute after bringing entities
    if (CallExecuteAfterBringingEntities) {
        ParallelFillCommunicator(rModelPart, r_data_communicator).Execute();
    }

    KRATOS_CATCH("");
}

void GatherModelPartUtility::GatherNodesFromOtherPartitions(
    ModelPart& rModelPart,
    const std::map<int, std::vector<std::size_t>>& rNodesToBring,
    const bool CallExecuteAfterBringingEntities,
    const int EchoLevel
    )
{
    KRATOS_TRY

    // Retrieving the model part and the communicator
    const auto& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();

    // Call auxiliary methods
    const std::size_t nodes_to_bring = r_data_communicator.SumAll(rNodesToBring.size());
    if (nodes_to_bring > 0) {
        GatherEntityFromOtherPartitions<Node>(rModelPart, rNodesToBring, EchoLevel);
    }

    // Execute after bringing entities
    if (CallExecuteAfterBringingEntities) {
        ParallelFillCommunicator(rModelPart, r_data_communicator).Execute();
    }

    KRATOS_CATCH("");
}

void GatherModelPartUtility::GatherElementsFromOtherPartitions(
    ModelPart& rModelPart,
    const std::map<int, std::vector<std::size_t>>& rElementsToBring,
    const bool CallExecuteAfterBringingEntities,
    const int EchoLevel
    )
{
    KRATOS_TRY

    // Retrieving the model part and the communicator
    const auto& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();

    // Call auxiliary methods
    const std::size_t elements_to_bring = r_data_communicator.SumAll(rElementsToBring.size());
    if (elements_to_bring > 0) {
        GatherEntityFromOtherPartitions<Element>(rModelPart, rElementsToBring, EchoLevel);
    }

    // Execute after bringing entities
    if (CallExecuteAfterBringingEntities) {
        ParallelFillCommunicator(rModelPart, r_data_communicator).Execute();
    }

    KRATOS_CATCH("");
}

void GatherModelPartUtility::GatherConditionsFromOtherPartitions(
    ModelPart& rModelPart,
    const std::map<int, std::vector<std::size_t>>& rConditionsToBring,
    const bool CallExecuteAfterBringingEntities,
    const int EchoLevel
    )
{
    KRATOS_TRY

    // Retrieving the model part and the communicator
    const auto& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();

    // Call auxiliary methods
    const std::size_t conditions_to_bring = r_data_communicator.SumAll(rConditionsToBring.size());
    if (conditions_to_bring > 0) {
        GatherEntityFromOtherPartitions<Condition>(rModelPart, rConditionsToBring, EchoLevel);
    }

    // Execute after bringing entities
    if (CallExecuteAfterBringingEntities) {
        ParallelFillCommunicator(rModelPart, r_data_communicator).Execute();
    }

    KRATOS_CATCH("");
}

template <class TObjectType>
void GatherModelPartUtility::GatherEntityFromOtherPartitions(
    ModelPart& rModelPart,
    const std::map<int, std::vector<std::size_t>>& rEntitiesToBring,
    const int EchoLevel
    )
{
    /* First make all partitions aware of which communications are needed (with the current information we only know the entities of we want to bring in current partition) */

    // Entity name
    std::string entity_name;
    if constexpr (std::is_same<TObjectType, Node>::value) {
        entity_name = "Node";
    } else if constexpr (std::is_same<TObjectType, Element>::value) {
        entity_name = "Element";
    } else if constexpr (std::is_same<TObjectType, Condition>::value) {
        entity_name = "Condition";
    } else {
        KRATOS_ERROR << "Entity type not supported" << std::endl;
    }

    /// Type alias for the container of entities.
    using ContainerType = PointerVectorSet<TObjectType,
        IndexedObject,
        std::less<typename IndexedObject::result_type>,
        std::equal_to<typename IndexedObject::result_type>,
        typename TObjectType::Pointer,
        std::vector<typename TObjectType::Pointer>
    >;
    
    // Allocating the temporary entities to bring container
    ContainerType entities_to_bring;
    std::size_t counter = 0;
    for (const auto& r_bring : rEntitiesToBring) {
        counter += r_bring.second.size();
    }
    entities_to_bring.reserve(counter);

    // Retrieve MPI data
    const auto& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();
    const int rank = r_data_communicator.Rank();
    const int world_size = r_data_communicator.Size();

    // First counting how many entities transfer for partition
    int tag_send = 1;
    std::vector<int> other_partition_indices;
    other_partition_indices.reserve(world_size - 1);
    for (int i_rank = 0; i_rank < world_size; ++i_rank) {
        if (i_rank != rank) other_partition_indices.push_back(i_rank);
    }
    std::map<int, std::vector<std::size_t>> send_entities;
    for (int i_rank = 0; i_rank < world_size; ++i_rank) {
        if (i_rank == rank) {
            for (auto index : other_partition_indices) {
                std::vector<std::size_t> send_vector;
                r_data_communicator.Recv(send_vector, index, tag_send);
                // Just adding in case not empty
                if (send_vector.size() > 0) {
                    send_entities[index] = send_vector;
                }
            }
        } else {
            auto it_find = rEntitiesToBring.find(i_rank);
            // Sending in case defined
            if (it_find != rEntitiesToBring.end()) {
                r_data_communicator.Send(it_find->second, i_rank, tag_send);
            } else { // Sending empty vector in case not defined
                std::vector<std::size_t> empty_vector;
                r_data_communicator.Send(empty_vector, i_rank, tag_send);
            }
        }
    }

    // Depending of the echo level, print info
    if (EchoLevel > 0) {
        std::stringstream buffer;
        buffer << "\nRank " << rank << " has to bring " << entity_name << "s from other partitions:" << std::endl;
        for (auto bring : rEntitiesToBring) {
            buffer << "\tFrom rank " << bring.first << ": " << bring.second.size() << " " << entity_name << "s\n\t" << bring.second << std::endl;
        }
        buffer << "\nRank " << rank << " has to send " << entity_name << "s from other partitions:" << std::endl;
        for (auto send : send_entities) {
            buffer << "\tTo rank " << send.first << ": " << send.second.size() << " " << entity_name << "s\n\t" << send.second << std::endl;
        }
        buffer << std::endl;
        KRATOS_INFO("GatherModelPartUtility") << buffer.str();
    }

    // Use serializer
    for (int i_rank = 0; i_rank < world_size; ++i_rank) {
        // Receiving
        if (i_rank == rank) {
            for (auto& r_bring : rEntitiesToBring) {
                const int origin_rank = r_bring.first;
                for (auto index : r_bring.second) {
                    std::string recv_buffer;
                    r_data_communicator.Recv(recv_buffer, origin_rank, static_cast<int>(index));
                    StreamSerializer serializer;
                    const auto p_serializer_buffer = dynamic_cast<std::stringstream*>(serializer.pGetBuffer());
                    p_serializer_buffer->write(recv_buffer.data(), recv_buffer.size());
                    if constexpr (std::is_same<TObjectType, Node>::value) {
                        Node::Pointer p_new_node;
                        serializer.load("bring_node_" + std::to_string(index), p_new_node);
                        KRATOS_DEBUG_ERROR_IF(rModelPart.HasNode(p_new_node->Id())) << "The node " << p_new_node->Id() << " from rank: " << origin_rank << " already exists in rank: " << rank << std::endl;
                        entities_to_bring.push_back(p_new_node);
                    } else if constexpr (std::is_same<TObjectType, Element>::value) {
                        Element::Pointer p_new_element;
                        serializer.load("bring_element_" + std::to_string(index), p_new_element);
                        KRATOS_DEBUG_ERROR_IF(rModelPart.HasElement(p_new_element->Id())) << "The element " << p_new_element->Id() << " from rank: " << origin_rank << " already exists in rank: " << rank << std::endl;
                        entities_to_bring.push_back(p_new_element);
                    } else if constexpr (std::is_same<TObjectType, Condition>::value) {
                        Condition::Pointer p_new_condition;
                        serializer.load("bring_condition_" + std::to_string(index), p_new_condition);
                        KRATOS_DEBUG_ERROR_IF(rModelPart.HasCondition(p_new_condition->Id())) << "The condition " << p_new_condition->Id() << " from rank: " << origin_rank << " already exists in rank: " << rank << std::endl;
                        entities_to_bring.push_back(p_new_condition);
                    } else {
                        KRATOS_ERROR << "Entity type not supported" << std::endl;
                    }
                }
            }
        } else { // Sending
            auto it_find = send_entities.find(i_rank);
            if (it_find != send_entities.end()) {
                for (auto index : it_find->second) {
                    StreamSerializer serializer;
                    if constexpr (std::is_same<TObjectType, Node>::value) {
                        KRATOS_DEBUG_ERROR_IF_NOT(rModelPart.HasNode(index)) << "Node with index " << index << " not found in model part" << std::endl;
                        auto p_send_node = rModelPart.pGetNode(index);
                        serializer.save("bring_node_" + std::to_string(index), p_send_node);
                    } else if constexpr (std::is_same<TObjectType, Element>::value) {
                        KRATOS_DEBUG_ERROR_IF_NOT(rModelPart.HasElement(index)) << "Element with index " << index << " not found in model part" << std::endl;
                        auto p_send_element = rModelPart.pGetElement(index);
                        serializer.save("bring_element_" + std::to_string(index), p_send_element);
                    } else if constexpr (std::is_same<TObjectType, Condition>::value) {
                        KRATOS_DEBUG_ERROR_IF_NOT(rModelPart.HasCondition(index)) << "Condition with index " << index << " not found in model part" << std::endl;
                        auto p_send_condition = rModelPart.pGetCondition(index);
                        serializer.save("bring_condition_" + std::to_string(index), p_send_condition);
                    } else {
                        KRATOS_ERROR << "Entity type not supported" << std::endl;
                    }
                    const auto p_serializer_buffer = dynamic_cast<std::stringstream*>(serializer.pGetBuffer());
                    const std::string& r_send_buffer = p_serializer_buffer->str();
                    r_data_communicator.Send(r_send_buffer, i_rank, static_cast<int>(index));
                }
            }
        }
    }

    // Add to model part
    if constexpr (std::is_same<TObjectType, Node>::value) {
        rModelPart.AddNodes(entities_to_bring.begin(), entities_to_bring.end());
    } else if constexpr (std::is_same<TObjectType, Element>::value) {
        rModelPart.AddElements(entities_to_bring.begin(), entities_to_bring.end());
    } else if constexpr (std::is_same<TObjectType, Condition>::value) {
        rModelPart.AddConditions(entities_to_bring.begin(), entities_to_bring.end());
    } else {
        KRATOS_ERROR << "Entity type not supported" << std::endl;
    }
}

template void GatherModelPartUtility::GatherEntityFromOtherPartitions<Node>(ModelPart& rModelPart, const std::map<int, std::vector<std::size_t>>& rEntitiesToBring, const int EchoLevel);
template void GatherModelPartUtility::GatherEntityFromOtherPartitions<Element>(ModelPart& rModelPart, const std::map<int, std::vector<std::size_t>>& rEntitiesToBring, const int EchoLevel);
template void GatherModelPartUtility::GatherEntityFromOtherPartitions<Condition>(ModelPart& rModelPart, const std::map<int, std::vector<std::size_t>>& rEntitiesToBring, const int EchoLevel);

std::string GatherModelPartUtility::Info() const
{
    std::stringstream buffer;
    buffer << "GatherModelPartUtility";
    return buffer.str();
}

void GatherModelPartUtility::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "GatherModelPartUtility" << std::endl;
}

void GatherModelPartUtility::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.
