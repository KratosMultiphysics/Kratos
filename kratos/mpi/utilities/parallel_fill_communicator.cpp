//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Michael Andre, https://github.com/msandre
//

#include "parallel_fill_communicator.h"

#include "includes/model_part.h"
#include "includes/data_communicator.h"
#include "processes/graph_coloring_process.h"
#include "mpi/includes/mpi_communicator.h"

namespace Kratos
{

ParallelFillCommunicator::ParallelFillCommunicator(ModelPart& r_model_part)
    : mrBaseModelPart(r_model_part)
{}

void ParallelFillCommunicator::Execute()
{
    KRATOS_TRY
    ComputeCommunicationPlan(mrBaseModelPart);
    KRATOS_CATCH("");
}

void ParallelFillCommunicator::PrintDebugInfo()
{
    PrintModelPartDebugInfo(mrBaseModelPart);
}

void ParallelFillCommunicator::PrintModelPartDebugInfo(const ModelPart& rModelPart)
{
    KRATOS_TRY

    std::cout.flush();
    const DataCommunicator& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();
    r_data_communicator.Barrier();

    int rank = r_data_communicator.Rank();

    //get number of processors
    int num_processors = r_data_communicator.Size();

    for (int i = 0; i < num_processors; i++)
    {
        if (rank == i)
        {
            std::stringstream message;

            message << " *************************************** " << std::endl;
            message << " proc = " << rank << "communication colors " << rModelPart.GetCommunicator().NeighbourIndices() << std::endl;

            //print ghost mesh
            message << " proc = " << rank << " ghost mesh" << std::endl;
            for (ModelPart::NodesContainerType::const_iterator it = rModelPart.GetCommunicator().GhostMesh().NodesBegin();
                    it != rModelPart.GetCommunicator().GhostMesh().NodesEnd();
                    ++it)
            {
                KRATOS_ERROR_IF(it->FastGetSolutionStepValue(PARTITION_INDEX)==rank)
                << "error partition index can not be = to rank for ghost nodes" << it->Id();
                message << it->Id() << " " ;
            }
            message << std::endl;

            //print local mesh
            message << " proc = " << rank << " local mesh" << std::endl;
            for (ModelPart::NodesContainerType::const_iterator it = rModelPart.GetCommunicator().LocalMesh().NodesBegin();
                    it != rModelPart.GetCommunicator().LocalMesh().NodesEnd();
                    ++it)
            {
                KRATOS_ERROR_IF(it->FastGetSolutionStepValue(PARTITION_INDEX)!=rank)
                << "error partition index can not be != from rank for local nodes" << it->Id();
                message << it->Id() << " " ;
            }
            message << std::endl;

            //print interface mesh
            message << " proc = " << rank << " interface mesh" << std::endl;
            for (ModelPart::NodesContainerType::const_iterator it = rModelPart.GetCommunicator().InterfaceMesh().NodesBegin();
                    it != rModelPart.GetCommunicator().InterfaceMesh().NodesEnd();
                    ++it)
            {
                message << it->Id() << " " ;
            }
            message << std::endl;

            //now print everything color by color
            int destination = 0;
            message << "NeighbourIndices " ;
            const vector<int>& neighbours_indices = rModelPart.GetCommunicator().NeighbourIndices();
            for (unsigned int i_color = 0; i_color < neighbours_indices.size(); i_color++)
                message << neighbours_indices[i_color] << " " ;
            message << std::endl;
            for (unsigned int i_color = 0; i_color < neighbours_indices.size(); i_color++)
            {
                message << "color = " << i_color << std::endl;
                if ((destination = neighbours_indices[i_color]) >= 0)
                {
                    message << "ghost mesh for color --> " << i_color << std::endl;
                    for (ModelPart::NodesContainerType::const_iterator it = rModelPart.GetCommunicator().GhostMesh(i_color).NodesBegin();
                            it != rModelPart.GetCommunicator().GhostMesh(i_color).NodesEnd();
                            ++it)
                    {
                        KRATOS_ERROR_IF(it->FastGetSolutionStepValue(PARTITION_INDEX)==rank)
                        << "error partition index can not be = to rank for ghost nodes" << it->Id();
                        message << it->Id() << " " ;
                    }

                    message << "finished printing ghost mesh for color --> " << i_color<< std::endl;

                    message << "local mesh for color --> " << i_color << std::endl;
                    for (ModelPart::NodesContainerType::const_iterator it = rModelPart.GetCommunicator().LocalMesh(i_color).NodesBegin();
                            it != rModelPart.GetCommunicator().LocalMesh(i_color).NodesEnd();
                            ++it)
                    {
                        KRATOS_ERROR_IF(it->FastGetSolutionStepValue(PARTITION_INDEX)!=rank)
                        << "error partition index can not be != from rank for local nodes" << it->Id();
                        message << it->Id() << " " ;
                    }
                    message << "finished printing local mesh for color --> " << i_color<< std::endl;

                    message << "interface mesh for color --> " << i_color << std::endl;
                    for (ModelPart::NodesContainerType::const_iterator it = rModelPart.GetCommunicator().InterfaceMesh(i_color).NodesBegin();
                            it != rModelPart.GetCommunicator().InterfaceMesh(i_color).NodesEnd();
                            ++it)
                    {
                        message << it->Id() << " " ;
                    }
                    message << "finished printing interface mesh for color --> " << i_color<< std::endl;
                }
                else
                {
                    if(rModelPart.GetCommunicator().GhostMesh(i_color).Nodes().size()!=0)
                    {
                        message << "rank = " << rank << " color = " << i_color << std::endl;
                        KRATOS_ERROR << "nodes found in ghost mesh when communication is not expected";
                    }
                    if(rModelPart.GetCommunicator().LocalMesh(i_color).Nodes().size()!=0)
                    {
                        message << "local mesh for color --> " << i_color << "*********************************" <<  std::endl;
                        for (ModelPart::NodesContainerType::const_iterator it = rModelPart.GetCommunicator().LocalMesh(i_color).NodesBegin();
                                it != rModelPart.GetCommunicator().LocalMesh(i_color).NodesEnd();
                                ++it)
                        {
                            KRATOS_ERROR_IF(it->FastGetSolutionStepValue(PARTITION_INDEX)!=rank)
                            << "error partition index can not be != from rank for local nodes" << it->Id();
                            message << it->Id() << " " << it->FastGetSolutionStepValue(PARTITION_INDEX) << std::endl ;
                        }
                        message << "finished printing local mesh for color --> " << i_color<< std::endl;
                        message << "nodes found in local mesh when communication is not expected" << std::endl;
                        KRATOS_ERROR << "nodes found in local mesh when communication is not expected";
                    }
                    KRATOS_ERROR_IF(rModelPart.GetCommunicator().InterfaceMesh(i_color).Nodes().size()!=0)
                    << "nodes found in interface mesh when communication is not expected";
                }
            }

            message << "finished printing proc -> " << rank << "*********************" << std::endl;
            message << std::endl;
            KRATOS_INFO_ALL_RANKS("ParallelFillCommunicator Debug Info") << message.str();

        }

        r_data_communicator.Barrier();
    }
    KRATOS_CATCH("");
}

std::string ParallelFillCommunicator::Info() const
{
    std::stringstream buffer;
    buffer << "ParallelFillCommunicator";
    return buffer.str();
}

void ParallelFillCommunicator::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "ParallelFillCommunicator" << std::endl;
}

void ParallelFillCommunicator::PrintData(std::ostream& rOStream) const
{
}


void ParallelFillCommunicator::ComputeCommunicationPlan(ModelPart& rModelPart)
{
    KRATOS_TRY;

    if (rModelPart.NumberOfNodes() > 0) {
        KRATOS_ERROR_IF_NOT(rModelPart.NodesBegin()->SolutionStepsDataHas(PARTITION_INDEX)) << "\"PARTITION_INDEX\" missing as solution step variable for nodes of ModelPart \"" << rModelPart.Name() << "\"!" << std::endl;
    }

    constexpr unsigned root_id = 0;

    Communicator::Pointer pnew_comm = Kratos::make_shared< MPICommunicator >(&rModelPart.GetNodalSolutionStepVariablesList(), DataCommunicator::GetDefault());
    rModelPart.SetCommunicator(pnew_comm);

    const auto& r_data_communicator = pnew_comm->GetDataCommunicator();

    // Check if the nodes have been assigned a partition index (i.e. some value different from 0). If not issue a warning
    int non_zero_partition_index_found = 0;
    for (const auto& r_node : rModelPart.Nodes()) {
        const int node_partition_index = r_node.FastGetSolutionStepValue(PARTITION_INDEX);
        if (node_partition_index != 0) {
            non_zero_partition_index_found = 1;
            break;
        }
    }

    non_zero_partition_index_found = r_data_communicator.SumAll(non_zero_partition_index_found);

    KRATOS_WARNING_IF("ParallelFillCommunicator", r_data_communicator.Size() > 1 && non_zero_partition_index_found == 0) << "All nodes have a PARTITION_INDEX index of 0! This could mean that PARTITION_INDEX was not assigned" << std::endl;

    // Get rank of current processor.
    const unsigned my_rank = r_data_communicator.Rank();

    // Get number of processors.
    const unsigned num_processors = r_data_communicator.Size();
    // Find all ghost nodes on this process and mark the corresponding neighbour process for communication.
    vector<bool> receive_from_neighbour(num_processors, false);
    for (const auto& rNode : rModelPart.Nodes())
    {
        const unsigned partition_index = rNode.FastGetSolutionStepValue(PARTITION_INDEX);
        KRATOS_ERROR_IF(partition_index >= num_processors) << "The partition index is out of range. Invalid model part." << std::endl;
        if(partition_index != my_rank)
            receive_from_neighbour[partition_index] = true;
    }

    // Make a list of my receive process ids.
    std::vector<unsigned> my_receive_neighbours;
    my_receive_neighbours.reserve(30);
    for (unsigned p_id = 0; p_id < num_processors; ++p_id)
    {
        if (receive_from_neighbour[p_id])
            my_receive_neighbours.push_back(p_id);
    }

    // Initialize arrays for all neighbour id lists on root process.
    std::vector<std::size_t> number_of_receive_neighbours;
    std::vector<std::vector<unsigned>> receive_neighbours;
    if (my_rank == root_id)
    {
        number_of_receive_neighbours.resize(num_processors);
        receive_neighbours.resize(num_processors);
    }
    {
        std::vector<std::size_t> send_buf{my_receive_neighbours.size()};
        r_data_communicator.Gather(send_buf, number_of_receive_neighbours, root_id);
    }
    if (my_rank == root_id)
    {
        for (unsigned p_id = 0; p_id < num_processors; ++p_id)
            receive_neighbours[p_id].resize(number_of_receive_neighbours[p_id]);
    }

    // Fill the neighbour id lists of the partitions on root.
    if (my_rank == root_id) // On root we directly copy the data without calling MPI.
        std::copy(my_receive_neighbours.begin(), my_receive_neighbours.end(), receive_neighbours[root_id].begin());
    // Gather the remaining id lists to root.
    for (unsigned p_id = 1; p_id < num_processors; ++p_id)
    {
        if (my_rank == root_id)
        {
            r_data_communicator.Recv(receive_neighbours[p_id], p_id, p_id);
        }
        else if (my_rank == p_id)
        {
            r_data_communicator.Send(my_receive_neighbours, root_id, p_id);
        }
    }

    // Create the colored graph for communication.
    DenseMatrix<int> domains_colored_graph;
    int max_color_found = -1;
    if (my_rank == root_id)
    {
        ///@TODO for large problems, this should use a compressed matrix.
        DenseMatrix<int> domains_graph = ZeroMatrix(num_processors, num_processors);
        for (unsigned index1 = 0; index1 < num_processors; ++index1)
            for (unsigned index2 : receive_neighbours[index1])
            {
                KRATOS_ERROR_IF(index1 == index2) << "Trying to communicate with the node itself." << std::endl;
                domains_graph(index1, index2) = 1;
                domains_graph(index2, index1) = 1;
            }

        // max_color is overwritten by the GraphColoringProcess.
        int max_color = 2 * num_processors; // Max. number of one-directional communications (this has no effect).
        GraphColoringProcess coloring_process(num_processors, domains_graph, domains_colored_graph, max_color);
        coloring_process.Execute();
        // Count max colors.
        for (unsigned p_id = 0; p_id < num_processors; ++p_id)
            for (int j = 0; j < max_color; ++j)
                if (domains_colored_graph(p_id, j) != -1 && max_color_found < j) max_color_found = j;

        max_color_found += 1;
    }

    // Broadcast max_color_found.
    r_data_communicator.Broadcast(max_color_found, root_id);

    // Now send the colors of the communication to the processors.
    std::vector<int> colors(max_color_found);
    if (my_rank == root_id) // On root we directly copy the data.
    {
        for (int j = 0; j < max_color_found; ++j)
            colors[j] = domains_colored_graph(root_id, j);
    }
    // Send the remaining color patterns to processes.
    std::vector<int> send_colors(max_color_found);
    for (unsigned p_id = 1; p_id < num_processors; ++p_id)
    {
        if (my_rank == root_id)
        {
            for (int j = 0; j < max_color_found; ++j)
            {
                send_colors[j] = domains_colored_graph(p_id, j);
            }
            r_data_communicator.Send(send_colors, p_id, p_id);
        }
        else if (my_rank == p_id)
        {
            r_data_communicator.Recv(colors, root_id, p_id);
        }
    }

    InitializeParallelCommunicationMeshes(rModelPart, colors, my_rank);

    KRATOS_CATCH("");
}

void ParallelFillCommunicator::InitializeParallelCommunicationMeshes(
    ModelPart& rModelPart,
    const std::vector<int>& rColors,
    unsigned MyRank)
{
    KRATOS_TRY;
    // Allocate space needed in the communicator.
    rModelPart.GetCommunicator().SetNumberOfColors(rColors.size());
    rModelPart.GetCommunicator().NeighbourIndices().resize(rColors.size());
    for (unsigned color = 0; color < rColors.size(); ++color)
    {
        rModelPart.GetCommunicator().LocalMesh(color).Nodes().clear();
        rModelPart.GetCommunicator().GhostMesh(color).Nodes().clear();
        rModelPart.GetCommunicator().InterfaceMesh(color).Nodes().clear();
    }

    // For each color fill the list of ghost and local nodes and the
    // interface
    // mesh.
    for (unsigned color = 0; color < rColors.size(); ++color)
    {
        rModelPart.GetCommunicator().NeighbourIndices()[color] = rColors[color];
        GenerateMeshes(rColors[color], MyRank, color, rModelPart);
    }

    // Fill the list of all of the nodes to be communicated.
    ModelPart::NodesContainerType& r_local_nodes = rModelPart.GetCommunicator().LocalMesh().Nodes();
    ModelPart::NodesContainerType& r_ghost_nodes = rModelPart.GetCommunicator().GhostMesh().Nodes();
    ModelPart::NodesContainerType& r_interface_nodes = rModelPart.GetCommunicator().InterfaceMesh().Nodes();
    r_local_nodes.clear();
    r_ghost_nodes.clear();
    r_interface_nodes.clear();

    KRATOS_ERROR_IF(r_local_nodes.size() != 0)
        << "Local size can't be zero." << std::endl;

    // Fill nodes for LocalMesh and GhostMesh.
    for (auto it_node = rModelPart.NodesBegin(); it_node != rModelPart.NodesEnd(); ++it_node)
    {
        const unsigned index = it_node->FastGetSolutionStepValue(PARTITION_INDEX);
        if (index == MyRank)
        {
            r_local_nodes.push_back(*(it_node.base()));
        }
        else
        {
            r_ghost_nodes.push_back(*(it_node.base()));
        }
    }

    // Fill nodes for the InterfaceMesh.
    for (ModelPart::MeshType& r_color_interface_mesh : rModelPart.GetCommunicator().InterfaceMeshes())
    {
        ModelPart::NodesContainerType& r_color_interface_nodes =
            r_color_interface_mesh.Nodes();
        for (auto it = r_color_interface_nodes.ptr_begin(); it != r_color_interface_nodes.ptr_end(); ++it)
        {
            r_interface_nodes.push_back(*it);
        }
    }
    r_interface_nodes.Unique();
    r_local_nodes.Unique();
    r_ghost_nodes.Unique();

    // Assign elements and conditions for LocalMesh.
    rModelPart.GetCommunicator().LocalMesh().Elements().clear();
    rModelPart.GetCommunicator().LocalMesh().Conditions().clear();
    rModelPart.GetCommunicator().LocalMesh().Elements() = rModelPart.Elements();
    rModelPart.GetCommunicator().LocalMesh().Conditions() = rModelPart.Conditions();

    // Call the sub model part.
    for (ModelPart& r_sub_model_part : rModelPart.SubModelParts())
    {
        ComputeCommunicationPlan(r_sub_model_part);
    }

    KRATOS_CATCH("");
}

void ParallelFillCommunicator::GenerateMeshes(int NeighbourPID, int MyPID, unsigned Color, ModelPart& rModelPart)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(NeighbourPID == MyPID)
        << "Neighbour process coincides with rank! this should not happen."
        << std::endl;

    if (NeighbourPID == -1) // Don't communicate with this process.
        return;

    const DataCommunicator& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();

    ModelPart::NodesContainerType& r_ghost_nodes =
        rModelPart.GetCommunicator().GhostMesh(Color).Nodes();
    r_ghost_nodes.clear();

    // Fill nodes for GhostMesh(Color).
    for (auto it = rModelPart.NodesBegin(); it != rModelPart.NodesEnd(); ++it)
    {
        const int index = it->FastGetSolutionStepValue(PARTITION_INDEX);
        if (index == NeighbourPID)
        {
            r_ghost_nodes.push_back(*(it.base()));
        }
    }
    unsigned num_ghost_nodes = r_ghost_nodes.size();
    r_ghost_nodes.Unique();
    KRATOS_ERROR_IF(num_ghost_nodes != r_ghost_nodes.size())
        << "The list of nodes to receive has repeated nodes." << std::endl;

    std::vector<int> ids_to_receive(r_ghost_nodes.size());

    { // Fill receive ids (ids of ghost nodes).
        int i = 0;
        for (const ModelPart::NodeType& rNode : r_ghost_nodes)
        {
            ids_to_receive[i++] = rNode.Id();
        }
    }

    std::vector<int> ids_to_send;
    { // Syncronize how many nodes need to be sent/received.
        int send_tag = Color;
        int receive_tag = Color;
        std::size_t recv_buf = r_data_communicator.SendRecv(ids_to_receive.size(), NeighbourPID, send_tag, NeighbourPID, receive_tag);
        ids_to_send.resize(recv_buf);

        // Send/receive node ids.
        r_data_communicator.SendRecv(ids_to_receive, NeighbourPID, send_tag, ids_to_send, NeighbourPID, receive_tag);
    }

    // Fill nodes for LocalMesh(Color).
    ModelPart::NodesContainerType& r_local_nodes =
        rModelPart.GetCommunicator().LocalMesh(Color).Nodes();
    r_local_nodes.clear();
    for (int id : ids_to_send)
    {
        KRATOS_DEBUG_ERROR_IF(rModelPart.Nodes().find(id) == rModelPart.Nodes().end()) << "Trying to add Node with Id #" << id << " to the local mesh, but the node does not exist in the ModelPart!" << std::endl;
        r_local_nodes.push_back(rModelPart.Nodes()(id));
    }

    for (const ModelPart::NodeType& r_node : r_local_nodes)
    {
        KRATOS_ERROR_IF(r_node.FastGetSolutionStepValue(PARTITION_INDEX) != MyPID) << "A node in the local mesh is trying to communicate to the wrong partition."
                                                                                    << std::endl;
    }

    r_local_nodes.Unique();
    KRATOS_ERROR_IF(r_local_nodes.size() != ids_to_send.size())
        << "Impossible situation. Something went wrong." << std::endl;

    // Fill InterfaceMesh(Color) with local and ghost nodes.
    ModelPart::NodesContainerType& r_interface_nodes =
        rModelPart.GetCommunicator().InterfaceMesh(Color).Nodes();
    r_interface_nodes.clear();

    for (auto it = r_ghost_nodes.begin(); it != r_ghost_nodes.end(); ++it)
    {
        r_interface_nodes.push_back(*(it.base()));
    }

    for (auto it = r_local_nodes.begin(); it != r_local_nodes.end(); it++)
    {
        r_interface_nodes.push_back(*(it.base()));
    }

    unsigned num_interface_nodes = r_interface_nodes.size();
    r_interface_nodes.Unique();
    KRATOS_ERROR_IF(num_interface_nodes != r_interface_nodes.size())
        << "Something went wrong in the interface nodes." << std::endl;

    KRATOS_CATCH("");
}

}
