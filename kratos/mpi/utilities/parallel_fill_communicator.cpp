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
//

// System includes

// External includes

// Project includes
#include "parallel_fill_communicator.h"
#include "includes/model_part.h"
#include "processes/graph_coloring_process.h"
#include "mpi/includes/mpi_communicator.h"

namespace Kratos
{

ParallelFillCommunicator::ParallelFillCommunicator(ModelPart& rModelPart)
    : FillCommunicator(rModelPart)
{}

ParallelFillCommunicator::ParallelFillCommunicator(
    ModelPart& rModelPart,
    const DataCommunicator& rDataCommunicator)
    : FillCommunicator(rModelPart, rDataCommunicator)
{
    KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "Trying to create an ParallelFillCommunicator with a non-distributed DataCommunicator!" << std::endl;
}

void ParallelFillCommunicator::Execute()
{
    KRATOS_TRY
    mPartitionIndexCheckPerformed = false;
    auto& r_base_model_part = GetBaseModelPart();
    ComputeCommunicationPlan(r_base_model_part);

    // Depending of the echo level, print the debug info
    if (this->GetEchoLevel() == FillCommunicatorEchoLevel::DEBUG_INFO) {
        PrintModelPartDebugInfo(r_base_model_part);
    }

    KRATOS_CATCH("");
}

void ParallelFillCommunicator::PrintModelPartDebugInfo(const ModelPart& rModelPart)
{
    KRATOS_TRY

    std::cout.flush();
    const auto& r_communicator = rModelPart.GetCommunicator();
    const auto& r_data_communicator = r_communicator.GetDataCommunicator();
    r_data_communicator.Barrier();

    int rank = r_data_communicator.Rank();

    // Get number of processors
    int num_processors = r_data_communicator.Size();

    for (int i = 0; i < num_processors; i++) {
        if (rank == i) {
            std::stringstream message;

            message << " *************************************** " << std::endl;
            message << " Proc = " << rank << "communication colors " << r_communicator.NeighbourIndices() << std::endl;

            // Print ghost mesh
            message << " Proc = " << rank << " ghost mesh" << std::endl;
            for (auto it = r_communicator.GhostMesh().NodesBegin(); it != r_communicator.GhostMesh().NodesEnd(); ++it) {
                KRATOS_ERROR_IF(it->FastGetSolutionStepValue(PARTITION_INDEX)==rank) << "Error partition index can not be = to rank for ghost nodes" << it->Id() << std::endl;
                message << it->Id() << " " ;
            }
            message << std::endl;

            // Print local mesh
            message << " proc = " << rank << " local mesh" << std::endl;
            for (auto it = r_communicator.LocalMesh().NodesBegin(); it != r_communicator.LocalMesh().NodesEnd(); ++it) {
                KRATOS_ERROR_IF(it->FastGetSolutionStepValue(PARTITION_INDEX)!=rank) << "Error partition index can not be != from rank for local nodes" << it->Id() << std::endl;
                message << it->Id() << " " ;
            }
            message << std::endl;

            // Print interface mesh
            message << " Proc = " << rank << " interface mesh" << std::endl;
            for (auto it = r_communicator.InterfaceMesh().NodesBegin(); it != r_communicator.InterfaceMesh().NodesEnd(); ++it) {
                message << it->Id() << " " ;
            }
            message << std::endl;

            //now print everything color by color
            int destination = 0;
            message << "NeighbourIndices " ;
            const auto& r_neighbours_indices = r_communicator.NeighbourIndices();
            for (unsigned int i_color = 0; i_color < r_neighbours_indices.size(); i_color++)
                message << r_neighbours_indices[i_color] << " " ;
            message << std::endl;
            for (unsigned int i_color = 0; i_color < r_neighbours_indices.size(); i_color++) {
                message << "Color = " << i_color << std::endl;
                if ((destination = r_neighbours_indices[i_color]) >= 0) {
                    message << "Ghost mesh for color --> " << i_color << std::endl;
                    for (auto it = r_communicator.GhostMesh(i_color).NodesBegin(); it != r_communicator.GhostMesh(i_color).NodesEnd(); ++it) {
                        KRATOS_ERROR_IF(it->FastGetSolutionStepValue(PARTITION_INDEX)==rank) << "Error partition index can not be = to rank for ghost nodes" << it->Id() << std::endl;
                        message << it->Id() << " " ;
                    }

                    message << "Finished printing ghost mesh for color --> " << i_color<< std::endl;

                    message << "Local mesh for color --> " << i_color << std::endl;
                    for (auto it = r_communicator.LocalMesh(i_color).NodesBegin(); it != r_communicator.LocalMesh(i_color).NodesEnd(); ++it) {
                        KRATOS_ERROR_IF(it->FastGetSolutionStepValue(PARTITION_INDEX)!=rank) << "Error partition index can not be != from rank for local nodes" << it->Id() << std::endl;
                        message << it->Id() << " " ;
                    }
                    message << "Finished printing local mesh for color --> " << i_color<< std::endl;

                    message << "Interface mesh for color --> " << i_color << std::endl;
                    for (auto it = r_communicator.InterfaceMesh(i_color).NodesBegin(); it != r_communicator.InterfaceMesh(i_color).NodesEnd(); ++it) {
                        message << it->Id() << " " ;
                    }
                    message << "Finished printing interface mesh for color --> " << i_color<< std::endl;
                } else {
                    if(r_communicator.GhostMesh(i_color).Nodes().size()!=0) {
                        message << "Rank = " << rank << " color = " << i_color << std::endl;
                        KRATOS_ERROR << "Nodes found in ghost mesh when communication is not expected" << std::endl;
                    }
                    if(r_communicator.LocalMesh(i_color).Nodes().size()!=0) {
                        message << "Local mesh for color --> " << i_color << "*********************************" <<  std::endl;
                        for (auto it = r_communicator.LocalMesh(i_color).NodesBegin(); it != r_communicator.LocalMesh(i_color).NodesEnd(); ++it) {
                            KRATOS_ERROR_IF(it->FastGetSolutionStepValue(PARTITION_INDEX)!=rank) << "Error partition index can not be != from rank for local nodes" << it->Id() << std::endl;
                            message << it->Id() << " " << it->FastGetSolutionStepValue(PARTITION_INDEX) << std::endl ;
                        }
                        message << "Finished printing local mesh for color --> " << i_color<< std::endl;
                        message << "Nodes found in local mesh when communication is not expected" << std::endl;
                        KRATOS_ERROR << "Nodes found in local mesh when communication is not expected" << std::endl;
                    }
                    KRATOS_ERROR_IF(r_communicator.InterfaceMesh(i_color).Nodes().size()!=0)
                    << "Nodes found in interface mesh when communication is not expected" << std::endl;
                }
            }

            message << "Finished printing proc -> " << rank << "*********************\n" << std::endl;
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

    constexpr unsigned root_id = 0;

    Communicator::Pointer pnew_comm = Kratos::make_shared< MPICommunicator >(&rModelPart.GetNodalSolutionStepVariablesList(), mrDataComm);
    rModelPart.SetCommunicator(pnew_comm);

    // Check if the nodes have been assigned a partition index (i.e. some value different from 0). If not issue a warning
    if (!mPartitionIndexCheckPerformed) { // needs to be protected bcs this function is called recursively for SubModelParts
        mPartitionIndexCheckPerformed = true;

        if (rModelPart.NumberOfNodes() > 0) {
            KRATOS_ERROR_IF_NOT(rModelPart.NodesBegin()->SolutionStepsDataHas(PARTITION_INDEX)) << "\"PARTITION_INDEX\" missing as solution step variable for nodes of ModelPart \"" << rModelPart.Name() << "\"!" << std::endl;
        }

        bool non_zero_partition_index_found = false;
        for (const auto& r_node : rModelPart.Nodes()) {
            const int node_partition_index = r_node.FastGetSolutionStepValue(PARTITION_INDEX);
            if (node_partition_index != 0) {
                non_zero_partition_index_found = true;
                break;
            }
        }

        non_zero_partition_index_found = mrDataComm.OrReduceAll(non_zero_partition_index_found);

        KRATOS_WARNING_IF("ParallelFillCommunicator", mrDataComm.Size() > 1 && !non_zero_partition_index_found) << "All nodes have a PARTITION_INDEX index of 0! This could mean that PARTITION_INDEX was not assigned" << std::endl;
    }

    // Get rank of current processor.
    const int my_rank = mrDataComm.Rank();

    // Get number of processors.
    const int num_processors = mrDataComm.Size();
    // Find all ghost nodes on this process and mark the corresponding neighbour process for communication.
    DenseVector<bool> receive_from_neighbour(num_processors, false);
    for (const auto& rNode : rModelPart.Nodes())
    {
        const int partition_index = rNode.FastGetSolutionStepValue(PARTITION_INDEX);
        KRATOS_ERROR_IF(partition_index >= num_processors) << "The partition index (" << partition_index << ") is out of range (size = " << num_processors << "). Invalid model part." << std::endl;
        if(partition_index != my_rank)
            receive_from_neighbour[partition_index] = true;
    }

    // Make a list of my receive process ids.
    std::vector<int> my_receive_neighbours;
    my_receive_neighbours.reserve(30);
    for (int p_id = 0; p_id < num_processors; ++p_id)
    {
        if (receive_from_neighbour[p_id])
            my_receive_neighbours.push_back(p_id);
    }

    // Initialize arrays for all neighbour id lists on root process.
    std::vector<std::size_t> number_of_receive_neighbours;
    std::vector<std::vector<int>> receive_neighbours;
    if (my_rank == root_id)
    {
        number_of_receive_neighbours.resize(num_processors);
        receive_neighbours.resize(num_processors);
    }
    {
        std::vector<std::size_t> send_buf{my_receive_neighbours.size()};
        mrDataComm.Gather(send_buf, number_of_receive_neighbours, root_id);
    }
    if (my_rank == root_id)
    {
        for (int p_id = 0; p_id < num_processors; ++p_id)
            receive_neighbours[p_id].resize(number_of_receive_neighbours[p_id]);
    }

    // Fill the neighbour id lists of the partitions on root.
    if (my_rank == root_id) // On root we directly copy the data without calling MPI.
        std::copy(my_receive_neighbours.begin(), my_receive_neighbours.end(), receive_neighbours[root_id].begin());
    // Gather the remaining id lists to root.
    for (int p_id = 1; p_id < num_processors; ++p_id)
    {
        if (my_rank == root_id)
        {
            mrDataComm.Recv(receive_neighbours[p_id], p_id, p_id);
        }
        else if (my_rank == p_id)
        {
            mrDataComm.Send(my_receive_neighbours, root_id, p_id);
        }
    }

    // Create the colored graph for communication.
    DenseMatrix<int> domains_colored_graph;
    int max_color_found = -1;
    if (my_rank == root_id)
    {
        ///@TODO for large problems, this should use a compressed matrix.
        DenseMatrix<int> domains_graph = ZeroMatrix(num_processors, num_processors);
        for (int index1 = 0; index1 < num_processors; ++index1)
            for (int index2 : receive_neighbours[index1])
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
        for (int p_id = 0; p_id < num_processors; ++p_id)
            for (int j = 0; j < max_color; ++j)
                if (domains_colored_graph(p_id, j) != -1 && max_color_found < j) max_color_found = j;

        max_color_found += 1;
    }

    // Broadcast max_color_found.
    mrDataComm.Broadcast(max_color_found, root_id);

    // Now send the colors of the communication to the processors.
    std::vector<int> colors(max_color_found);
    if (my_rank == root_id) // On root we directly copy the data.
    {
        for (int j = 0; j < max_color_found; ++j)
            colors[j] = domains_colored_graph(root_id, j);
    }
    // Send the remaining color patterns to processes.
    std::vector<int> send_colors(max_color_found);
    for (int p_id = 1; p_id < num_processors; ++p_id)
    {
        if (my_rank == root_id)
        {
            for (int j = 0; j < max_color_found; ++j)
            {
                send_colors[j] = domains_colored_graph(p_id, j);
            }
            mrDataComm.Send(send_colors, p_id, p_id);
        }
        else if (my_rank == p_id)
        {
            mrDataComm.Recv(colors, root_id, p_id);
        }
    }

    InitializeParallelCommunicationMeshes(rModelPart, colors, my_rank);

    KRATOS_CATCH("");
}

void ParallelFillCommunicator::InitializeParallelCommunicationMeshes(
    ModelPart& rModelPart,
    const std::vector<int>& rColors,
    const int MyRank
    )
{
    KRATOS_TRY;
    // Allocate space needed in the communicator.
    rModelPart.GetCommunicator().SetNumberOfColors(rColors.size());
    rModelPart.GetCommunicator().NeighbourIndices().resize(rColors.size());
    for (unsigned color = 0; color < rColors.size(); ++color) {
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
        const int index = it_node->FastGetSolutionStepValue(PARTITION_INDEX);
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
    rModelPart.GetCommunicator().LocalMesh().MasterSlaveConstraints().clear();
    rModelPart.GetCommunicator().LocalMesh().Elements() = rModelPart.Elements();
    rModelPart.GetCommunicator().LocalMesh().Conditions() = rModelPart.Conditions();
    rModelPart.GetCommunicator().LocalMesh().MasterSlaveConstraints() = rModelPart.MasterSlaveConstraints();

    // Call the sub model part. We don't need to fully compute info in the sub model part, as the info is already computed in the parent model part. This reduces also number of communication calls.
    for (ModelPart& r_sub_model_part : rModelPart.SubModelParts()) {
        ComputeCommunicationPlanSubModelPart(r_sub_model_part, rColors);
    }

    KRATOS_CATCH("");
}

void ParallelFillCommunicator::ComputeCommunicationPlanSubModelPart(
    ModelPart& rSubModelPart,
    const std::vector<int>& rColors
    )
{
    KRATOS_TRY;

    // Set communicator for the sub model part.
    Communicator::Pointer pnew_comm = Kratos::make_shared< MPICommunicator >(&rSubModelPart.GetNodalSolutionStepVariablesList(), mrDataComm);
    rSubModelPart.SetCommunicator(pnew_comm);

    // Get the parent model part.
    ModelPart& r_parent_model_part = rSubModelPart.GetParentModelPart();

    // Get communicator
    auto& r_communicator = rSubModelPart.GetCommunicator();
    auto& r_parent_communicator = r_parent_model_part.GetCommunicator();

    // Create a local mesh for the sub model part.
    Communicator::MeshType::Pointer p_local_mesh = Kratos::make_shared<Communicator::MeshType>();
    r_communicator.SetLocalMesh(p_local_mesh);

    // Fill the list of all of the nodes to be communicated.
    auto& r_local_nodes = p_local_mesh->Nodes();
    auto& r_ghost_nodes = r_communicator.GhostMesh().Nodes();
    auto& r_interface_nodes = r_communicator.InterfaceMesh().Nodes();
    r_ghost_nodes.clear();
    r_interface_nodes.clear();

    // Parent model part nodes
    const auto& r_parent_local_nodes = r_parent_communicator.LocalMesh().Nodes();
    const auto& r_parent_ghost_nodes = r_parent_communicator.GhostMesh().Nodes();
    const auto& r_parent_interface_nodes = r_parent_communicator.InterfaceMesh().Nodes();

    // Interface nodes we use the information from the parent model part.
    for (auto it_node = rSubModelPart.NodesBegin(); it_node != rSubModelPart.NodesEnd(); ++it_node) {
        const std::size_t index = it_node->Id();
        if (r_parent_local_nodes.find(index) != r_parent_local_nodes.end()) {
            r_local_nodes.push_back(*(it_node.base()));
        }
        if (r_parent_ghost_nodes.find(index) != r_parent_ghost_nodes.end()) {
            r_ghost_nodes.push_back(*(it_node.base()));
        }
        if (r_parent_interface_nodes.find(index) != r_parent_interface_nodes.end()) {
            r_interface_nodes.push_back(*(it_node.base()));
        }
    }

    // // Calling Unique() on the nodes container will remove duplicates.
    // // NOTE: Not required as already computed in base model part.
    // r_interface_nodes.Unique();
    // r_local_nodes.Unique();
    // r_ghost_nodes.Unique();

    // Allocate space needed in the communicator.
    r_communicator.SetNumberOfColors(rColors.size());
    r_communicator.NeighbourIndices().resize(rColors.size());
    for (unsigned color = 0; color < rColors.size(); ++color) {
        r_communicator.LocalMesh(color).Nodes().clear();
        r_communicator.GhostMesh(color).Nodes().clear();
        r_communicator.InterfaceMesh(color).Nodes().clear();
    }

    // For each color fill the list of ghost and local nodes and the  interface mesh.
    for (unsigned color = 0; color < rColors.size(); ++color) {
        // Set the neighbour index.
        r_communicator.NeighbourIndices()[color] = rColors[color];

        // Get nodes for GhostMesh(Color) and LocalMesh(Color).
        auto& r_ghost_nodes = r_communicator.GhostMesh(color).Nodes();
        r_ghost_nodes.clear();
        auto& r_local_nodes = r_communicator.LocalMesh(color).Nodes();
        r_local_nodes.clear();

        // Get the parent model part meshed nodes.
        const auto& r_parent_ghost_nodes = r_parent_communicator.GhostMesh(color).Nodes();
        const auto& r_parent_local_nodes = r_parent_communicator.LocalMesh(color).Nodes();

        // Fill nodes for LocalMesh(Color) and GhostMesh(Color).
        for (auto it_node = rSubModelPart.NodesBegin(); it_node != rSubModelPart.NodesEnd(); ++it_node) {
            const std::size_t index = it_node->Id();
            if (r_parent_ghost_nodes.find(index) != r_parent_ghost_nodes.end()) {
                r_ghost_nodes.push_back(*(it_node.base()));
            }
            if (r_parent_local_nodes.find(index) != r_parent_local_nodes.end()) {
                r_local_nodes.push_back(*(it_node.base()));
            }
        }

        // Fill InterfaceMesh(Color) with local and ghost nodes.
        auto& r_interface_nodes = r_communicator.InterfaceMesh(color).Nodes();
        r_interface_nodes.clear();
        r_interface_nodes.reserve(r_local_nodes.size() + r_ghost_nodes.size());

        for (auto it = r_ghost_nodes.begin(); it != r_ghost_nodes.end(); ++it) {
            r_interface_nodes.push_back(*(it.base()));
        }

        for (auto it = r_local_nodes.begin(); it != r_local_nodes.end(); it++) {
            r_interface_nodes.push_back(*(it.base()));
        }
    }

    // Assign elements and conditions for LocalMesh.
    auto& r_submodel_part_mesh = rSubModelPart.GetMesh();
    if (r_submodel_part_mesh.pElements()) {
        p_local_mesh->SetElements(r_submodel_part_mesh.pElements());
    }
    if (r_submodel_part_mesh.pConditions()) {
        p_local_mesh->SetConditions(r_submodel_part_mesh.pConditions());
    }
    p_local_mesh->MasterSlaveConstraints() = rSubModelPart.MasterSlaveConstraints();

    // Call the sub model part.
    for (ModelPart& r_sub_model_part : rSubModelPart.SubModelParts()) {
        ComputeCommunicationPlanSubModelPart(r_sub_model_part, rColors);
    }

    KRATOS_CATCH("");
}

void ParallelFillCommunicator::GenerateMeshes(
    const int NeighbourPID, 
    const int MyPID, 
    const unsigned int Color, 
    ModelPart& rModelPart
    )
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
    { // Synchronize how many nodes need to be sent/received.
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
    r_interface_nodes.reserve(r_local_nodes.size() + r_ghost_nodes.size());

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
