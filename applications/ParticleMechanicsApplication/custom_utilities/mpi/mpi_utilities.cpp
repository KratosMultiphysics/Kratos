//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Meßmer
//

// Kratos includes
#include "mpi_utilities.h"
#include "particle_mechanics_application_variables.h"

namespace Kratos {

    void MPM_MPI_Utilities::TransferElements( ModelPart& rModelPart, std::vector<ElementsContainerType>& rSendElements, std::vector<ElementsContainerType>& rRecvElements)
    {
        const unsigned int rank = rModelPart.GetCommunicator().MyPID();
        const unsigned int size = rModelPart.GetCommunicator().TotalProcesses();

        KRATOS_ERROR_IF(rSendElements.size() != size)
            << "MPM_MPI_Utilities::TransferElements: " << "Size of element container: " << rSendElements.size()
            << " does not match the number of processes: " << size << "!" << std::endl;

        KRATOS_WARNING_IF("MPM_MPI_Utilities::TransferElements", rSendElements[rank].size() != 0)
            << "Rank: " << rank << " holds a sending element container with one or more elements with itself "
            << "(Rank: " << rank << ") as destination." << std::endl
            << "Theses elements will be cleared." << std::endl;

        // Exchange elements
        rModelPart.GetCommunicator().TransferObjects(rSendElements, rRecvElements);

        // Remove sent elements from current ModelPart
        for( unsigned int i = 0; i < rSendElements.size(); ++i){
            for(auto it = rSendElements[i].begin(); it != rSendElements[i].end(); ++it){
                it->GetGeometry().clear();
                it->Reset(ACTIVE);
                rModelPart.RemoveElementFromAllLevels(it->Id());
                rModelPart.GetCommunicator().LocalMesh().RemoveElement(it->Id());
            }
        }

        rSendElements.clear();
    }

    void MPM_MPI_Utilities::TransferConditions(ModelPart& rModelPart, std::vector<ConditionsContainerType>& rSendConditions, std::vector<ConditionsContainerType>& rRecvConditions)
    {
        const unsigned int rank = rModelPart.GetCommunicator().MyPID();
        const unsigned int size = rModelPart.GetCommunicator().TotalProcesses();

        KRATOS_ERROR_IF(rSendConditions.size() != size)
            << "MPM_MPI_Utilities::TransferConditions: " << "Size of condition container: " << rSendConditions.size()
            << " does not match the number of processes: " << size << "!" << std::endl;

        KRATOS_WARNING_IF("MPM_MPI_Utilities::TransferConditions", rSendConditions[rank].size() != 0)
            << "Rank: " << rank << " holds a sending condition container with one or more conditions with itself "
            << "(Rank: " << rank << ") as destination." << std::endl
            << "Theses conditions will be cleared." << std::endl;

        // Exchange conditions
        rModelPart.GetCommunicator().TransferObjects(rSendConditions, rRecvConditions);

        // Remove sent conditions from current ModelPart
        for( unsigned int i = 0; i < rSendConditions.size(); ++i){
            for(auto it = rSendConditions[i].begin(); it != rSendConditions[i].end(); ++it){
                it->GetGeometry().clear();
                it->Reset(ACTIVE);
                rModelPart.RemoveConditionFromAllLevels(it->Id());
                rModelPart.GetCommunicator().LocalMesh().RemoveCondition(it->Id());
            }
        }
        rSendConditions.clear();
    }

    void MPM_MPI_Utilities::SynchronizeActiveDofsAtInterface(ModelPart& rModelPart){
        const unsigned int size = rModelPart.GetCommunicator().TotalProcesses();

        // Get neigbours of each rank
        auto neighbours = rModelPart.GetCommunicator().NeighbourIndices();
        std::vector<int> neighbours_custom;
        for( auto neighbour : neighbours){
            if(neighbour >= 0)
                neighbours_custom.push_back(neighbour);
        }

        // Fill send nodes containers
        std::vector<NodesContainerType> send_nodes_container(size);
        NodesContainerType send_nodes;

        // Send active ghost nodes to their home partition.
        const unsigned int number_of_ghost_nodes = rModelPart.GetCommunicator().GhostMesh().NumberOfNodes();
        const auto ghost_nodes_it_begin = rModelPart.GetCommunicator().GhostMesh().NodesBegin();
        for( unsigned int i = 0; i < number_of_ghost_nodes; ++i){
            auto node_it = ghost_nodes_it_begin + i;
            if( node_it->Is(ACTIVE) ){
                const unsigned int partition = node_it->FastGetSolutionStepValue(PARTITION_INDEX);
                send_nodes_container[partition].push_back(*node_it.base());
            }
        }

        // Transfer nodes among proc's
        std::vector<NodesContainerType> recv_nodes_container(size);
        rModelPart.GetCommunicator().TransferObjects(send_nodes_container, recv_nodes_container);

        // Set recieved nodes active
        // Attention: Right now this function only requires the nodal id;
        // TODO: Only send nodal id. Keep this for a while, in case more information is required.
        auto& r_interface_nodes = rModelPart.GetCommunicator().InterfaceMesh().Nodes();
        for( auto neighbour : neighbours_custom){
            for( auto it = recv_nodes_container[neighbour].begin(); it != recv_nodes_container[neighbour].end(); ++it){
                const unsigned int node_id = it->Id();
                auto found_node_it = r_interface_nodes.find(node_id);
                if( found_node_it != r_interface_nodes.end() ){
                    found_node_it->Set(ACTIVE, true);
                }
            }
        }
    }

    // ##############################################################################################
    // !!!This function is no longer needed. However, function is kept for some while just in case.!!!
    // ##############################################################################################
    void MPM_MPI_Utilities::SynchronizeNodalDisplacementAtInterface(ModelPart& rModelPart)
    {
        const unsigned int rank = rModelPart.GetCommunicator().MyPID();
        const unsigned int size = rModelPart.GetCommunicator().TotalProcesses();

        // Get neigbours of each rank
        auto neighbours = rModelPart.GetCommunicator().NeighbourIndices();
        std::vector<int> neighbours_custom;
        for( auto neighbour : neighbours){
            if(neighbour >= 0)
                neighbours_custom.push_back(neighbour);
        }

        // Fill send nodes containers
        std::vector<NodesContainerType> send_nodes_container(size);
        NodesContainerType send_nodes;

        auto& r_ghost_nodes = rModelPart.GetCommunicator().GhostMesh().Nodes();
        const auto ghost_nodes_it_begin = rModelPart.GetCommunicator().GhostMesh().NodesBegin();
        const unsigned int number_ghost_nodes = rModelPart.GetCommunicator().GhostMesh().NumberOfNodes();

        // for( int i = 0; i < number_ghost_nodes; ++i){
        //     auto ghost_node_it = ghost_nodes_it_begin + i;
        //     ghost_node_it->SetValue(IS_GHOST_NODE, true);
        // }

        auto node_it_begin = rModelPart.GetCommunicator().InterfaceMesh().NodesBegin();
        const unsigned int number_of_interface_nodes = rModelPart.GetCommunicator().InterfaceMesh().NumberOfNodes();
        for( auto neighbour : neighbours_custom){
            for( unsigned int i = 0; i < number_of_interface_nodes; ++i){
                auto node_it = node_it_begin + i;
                int node_id = node_it->Id();
                auto disp = node_it->FastGetSolutionStepValue(DISPLACEMENT);
                if( norm_2(disp) != 0 ) //&&!node_it->GetValue(IS_GHOST_NODE)) //Only send and non-ghost nodes
                    send_nodes_container[neighbour].push_back(*node_it.base());
            }
        }

        std::vector<NodesContainerType> recv_nodes_container(size);
        rModelPart.GetCommunicator().TransferObjects(send_nodes_container, recv_nodes_container);
        auto& interface_nodes = rModelPart.GetCommunicator().InterfaceMesh().Nodes();
        // Update ghost nodes
        for( auto neighbour : neighbours_custom){
            for( auto it = recv_nodes_container[neighbour].begin(); it != recv_nodes_container[neighbour].end(); ++it){
                int node_id = it->Id();

                auto displacement = it->FastGetSolutionStepValue(DISPLACEMENT);
                auto ghost_node = interface_nodes.find(node_id);
                if( ghost_node != interface_nodes.end() )
                    ghost_node->FastGetSolutionStepValue(DISPLACEMENT) = displacement;
            }
        }
        rModelPart.GetCommunicator().GetDataCommunicator().Barrier();
    }

    void MPM_MPI_Utilities::SetMPICommunicator( ModelPart& SourceModelPart, ModelPart& DestModelPart){
        DestModelPart.SetCommunicator(SourceModelPart.pGetCommunicator());
    }

    void MPM_MPI_Utilities::ClearLocalElementsFromCommunicator( ModelPart& rModelPart){
        rModelPart.GetCommunicator().LocalMesh().Elements().clear();
        // Remove non-grid conditions
        const auto condition_it_begin = rModelPart.ConditionsBegin();
        const unsigned int number_of_conditions = rModelPart.NumberOfConditions();
        for( unsigned int i = 0; i < number_of_conditions; ++i){
            auto condition_it = condition_it_begin + i;
            if( condition_it->Is(BOUNDARY) ){
                rModelPart.GetCommunicator().LocalMesh().RemoveCondition(&*condition_it);
            }
        }
    }

} // end namespace Kratos