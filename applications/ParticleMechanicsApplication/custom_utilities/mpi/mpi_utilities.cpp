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

namespace Kratos {

    void MPM_MPI_Utilities::TransferElements( ModelPart& rModelPart, std::vector<ElementsContainerType>& rSendElements)
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

        // Recieving elements container
        std::vector<ElementsContainerType> recv_elements(size);

        // Exchange elements
        rModelPart.GetCommunicator().TransferObjects(rSendElements, recv_elements);

        // Remove sent elements from current ModelPart
        for( unsigned int i = 0; i < rSendElements.size(); ++i){
            for(auto it = rSendElements[i].begin(); it != rSendElements[i].end(); ++it){
                it->GetGeometry().clear();
                it->Reset(ACTIVE);
                rModelPart.RemoveElementFromAllLevels(it->Id());
            }
        }

        // Add recieved elements to current ModelPart
        for( unsigned int i = 0; i < recv_elements.size(); ++i){
            for(auto it = recv_elements[i].begin(); it != recv_elements[i].end(); ++it){
                rModelPart.AddElement(*it.base());
            }
        }
        rSendElements.clear();
    }

    void MPM_MPI_Utilities::TransferConditions(ModelPart& rModelPart, std::vector<ConditionsContainerType>& rSendConditions)
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

        // Recieving conditions container
        std::vector<ConditionsContainerType> recv_conditions(size);

        // Exchange conditions
        rModelPart.GetCommunicator().TransferObjects(rSendConditions, recv_conditions);

        // Remove sent conditions from current ModelPart
        for( unsigned int i = 0; i < rSendConditions.size(); ++i){
            for(auto it = rSendConditions[i].begin(); it != rSendConditions[i].end(); ++it){
                it->GetGeometry().clear();
                it->Reset(ACTIVE);
                rModelPart.RemoveConditionFromAllLevels(it->Id());
            }
        }

        // Add recieved conditions to current ModelPart
        for( unsigned int i = 0; i < recv_conditions.size(); ++i){
            for(auto it = recv_conditions[i].begin(); it != recv_conditions[i].end(); ++it){
                rModelPart.AddCondition(*it.base());
            }
        }
        rSendConditions.clear();
    }

} // end namespace Kratos