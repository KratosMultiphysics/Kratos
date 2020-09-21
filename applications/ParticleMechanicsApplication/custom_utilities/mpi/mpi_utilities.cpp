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

    void MPM_MPI_Utilities::TransferElements( ModelPart& rMPMModelPart, std::vector<ElementsContainerType>& rSendElements)
    {
        const unsigned int rank = rMPMModelPart.GetCommunicator().MyPID();
        const unsigned int size = rMPMModelPart.GetCommunicator().TotalProcesses();

        std::vector<ElementsContainerType> recv_elements(size);

        rMPMModelPart.GetCommunicator().TransferObjects(rSendElements, recv_elements);

        // Remove sent elements from current ModelPart
        for( unsigned int i = 0; i < rSendElements.size(); ++i){
            for(ElementsContainerType::iterator it = rSendElements[i].begin(); it != rSendElements[i].end(); ++it){
                it->GetGeometry().clear();
                it->Reset(ACTIVE);
                rMPMModelPart.RemoveElementFromAllLevels(it->Id());
            }
        }

        // Add recieved elements to current ModelPart
        for( unsigned int i = 0; i < recv_elements.size(); ++i){
            if( rank != i){
                for(ElementsContainerType::iterator it = recv_elements[i].begin(); it != recv_elements[i].end(); ++it){
                    rMPMModelPart.AddElement(*it.base());
                }
            }
        }
    }

    void MPM_MPI_Utilities::TransferConditions(ModelPart& rMPMModelPart, std::vector<ConditionsContainerType>& rSendConditions)
    {
        const unsigned int rank = rMPMModelPart.GetCommunicator().MyPID();
        const unsigned int size = rMPMModelPart.GetCommunicator().TotalProcesses();

        std::vector<ConditionsContainerType> RecvConditions(size);

        rMPMModelPart.GetCommunicator().TransferObjects(rSendConditions, RecvConditions);

        // Remove sent elements from current ModelPart
        for( unsigned int i = 0; i < rSendConditions.size(); ++i){
            for(ConditionsContainerType::iterator it = rSendConditions[i].begin(); it != rSendConditions[i].end(); ++it){
                rMPMModelPart.RemoveConditionFromAllLevels(it->Id());
            }
        }

        // Add recieved elements to current ModelPart
        for( unsigned int i = 0; i < RecvConditions.size(); ++i){
            if( rank != i){
                for(ConditionsContainerType::iterator it = RecvConditions[i].begin(); it != RecvConditions[i].end(); ++it){
                    rMPMModelPart.AddCondition(*it.base());
                }
            }
        }

    }

} // end namespace Kratos