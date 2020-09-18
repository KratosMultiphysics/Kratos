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

    void MPM_MPI_Utilities::TransferElements( ModelPart& rMPMModelPart, std::vector<ElementsContainerType>& SendElements)
    {
        const unsigned int rank = rMPMModelPart.GetCommunicator().MyPID();
        const unsigned int size = rMPMModelPart.GetCommunicator().TotalProcesses();

        std::vector<ElementsContainerType> RecvElements(size);

        rMPMModelPart.GetCommunicator().TransferObjects(SendElements, RecvElements);

        // Remove sent elements from current ModelPart
        for( unsigned int i = 0; i < SendElements.size(); ++i){
            for(ElementsContainerType::iterator it = SendElements[i].begin(); it != SendElements[i].end(); ++it){
                it->GetGeometry().clear();
                it->Reset(ACTIVE);
                rMPMModelPart.RemoveElementFromAllLevels(it->Id());
            }
        }

        // Add recieved elements to current ModelPart
        for( unsigned int i = 0; i < RecvElements.size(); ++i){
            if( rank != i){
                for(ElementsContainerType::iterator it = RecvElements[i].begin(); it != RecvElements[i].end(); ++it){
                    rMPMModelPart.AddElement(*it.base());
                }
            }
        }
    }

} // end namespace Kratos