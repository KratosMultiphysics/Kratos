//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

// Project includes
#include "mpm_mpi_search.h"

namespace Kratos {

void MPM_MPI_SEARCH::SearchElements(ModelPart& rMPMModelPart, ModelPart& rBackgroundGridModelPart,
                        std::vector<Element::Pointer>& rMissingElements,
                        std::vector<Condition::Pointer>& rMissingConditions,
                        const std::size_t MaxNumberOfResults, const double Tolerance)
{
    ElementsArrayType Elements;
    bool sent = false;
    for( auto& el : rMissingElements){
        Elements.push_back(el);
        sent = true;
        break;
    }
    std::vector<ElementsArrayType> SendElements(2);
    if( rMPMModelPart.GetCommunicator().MyPID() == 0){
    SendElements[1] = Elements;
    }
    else{
    SendElements[0] = Elements;
    }
    KRATOS_WATCH( SendElements[1].GetContainer() );
    KRATOS_WATCH( SendElements[0].GetContainer() );
    std::vector<ElementsArrayType> RecElements(2);
    // First test to send an Element
    rMPMModelPart.GetCommunicator().TransferObjects(SendElements,RecElements);
    KRATOS_WATCH( RecElements[0].GetContainer() );
    KRATOS_WATCH( RecElements[1].GetContainer() );
}

} // namespace Kratos