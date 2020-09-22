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

#ifndef KRATOS_MPM_MPI_UTILITIES_INCLUDE_H
#define KRATOS_MPM_MPI_UTILITIES_INCLUDE_H

// system includes

// external includes

// kratos includes
#include "includes/model_part.h"

namespace Kratos {
namespace MPM_MPI_Utilities{

    // Typedefs
    typedef ModelPart::ElementsContainerType ElementsContainerType;
    typedef ModelPart::ConditionsContainerType ConditionsContainerType;

    /**
     * @brief Interface to exchange elements.
     * @details Interface to exchange elements betweem mpi-processes and remove/add elements from/to respective ModelPart.
     * @param rSendElements list of objects to be send.      SendObjects[i] -> Objects to   process i
     **/
    void TransferElements(ModelPart& rModelPart,
                          std::vector<ElementsContainerType>& rSendElements);

    /**
     * @brief Interface to exchange conditions.
     * @details Interface to exchange conditions betweem mpi-processes and remove/add conditions from/to respective ModelPart.
     * @param rSendCondition list of conditions to be send.      SendObjects[i] -> Objects to   process i
     **/
    void TransferConditions(ModelPart& rModelPart,
                            std::vector<ConditionsContainerType>& rSendCondition);

} // end MPM_MPI_Utilities
} // end namespace Kratos

#endif // end KRATOS_MPM_MPI_UTILITIES_INCLUDE_H
