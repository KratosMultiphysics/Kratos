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

    // Interface to exchange elements and remove/add elements to respective ModelPart
    void TransferElements(ModelPart& rMPMModelPart,
                          std::vector<ElementsContainerType>& rSendElements);

    // Interface to exchange conditions and remove/add conditions from/to respective ModelPart
    void TransferConditions(ModelPart& rMPMModelPart,
                            std::vector<ConditionsContainerType>& rSendCondition);

} // end MPM_MPI_Utilities
} // end namespace Kratos

#endif // end KRATOS_MPM_MPI_UTILITIES_INCLUDE_H
