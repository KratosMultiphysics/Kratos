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
#include "includes/data_communicator.h"




namespace Kratos {

///@addtogroup ParticleMechanicsApplication
///@{

///@name Kratos Classes
///@{



/**
 * @class MPM_MPI_Utilities
 * @ingroup ParticleMechanicsApplication
 * @brief Provides place to add mpi related utility functions.
 */
class KRATOS_API(PARTICLE_MECHANICS_APPLICATION) MPM_MPI_Utilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MPM_MPI_Utilities
    KRATOS_CLASS_POINTER_DEFINITION(MPM_MPI_Utilities);

    typedef ModelPart::ElementsContainerType ElementsContainerType;
    typedef ModelPart::ConditionsContainerType ConditionsContainerType;
    typedef ModelPart::NodesContainerType NodesContainerType;

    ///@}
    ///@name Static Operations
    ///@{

    /**
     * @brief Interface to exchange elements.
     * @details Interface to exchange elements betweem mpi-processes and removes elements from sending ModelPart.
     * @param rSendElements list of objects to be send.      SendObjects[i] -> Objects to   process i
     * @param rRecvElements list of objects to be recieved.  RecvObjects[i] -> Objects from process i
     **/
    static void TransferElements(ModelPart& rModelPart,
                          std::vector<ElementsContainerType>& rSendElements,
                          std::vector<ElementsContainerType>& rRecvElements);

    /**
     * @brief Interface to exchange conditions.
     * @details Interface to exchange conditions betweem mpi-processes and removes conditions from sending ModelPart.
     * @param rSendCondition list of conditions to be send.      SendObjects[i] -> Objects to   process i
     * @param rRecvCondition list of objects to be recieved.     RecvObjects[i] -> Objects from process i
     **/
    static void TransferConditions(ModelPart& rModelPart,
                            std::vector<ConditionsContainerType>& rSendCondition,
                            std::vector<ConditionsContainerType>& rRecvCondition);

    /**
     * @brief Synchronize nodal displacement at interface mesh;
     **/
    static void SynchronizeNodalDisplacementAtInterface(ModelPart& rModelPart);

    /**
     * @brief Set DestModelPart communicator to SourceModelPart communicator.
     **/
    static void SetMPICommunicator( ModelPart& SourceModelPart, ModelPart& DestModelPart);

    ///@}

}; // end class MPM_MPI_Utilities
///@} classes
///@} addToGroup
} // end namespace Kratos

#endif // end KRATOS_MPM_MPI_UTILITIES_INCLUDE_H
