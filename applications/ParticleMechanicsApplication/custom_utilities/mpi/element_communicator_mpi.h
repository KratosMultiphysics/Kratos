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

#ifndef KRATOS_ELEMENT_COMMUNICATOR_MPI_INCLUDE_H
#define KRATOS_ELEMENT_COMMUNICATOR_MPI_INCLUDE_H

// Project inlcudes
#include "includes/model_part.h"
#include "mpi.h"

namespace Kratos {

namespace ElementCommunicatorMPI{

    void MPI_Search(ModelPart& rMPMModelPart, ModelPart& rBackgroundGridModelPart,
                    std::vector<typename Condition::Pointer>& rMissingConditions,
                    const std::size_t MaxNumberOfResults, const double Tolerance);
    void SentCondition(ModelPart& rMPMModelPart, Condition::Pointer cond, bool empty,
                    MPI_Datatype& message_type, int destination, int tag);
    int RecieveCondition(ModelPart& rMPMModelPart, ModelPart& rBackgroundGridModelPart,
                    std::vector<typename Condition::Pointer>& rMissingCondition,
                    bool& empty, MPI_Datatype& message_type, int source, int tag);
    void CreateMPIDataType(MPI_Datatype& type, bool condition);
} // namespace ElementCommunicatorMPI
} // namespace Kratos

#endif // KRATOS_MPM_MPI_UTILITY_IMPORT_H