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

#ifndef KRATOS_MPM_MPI_UTILITY_IMPORT_H
#define KRATOS_MPM_MPI_UTILITY_IMPORT_H

// Project inlcudes
#include "includes/model_part.h"

namespace Kratos {

namespace MPM_MPI_Utility{
    void CopyMPICommunicator( ModelPart& SourceModelPart, ModelPart& DestModelPart);

} // namespace MPM_MPI_Utility
} // namespace Kratos

#endif // KRATOS_MPM_MPI_UTILITY_IMPORT_H