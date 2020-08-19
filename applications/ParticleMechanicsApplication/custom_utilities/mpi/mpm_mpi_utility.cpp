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
#include "mpm_mpi_utility.h"

namespace Kratos {

    void MPM_MPI_Utility::CopyMPICommunicator( ModelPart& SourceModelPart, ModelPart& DestModelPart){
        DestModelPart.SetCommunicator(SourceModelPart.pGetCommunicator());
    }

} // namespace Kratos