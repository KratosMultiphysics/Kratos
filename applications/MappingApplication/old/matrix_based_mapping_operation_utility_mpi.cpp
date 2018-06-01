//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "matrix_based_mapping_operation_utility_mpi.h"

namespace Kratos
{
    /***********************************************************************************/
    /* PUBLIC Methods */
    /***********************************************************************************/
    MatrixBasedMappingOperationUtilityMPI::MatrixBasedMappingOperationUtilityMPI(
        ModelPartPointerType pInterfaceModelPart)
        : MappingOperationUtility(pInterfaceModelPart)
    {
        KRATOS_ERROR_IF_NOT(SparseSpaceType::IsDistributed())
            << "Wrong SparseSpace used!" << std::endl;





    }

    /***********************************************************************************/
    /* PROTECTED Methods */
    /***********************************************************************************/


    /***********************************************************************************/
    /* PRIVATE Methods */
    /***********************************************************************************/


}  // namespace Kratos.
