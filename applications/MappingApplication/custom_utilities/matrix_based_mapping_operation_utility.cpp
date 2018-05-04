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
#include "matrix_based_mapping_operation_utility.h"
#include "custom_utilities/mapper_typedefs.h"

namespace Kratos
{
    /***********************************************************************************/
    /* PUBLIC Methods */
    /***********************************************************************************/
    template<>
    MatrixBasedMappingOperationUtility<MapperDefinitions::SparseSpaceType,
        MapperDefinitions::DenseSpaceType>::MatrixBasedMappingOperationUtility()
        : MappingOperationUtility<MapperDefinitions::SparseSpaceType,
          MapperDefinitions::DenseSpaceType>()
    {
        KRATOS_WATCH("Non-MPI-Consructor")
    }


    /***********************************************************************************/
    /* PROTECTED Methods */
    /***********************************************************************************/


    /***********************************************************************************/
    /* PRIVATE Methods */
    /***********************************************************************************/

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class MatrixBasedMappingOperationUtility< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;


}  // namespace Kratos.
