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
#include "custom_utilities/interface_vector_container.h"
#include "mapper_mpi_define.h"
#include "custom_utilities/mapper_utilities.h"

namespace Kratos
{
typedef typename MPIMapperDefinitions::SparseSpaceType SparseSpaceType;
typedef typename MPIMapperDefinitions::DenseSpaceType  DenseSpaceType;

typedef InterfaceVectorContainer<SparseSpaceType, DenseSpaceType> VectorContainerType;

/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
template<>
void VectorContainerType::UpdateSystemVectorFromModelPart(const Variable<double>& rVariable,
                                                          const Kratos::Flags& rMappingOptions)
{
    constexpr bool in_parallel = false; // accessing the trilinos vectors is not threadsafe in the default configuration!
    MapperUtilities::UpdateSystemVectorFromModelPart((*mpInterfaceVector)[0], mrModelPart, rVariable, rMappingOptions, in_parallel);
}

template<>
void VectorContainerType::UpdateModelPartFromSystemVector(const Variable<double>& rVariable,
                                                          const Kratos::Flags& rMappingOptions)
{
    constexpr bool in_parallel = false; // accessing the trilinos vectors is not threadsafe in the default configuration!
    MapperUtilities::UpdateModelPartFromSystemVector((*mpInterfaceVector)[0], mrModelPart, rVariable, rMappingOptions, in_parallel);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class InterfaceVectorContainer< SparseSpaceType, DenseSpaceType >;

}  // namespace Kratos.
