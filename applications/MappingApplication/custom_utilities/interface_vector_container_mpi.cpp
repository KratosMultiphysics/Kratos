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
#include "interface_vector_container.h"
#include "custom_utilities/mapper_typedefs.h"
#include "custom_utilities/mapper_utilities.h"

namespace Kratos
{
typedef typename MapperDefinitions::MPISparseSpaceType MPISparseSpaceType;
typedef typename MapperDefinitions::DenseSpaceType DenseSpaceType;

typedef InterfaceVectorContainer<MPISparseSpaceType, DenseSpaceType> MPIVectorContainerType;

/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
template<>
void MPIVectorContainerType::UpdateSystemVectorFromModelPart(const DoubleVariableType& rVariable,
                                                          const Kratos::Flags& rMappingOptions)
{
    MapperUtilities::UpdateSystemVectorFromModelPart((*mpInterfaceVector)[0], mrModelPart, rVariable, rMappingOptions);
}

template<>
void MPIVectorContainerType::UpdateSystemVectorFromModelPart(const ComponentVariableType& rVariable,
                                                          const Kratos::Flags& rMappingOptions)
{
    MapperUtilities::UpdateSystemVectorFromModelPart((*mpInterfaceVector)[0], mrModelPart, rVariable, rMappingOptions);
}

template<>
void MPIVectorContainerType::UpdateModelPartFromSystemVector(const DoubleVariableType& rVariable,
                                                          const Kratos::Flags& rMappingOptions)
{
    MapperUtilities::UpdateModelPartFromSystemVector((*mpInterfaceVector)[0], mrModelPart, rVariable, rMappingOptions);
}

template<>
void MPIVectorContainerType::UpdateModelPartFromSystemVector(const ComponentVariableType& rVariable,
                                                          const Kratos::Flags& rMappingOptions)
{
    MapperUtilities::UpdateModelPartFromSystemVector((*mpInterfaceVector)[0], mrModelPart, rVariable, rMappingOptions);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class InterfaceVectorContainer< MPISparseSpaceType, DenseSpaceType >;

}  // namespace Kratos.
