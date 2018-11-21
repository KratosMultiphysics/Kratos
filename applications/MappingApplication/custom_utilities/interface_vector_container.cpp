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
typedef typename MapperDefinitions::SparseSpaceType SparseSpaceType;
typedef typename MapperDefinitSettingsions::DenseSpaceType DenseSpaceType;

typedef InterfaceVectorContainer<SparseSpaceType, DenseSpaceType> UtilityType;

typedef std::size_t IndexType;
typedef std::size_t SizeType;

/***********************************************************************************/
/* Functions for internal use in this file */
/***********************************************************************************/


/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
template<>
UtilityType::InterfaceVectorContainer(ModelPart& rModelPart)
        : mrModelPart(rModelPart)
{

}

template<>
void UtilityType::UpdateSystemVectorFromModelPart(const DoubleVariableType& rVariable,
                                                  const Kratos::Flags& rMappingOptions)
{


}

template<>
void UtilityType::UpdateSystemVectorFromModelPart(const ComponentVariableType& rVariable,
                                                  const Kratos::Flags& rMappingOptions)
{


}

template<>
void UtilityType::UpdateSystemVectorFromModelPart(const Array3VariableType& rVariable,
                                                  const Kratos::Flags& rMappingOptions)
{
    KRATOS_ERROR << "MultiVectors are not yet implemented" << std::endl;
}

template<>
void UtilityType::UpdateModelPartFromSystemVector(const DoubleVariableType& rVariable,
                                                  const Kratos::Flags& rMappingOptions)
{


}

template<>
void UtilityType::UpdateModelPartFromSystemVector(const ComponentVariableType& rVariable,
                                                  const Kratos::Flags& rMappingOptions)
{


}

template<>
void UtilityType::UpdateModelPartFromSystemVector(const Array3VariableType& rVariable,
                                                  const Kratos::Flags& rMappingOptions)
{
    KRATOS_ERROR << "MultiVectors are not yet implemented" << std::endl;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class InterfaceVectorContainer< SparseSpaceType, DenseSpaceType >;


}  // namespace Kratos.
