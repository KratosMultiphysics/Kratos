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
{

template< class TVectorType, class TVarType >
void TUpdateSystemVectorFromModelPart(TVectorType& rVector,
                        ModelPart& rModelPart,
                        const TVarType& rVariable,
                        const Kratos::Flags& rMappingOptions)
{
    // Here we construct a function pointer to not have the if all the time inside the loop
    const auto fill_fct = MapperUtilities::GetFillFunction<TVarType>(rMappingOptions);

    const int num_local_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = rModelPart.GetCommunicator().LocalMesh().NodesBegin();

    #pragma omp parallel for
    for (int i=0; i<num_local_nodes; i++) {
        fill_fct(*(nodes_begin + i), rVariable, rVector[i]);
    }
}

template< class TVectorType, class TVarType >
void TUpdateModelPartFromSystemVector(TVectorType& rVector,
            ModelPart& rModelPart,
            const TVarType& rVariable,
            const Kratos::Flags& rMappingOptions)
{
    const double factor = rMappingOptions.Is(MapperFlags::SWAP_SIGN) ? -1.0 : 1.0;

    // Here we construct a function pointer to not have the if all the time inside the loop
    const auto update_fct = std::bind(MapperUtilities::GetUpdateFunction<TVarType>(rMappingOptions),
                                        std::placeholders::_1,
                                        std::placeholders::_2,
                                        std::placeholders::_3,
                                        factor);
    const int num_local_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
    const auto nodes_begin = rModelPart.GetCommunicator().LocalMesh().NodesBegin();

    #pragma omp parallel for
    for (int i=0; i<num_local_nodes; i++) {
        update_fct(*(nodes_begin + i), rVariable, rVector[i]);
    }
}

}

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
    TUpdateSystemVectorFromModelPart(mpInterfaceVector, mrModelPart, rVariable, rMappingOptions);
}

template<>
void UtilityType::UpdateSystemVectorFromModelPart(const ComponentVariableType& rVariable,
                                                  const Kratos::Flags& rMappingOptions)
{
    TUpdateSystemVectorFromModelPart(mpInterfaceVector, mrModelPart, rVariable, rMappingOptions);
}

template<>
void UtilityType::UpdateModelPartFromSystemVector(const DoubleVariableType& rVariable,
                                                  const Kratos::Flags& rMappingOptions)
{
    TUpdateModelPartFromSystemVector(mpInterfaceVector, mrModelPart, rVariable, rMappingOptions);
}

template<>
void UtilityType::UpdateModelPartFromSystemVector(const ComponentVariableType& rVariable,
                                                  const Kratos::Flags& rMappingOptions)
{
    TUpdateModelPartFromSystemVector(mpInterfaceVector, mrModelPart, rVariable, rMappingOptions);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class InterfaceVectorContainer< SparseSpaceType, DenseSpaceType >;

}  // namespace Kratos.
