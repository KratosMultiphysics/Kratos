//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
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
    constexpr bool in_parallel = false; // accessing the Trilinos vectors is not threadsafe in the default configuration!

    auto& r_vector = (*mpInterfaceVector)[0];

    switch (mInterfaceEntityType) {
        case InterfaceEntityType::NODES:
            MapperUtilities::UpdateSystemVectorFromModelPartNodes(
                r_vector, mrModelPart, rVariable, rMappingOptions, in_parallel);
            break;

        case InterfaceEntityType::ELEMENTS:
            MapperUtilities::UpdateSystemVectorFromModelPartElements(
                r_vector, mrModelPart, rVariable, rMappingOptions, in_parallel);
            break;

        case InterfaceEntityType::CONDITIONS:
            MapperUtilities::UpdateSystemVectorFromModelPartConditions(
                r_vector, mrModelPart, rVariable, rMappingOptions, in_parallel);
            break;

        case InterfaceEntityType::GEOMETRIES:
            MapperUtilities::UpdateSystemVectorFromModelPartGeometries(
                r_vector, mrModelPart, rVariable, rMappingOptions, in_parallel);
            break;
    }
}

template<>
void VectorContainerType::UpdateModelPartFromSystemVector(const Variable<double>& rVariable,
                                                          const Kratos::Flags& rMappingOptions)
{
    constexpr bool in_parallel = false; // accessing the Trilinos vectors is not threadsafe in the default configuration!

    const auto& r_vector = (*mpInterfaceVector)[0];

    switch (mInterfaceEntityType) {
        case InterfaceEntityType::NODES:
            MapperUtilities::UpdateModelPartNodesFromSystemVector(
                r_vector, mrModelPart, rVariable, rMappingOptions, in_parallel);
            break;

        case InterfaceEntityType::ELEMENTS:
            MapperUtilities::UpdateModelPartElementsFromSystemVector(
                r_vector, mrModelPart, rVariable, rMappingOptions, in_parallel);
            break;

        case InterfaceEntityType::CONDITIONS:
            MapperUtilities::UpdateModelPartConditionsFromSystemVector(
                r_vector, mrModelPart, rVariable, rMappingOptions, in_parallel);
            break;

        case InterfaceEntityType::GEOMETRIES:
            MapperUtilities::UpdateModelPartGeometriesFromSystemVector(
                r_vector, mrModelPart, rVariable, rMappingOptions, in_parallel);
            break;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class InterfaceVectorContainer<SparseSpaceType, DenseSpaceType>;

}  // namespace Kratos
