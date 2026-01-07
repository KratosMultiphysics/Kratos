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
#include "interface_vector_container.h"
#include "mappers/mapper_define.h"
#include "custom_utilities/mapper_utilities.h"

namespace Kratos
{
typedef typename MapperDefinitions::SparseSpaceType SparseSpaceType;
typedef typename MapperDefinitions::DenseSpaceType  DenseSpaceType;

typedef InterfaceVectorContainer<SparseSpaceType, DenseSpaceType> VectorContainerType;

/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
template<>
void VectorContainerType::UpdateSystemVectorFromModelPart(const Variable<double>& rVariable,
                                                          const Kratos::Flags& rMappingOptions)
{
    switch (mInterfaceEntityType) {
        case InterfaceEntityType::NODES:
            MapperUtilities::UpdateSystemVectorFromModelPartNodes(*mpInterfaceVector, mrModelPart, rVariable, rMappingOptions);
            break;

        case InterfaceEntityType::ELEMENTS:
            MapperUtilities::UpdateSystemVectorFromModelPartElements(*mpInterfaceVector, mrModelPart, rVariable, rMappingOptions);
            break;

        case InterfaceEntityType::CONDITIONS:
            MapperUtilities::UpdateSystemVectorFromModelPartConditions(*mpInterfaceVector, mrModelPart, rVariable, rMappingOptions);
            break;

        case InterfaceEntityType::GEOMETRIES:
            MapperUtilities::UpdateSystemVectorFromModelPartGeometries(*mpInterfaceVector, mrModelPart, rVariable, rMappingOptions);
            break;
    }
}

template<>
void VectorContainerType::UpdateModelPartFromSystemVector(const Variable<double>& rVariable,
                                                          const Kratos::Flags& rMappingOptions)
{
    switch (mInterfaceEntityType) {
        case InterfaceEntityType::NODES:
            MapperUtilities::UpdateModelPartNodesFromSystemVector(*mpInterfaceVector, mrModelPart, rVariable, rMappingOptions);
            break;

        case InterfaceEntityType::ELEMENTS:
            MapperUtilities::UpdateModelPartElementsFromSystemVector(*mpInterfaceVector, mrModelPart, rVariable, rMappingOptions);
            break;

        case InterfaceEntityType::CONDITIONS:
            MapperUtilities::UpdateModelPartConditionsFromSystemVector(*mpInterfaceVector, mrModelPart, rVariable, rMappingOptions);
            break;

        case InterfaceEntityType::GEOMETRIES:
            MapperUtilities::UpdateModelPartGeometriesFromSystemVector(*mpInterfaceVector, mrModelPart, rVariable, rMappingOptions);
            break;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class InterfaceVectorContainer< SparseSpaceType, DenseSpaceType >;

}  // namespace Kratos.
