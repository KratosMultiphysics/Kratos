//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#if !defined(KRATOS_MAPPER_BACKEND_H_INCLUDED)
#define KRATOS_MAPPER_BACKEND_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_searching/interface_communicator.h"
#include "custom_utilities/mapping_matrix_utilities.h"
#include "custom_utilities/interface_vector_container.h"

namespace Kratos {

template<class TSparseSpace, class TDenseSpace>
struct MapperBackend
{
    using InterfaceCommunicatorType = InterfaceCommunicator;
    using MappingMatrixUtilitiesType = MappingMatrixUtilities<TSparseSpace, TDenseSpace>;
    using InterfaceVectorContainerType = InterfaceVectorContainer<TSparseSpace, TDenseSpace>;
};

}  // namespace Kratos.

#endif // KRATOS_MAPPER_BACKEND_H_INCLUDED defined
