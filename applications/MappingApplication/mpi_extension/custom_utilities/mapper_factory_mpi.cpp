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
#include "custom_utilities/mapper_factory.h"
#include "mapper_mpi_define.h"


namespace Kratos {

template<> KRATOS_API(MAPPING_APPLICATION) std::unordered_map<std::string, typename Mapper<MPIMapperDefinitions::SparseSpaceType,
    MPIMapperDefinitions::DenseSpaceType>::Pointer>& MapperFactory::GetRegisteredMappersList<MPIMapperDefinitions::SparseSpaceType,
    MPIMapperDefinitions::DenseSpaceType>();
template<>
std::unordered_map<std::string, typename Mapper<MPIMapperDefinitions::SparseSpaceType,
    MPIMapperDefinitions::DenseSpaceType>::Pointer>& MapperFactory::GetRegisteredMappersList<MPIMapperDefinitions::SparseSpaceType,
    MPIMapperDefinitions::DenseSpaceType>()
{
    static std::unordered_map<std::string, typename Mapper<MPIMapperDefinitions::SparseSpaceType, MPIMapperDefinitions::DenseSpaceType>::Pointer> registered_mappers;

    return registered_mappers;
}

}  // namespace Kratos.

