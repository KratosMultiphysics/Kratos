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
#include "factories/mapper_factory.h"
#include "mapper_mpi_define.h"


namespace Kratos {

typedef typename MPIMapperDefinitions::SparseSpaceType SparseSpaceType;
typedef typename MPIMapperDefinitions::DenseSpaceType  DenseSpaceType;

template<>
std::unordered_map<std::string, typename Mapper<SparseSpaceType, DenseSpaceType>::Pointer>& MapperFactory<SparseSpaceType,
    DenseSpaceType>::GetRegisteredMappersList()
{
    static std::unordered_map<std::string, typename Mapper<SparseSpaceType, DenseSpaceType>::Pointer> registered_mappers;

    return registered_mappers;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class MapperFactory< SparseSpaceType, DenseSpaceType >;

}  // namespace Kratos.

