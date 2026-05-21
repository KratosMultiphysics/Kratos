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

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class MapperFactory< SparseSpaceType, DenseSpaceType >;

}  // namespace Kratos.

