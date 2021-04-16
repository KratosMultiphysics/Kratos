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
#include "interpolative_mapper_base.h"
#include "custom_utilities/mapper_typedefs.h"
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "custom_searching/mpi/interface_communicator_mpi.h"
#endif

namespace Kratos {

template<>
void InterpolativeMapperBase<MapperDefinitions::SparseSpaceType,
    MapperDefinitions::DenseSpaceType>::InitializeInterfaceCommunicator()
{
    mpIntefaceCommunicator = Kratos::make_unique<InterfaceCommunicator>(mrModelPartOrigin,
                                                                        mMapperLocalSystems,
                                                                        mMapperSettings);
}

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template<>
void InterpolativeMapperBase<MapperDefinitions::MPISparseSpaceType,
    MapperDefinitions::DenseSpaceType>::InitializeInterfaceCommunicator()
{
    mpIntefaceCommunicator = Kratos::make_unique<InterfaceCommunicatorMPI>(mrModelPartOrigin,
                                                                           mMapperLocalSystems,
                                                                           mMapperSettings);
}
#endif


///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class InterpolativeMapperBase< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template class InterpolativeMapperBase< MapperDefinitions::MPISparseSpaceType, MapperDefinitions::DenseSpaceType >;
#endif

}  // namespace Kratos.
