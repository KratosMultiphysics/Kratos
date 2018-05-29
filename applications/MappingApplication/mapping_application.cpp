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
#include "mapping_application.h"
#include "mapping_application_variables.h"

#include "custom_utilities/mapper_typedefs.h"

#include "geometries/point_3d.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/prism_3d_6.h"
#include "geometries/hexahedra_3d_8.h"

#include "custom_utilities/mapper_factory.h"

#include "custom_mappers/nearest_neighbor_mapper.h"
#include "custom_mappers/nearest_element_mapper.h"

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "mpi.h"
#define KRATOS_REGISTER_MAPPER(MapperType, MapperName)                                            \
{                                                                                                 \
ModelPart dummy_model_part("dummy");                                                              \
/* Registering the non-MPI mapper */                                                              \
MapperFactory::Register<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>    \
    (MapperName, Kratos::make_shared<MapperType<                                                  \
    MapperDefinitions::SparseSpaceType,MapperDefinitions::DenseSpaceType>>                        \
    (dummy_model_part, dummy_model_part));                                                        \
/* Registering the MPI mapper */                                                                  \
MapperFactory::Register<MapperDefinitions::MPISparseSpaceType, MapperDefinitions::DenseSpaceType> \
    (MapperName, Kratos::make_shared<MapperType<                                                  \
    MapperDefinitions::MPISparseSpaceType,MapperDefinitions::DenseSpaceType>>                     \
    (dummy_model_part, dummy_model_part));                                                        \
}
#else
#define KRATOS_REGISTER_MAPPER(MapperType, MapperName)                                            \
{                                                                                                 \
ModelPart dummy_model_part("dummy");                                                              \
MapperFactory::Register<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>    \
    (MapperName, Kratos::make_shared<MapperType<                                                  \
    MapperDefinitions::SparseSpaceType,MapperDefinitions::DenseSpaceType>>                        \
    (dummy_model_part, dummy_model_part));                                                        \
}
#endif

namespace Kratos
{

KratosMappingApplication::KratosMappingApplication() :
    KratosApplication("MappingApplication"),
    mInterfaceObject(array_1d<double, 3>(0.0)),
    mInterfaceNode(),
    mInterfaceGeometryObject()
{}

void KratosMappingApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();

    std::stringstream banner;
    banner << "    KRATOS ______  ___                      _____  "                          << std::endl;
    banner << "           ___   |/  /_____ ___________________(_)_____________ _  "          << std::endl;
    banner << "           __  /|_/ /_  __ `/__  __ \\__  __ \\_  /__  __ \\_  __ `/  "       << std::endl;
    banner << "           _  /  / / / /_/ /__  /_/ /_  /_/ /  / _  / / /  /_/ /  "           << std::endl;
    banner << "           /_/  /_/  \\__,_/ _  .___/_  .___//_/  /_/ /_/_\\__, /  "          << std::endl;
    banner << "                            /_/     /_/                 /____/ Application"   << std::endl;

    banner << "Initializing KratosMappingApplication... " << std::endl;

    int rank = 0;

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized) MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (rank == 0) std::cout << banner.str();

    // registering the mappers using the registration-macro
    // note that this takes care to register the mappers also in MPI
    KRATOS_REGISTER_MAPPER(NearestNeighborMapper, "nearest_neighbor");
    KRATOS_REGISTER_MAPPER(NearestElementMapper,  "nearest_element");

    KRATOS_REGISTER_VARIABLE( INTERFACE_EQUATION_ID )

    // Needed to exchange Information abt the found neighbors (i.e. only for debugging)
    KRATOS_REGISTER_VARIABLE( NEIGHBOR_RANK )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NEIGHBOR_COORDINATES )
}
}  // namespace Kratos.