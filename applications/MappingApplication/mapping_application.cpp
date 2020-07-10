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
#include "containers/model.h"

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

// Macro to register the mapper WITHOUT MPI
#define KRATOS_REGISTER_MAPPER_SERIAL(MapperType, MapperName)                                         \
    {                                                                                                 \
    Model current_model;                                                                              \
    ModelPart& dummy_model_part = current_model.CreateModelPart("dummy");                             \
    MapperFactory::Register<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>    \
        (MapperName, Kratos::make_shared<MapperType<                                                  \
        MapperDefinitions::SparseSpaceType,MapperDefinitions::DenseSpaceType>>                        \
        (dummy_model_part, dummy_model_part));                                                        \
    }

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
    // Macro to register the mapper WITH MPI
    #define KRATOS_REGISTER_MAPPER_MPI(MapperType, MapperName)                                            \
        {                                                                                                 \
        Model current_model;                                                                              \
        ModelPart& dummy_model_part = current_model.CreateModelPart("dummy");                             \
        MapperFactory::Register<MapperDefinitions::MPISparseSpaceType, MapperDefinitions::DenseSpaceType> \
            (MapperName, Kratos::make_shared<MapperType<                                                  \
            MapperDefinitions::MPISparseSpaceType,MapperDefinitions::DenseSpaceType>>                     \
            (dummy_model_part, dummy_model_part));                                                        \
        }
    // Macro to register the mapper WITH and WITHOUT MPI
    #define KRATOS_REGISTER_MAPPER(MapperType, MapperName)                                            \
        KRATOS_REGISTER_MAPPER_SERIAL(MapperType, MapperName)                                         \
        KRATOS_REGISTER_MAPPER_MPI(MapperType, MapperName)
    #else
    // Macro to register the mapper WITH and WITHOUT MPI
    #define KRATOS_REGISTER_MAPPER(MapperType, MapperName)                                            \
        KRATOS_REGISTER_MAPPER_SERIAL(MapperType, MapperName)
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
    KRATOS_INFO("") << "    KRATOS ______  ___                      _____\n"
                    << "           ___   |/  /_____ ___________________(_)_____________ _\n"
                    << "           __  /|_/ /_  __ `/__  __ \\__  __ \\_  /__  __ \\_  __ `/\n"
                    << "           _  /  / / / /_/ /__  /_/ /_  /_/ /  / _  / / /  /_/ /\n"
                    << "           /_/  /_/  \\__,_/ _  .___/_  .___//_/  /_/ /_/_\\__, /\n"
                    << "                            /_/     /_/                 /____/\n"
                    << "Initializing KratosMappingApplication..." << std::endl;

    // registering the mappers using the registration-macro
    KRATOS_REGISTER_MAPPER(NearestNeighborMapper, "nearest_neighbor");
    KRATOS_REGISTER_MAPPER(NearestElementMapper,  "nearest_element");

    KRATOS_REGISTER_VARIABLE( INTERFACE_EQUATION_ID )
    KRATOS_REGISTER_VARIABLE( PAIRING_STATUS )
    KRATOS_REGISTER_VARIABLE( CURRENT_COORDINATES )
}
}  // namespace Kratos.
