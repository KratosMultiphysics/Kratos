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
#include "containers/model.h"

#include "mapping_application.h"
#include "mapping_application_variables.h"
#include "mappers/mapper_define.h"
#include "custom_utilities/mapper_backend.h"

#include "geometries/point_3d.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/prism_3d_6.h"
#include "geometries/hexahedra_3d_8.h"

#include "factories/mapper_factory.h"

#include "custom_mappers/nearest_neighbor_mapper.h"
#include "custom_mappers/nearest_element_mapper.h"
#include "custom_mappers/barycentric_mapper.h"
#include "custom_mappers/projection_3D_2D_mapper.h"
#include "custom_mappers/coupling_geometry_mapper.h"
#include "custom_mappers/nearest_neighbor_mapper_iga.h"

// Macros for registering mappers
// will be removed once using the core factories
#define KRATOS_REGISTER_MAPPER(MapperType, MapperName)                                                \
    {                                                                                                 \
    Model current_model;                                                                              \
    ModelPart& dummy_model_part = current_model.CreateModelPart("dummy");                             \
    MapperFactory<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>::Register    \
        (MapperName, Kratos::make_shared<MapperType<                                                  \
        MapperDefinitions::SparseSpaceType,MapperDefinitions::DenseSpaceType>>                        \
        (dummy_model_part, dummy_model_part));                                                        \
    }

#define KRATOS_REGISTER_MAPPER_WITH_BACKEND(MapperType, MapperName)                                   \
    {                                                                                                 \
    Model current_model;                                                                              \
    ModelPart& dummy_model_part = current_model.CreateModelPart("dummy");                             \
    MapperFactory<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>::Register    \
        (MapperName, Kratos::make_shared<MapperType<                                                  \
        MapperDefinitions::SparseSpaceType,MapperDefinitions::DenseSpaceType,MapperBackend<           \
                MapperDefinitions::SparseSpaceType,MapperDefinitions::DenseSpaceType>>>               \
        (dummy_model_part, dummy_model_part));                                                        \
    }

Kratos::KratosApplication* CreateApplication()
{
    return new Kratos::KratosMappingApplication();
}

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

    KRATOS_REGISTER_MAPPER_WITH_BACKEND(NearestNeighborMapper, "nearest_neighbor");
    KRATOS_REGISTER_MAPPER_WITH_BACKEND(NearestElementMapper,  "nearest_element");
    KRATOS_REGISTER_MAPPER_WITH_BACKEND(BarycentricMapper,     "barycentric");
    KRATOS_REGISTER_MAPPER_WITH_BACKEND(Projection3D2DMapper,  "projection_3D_2D");

    KRATOS_REGISTER_MAPPER_WITH_BACKEND(NearestNeighborMapperIGA, "nearest_neighbor_iga");

    KRATOS_REGISTER_MAPPER(CouplingGeometryMapper,  "coupling_geometry");

    KRATOS_REGISTER_MODELER("MappingGeometriesModeler", mMappingGeometriesModeler);
    KRATOS_REGISTER_MODELER("IgaFEMMappingGeometriesModeler", mIgaFEMMappingGeometriesModeler);

    KRATOS_REGISTER_VARIABLE( INTERFACE_EQUATION_ID )
    KRATOS_REGISTER_VARIABLE( PAIRING_STATUS )
    KRATOS_REGISTER_VARIABLE( CURRENT_COORDINATES )
    KRATOS_REGISTER_VARIABLE( IS_PROJECTED_LOCAL_SYSTEM)
    KRATOS_REGISTER_VARIABLE( IS_DUAL_MORTAR)
}
}  // namespace Kratos.
