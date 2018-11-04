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

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/prism_3d_6.h"
#include "geometries/hexahedra_3d_8.h"

#include "custom_utilities/mapper_factory.h"

#include "custom_mappers/nearest_neighbor_mapper.h"
#include "custom_mappers/nearest_element_mapper.h"

#ifdef KRATOS_USING_MPI
#include "mpi.h"
#endif

namespace Kratos
{

KratosMappingApplication::KratosMappingApplication() :
    KratosApplication("MappingApplication"),
    mInterfaceObject(0.0, 0.0, 0.0),
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
    if (mpi_initialized)   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (rank == 0) std::cout << banner.str();

    Model dummy_model;

    ModelPart& dummy_model_part = dummy_model.CreateModelPart("dummy");

    MapperFactory::Register("nearest_neighbor", Kratos::make_shared<NearestNeighborMapper>(dummy_model_part, dummy_model_part));
    MapperFactory::Register("nearest_element",  Kratos::make_shared<NearestElementMapper>(dummy_model_part, dummy_model_part));

    // Needed to exchange Information abt the found neighbors (i.e. only for debugging)
    KRATOS_REGISTER_VARIABLE( NEIGHBOR_RANK )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NEIGHBOR_COORDINATES )

}
}  // namespace Kratos.