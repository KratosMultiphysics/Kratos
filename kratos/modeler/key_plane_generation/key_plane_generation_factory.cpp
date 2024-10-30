//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

// External includes

// Project includes
#include "key_plane_generation_factory.h"
#include "key_plane_generation_by_bounding_box.h"
#include "key_plane_generation_by_outer_shell.h"
#include "key_plane_generation_with_refinement.h"

namespace Kratos {

VoxelMesherKeyPlaneGeneration::Pointer KeyPlaneGenerationFactory::Create(
    VoxelMeshGeneratorModeler& rModeler,
    Parameters GenerationSettings) const
{
    KRATOS_TRY;

    if (GenerationSettings.Has("type") && !GenerationSettings["type"].IsString()) {
        KRATOS_ERROR << "Voxel mesh generation modeler key plane factory error: \"type\" must be a string" << std::endl;
    }
    KRATOS_ERROR_IF_NOT(GenerationSettings.Has("Parameters"))
        << "Voxel mesh generation modeler key plane factory error: \"Parameters\" argument is required" << std::endl;

    const std::string type = GenerationSettings.Has("type") ? GenerationSettings["type"].GetString() : "bounding_box";
    KRATOS_ERROR_IF_NOT(Has(type))
        << "Trying to construct key plane generation strategy with unregistered type \"" << type << "\"" << std::endl
        << "The list of available options are: " << std::endl
        << KratosComponents<RegisteredType>() << std::endl;

    return KratosComponents<RegisteredType>::Get(type).Create(rModeler, GenerationSettings["Parameters"]);

    KRATOS_CATCH("");
}

template class KratosComponents<KeyPlaneGenerationFactory::RegisteredType>;


void RegisterVoxelMesherKeyPlaneGeneration()
{
    KeyPlaneGenerationFactory::Register<KeyPlaneGenerationByBoundingBox>("bounding_box");
    KeyPlaneGenerationFactory::Register<KeyPlaneGenerationByOuterShell>("outer_shell");
    KeyPlaneGenerationFactory::Register<KeyPlaneGenerationWithRefinement>("outer_shell_with_refinement");
}

}