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

#pragma once

// System includes

// External includes

// Project includes
#include "modeler/internals/voxel_mesher_component_factory.h"
#include "key_plane_generation.h"

namespace Kratos {

class KeyPlaneGenerationFactory: public VoxelMesherComponentFactory<VoxelMesherKeyPlaneGeneration> {
public:
    VoxelMesherKeyPlaneGeneration::Pointer Create(VoxelMeshGeneratorModeler& rModeler, Parameters ColoringSettings) const;

};

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<KeyPlaneGenerationFactory::RegisteredType>;

void RegisterVoxelMesherKeyPlaneGeneration();

}