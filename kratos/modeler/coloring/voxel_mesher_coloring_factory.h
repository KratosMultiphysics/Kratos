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
#include "voxel_mesher_coloring.h"

namespace Kratos {

using VoxelMesherColoringFactory = VoxelMesherComponentFactory<VoxelMesherColoring>;

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<VoxelMesherColoringFactory::RegisteredType>;

void RegisterVoxelMesherColoring();

}