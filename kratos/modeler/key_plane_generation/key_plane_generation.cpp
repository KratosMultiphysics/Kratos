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
#include "key_plane_generation.h"

namespace Kratos {

VoxelMesherKeyPlaneGeneration::VoxelMesherKeyPlaneGeneration(VoxelMeshGeneratorModeler& rModeler, Parameters GenerationParameters):
    mrModeler(rModeler),
    mParameters(GenerationParameters)
{}


void VoxelMesherKeyPlaneGeneration::ValidateParameters()
{
    mParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
}


Parameters VoxelMesherKeyPlaneGeneration::GetParameters() const
{
    return mParameters;
}


Parameters VoxelMesherKeyPlaneGeneration::GetDefaultParameters() const
{
    return Parameters(R"({
        "voxel_sizes": [],
        "min_point": [],
        "max_point": []
    })");
}


void VoxelMesherKeyPlaneGeneration::AddKeyPlane(std::size_t Direction, double Coordinate)
{
    mrModeler.mKeyPlanes[Direction].push_back(Coordinate);
}


const ModelPart& VoxelMesherKeyPlaneGeneration::GetInputModelPart() const {
    return *(mrModeler.mpInputModelPart);
}

const ModelPart& VoxelMesherKeyPlaneGeneration::GetModelPart(const std::string& rMyModelPart) const {
    return mrModeler.CreateAndGetModelPart(rMyModelPart);
}

}