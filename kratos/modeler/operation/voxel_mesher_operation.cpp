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
#include "voxel_mesher_operation.h"

#include "modeler/voxel_mesh_generator_modeler.h"

namespace Kratos {

VoxelMesherOperation::VoxelMesherOperation(VoxelMeshGeneratorModeler& rModeler, Parameters OperationParameters):
    mrModeler(rModeler),
    mParameters(OperationParameters)
{}


void VoxelMesherOperation::ValidateParameters()
{
    mParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
}


Parameters VoxelMesherOperation::GetParameters() const
{
    return mParameters;
}


Parameters VoxelMesherOperation::GetDefaultParameters() const
{
    return Parameters(R"({})");
}


ModelPart& VoxelMesherOperation::GetModelPart(const std::string& rName) const {
    return mrModeler.mpModel->GetModelPart(rName);
}


ModelPart& VoxelMesherOperation::CreateAndGetModelPart(std::string const& rFullName) const {
    return mrModeler.CreateAndGetModelPart(rFullName);
}


const VoxelMesherOperation::CartesianMeshColors& VoxelMesherOperation::GetMeshColors() const
{
    return mrModeler.mColors;
}


std::size_t VoxelMesherOperation::CalculateCenterOfElementPosition(double Coordinate, int ThisDimension) const
{
    return mrModeler.mColors.CalculateCenterOfElementPosition(Coordinate, ThisDimension);
}

}