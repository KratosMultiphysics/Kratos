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

std::size_t VoxelMesherOperation::GetNodeIndex(
    std::size_t I, 
    std::size_t J, 
    std::size_t K) const
{
    return mrModeler.mMeshingData.GetNodeIndex(I,J,K);
}

VoxelMeshGeneratorModeler::CartesianNodalData& VoxelMesherOperation::GetNodalData(
    std::size_t I, 
    std::size_t J, 
    std::size_t K) const
{
    return mrModeler.mMeshingData.GetNodalData(I, J, K);
}


VoxelMesherOperation::CartesianMeshColors& VoxelMesherOperation::GetMeshColors() const
{
    return mrModeler.mColors;
}

const array_1d<std::size_t, 3>& VoxelMesherOperation::GetNumberOfDivisions() const {
    return mrModeler.mMeshingData.GetNumberOfDivisions();
}


std::size_t VoxelMesherOperation::CalculateCenterOfElementPosition(double Coordinate, int ThisDimension) const
{
    return mrModeler.mColors.CalculateCenterOfElementPosition(Coordinate, ThisDimension);
}

}