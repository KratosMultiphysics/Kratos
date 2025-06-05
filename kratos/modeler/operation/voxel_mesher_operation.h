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
#include "includes/kratos_parameters.h"

#include "modeler/voxel_mesh_generator_modeler.h"
#include "modeler/internals/cartesian_mesh_colors.h"

namespace Kratos {

class KRATOS_API(KRATOS_CORE) VoxelMesherOperation {

    VoxelMeshGeneratorModeler& mrModeler;
    Parameters mParameters;

public:

    using CartesianMeshColors = Internals::CartesianMeshColors;

    KRATOS_CLASS_POINTER_DEFINITION(VoxelMesherOperation);

    VoxelMesherOperation(VoxelMeshGeneratorModeler& rModeler, Parameters OperationParameters);

    virtual ~VoxelMesherOperation() = default;

    virtual void ValidateParameters();

    virtual void Execute() = 0;

protected:

    Parameters GetParameters() const;

    virtual Parameters GetDefaultParameters() const;

    // internal interface to VoxelMeshGeneratorModeler

    ModelPart& GetModelPart(const std::string& rName) const;

    ModelPart& CreateAndGetModelPart(const std::string& rFullName) const;

    CartesianMeshColors& GetMeshColors() const;

    const array_1d<std::size_t, 3> &GetNumberOfDivisions() const;

    std::vector<double>& GetKeyPlanes(std::size_t Direction) const {
        return mrModeler.mKeyPlanes[Direction];
    }

    std::size_t GetNodeIndex(
        std::size_t I, 
        std::size_t J, 
        std::size_t K) const;

    VoxelMeshGeneratorModeler::CartesianNodalData& GetNodalData(
        std::size_t I, 
        std::size_t J, 
        std::size_t K) const; 

    std::size_t CalculateCenterOfElementPosition(double Coordinate, int ThisDimension) const;

};

}
