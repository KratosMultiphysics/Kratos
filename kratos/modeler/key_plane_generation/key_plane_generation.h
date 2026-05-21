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

namespace Kratos {

class KRATOS_API(KRATOS_CORE) VoxelMesherKeyPlaneGeneration {

    VoxelMeshGeneratorModeler& mrModeler;
    Parameters mParameters;

public:

    KRATOS_CLASS_POINTER_DEFINITION(VoxelMesherKeyPlaneGeneration);

    VoxelMesherKeyPlaneGeneration(VoxelMeshGeneratorModeler& rModeler, Parameters GenerationParameters);

    virtual ~VoxelMesherKeyPlaneGeneration() = default;

    virtual void ValidateParameters();

    virtual void Generate() = 0;

protected:

    Parameters GetParameters() const;

    virtual Parameters GetDefaultParameters() const;

    // internal interface to VoxelMeshGeneratorModeler

    void AddKeyPlane(std::size_t Direction, double Coordinate);

    const ModelPart& GetInputModelPart() const;

    const ModelPart& GetModelPart(const std::string& rMyName) const;

};

}
