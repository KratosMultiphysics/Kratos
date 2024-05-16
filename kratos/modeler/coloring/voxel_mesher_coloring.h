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
#include "geometries/geometry_data.h"
#include "includes/kratos_parameters.h"

#include "modeler/voxel_mesh_generator_modeler.h"
#include "modeler/internals/cartesian_mesh_colors.h"
#include "modeler/internals/skin_intersection.h"

namespace Kratos {

class VoxelMesherColoring {

    VoxelMeshGeneratorModeler& mrModeler;
    Parameters mParameters;

public:

    KRATOS_CLASS_POINTER_DEFINITION(VoxelMesherColoring);

    using CartesianMeshColors = Internals::CartesianMeshColors;
    using SkinIntersection = Internals::SkinIntersection;

    VoxelMesherColoring(VoxelMeshGeneratorModeler& rModeler, Parameters ColoringParameters);

    virtual ~VoxelMesherColoring() = default;

    virtual void ValidateParameters();

    virtual void Apply() const = 0;

protected:

    Parameters GetParameters() const;

    virtual Parameters GetDefaultParameters() const;

    void CheckGeometryType(const GeometryData::KratosGeometryType &rType) const;

    // internal interface to VoxelMeshGeneratorModeler

    ModelPart& GetModelPart(const std::string& rName) const;

    const ModelPart& GetInputModelPart() const;

    CartesianMeshColors& GetMeshColors() const;

    SkinIntersection& GetIntersections(int Color) const;

    const std::vector<double>& GetKeyPlanes(std::size_t Direction) const;
};

}
