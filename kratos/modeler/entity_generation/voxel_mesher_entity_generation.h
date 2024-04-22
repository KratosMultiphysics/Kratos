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
#include "modeler/internals/skin_intersection.h"

namespace Kratos {

class VoxelMesherEntityGeneration {

    VoxelMeshGeneratorModeler& mrModeler;
    Parameters mParameters;

public:
    using CartesianMeshColors = Internals::CartesianMeshColors;
    using SkinIntersection = Internals::SkinIntersection;

    KRATOS_CLASS_POINTER_DEFINITION(VoxelMesherEntityGeneration);

    VoxelMesherEntityGeneration(VoxelMeshGeneratorModeler& rModeler, Parameters GenerationParameters);

    virtual ~VoxelMesherEntityGeneration() = default;

    virtual void ValidateParameters();

    virtual void Generate();

    virtual void Generate(ModelPart& rModelPart, Parameters parameters) = 0;

protected:

    Parameters GetParameters() const;

    virtual Parameters GetDefaultParameters() const;

    // internal interface to VoxelMeshGeneratorModeler

    ModelPart& GetModelPart(const std::string& rName) const;

    ModelPart& CreateAndGetModelPart(const std::string& rFullName) const;

    const CartesianMeshColors& GetMeshColors() const;

    Properties::Pointer CreateAndGetProperty(ModelPart& rModelPart, std::size_t PropertyId) const;

    const array_1d<std::size_t, 3> &GetNumberOfDivisions() const;

    const SkinIntersection& GetIntersections(int Color) const;

    SkinIntersection& GetIntersections(int Color);

    bool IntersectionsGenerated(int Color) const;

    void AddNodesToModelPart(ModelPart& rModelPart, ModelPart::NodesContainerType& rNewNodes) const;

    void SetStartIds(ModelPart& rModelPart);

    std::size_t GetStartElementId() const;

    std::size_t GetStartConditionId() const;

    array_1d<std::size_t, 3> GetNumberOfCells() const;

    Node::Pointer GenerateOrRetrieveNode(
        ModelPart& rTheVolumeModelPart,
        ModelPart::NodesContainerType& rThisNodes,
        const std::size_t I,
        const std::size_t J,
        const std::size_t K);

    Node::Pointer GenerateOrRetrieveQuadraticNode(
        ModelPart& rTheVolumeModelPart,
        ModelPart::NodesContainerType& rThisNodes,
        const std::size_t I,
        const std::size_t J,
        const std::size_t K);

    Node::Pointer GenerateNode(ModelPart& rModelPart, const Point& rCoordinates);

    void GetLinearCellNodes(
        Element::NodesArrayType& rCellNodes,
        ModelPart& rTheVolumeModelPart,
        ModelPart::NodesContainerType& rThisNodes,
        const std::size_t I,
        const std::size_t J,
        const std::size_t K);

    void GetQuadraticCellNodes(
        Element::NodesArrayType& rCellNodes,
        ModelPart& rTheVolumeModelPart,
        ModelPart::NodesContainerType& rThisNodes,
        const std::size_t I,
        const std::size_t J,
        const std::size_t K);

    void InitializeQuadraticData();
};

}
