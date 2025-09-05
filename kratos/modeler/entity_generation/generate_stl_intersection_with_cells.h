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
#include "voxel_mesher_entity_generation.h"

namespace Kratos {

class GenerateStlIntersectionWithCells: public VoxelMesherEntityGeneration {
public:
    using WorkGeometryType = Internals::SkinIntersection::SplitGeometryType;

    GenerateStlIntersectionWithCells(VoxelMeshGeneratorModeler& rModeler, Parameters GenerationParameters):
        VoxelMesherEntityGeneration(rModeler, GenerationParameters)
    {}

    ~GenerateStlIntersectionWithCells() override = default;

    Parameters GetDefaultParameters() const override;

    void Generate(ModelPart& rModelPart, Parameters parameters) override;

private:

    PointerVector<Node> GenerateCutGeometryNodes(
        ModelPart& rTheCutModelPart,
        const std::vector<typename WorkGeometryType::Pointer>& rCutGeometries
        );

    void CheckGeometryType(const GeometryData::KratosGeometryType &rType) const;

    void GenerateIntersection(const ModelPart& rReferenceModelPart, const std::string& rInputEntities, int Color);

    void GenerateIntersection(
        const Geometry<Node>& rGeometry,
        int Color,
        const CartesianMeshColors& rMeshColors,
        SkinIntersection& rSkinIntersection
        );
};

}