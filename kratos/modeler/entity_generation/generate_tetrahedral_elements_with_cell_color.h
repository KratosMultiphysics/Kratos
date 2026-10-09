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

class GenerateTetrahedralElementsWithCellColor: public VoxelMesherEntityGeneration {
public:
    GenerateTetrahedralElementsWithCellColor(VoxelMeshGeneratorModeler& rModeler, Parameters GenerationParameters):
        VoxelMesherEntityGeneration(rModeler, GenerationParameters)
    {}

    ~GenerateTetrahedralElementsWithCellColor() override = default;

    Parameters GetDefaultParameters() const override;

    void Generate(ModelPart& rModelPart, Parameters parameters) override;

private:

    void CreateTetrahedraInCell(
        std::vector<ModelPart::ElementType::Pointer>& rElements,
        Element::NodesArrayType& rCellNodes,
        const std::size_t StartId,
        Properties::Pointer& pProperties,
        const Element& rPrototype
        );
};

}