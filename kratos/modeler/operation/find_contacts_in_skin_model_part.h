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
#include "includes/model_part.h"
#include "voxel_mesher_operation.h"

namespace Kratos {

class FindContactsInSkinModelPart: public VoxelMesherOperation {

    struct ContactContainer
    {
        std::unordered_map<int,std::vector<std::size_t>> mNodeMap;
        std::unordered_map<int,std::vector<std::size_t>> mConditionsMap;
        std::unordered_map<int,ModelPart*> mModelPartMap;
        std::vector<int> mColorVector;
    };

public:
    FindContactsInSkinModelPart(VoxelMeshGeneratorModeler& rModeler, Parameters OperationParameters):
        VoxelMesherOperation(rModeler, OperationParameters)
    {}

    ~FindContactsInSkinModelPart() override = default;

    Parameters GetDefaultParameters() const override;

    void ValidateParameters() override;

    void Apply() const override;

    void FindConditionContact(
        Condition& rCondition,
        const int CellColor,
        ContactContainer& rContactContainer) const;

    void GetClosestContactColorFromNeighbours(
        const std::vector<int>& rContactColors,
        const std::vector<int>& rNeighbourColors,
        const std::vector<array_1d<std::size_t, 3>>& rNeighbourIndexes,
        const array_1d<double,3>& rConditionCenter,
        const array_1d<double,3>& rConditionNormal,
        const CartesianMeshColors& rColors,
        int& rPairingStatus,
        double& rMinDistance,
        int& rTemporaryColor) const;

    void GetCellNeihgbourColors(
        const array_1d<std::size_t,3>& rCellIndexes,
        std::vector<int>& rNeighbourColors,
        const CartesianMeshColors& rColors) const;

    void GetCellNeihgbourColorsConnectedByFace(
        const array_1d<std::size_t,3>& rCellIndexes,
        const CartesianMeshColors& rColors,
        std::vector<int>& rNeighbourColors,
        std::vector<array_1d<std::size_t, 3>>& rNeighbourIndexes) const;


    bool IsElementInsideBounds(const CartesianMeshColors& rColors, int i, int j, int k) const;

};

}