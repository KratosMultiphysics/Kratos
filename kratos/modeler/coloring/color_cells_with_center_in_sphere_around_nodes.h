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
#include "voxel_mesher_coloring.h"

namespace Kratos {

class ColorCellsWithCenterInSphereArounNodes: public VoxelMesherColoring {

public:
    ColorCellsWithCenterInSphereArounNodes(VoxelMeshGeneratorModeler& rModeler, Parameters ColoringParameters);

    ~ColorCellsWithCenterInSphereArounNodes() override = default;

    void Apply() const override;

};

}