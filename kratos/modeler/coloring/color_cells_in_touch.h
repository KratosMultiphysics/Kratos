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
#include "geometries/geometry.h"
#include "voxel_mesher_coloring.h"

namespace Kratos {

class ColorCellsInTouch: public VoxelMesherColoring {

public:
    ColorCellsInTouch(VoxelMeshGeneratorModeler& rModeler, Parameters ColoringParameters);

    ~ColorCellsInTouch() override = default;

    void Apply() const override;

private:

    void ApplyColorToCellsInTouchWithGeometry(
        const Geometry<Node>& rGeometry,
        int InsideColor,
        int OutsideColor) const;

};

}