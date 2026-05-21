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

class ColorConnectedCellsInTouch: public VoxelMesherColoring {

public:
    ColorConnectedCellsInTouch(VoxelMeshGeneratorModeler& rModeler, Parameters ColoringParameters);

    ~ColorConnectedCellsInTouch() override = default;

    void Apply() const override;

private:

    void ApplyColorToConnectedCellsInTouchWithGeometry(
        const Geometry<Node>& rGeometry,
        const int InsideColor,
        const int CellColor
        ) const;

    void ColorConnectedCellsToThisCell(
        const std::size_t I,
        const std::size_t J,
        const std::size_t K,
        const int InsideColor,
        const int CellColor,
        CartesianMeshColors& rColors) const;

};

}
