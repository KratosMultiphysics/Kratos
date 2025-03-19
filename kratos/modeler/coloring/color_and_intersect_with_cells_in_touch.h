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

class ColorAndIntersectWithCellsInTouch: public VoxelMesherColoring {

public:
    ColorAndIntersectWithCellsInTouch(VoxelMeshGeneratorModeler& rModeler, Parameters ColoringParameters);

    ~ColorAndIntersectWithCellsInTouch() override = default;

    void Apply() const override;

    Parameters GetDefaultParameters() const override;

private:

    void ApplyColorToCellsInTouchWithGeometryAndComputeIntersection(
        const Geometry<Node>& rGeometry,
        int InsideColor,
        Internals::SkinIntersection& rSkinIntersection) const;

};

}