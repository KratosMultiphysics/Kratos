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

class ColorOuterFacesOfCellsWithColors: public VoxelMesherColoring {

public:
    ColorOuterFacesOfCellsWithColors(VoxelMeshGeneratorModeler& rModeler, Parameters ColoringParameters);

    ~ColorOuterFacesOfCellsWithColors() override = default;

    void Apply() const override;

private:

    void ApplyColorIfOuterFace(
        const int InterfaceColor,
        const int CellColor,
        const array_1d<std::size_t, 3>& rCellIndices) const;

};

}