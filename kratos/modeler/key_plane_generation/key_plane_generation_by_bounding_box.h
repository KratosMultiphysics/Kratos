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
#include "key_plane_generation.h"

namespace Kratos {

/// Creates the key planes positions in x,y,z directions
/** It takes the min_point and max_point of the bounding box and
 *  divides each direction into uniform divisions with length as
 *  close as possible to the given size.
 */
class KeyPlaneGenerationByBoundingBox: public VoxelMesherKeyPlaneGeneration {
public:
    KeyPlaneGenerationByBoundingBox(VoxelMeshGeneratorModeler& rModeler, Parameters GenerationParameters):
        VoxelMesherKeyPlaneGeneration(rModeler, GenerationParameters)
    {}

    ~KeyPlaneGenerationByBoundingBox() override = default;

    void ValidateParameters() override;

    void Generate() override;

protected:

    virtual void Generate(Parameters parameters);

};

}
