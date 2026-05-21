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
#include "key_plane_generation_by_bounding_box.h"

namespace Kratos {

class KeyPlaneGenerationByOuterShell: public KeyPlaneGenerationByBoundingBox {
public:
    KeyPlaneGenerationByOuterShell(VoxelMeshGeneratorModeler& rModeler, Parameters GenerationParameters):
        KeyPlaneGenerationByBoundingBox(rModeler, GenerationParameters)
    {}

    ~KeyPlaneGenerationByOuterShell() override = default;

    void ValidateParameters() override
    {
        KRATOS_TRY;
        Parameters parameters = this->GetParameters();
        parameters.ValidateAndAssignDefaults(this->GetDefaultParameters());

        KRATOS_ERROR_IF(parameters["voxel_sizes"].size() != 3) << "voxel_sizes should be defined for this generator as an array of 3 sizes for x,y,z directions" << std::endl;

        KRATOS_CATCH("");
    }

    void Generate() override
    {
        Parameters parameters = this->GetParameters();
        Parameters bounding_box_parameters = parameters.Clone();

        auto dx = parameters["voxel_sizes"].GetVector();
        // Now we find the min and max value for the node coordinates
        Vector min_v(3,std::numeric_limits<double>::max()), max_v(3, -std::numeric_limits<double>::max());
        for(auto& node: GetInputModelPart().Nodes()){
            double x = node.X(), y = node.Y(), z = node.Z();
            min_v[0] = std::min(min_v[0], x); max_v[0] = std::max(max_v[0], x);
            min_v[1] = std::min(min_v[1], y); max_v[1] = std::max(max_v[1], y);
            min_v[2] = std::min(min_v[2], z); max_v[2] = std::max(max_v[2], z);
        }
        // Now we fill the PP and complete it
        double epsilon = parameters["margin"].GetDouble();
        bounding_box_parameters.AddVector("min_point", min_v - ((1.0 + epsilon)*dx));
        bounding_box_parameters.AddVector("max_point", max_v + ((1.0 + epsilon)*dx));

        KeyPlaneGenerationByBoundingBox::Generate(bounding_box_parameters);
    }

protected:

    Parameters GetDefaultParameters() const override
    {
        return Parameters(R"({
            "voxel_sizes": [],
            "margin": 0.01
        })");
    }

};

}
