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

// System includes

// External includes

// Project includes
#include "key_plane_generation_by_bounding_box.h"

namespace Kratos {

void KeyPlaneGenerationByBoundingBox::ValidateParameters()
{
    KRATOS_TRY;
    Parameters parameters = this->GetParameters();
    parameters.ValidateAndAssignDefaults(this->GetDefaultParameters());

    KRATOS_ERROR_IF(parameters["voxel_sizes"].size() != 3) << "voxel_sizes should be defined for this generator as an array of 3 sizes for x,y,z directions" << std::endl;
    KRATOS_ERROR_IF(parameters["min_point"].size() != 3) << "min_point should be defined for this generator as an array of 3 coordinates of the min point" << std::endl;
    KRATOS_ERROR_IF(parameters["max_point"].size() != 3) << "max_point should be defined for this generator as an array of 3 coordinates of the max point" << std::endl;

    KRATOS_CATCH("");
}

void KeyPlaneGenerationByBoundingBox::Generate()
{
    Parameters parameters = GetParameters();
    Generate(parameters);
}

void KeyPlaneGenerationByBoundingBox::Generate(Parameters parameters)
{
    array_1d<double,3> voxel_sizes = parameters["voxel_sizes"].GetVector();
    array_1d<double,3> min_point = parameters["min_point"].GetVector();
    array_1d<double,3> max_point = parameters["max_point"].GetVector();

    for(std::size_t i_direction = 0 ; i_direction < 3 ; i_direction++){
        double input_voxel_size = voxel_sizes[i_direction];
        const double min_coordinate = min_point[i_direction];
        const double max_coordinate = max_point[i_direction];

        KRATOS_ERROR_IF(input_voxel_size == 0.00) << "voxel_sizes in direction " << i_direction << " cannot be 0.00";

        const double length = (max_coordinate - min_coordinate);

        KRATOS_ERROR_IF_NOT(length>0.0) << "Negative or zero length of voxelization bounding box in " << i_direction << " direction" << std::endl;

        std::size_t number_of_divisions = static_cast<std::size_t>(std::round(length / input_voxel_size));
        number_of_divisions= (number_of_divisions == 0) ? 1 : number_of_divisions;

        const double voxel_size = length / number_of_divisions;

        for( std::size_t i = 0 ; i < number_of_divisions + 1 ; i++){
            AddKeyPlane(i_direction, min_coordinate + i*voxel_size);
        }
    }
}

}
