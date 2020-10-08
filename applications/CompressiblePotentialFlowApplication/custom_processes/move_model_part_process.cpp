//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Núñez,
//


#include "move_model_part_process.h"
#include "utilities/geometrical_transformation_utilities.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
// Constructor for MoveModelPartProcess Process
MoveModelPartProcess::MoveModelPartProcess(ModelPart& rModelPart,
                    Parameters ThisParameters
                ):
    Process(),
    mrModelPart(rModelPart)
{
    Parameters default_parameters = Parameters(R"(
    {
        "origin"                        : [0.0,0.0,0.0],
        "rotation_point"                : [0.0,0.0,0.0],
        "rotation_axis"                 : [0.0,0.0,1.0],
        "rotation_angle"                : 0.0,
        "sizing_multiplier"             : 1.0

    })" );
    bool assign_rotation_point = false;
    if (ThisParameters.Has("rotation_point")) {
        assign_rotation_point = true;
    }
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);


    mOrigin = ThisParameters["origin"].GetVector();
    if (assign_rotation_point){
        mRotationPoint = ThisParameters["rotation_point"].GetVector();
    }
    else{
        mRotationPoint = mOrigin;
    }
    mRotationAngle = ThisParameters["rotation_angle"].GetDouble();
    mRotationAxis = ThisParameters["rotation_axis"].GetVector();
    mSizingMultiplier = ThisParameters["sizing_multiplier"].GetDouble();
}

void MoveModelPartProcess::Execute()
{
    KRATOS_TRY;

    Matrix translation_matrix = ZeroMatrix(4,4);
    GeometricalTransformationUtilities::CalculateTranslationMatrix(1.0, translation_matrix, mOrigin);

    Matrix rotation_matrix = ZeroMatrix(4,4);
    GeometricalTransformationUtilities::CalculateRotationMatrix(mRotationAngle, rotation_matrix, mRotationAxis, mRotationPoint);

    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode)
    {
        auto &r_coordinates = rNode.Coordinates();
        for (std::size_t i_dim = 0; i_dim < r_coordinates.size(); i_dim++){
            if (mSizingMultiplier > std::numeric_limits<double>::epsilon()){
                r_coordinates[i_dim] *= mSizingMultiplier;
            }
            if (norm_2(mOrigin) > std::numeric_limits<double>::epsilon()) {
                r_coordinates[i_dim] += translation_matrix(i_dim, 3);
            }
        }

        if (std::abs(mRotationAngle) > 0.0){
            Vector aux_coordinates_copy = r_coordinates;
            aux_coordinates_copy.resize(4, true);
            aux_coordinates_copy[3] = 1.0;
            Vector rotated_coordinates = prod(rotation_matrix, aux_coordinates_copy);

            for (std::size_t i_dim = 0; i_dim < 3; i_dim++){
                r_coordinates[i_dim] = rotated_coordinates[i_dim];
            }
        }
    });

    KRATOS_CATCH("");
}
}// Namespace Kratos
