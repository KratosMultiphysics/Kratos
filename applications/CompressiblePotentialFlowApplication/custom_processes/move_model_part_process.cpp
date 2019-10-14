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
    mSizingMultiplier = ThisParameters["sizing_multiplier"].GetDouble();
}

void MoveModelPartProcess::Execute()
{
    KRATOS_TRY;

    #pragma omp parallel for
    for(int i = 0; i <  static_cast<int>(mrModelPart.NumberOfNodes()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        auto &r_coordinates = it_node->Coordinates();

        for (std::size_t i_dim = 0; i_dim<3;i_dim++){
            r_coordinates[i_dim] = mSizingMultiplier*r_coordinates[i_dim]+mOrigin[i_dim];
        }

        if (mRotationAngle != 0.0){
            array_1d<double, 3> old_coordinates = r_coordinates;
            // X-Y plane rotation
            r_coordinates[0] = mRotationPoint[0]+cos(mRotationAngle)*(old_coordinates[0]-mRotationPoint[0])-
                            sin(mRotationAngle)*(old_coordinates[1]-mRotationPoint[1]);
            r_coordinates[1] = mRotationPoint[1]+sin(mRotationAngle)*(old_coordinates[0]-mRotationPoint[0])-
                            cos(mRotationAngle)*(old_coordinates[1]-mRotationPoint[1]);
        }
    }

    KRATOS_CATCH("");
}
}// Namespace Kratos
