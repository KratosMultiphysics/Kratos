//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Sebastian Ares de Parga Regalado
//

// System includes
#include <string>
#include <iostream>

// External includes

// Application includes
#include "rotating_frame_utility.h"

namespace Kratos 
{

void RotatingFrameUtility::ApplyVelocityToRotatingObject(
    ModelPart& rModelPart,
    const array_1d<double, 3>& rAxisOfRotation, 
    const double& rOmega, 
    const array_1d<double, 3>& rCenterOfRotation)
{
    KRATOS_TRY;

    // Creating the angular velocity vector
    array_1d<double, 3> angular_velocity_vector = rOmega * rAxisOfRotation;

    thread_local array_1d<double, 3> position_vector;
    thread_local array_1d<double, 3> velocity_vector;

    // Apply the velocity calculations to each node in parallel
    block_for_each(rModelPart.Nodes(), [&](Node& rNode) {
        // Getting the current coordinates of the node
        auto& r_point = rNode.Coordinates();

        // Calculating the position vector (relative to the rotation center)
        position_vector = r_point - rCenterOfRotation;

        // Computing the velocity due to rotation (v = omega cross r)
        velocity_vector = MathUtils<double>::CrossProduct(angular_velocity_vector, position_vector);

        // Setting the node's velocity
        // rNode.FastGetSolutionStepValue(VELOCITY_X) = velocity_vector[0];
        // rNode.FastGetSolutionStepValue(VELOCITY_Y) = velocity_vector[1];
        // rNode.FastGetSolutionStepValue(VELOCITY_Z) = velocity_vector[2];
        rNode.FastGetSolutionStepValue(VELOCITY) = velocity_vector;

        // Fix the velocity components
        rNode.Fix(VELOCITY_X);
        rNode.Fix(VELOCITY_Y);
        rNode.Fix(VELOCITY_Z);
    });

    KRATOS_CATCH("");
}

void RotatingFrameUtility::ApplyRotationAndMeshDisplacement(
    ModelPart& rRotatingFrameModelPart,
    const array_1d<double, 3>& rAxisOfRotation, 
    const double& rTheta, 
    const array_1d<double, 3>& rCenterOfRotation)
{
    KRATOS_TRY;

    // Quaternion constants
    const double a = std::cos(rTheta / 2);
    const double b = -rAxisOfRotation[0] * std::sin(rTheta / 2);
    const double c = -rAxisOfRotation[1] * std::sin(rTheta / 2);
    const double d = -rAxisOfRotation[2] * std::sin(rTheta / 2);

    // Creating a quaternion rotation matrix
    BoundedMatrix<double, 3, 3> rot_matrix;
    rot_matrix(0, 0) = a*a+b*b-c*c-d*d;
    rot_matrix(0, 1) = 2*(b*c-a*d);
    rot_matrix(0, 2) = 2*(b*d+a*c);
    rot_matrix(1, 0) = 2*(b*c+a*d);
    rot_matrix(1, 1) = a*a+c*c-b*b-d*d;
    rot_matrix(1, 2) = 2*(c*d-a*b);
    rot_matrix(2, 0) = 2*(b*d-a*c);
    rot_matrix(2, 1) = 2*(c*d+a*b);
    rot_matrix(2, 2) = a*a+d*d-b*b-c*c;
    
    thread_local array_1d<double, 3> rotated_point;

    // Apply the rotation and mesh displacement calculations to each node in parallel
    block_for_each(rRotatingFrameModelPart.Nodes(), [&](Node& rNode) {
        // Getting the initial coordinates of the node
        auto& r_point = rNode.GetInitialPosition().Coordinates();

        // Shifting the rotation center to the origin
        auto centered_point = r_point - rCenterOfRotation;

        // Applying the rotation
        rotated_point = prod(centered_point, rot_matrix);

        // Shifting the point back and updating the node coordinates
        rotated_point += rCenterOfRotation;
        rNode.FastGetSolutionStepValue(MESH_DISPLACEMENT_X) = rotated_point[0] - rNode.X0();
        rNode.FastGetSolutionStepValue(MESH_DISPLACEMENT_Y) = rotated_point[1] - rNode.Y0();
        rNode.FastGetSolutionStepValue(MESH_DISPLACEMENT_Z) = rotated_point[2] - rNode.Z0();
        rNode.Fix(MESH_DISPLACEMENT_X);
        rNode.Fix(MESH_DISPLACEMENT_Y);
        rNode.Fix(MESH_DISPLACEMENT_Z);
        rNode.X() = rotated_point[0];
        rNode.Y() = rotated_point[1];
        rNode.Z() = rotated_point[2];
    });

    KRATOS_CATCH("");
}

}  // namespace Kratos.

