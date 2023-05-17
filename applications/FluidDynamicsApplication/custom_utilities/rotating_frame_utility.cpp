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
    ModelPart& rRotatingObjectModelPart,
    const double Omega,
    const array_1d<double, 3>& rAxisOfRotation)
{
    KRATOS_TRY;

    // Get the array of nodes in the rotating object model part
    auto& nodes_array = rRotatingObjectModelPart.Nodes();

    // Declare thread-local variables outside the block_for_each
    thread_local double radius;
    thread_local double velocity_magnitude;
    thread_local array_1d<double, 3> velocity_direction;
    thread_local array_1d<double, 3> velocity_vector;
    thread_local array_1d<double, 3> velocity;

    // Apply the velocity calculations to each node in parallel
    block_for_each(nodes_array, [&](Node& rNode) {
        radius = std::sqrt(rNode.X() * rNode.X() + rNode.Y() * rNode.Y() + rNode.Z() * rNode.Z());

        // Apply component by component the corresponding velocity.
        velocity_magnitude = Omega * radius;
        velocity_direction = MathUtils<double>::CrossProduct(rAxisOfRotation, rNode.Coordinates());
        velocity_vector = velocity_direction / norm_2(velocity_direction);
        velocity = velocity_magnitude * velocity_vector;

        // Set the solution step values for the velocity components using FastGetSolutionStepValue
        rNode.FastGetSolutionStepValue(VELOCITY_X) = velocity[0];
        rNode.FastGetSolutionStepValue(VELOCITY_Y) = velocity[1];
        rNode.FastGetSolutionStepValue(VELOCITY_Z) = velocity[2];

        // Fix the velocity components to ensure they remain constant during the simulation
        rNode.Fix(VELOCITY_X);
        rNode.Fix(VELOCITY_Y);
        rNode.Fix(VELOCITY_Z);
    });

    KRATOS_CATCH("");
}

void RotatingFrameUtility::ApplyRotationAndMeshDisplacement(
    ModelPart& rRotatingFrameModelPart,
    const array_1d<double, 3>& rAxis, 
    const double& rTheta, 
    const array_1d<double, 3>& rCenter)
{
    KRATOS_TRY;

    // Quaternion constants
    const double a = std::cos(rTheta / 2);
    const double b = -rAxis[0] * std::sin(rTheta / 2);
    const double c = -rAxis[1] * std::sin(rTheta / 2);
    const double d = -rAxis[2] * std::sin(rTheta / 2);

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

    // Apply the rotation and mesh displacement calculations to each node in parallel
    block_for_each(rRotatingFrameModelPart.Nodes(), [&](Node& rNode) {
        // Getting the initial coordinates of the node
        auto& point = rNode.GetInitialPosition().Coordinates();
        // point[0] = node.X0();
        // point[0] = node.X0();
        // point[0] = node.X0();

        // Shifting the rotation center to the origin
        point -= rCenter;

        // Applying the rotation
        array_1d<double, 3> rotated_point = prod(rot_matrix, point);

        // Shifting the point back and updating the node coordinates
        rotated_point += rCenter;
        rNode.FastGetSolutionStepValue(DISPLACEMENT_X) = rotated_point[0] - rNode.X0();
        rNode.FastGetSolutionStepValue(DISPLACEMENT_Y) = rotated_point[1] - rNode.Y0();
        rNode.FastGetSolutionStepValue(DISPLACEMENT_Z) = rotated_point[2] - rNode.Z0();
        rNode.Fix(DISPLACEMENT_X);
        rNode.Fix(DISPLACEMENT_Y);
        rNode.Fix(DISPLACEMENT_Z);
        rNode.X() = rotated_point[0];
        rNode.Y() = rotated_point[1];
        rNode.Z() = rotated_point[2];
    });

    KRATOS_CATCH("");
}


// void RotatingFrameUtility::ApplyRotationAndMeshDisplacement(
//     ModelPart& rRotatingFrameModelPart,
//     const Matrix& rTransformationMatrix)
// {
//     KRATOS_TRY;

//     // Get the array of nodes in the rotating frame model part
//     auto& nodes_array = rRotatingFrameModelPart.Nodes();

//     // Declare thread-local variables outside the block_for_each
//     thread_local array_1d<double, 4> homogeneous_coordinates;
//     thread_local array_1d<double, 4> rotated_coordinates;
//     thread_local array_1d<double, 3> rotated_coor;
//     thread_local double mesh_displacement_x;
//     thread_local double mesh_displacement_y;
//     thread_local double mesh_displacement_z;

//     // Apply the rotation and mesh displacement calculations to each node in parallel
//     block_for_each(nodes_array, [&](Node& rNode) {

//         // Convert point to homogeneous coordinates
//         // homogeneous_coordinates.resize(4);
//         homogeneous_coordinates[0] = rNode.X0();
//         homogeneous_coordinates[1] = rNode.Y0();
//         homogeneous_coordinates[2] = rNode.Z0();
//         homogeneous_coordinates[3] = 1.0;

//         // Apply transformation to homogeneous coordinates
//         // rotated_coordinates.resize(4);
//         noalias(rotated_coordinates) = prod(rTransformationMatrix, homogeneous_coordinates);

//         // Convert back to 3D coordinates
//         // rotated_coor.resize(3);
//         rotated_coor[0] = rotated_coordinates[0] / rotated_coordinates[3];
//         rotated_coor[1] = rotated_coordinates[1] / rotated_coordinates[3];
//         rotated_coor[2] = rotated_coordinates[2] / rotated_coordinates[3];

//         // Calculate mesh displacement
//         mesh_displacement_x = rotated_coor[0] - rNode.X0();
//         mesh_displacement_y = rotated_coor[1] - rNode.Y0();
//         mesh_displacement_z = rotated_coor[2] - rNode.Z0();

//         // Set the solution step values for mesh displacement using FastGetSolutionStepValue
//         rNode.FastGetSolutionStepValue(MESH_DISPLACEMENT_X) = mesh_displacement_x;
//         rNode.FastGetSolutionStepValue(MESH_DISPLACEMENT_Y) = mesh_displacement_y;
//         rNode.FastGetSolutionStepValue(MESH_DISPLACEMENT_Z) = mesh_displacement_z;

//         // Fix the mesh displacement components to ensure they remain constant during the simulation
//         rNode.Fix(MESH_DISPLACEMENT_X);
//         rNode.Fix(MESH_DISPLACEMENT_Y);
//         rNode.Fix(MESH_DISPLACEMENT_Z);

//         // Update the coordinates of the node using FastGetSolutionStepValue
//         rNode.X() = rotated_coor[0];
//         rNode.Y() = rotated_coor[1];
//         rNode.Z() = rotated_coor[2];
//     });

//     KRATOS_CATCH("");
// }

}  // namespace Kratos.

