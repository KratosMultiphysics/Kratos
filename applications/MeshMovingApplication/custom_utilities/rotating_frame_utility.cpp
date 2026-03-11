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

// External includes

// Project includes
#include "custom_utilities/rotating_frame_utility.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{

void RotatingFrameUtility::AssignRotationalVelocity(
    ModelPart& rModelPart,
    const array_1d<double,3>& rAxisOfRotation,
    const double Omega,
    const array_1d<double,3>& rCenterOfRotation)
{
    KRATOS_TRY

    const array_1d<double,3> angular_velocity_vector = Omega * rAxisOfRotation;

    block_for_each(rModelPart.Nodes(), [&](Node& rNode) {
        const auto& r_point = rNode.Coordinates();
        const array_1d<double,3> position_vector = r_point - rCenterOfRotation;
        const array_1d<double,3> velocity_vector =
            MathUtils<double>::CrossProduct(angular_velocity_vector, position_vector);

        noalias(rNode.FastGetSolutionStepValue(VELOCITY)) = velocity_vector;
    });

    KRATOS_CATCH("");
}

void RotatingFrameUtility::ApplyRotationAndMeshDisplacement(
    ModelPart& rModelPart,
    const array_1d<double,3>& rAxisOfRotation,
    const double Theta,
    const array_1d<double,3>& rCenterOfRotation)
{
    KRATOS_TRY

    const double a = std::cos(Theta / 2.0);
    const double b = -rAxisOfRotation[0] * std::sin(Theta / 2.0);
    const double c = -rAxisOfRotation[1] * std::sin(Theta / 2.0);
    const double d = -rAxisOfRotation[2] * std::sin(Theta / 2.0);

    BoundedMatrix<double,3,3> rot_matrix;
    rot_matrix(0,0) = a*a + b*b - c*c - d*d;
    rot_matrix(0,1) = 2.0 * (b*c - a*d);
    rot_matrix(0,2) = 2.0 * (b*d + a*c);
    rot_matrix(1,0) = 2.0 * (b*c + a*d);
    rot_matrix(1,1) = a*a + c*c - b*b - d*d;
    rot_matrix(1,2) = 2.0 * (c*d - a*b);
    rot_matrix(2,0) = 2.0 * (b*d - a*c);
    rot_matrix(2,1) = 2.0 * (c*d + a*b);
    rot_matrix(2,2) = a*a + d*d - b*b - c*c;

    block_for_each(rModelPart.Nodes(), [&](Node& rNode) {
        const auto& r_initial_position = rNode.GetInitialPosition().Coordinates();
        const array_1d<double,3> centered_point = r_initial_position - rCenterOfRotation;

        array_1d<double,3> rotated_point = prod(centered_point, rot_matrix);
        rotated_point += rCenterOfRotation;

        noalias(rNode.FastGetSolutionStepValue(MESH_DISPLACEMENT)) =
            rotated_point - rNode.GetInitialPosition();

        noalias(rNode.Coordinates()) = rotated_point;
    });

    KRATOS_CATCH("");
}

} // namespace Kratos