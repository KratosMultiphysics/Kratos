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

#pragma once

// System includes
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/mesh_moving_variables.h"
#include "utilities/math_utils.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

///@addtogroup MeshMovingApplication
///@{

/**
 * @class RotatingFrameUtility
 * @ingroup MeshMovingApplication
 * @brief Utility providing rigid-body rotational kinematics for model parts.
 *
 * @details
 * This utility provides helper functions to:
 *
 *  - Compute and assign rigid-body rotational velocity fields
 *  - Apply a rigid rotation to a model part based on reference coordinates
 *  - Compute the corresponding mesh displacement induced by the rotation
 *
 * The rotation is defined by:
 *  - an axis of rotation
 *  - a rotation angle or angular velocity
 *  - a center of rotation
 *
 * The utility operates purely at the kinematic level and **does not enforce
 * solver policies such as DOF fixity**. Such enforcement should be handled
 * at the process level if required.
 *
 * The implementation is designed to operate efficiently on model part nodes
 * and is compatible with Kratos parallel utilities.
 */

class KRATOS_API(MESH_MOVING_APPLICATION) RotatingFrameUtility
{
public:

    /**
     * @brief Assign rigid-body rotational velocity to nodes.
     *
     * Computes the velocity induced by rigid-body rotation:
     *
     *      v = ω × r
     *
     * where:
     *  - ω is the angular velocity vector
     *  - r is the position relative to the center of rotation
     *
     * The computed velocity is assigned to the VELOCITY variable of the nodes
     * in the provided model part.
     *
     * @param rModelPart Model part whose nodes receive the rotational velocity
     * @param rAxisOfRotation Axis of rotation
     * @param Omega Angular velocity magnitude
     * @param rCenterOfRotation Center of rotation
     */
    static void AssignRotationalVelocity(
        ModelPart& rModelPart,
        const array_1d<double,3>& rAxisOfRotation,
        const double Omega,
        const array_1d<double,3>& rCenterOfRotation);

    /**
     * @brief Apply rigid rotation and compute mesh displacement.
     *
     * Applies a rigid-body rotation to the nodes of the provided model part
     * using their reference coordinates. The resulting mesh displacement is
     * computed as:
     *
     *      MESH_DISPLACEMENT = rotated_position − initial_position
     *
     * The node coordinates are then updated consistently.
     *
     * @param rModelPart Model part whose nodes will be rotated
     * @param rAxisOfRotation Axis of rotation
     * @param Theta Rotation angle
     * @param rCenterOfRotation Center of rotation
     */
    static void ApplyRotationAndMeshDisplacement(
        ModelPart& rModelPart,
        const array_1d<double,3>& rAxisOfRotation,
        const double Theta,
        const array_1d<double,3>& rCenterOfRotation);

};

///@}

} // namespace Kratos