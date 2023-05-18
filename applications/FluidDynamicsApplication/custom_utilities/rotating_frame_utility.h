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

#if !defined(KRATOS_ROTATING_FRAME_UTILITY_H)
#define KRATOS_ROTATING_FRAME_UTILITY_H

// System includes
#include <string>
#include <iostream>

// External includes
#include "includes/mesh_moving_variables.h"

// Include kratos definitions
#include "includes/define.h"

// Application includes
#include "utilities/variable_utils.h"

// Project includes
#include "utilities/builtin_timer.h"
#include "utilities/math_utils.h"

namespace Kratos {
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

/**
 * @class RotateFrameUtility
 * @ingroup Kratos Core
 * @brief Utility for rotating a frame with a rotating object.
 * @details A utility that rotates a frame (model part) with a rotating object (solid)
 * based on the specified angular velocity and axis of rotation. It enables the
 * transformation of coordinates and displacements between the rotating and
 * non-rotating frames, ensuring consistent updates of the frame's entities
 * during simulation. Designed for use within the Kratos framework for various
 * applications such as sliding mesh problems and rotating machinery simulations.
 * @note Thread-safe and supports parallel execution using Kratos block for each.
 * @author Sebastian Ares de Parga Regalado
*/

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) RotatingFrameUtility
{
    public:
      ///@name Type Definitions
      ///@{

      ///@}
      ///@name Static Operations
      ///@{

      /**
       * @brief Apply velocity to rotating object.
       * 
       * This function applies velocity to a rotating object model part.
       * The velocity is computed based on the angular velocity and the
       * axis of rotation. The function calculates the radius of each node,
       * and then applies the corresponding velocity component-wise. The
       * velocity components are fixed to remain constant during the simulation,
       * and the solution step values for the velocity components are set
       * accordingly.
       * 
       * @param rRotatingObjectModelPart The rotating object model part.
       * @param Omega The angular velocity.
       * @param rAxisOfRotation The axis of rotation.
       */
      static void ApplyVelocityToRotatingObject(
          ModelPart& rRotatingObjectModelPart,
          const double Omega,
          const array_1d<double, 3>& rAxisOfRotation);

      /**
       * @brief Apply rotation and mesh displacement to a rotating frame.
       * This function applies rotation and mesh displacement to a rotating frame model 
       * part based on the provided rotation axis, angle and rotation center. 
       * It constructs a rotation matrix based on the rotation axis and angle. 
       * The rotation is applied to each node, which is then displaced 
       * by the difference between the rotated coordinates and the initial coordinates. 
       * The solution step values for the mesh displacement components are set accordingly, 
       * and the mesh displacement components are fixed to remain constant during the simulation. 
       * Finally, the node coordinates are updated with the rotated coordinates.
       * @param rRotatingFrameModelPart The rotating frame model part to which the rotation and mesh displacement will be applied.
       * @param rAxisOfRotation The rotation axis vector.
       * @param rTheta The rotation angle.
       * @param rCenterOfRotation The center of rotation.
       */
      static void ApplyRotationAndMeshDisplacement(
          ModelPart& rRotatingFrameModelPart,
          const array_1d<double, 3>& rAxisOfRotation, 
          const double& rTheta, 
          const array_1d<double, 3>& rCenterOfRotation);
          
}; // Class RotatingFrameUtility

///@}

}  // namespace Kratos.

#endif // KRATOS_ROTATING_FRAME_UTILITY_H