//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//	                Kratos default license: kratos/license.txt
//
//  Main Authors:   Máté Kelemen
//

#pragma once

// Project includes
#include "utilities/quaternion.h"

namespace Kratos
{


namespace Detail {


} // namespace Detail

///@name Kratos Classes
///@{


/** @brief Class for applying affine transformations
 *  @details for now, only a rotation followed by a translation on an array of size 3 is implemented.
 *           The net transformation is equivalent to:
 *           1) Translation to the reference frame (offset the origin)
 *           2) Specified rotation
 *           3) Reverse translation from the reference frame (undo origin offset)
 *           4) Specified translation
 *  @note angles in radians.
 */
class KRATOS_API(MESH_MOVING_APPLICATION) AffineTransform
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AffineTransform);

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Default constructor sets up an identity transform.
    AffineTransform();

    /** @brief Construct from axis and angle.
     *  @param rAxis axis of rotation (direction vector).
     *  @param angle angle of rotation (in radians).
     *  @param rReferencePoint a point on the axis of rotation.
     *  @param rTranslationVector translation vector.
     */
    AffineTransform(const array_1d<double,3>& rAxis,
                    const double angle,
                    const array_1d<double,3>& rReferencePoint,
                    const array_1d<double,3>& rTranslationVector);

    /** @brief Construct from euler angles.
     *  @param rEulerAngles euler angles (radians).
     *  @param rReferencePoint origin of rotation.
     *  @param rTranslationVector translation vector.
     *  @note The euler angles follow the convention specified by @ref Quaternion (Z, -X', Z").
     */
    AffineTransform(const array_1d<double,3>& rEulerAngles,
                    const array_1d<double,3>& rReferencePoint,
                    const array_1d<double,3>& rTranslationVector);

    /** @brief Construct from a quaternion.
     *  @param rQuaternion quaternion defining the rotation angle and axis.
     *  @param rReferencePoint a point on the axis of rotation.
     *  @param rTranslationVector translation vector.
     */
    AffineTransform(const Quaternion<double>& rQuaternion,
                    const array_1d<double,3>& rReferencePoint,
                    const array_1d<double,3>& rTranslationVector);

    ///@}
    ///@name Operations
    ///@{

    /** @brief Set rotation from an axis and an angle.
     *  @param rAxis axis of rotation.
     *  @param angle angle of rotation (radians).
     *  @param rReferencePoint a point on the axis of rotation.
     */
    void SetRotation(const array_1d<double,3>& rAxis,
                     const double angle,
                     const array_1d<double,3>& rReferencePoint);

    /** @brief Set rotation from euler angles.
     *  @param rEulerAngles euler angles (radians).
     *  @param rReferencePoint a point on the axis of rotation.
     *  @note The euler angles follow the convention specified by @ref Quaternion (Z, -X', Z").
     */
    void SetRotation(const array_1d<double,3>& rEulerAngles,
                     const array_1d<double,3>& rReferencePoint);

    /** @brief Set rotation from a quaternion.
     *  @param rQuaternion quaternion defining the axis and angle of rotation.
     *  @param rReferencePoint a point on the axis of rotation.
     */
    void SetRotation(const Quaternion<double>& rQuaternion,
                     const array_1d<double,3>& rReferencePoint);

    /** @brief Set translation part of the transformation.
     *  @param rTranslationVector translation vector.
     */
    void SetTranslation(const array_1d<double,3>& rTranslationVector);

    /// @brief Return the transformed version of the input vector.
    array_1d<double,3> Apply(const array_1d<double,3>& rPoint) const
    {
        return prod(mRotationMatrix, rPoint - mReferencePoint) + mReferencePoint + mTranslationVector;
    }

    ///@}

protected:
    ///@name Member Variables
    ///@{

    array_1d<double,3> mReferencePoint;

    array_1d<double,3> mTranslationVector;

    Matrix mRotationMatrix;

    ///@}
}; // class AffineTransform


///@}

} // namespace Kratos

