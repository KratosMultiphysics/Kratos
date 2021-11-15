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

#ifndef KRATOS_MESH_MOVING_LINEAR_TRANSFORMATION_INCLUDED
#define KRATOS_MESH_MOVING_LINEAR_TRANSFORMATION_INCLUDED

// Project includes
#include "utilities/quaternion.h"

namespace Kratos
{

///@name Kratos Classes
///@{


/** Class for applying linear transformations
 *  @details for now, only a rotation followed by a translation on an array of size 3 is implemented.
 *  The net transformation is equivalent to:
 *  1) Translation to the reference frame (offset the origin)
 *  2) Specified rotation
 *  3) Reverse translation from the reference frame (undo origin offset)
 *  4) Specified translation
 *  @note angles in radians
 */
class KRATOS_API(MESH_MOVING_APPLICATION) LinearTransform
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(LinearTransform);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor -> identity transform
    LinearTransform();
    
    /** Construct via axis & angle
     *  @param rAxis axis of rotation (direction vector)
     *  @param angle angle of rotation (in radians)
     *  @param rReferencePoint a point on the axis of rotation
     *  @param rTranslationVector translation vector
     */
    LinearTransform(const array_1d<double,3>& rAxis,
                    const double angle,
                    const array_1d<double,3>& rReferencePoint,
                    const array_1d<double,3>& rTranslationVector);

    /** Construct via euler angles
     *  @param rEulerAngles euler angles (radians)
     *  @param rReferencePoint origin of rotation
     *  @param rTranslationVector translation vector
     *  @note The euler angles follow the convention specified by @ref{Quaternion} (Z, -X', Z")
     */
    LinearTransform(const array_1d<double,3>& rEulerAngles,
                    const array_1d<double,3>& rReferencePoint,
                    const array_1d<double,3>& rTranslationVector);

    /** Construct via quaternion
     *  @param rQuaternion quaternion defining the rotation angle and axis
     *  @param rReferencePoint a point on the axis of rotation
     *  @param rTranslationVector translation vector
     */
    LinearTransform(const Quaternion<double>& rQuaternion,
                    const array_1d<double,3>& rReferencePoint,
                    const array_1d<double,3>& rTranslationVector);

    LinearTransform(const LinearTransform& rOther) = default;

    LinearTransform(LinearTransform&& rOther) = default;

    LinearTransform& operator=(const LinearTransform& rOther) = default;

    ///@}
    ///@name Operations
    ///@{

    /** Set rotation via axis & angle
     *  @param rAxis axis of rotation
     *  @param angle angle of rotation (radians)
     *  @param rReferencePoint a point on the axis of rotation
     */
    void SetRotation(const array_1d<double,3>& rAxis,
                     const double angle,
                     const array_1d<double,3>& rReferencePoint);

    /** Set rotation via euler angles
     *  @param rEulerAngles euler angles (radians)
     *  @param rReferencePoint a point on the axis of rotation
     *  @note The euler angles follow the convention specified by @ref{Quaternion} (Z, -X', Z")
     */
    void SetRotation(const array_1d<double,3>& rEulerAngles,
                     const array_1d<double,3>& rReferencePoint);

    /** Set rotation via quaternion
     *  @param rQuaternion quaternion specifying the axis and angle
     *  @param rReferencePoint a point on the axis of rotation
     */
    void SetRotation(const Quaternion<double>& rQuaternion,
                     const array_1d<double,3>& rReferencePoint);

    /** Set translation part of the transformation
     *  @param rTranslationVector translation vector
     */
    void SetTranslation(const array_1d<double,3>& rTranslationVector);

    /// Return the transformed version of the input vector
    array_1d<double,3> Apply(const array_1d<double,3>& rPoint) const;

    ///@}

protected:
    ///@name Member Variables
    ///@{

    array_1d<double,3> mReferencePoint;

    array_1d<double,3> mTranslationVector;

    Matrix mRotationMatrix;

    ///@}
}; // class LinearTransform


///@}

} // namespace Kratos

#endif // KRATOS_MESH_MOVING_LINEAR_TRANSFORMATION_INCLUDED