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

// Internal includes
#include "linear_transform.h"

namespace Kratos
{


LinearTransform::LinearTransform()
{
    KRATOS_TRY

    this->SetRotation(
        array_1d<double,3>({0.0, 0.0, 1.0}),
        0.0,
        ZeroVector(3)
    );

    this->SetTranslation(ZeroVector(3));

    KRATOS_CATCH("");
}


LinearTransform::LinearTransform(const array_1d<double,3>& rAxis,
                                 const double angle,
                                 const array_1d<double,3>& rReferencePoint,
                                 const array_1d<double,3>& rTranslationVector)
{
    KRATOS_TRY

    this->SetRotation(
        rAxis,
        angle,
        rReferencePoint);

    this->SetTranslation(rTranslationVector);

    KRATOS_CATCH("");
}


LinearTransform::LinearTransform(const array_1d<double,3>& rEulerAngles,
                                 const array_1d<double,3>& rReferencePoint,
                                 const array_1d<double,3>& rTranslationVector)
{
    KRATOS_TRY

    this->SetRotation(
        rEulerAngles,
        rReferencePoint);

    this->SetTranslation(rTranslationVector);

    KRATOS_CATCH("");
}


LinearTransform::LinearTransform(const Quaternion<double>& rQuaternion,
                                 const array_1d<double,3>& rReferencePoint,
                                 const array_1d<double,3>& rTranslationVector)
{
    KRATOS_TRY

    this->SetRotation(
        rQuaternion,
        rReferencePoint);

    this->SetTranslation(rTranslationVector);

    KRATOS_CATCH("");
}


void LinearTransform::SetRotation(const array_1d<double,3>& rAxis,
                                  const double angle,
                                  const array_1d<double,3>& rReferencePoint)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(std::abs(norm_2(rAxis)) < 1e-15)
    << "Axis of rotation must not be a null vector!";

    auto quaternion = Quaternion<double>::FromAxisAngle(
        rAxis[0],
        rAxis[1],
        rAxis[2],
        angle);

    this->SetRotation(quaternion, rReferencePoint);

    KRATOS_CATCH("");
}


void LinearTransform::SetRotation(const array_1d<double,3>& rEulerAngles,
                                  const array_1d<double,3>& rReferencePoint)
{
    KRATOS_TRY

    auto quaternion = Quaternion<double>::FromEulerAngles(rEulerAngles);

    this->SetRotation(quaternion, rReferencePoint);

    KRATOS_CATCH("");
}


void LinearTransform::SetRotation(const Quaternion<double>& rQuaternion,
                                  const array_1d<double,3>& rReferencePoint)
{
    KRATOS_TRY

    mRotationMatrix.resize(3, 3);
    mReferencePoint = rReferencePoint;

    rQuaternion.ToRotationMatrix(mRotationMatrix);

    KRATOS_CATCH("");
}


void LinearTransform::SetTranslation(const array_1d<double,3>& rTranslationVector)
{
    mTranslationVector = rTranslationVector;
}


array_1d<double,3> LinearTransform::Apply(const array_1d<double,3>& rPoint) const
{
    KRATOS_TRY

    return prod(mRotationMatrix, rPoint - mReferencePoint) + mReferencePoint + mTranslationVector;

    KRATOS_CATCH("");
}


} // namespace Kratos