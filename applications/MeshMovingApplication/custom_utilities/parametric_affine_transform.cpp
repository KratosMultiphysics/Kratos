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

// Project includes
#include "parametric_affine_transform.h"

namespace Kratos {


ParametricAffineTransform::ParametricAffineTransform(const Parameters axis,
                                                     const Parameters angle,
                                                     const Parameters referencePoint,
                                                     const Parameters translationVector)
    : mReferencePointFunction(referencePoint),
      mTranslationVectorFunction(translationVector)
{
    KRATOS_TRY

    VectorFunction<3> axis_function(axis);
    FunctionType angle_function(ExtractFunctionBody(angle));

    mQuaternionFunction = [axis_function, angle_function](const double x,
                                                          const double y,
                                                          const double z,
                                                          const double t,
                                                          const double X,
                                                          const double Y,
                                                          const double Z) mutable
    {
        KRATOS_TRY

        const auto axis = axis_function(x, y, z, t, X, Y, Z);
        const auto angle = angle_function.CallFunction(x, y, z, t, X, Y, Z);

        return Quaternion<double>::FromAxisAngle(
            axis[0],
            axis[1],
            axis[2],
            angle);

        KRATOS_CATCH("");
    };

    KRATOS_CATCH("");
}


ParametricAffineTransform::ParametricAffineTransform(const Parameters eulerAngles,
                                                     const Parameters referencePoint,
                                                     const Parameters translationVector)
    : mReferencePointFunction(referencePoint),
      mTranslationVectorFunction(translationVector)
{
    KRATOS_TRY

    VectorFunction<3> angle_function(eulerAngles);

    mQuaternionFunction = [angle_function](const double x,
                                           const double y,
                                           const double z,
                                           const double t,
                                           const double X,
                                           const double Y,
                                           const double Z) mutable
    {
        KRATOS_TRY
        const auto euler_angles = angle_function(x, y, z, t, X, Y, Z);
        return Quaternion<double>::FromEulerAngles(euler_angles);
        KRATOS_CATCH("");
    };

    KRATOS_CATCH("");
}


array_1d<double,3> ParametricAffineTransform::Apply(const array_1d<double,3>& rPoint,
                                                    const double t,
                                                    const double X,
                                                    const double Y,
                                                    const double Z)
{
    KRATOS_TRY

    // Evaluate current transformation parameters
    const Quaternion<double> quaternion = mQuaternionFunction(
        rPoint[0],
        rPoint[1],
        rPoint[2],
        t,
        X,
        Y,
        Z);

    const array_1d<double,3> reference_point = mReferencePointFunction(
        rPoint[0],
        rPoint[1],
        rPoint[2],
        t,
        X,
        Y,
        Z);

    const array_1d<double,3> translation_vector = mTranslationVectorFunction(
        rPoint[0],
        rPoint[1],
        rPoint[2],
        t,
        X,
        Y,
        Z);

    // Update rotation if necessary
    bool needs_updating = false;

    for (std::size_t i=0; i<4; ++i) {
        if (quaternion[i] != mQuaternion[i]) {
            mQuaternion = quaternion;
            needs_updating = true;
            break;
        }
    }

    for (std::size_t i=0; i<3; ++i) {
        if (reference_point[i] != mReferencePoint[i]) {
            needs_updating = true;
            break;
        }
    }

    if (needs_updating) {
        this->SetRotation(quaternion, reference_point);
    }

    // Update translation (it's more performant to just skip checking)
    this->SetTranslation(translation_vector);

    // Perform the transformation
    return AffineTransform::Apply(rPoint);

    KRATOS_CATCH("");
}


std::string ParametricAffineTransform::ExtractFunctionBody(const Parameters parameters)
{
    KRATOS_TRY

    std::string body;

    if (parameters.IsString()) {
        body = parameters.GetString();
    } else if (parameters.IsNumber()) {
        body = std::to_string(parameters.GetDouble());
    } else {
        KRATOS_ERROR << "Cannot extract function body from Parameters: " << parameters;
    }

    return body;

    KRATOS_CATCH("");
}


} // namespace Kratos