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

#ifndef KRATOS_MESH_MOVING_PARAMETRIC_LINEAR_TRANSFORMATION_INCLUDED
#define KRATOS_MESH_MOVING_PARAMETRIC_LINEAR_TRANSFORMATION_INCLUDED

// Project includes
#include "linear_transform.h"
#include "includes/kratos_parameters.h"
#include "utilities/function_parser_utility.h"

namespace Kratos
{

///@name Kratos Classes
///@{


/** Class for applying parametrically defined linear transformations
 *  @details see @ref{LinearTransform} for transform conventions. This class
 *  derives from @ref{LinearTransform}, but every parameter is parsed from
 *  an input @ref{Parameters} to a @ref{GenericFunctionUtility} that is
 *  evaluated on each application of the transform.
 */
class KRATOS_API(MESH_MOVING_APPLICATION) ParametricLinearTransform : protected LinearTransform
{
public:
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ParametricLinearTransform);

    using FunctionType = GenericFunctionUtility;

    ///@}
    ///@name Life Cycle
    ///@{
    
    /** Construct via axis & angle
     *  @param rAxis axis of rotation (array of size 3)
     *  @param rAngle angle of rotation (in radians, scalar)
     *  @param rReferencePoint a point on the axis of rotation (array of size 3)
     *  @param rTranslationVector translation vector (array of size 3)
     */
    ParametricLinearTransform(const Parameters axis,
                              const Parameters angle,
                              const Parameters referencePoint,
                              const Parameters translationVector);

    /** Construct via euler angles
     *  @param rEulerAngles euler angles (radians, array of size 3)
     *  @param rReferencePoint origin of rotation (array of size 3)
     *  @param rTranslationVector translation vector (array of size 3)
     *  @note The euler angles follow the convention specified by @ref{Quaternion} (Z, -X', Z")
     */
    ParametricLinearTransform(const Parameters eulerAngles,
                              const Parameters referencePoint,
                              const Parameters translationVector);

    ///@}
    ///@name Operations
    ///@{

    /** Evaluate the transformation at the given location and apply it to a position vector, returning its transformed version
     *  @note the transformation is not updated unless at least one transformation-related
     *  quantity has changed since the last call to ParametricLinearTransform::Apply.
     *  @param rPoint point to be transformed; it also acts as the current position input to @ref{GenericFunctionUtility}
     *  @param t time, see @ref{GenericFunctionUtility}
     *  @param X initial x-coordinate, see @ref{GenericFunctionUtility}
     *  @param Y initial y-coordinate, see @ref{GenericFunctionUtility}
     *  @param Z initial z-coordinate, see @ref{GenericFunctionUtility}
     */
    array_1d<double,3> Apply(const array_1d<double,3>& rPoint,
                             const double t,
                             const double X = 0.0,
                             const double Y = 0.0,
                             const double Z = 0.0);

    ///@}

private:

    ///@name Private type definitions
    ///@{

    // Necessary to parse constants
    static std::string ExtractFunctionBody(const Parameters parameters);

    /// Class for storing and evaluating an array of @ref{GenericFunctionUtility}
    template <std::size_t ArraySize>
    class VectorFunction : public std::array<FunctionType::Pointer,ArraySize>
    {
    public:
        VectorFunction(const Parameters rFunctions)
        {
            KRATOS_TRY

            KRATOS_ERROR_IF_NOT(rFunctions.IsArray()) << "VectorFunction can only be constructed from a parameter array";

            // Assign components - bounds checks are performed in Parameters::GetArrayItem
            for (std::size_t i=0; i<ArraySize; ++i) {
                const std::string function_string = ExtractFunctionBody(rFunctions.GetArrayItem(i));
                this->at(i) = std::make_shared<GenericFunctionUtility>(function_string);
            }

            KRATOS_CATCH("");
        }

        VectorFunction(const VectorFunction<ArraySize>& rOther) = default;

        array_1d<double,ArraySize> operator()(const double x,
                                              const double y,
                                              const double z,
                                              const double t,
                                              const double X,
                                              const double Y,
                                              const double Z)
    {
        KRATOS_TRY

        array_1d<double,ArraySize> output;

        for (std::size_t i=0; i<ArraySize; ++i) {
            output[i] = this->at(i)->CallFunction(x, y, z, t, X, Y, Z);
        }

        return output;

        KRATOS_CATCH("");
    }
    }; // class VectorFunction

    ///@}
    ///@name Member Variables
    ///@{

    VectorFunction<3> mReferencePointFunction;

    VectorFunction<3> mTranslationVectorFunction;

    // The quaternion getter has to be defined at runtime, depending on how the rotation was defined
    std::function<Quaternion<double>(const double, const double, const double, const double, const double, const double, const double)> mQuaternionFunction;

    // Store the last quaternion to compare to new ones
    Quaternion<double> mQuaternion;

    ///@}
}; // class ParametricLinearTransform


///@}

} // namespace Kratos

#endif