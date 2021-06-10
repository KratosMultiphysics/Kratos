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

#ifndef KRATOS_IMPOSE_MESH_MOTION_PROCESS_INCLUDED
#define KRATOS_IMPOSE_MESH_MOTION_PROCESS_INCLUDED

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "utilities/quaternion.h"

namespace Kratos
{

///@name Kratos Classes
///@{


/** Impose a rotation followed by translation on a ModelPart
 *  @details The transformation is equivalent to:
 *  1) Translation to the reference frame (offset the origin)
 *  2) Specified rotation
 *  3) Reverse translation from the reference frame (undo origin offset)
 *  4) Specified translation
 *  @note angles in radians
 */
class KRATOS_API(MESH_MOVING_APPLICATION) ImposeMeshMotionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SetHMapProcess
    KRATOS_CLASS_POINTER_DEFINITION(ImposeMeshMotionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor
     *  @details The rotation can be defined by either "euler_angles"
     *  or a "rotation_axis" and "rotation_angle" pair. Default parameters:
     *  {
     *      "model_part_name"       : "",
     *      "interval"              : [0.0, 1e30],
     *      "rotation_definition"   : "rotation_axis",
     *      "euler_angles"          : [0.0, 0.0, 0.0],
     *      "rotation_axis"         : [0.0, 0.0, 1.0],
     *      "reference_point"       : [0.0, 0.0, 0.0]
     *      "rotation_angle"        : 0,
     *      "translation_vector"    : [0.0, 0.0, 0.0],
     *  }
     *  @note The euler angles follow the convention specified by @ref{Quaternion} (Z, -X', Z")
     */
    ImposeMeshMotionProcess(ModelPart& rModelPart, Parameters parameters);

    ///@}
    ///@name Operations
    ///@{

    void LoadFromParameters(Parameters parameters);

    void LoadFromAxisAndAngle(const array_1d<double,3>& rRotationAxis,
                              double rotationAngle,
                              const array_1d<double,3>& rReferencePoint,
                              const array_1d<double,3>& rTranslationVector);

    void LoadFromEulerAngles(const array_1d<double,3>& rEulerAngles,
                             const array_1d<double,3>& rReferencePoint,
                             const array_1d<double,3>& rTranslationVector);

    void LoadFromQuaternion(const Quaternion<double>& rQuaternion,
                            const array_1d<double,3>& rReferencePoint,
                            const array_1d<double,3>& rTranslationVector);

    virtual void ExecuteInitializeSolutionStep() override;

    virtual const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Class name as string
    virtual std::string Info() const override;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rStream) const override;

    ///@}

private:
    ///@name Private operations
    ///@{

    void Transform(array_1d<double,3>& rPoint) const
    {
        rPoint = prod(mRotationMatrix, rPoint - mReferencePoint);
        rPoint += mReferencePoint + mTranslationVector;
    }

    ///@}
    ///@name Member Variables
    ///@{
        
    ModelPart& mrModelPart;

    Matrix mRotationMatrix;

    array_1d<double,3> mReferencePoint;

    array_1d<double,3> mTranslationVector;

    array_1d<double,2> mInterval;

    ///@}
};

///@}

///@}
///@name Input and output
///@{

/// Dump info to output stream
std::ostream& operator<<(std::ostream& rStream,
                         const ImposeMeshMotionProcess& rThis);

///@}

} // namespace Kratos

#endif // KRATOS_IMPOSE_MESH_MOTION_PROCESS_INCLUDED