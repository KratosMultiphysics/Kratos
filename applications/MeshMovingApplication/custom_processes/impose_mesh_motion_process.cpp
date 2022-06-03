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

// Internal inlcudes
#include "impose_mesh_motion_process.h"

// Project includes
#include "includes/checks.h"
#include "utilities/parallel_utilities.h"
#include "includes/mesh_moving_variables.h"


namespace Kratos
{


ImposeMeshMotionProcess::ImposeMeshMotionProcess(ModelPart& rModelPart, Parameters parameters) :
    Process(),
    mrModelPart(rModelPart)
{
    KRATOS_TRY;

    this->LoadFromParameters(parameters);

    KRATOS_CATCH("");
}


void ImposeMeshMotionProcess::LoadFromParameters(Parameters parameters)
{
    KRATOS_TRY;

    parameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
    Vector vector_parameter;

    // Parse interval
    vector_parameter = parameters.GetValue("interval").GetVector();
    KRATOS_CHECK(vector_parameter.size() == 2);
    mInterval[0] = vector_parameter[0];
    mInterval[1] = vector_parameter[1];

    // Parse translation
    array_1d<double,3> translation_vector;
    vector_parameter = parameters.GetValue("translation_vector").GetVector();
    if (vector_parameter.size() != 3) {
        KRATOS_ERROR << "'translation_vector' must be of size 3, but got" << vector_parameter.size();
    }

    translation_vector[0] = vector_parameter[0];
    translation_vector[1] = vector_parameter[1];
    translation_vector[2] = vector_parameter[2];

    // Parse reference point
    array_1d<double,3> reference_point;
    vector_parameter = parameters.GetValue("reference_point").GetVector();
    if (vector_parameter.size() != 3) {
        KRATOS_ERROR << "'reference_point' must be of size 3, but got" << vector_parameter.size();
    }

    reference_point[0] = vector_parameter[0];
    reference_point[1] = vector_parameter[1];
    reference_point[2] = vector_parameter[2];

    // Parse rotation depending on how it is defined
    std::string rotation_definition = parameters.GetValue("rotation_definition").GetString();

    // Parse rotation from axis and angle of rotation
    if (rotation_definition == "rotation_axis") {
        array_1d<double,3> rotation_axis;
        vector_parameter = parameters.GetValue("rotation_axis").GetVector();
        if (vector_parameter.size() != 3) {
            KRATOS_ERROR << "'rotation_axis' must be of size 3, but got" << vector_parameter.size();
        }

        rotation_axis[0] = vector_parameter[0];
        rotation_axis[1] = vector_parameter[1];
        rotation_axis[2] = vector_parameter[2];

        double rotation_angle = parameters.GetValue("rotation_angle").GetDouble();

        this->LoadFromAxisAndAngle(
            rotation_axis,
            rotation_angle,
            reference_point,
            translation_vector
        );
    } // if rotation_definition == "rotation_axis"

    // Parse rotation from euler angles
    else if (rotation_definition == "euler_angles") {
        array_1d<double,3> euler_angles;
        vector_parameter = parameters.GetValue("euler_angles").GetVector();
        if (vector_parameter.size() != 3) {
            KRATOS_ERROR << "'euler_angles' must be of size 3, but got" << vector_parameter.size();
        }

        euler_angles[0] = vector_parameter[0];
        euler_angles[1] = vector_parameter[1];
        euler_angles[2] = vector_parameter[2];

        this->LoadFromEulerAngles(
            euler_angles,
            reference_point,
            translation_vector
        );
    } // if rotation_definition == "euler_angles"

    // Rotation definition not implemented
    else {
        KRATOS_ERROR << "Invalid parameter for 'rotation_definition': " << rotation_definition
        << " (expecting 'rotation_axis' or 'euler_angles')";
    } // rotation_definition else

    KRATOS_CATCH("");
}


void ImposeMeshMotionProcess::LoadFromAxisAndAngle(const array_1d<double,3>& rRotationAxis,
                                                   double rotationAngle,
                                                   const array_1d<double,3>& rReferencePoint,
                                                   const array_1d<double,3>& rTranslationVector)
{
    KRATOS_TRY;

    // Check arguments
    if (std::abs(norm_2(rRotationAxis)) < 1e-15) {
        KRATOS_ERROR << "Axis of rotation must not be a null vector!";
    }

    auto quaternion = Quaternion<double>::FromAxisAngle(
        rRotationAxis[0],
        rRotationAxis[1],
        rRotationAxis[2],
        rotationAngle
    );

    this->LoadFromQuaternion(
        quaternion,
        rReferencePoint,
        rTranslationVector
    );

    KRATOS_CATCH("");
}


void ImposeMeshMotionProcess::LoadFromEulerAngles(const array_1d<double,3>& rEulerAngles,
                                                  const array_1d<double,3>& rReferencePoint,
                                                  const array_1d<double,3>& rTranslationVector)
{
    KRATOS_TRY;

    auto quaternion = Quaternion<double>::FromEulerAngles(rEulerAngles);

    this->LoadFromQuaternion(
        quaternion,
        rReferencePoint,
        rTranslationVector
    );

    KRATOS_CATCH("");
}


void ImposeMeshMotionProcess::LoadFromQuaternion(const Quaternion<double>& rQuaternion,
                                                 const array_1d<double,3>& rReferencePoint,
                                                 const array_1d<double,3>& rTranslationVector)
{
    KRATOS_TRY;

    mRotationMatrix.resize(3, 3);
    mReferencePoint = rReferencePoint;
    mTranslationVector = rTranslationVector;

    rQuaternion.ToRotationMatrix(mRotationMatrix);

    KRATOS_CATCH("");
}


void ImposeMeshMotionProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY;

    block_for_each(mrModelPart.Nodes(),
        [this](Node<3>& rNode) {
            array_1d<double,3> transformed_point = rNode;
            this->Transform(transformed_point);
            rNode.SetValue(MESH_DISPLACEMENT, transformed_point - rNode);
        }
    );

    KRATOS_CATCH("");
}


const Parameters ImposeMeshMotionProcess::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "model_part_name"       : "",
        "interval"              : [0.0, 1e30],
        "rotation_definition"   : "euler_angles",
        "euler_angles"          : [0.0, 0.0, 0.0],
        "rotation_axis"         : [0.0, 0.0, 1.0],
        "reference_point"       : [0.0, 0.0, 0.0],
        "rotation_angle"        : 0,
        "translation_vector"    : [0.0, 0.0, 0.0]
    })");
}


std::string ImposeMeshMotionProcess::Info() const
{
    return "ImposeMeshMotionProcess";
}


void ImposeMeshMotionProcess::PrintInfo(std::ostream& rStream) const{
    rStream << this->Info();
}


std::ostream& operator<<(std::ostream& rStream,
                         ImposeMeshMotionProcess& rThis)
{
    rThis.PrintInfo(rStream);
    rStream << std::endl;
    rThis.PrintData(rStream);
    return rStream;
}


} // namespace Kratos