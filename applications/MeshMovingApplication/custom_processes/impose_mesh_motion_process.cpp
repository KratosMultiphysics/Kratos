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
#include "custom_utilities/parametric_linear_transform.h"
#include "includes/checks.h"
#include "utilities/parallel_utilities.h"
#include "includes/mesh_moving_variables.h"


namespace Kratos
{


ImposeMeshMotionProcess::ImposeMeshMotionProcess(Model& rModel, Parameters parameters) :
    ImposeMeshMotionProcess(rModel.GetModelPart(parameters.GetValue("model_part_name").GetString()), parameters)
{
    
}


ImposeMeshMotionProcess::ImposeMeshMotionProcess(ModelPart& rModelPart, Parameters parameters) :
    Process(),
    mrModelPart(rModelPart),
    mIntervalUtility(parameters) // initialized here because IntervalUtility has no default constructor
{
    KRATOS_TRY;

    // Check model part name and validate parameters
    KRATOS_ERROR_IF(rModelPart.FullName() != parameters.GetValue("model_part_name").GetString())
    << "Expecting model part named '" << parameters.GetValue("model_part_name").GetString() << "' but got '" << rModelPart.FullName() << "' instead";

    // Process default parameters
    // 'rotation_angle' can either be numeric or an expression (string),
    // but the Parameters validator can't handle both.
    // -> default parameters must be tweaked to match the input type
    Parameters default_parameters = this->GetDefaultParameters();
    if (parameters.Has("rotation_angle") && parameters.GetValue("rotation_angle").IsString()) {
        default_parameters.RemoveValue("rotation_angle");
        default_parameters.AddValue("rotation_angle", Parameters("\"0\""));
    }

    parameters.ValidateAndAssignDefaults(default_parameters);

    // Parse general parameters
    mIntervalUtility = IntervalUtility(parameters);
    std::string rotation_definition = parameters.GetValue("rotation_definition").GetString();

    // Check whether all numeric parameters are constants
    //  - if yes -> LinearTransform
    //  - otherwise -> ParametricLinearTransform
    bool require_parametric_transform = false;

    const auto euler_angles_parameters = parameters.GetValue("euler_angles");
    const auto rotation_axis_parameters = parameters.GetValue("rotation_axis");
    const auto reference_point_parameters = parameters.GetValue("reference_point");
    const auto rotation_angle_parameters = parameters.GetValue("rotation_angle");
    const auto translation_vector_parameters = parameters.GetValue("translation_vector");

    if (rotation_angle_parameters.IsString()) {
        require_parametric_transform = true;
    }
    else {
        for (std::size_t i=0; i<3; ++i)
        {
            if (euler_angles_parameters.GetArrayItem(i).IsString()
                || rotation_axis_parameters.GetArrayItem(i).IsString()
                || reference_point_parameters.GetArrayItem(i).IsString()
                || translation_vector_parameters.GetArrayItem(i).IsString()) {
                    require_parametric_transform = true;
                    break;
            } // if any parameter is a string
        } // for i in range(3)
    } // if rotation_angle is not a string

    if (require_parametric_transform) {
        this->ParseAndSetParametricTransform(
            rotation_definition,
            euler_angles_parameters,
            rotation_axis_parameters,
            rotation_angle_parameters,
            reference_point_parameters,
            translation_vector_parameters);
    } // if transform is parametric
    else {
        this->ParseAndSetConstantTransform(
            rotation_definition,
            euler_angles_parameters,
            rotation_axis_parameters,
            rotation_angle_parameters,
            reference_point_parameters,
            translation_vector_parameters);
    } // transform is constant

    KRATOS_CATCH("");
}


void ImposeMeshMotionProcess::ParseAndSetConstantTransform(const std::string& rRotationDefinition,
                                                           const Parameters& rEulerAngles,
                                                           const Parameters& rRotationAxis,
                                                           const Parameters& rRotationAngle,
                                                           const Parameters& rReferencePoint,
                                                           const Parameters& rTranslationVector)
{
    KRATOS_TRY

    LinearTransform::Pointer p_transform;

    // Parse parameters - translation vector
    Vector vector_parameter = rTranslationVector.GetVector();
    if (vector_parameter.size() != 3) {
        KRATOS_ERROR << "'translation_vector' must be of size 3, but got" << vector_parameter.size();
    }
    array_1d<double,3> translation_vector;
    translation_vector[0] = vector_parameter[0];
    translation_vector[1] = vector_parameter[1];
    translation_vector[2] = vector_parameter[2];

    vector_parameter = rReferencePoint.GetVector();
    if (vector_parameter.size() != 3) {
        KRATOS_ERROR << "'reference_point' must be of size 3, but got" << vector_parameter.size();
    }
    array_1d<double,3> reference_point;
    reference_point[0] = vector_parameter[0];
    reference_point[1] = vector_parameter[1];
    reference_point[2] = vector_parameter[2];

    if (rRotationDefinition == "rotation_axis") {
        vector_parameter = rRotationAxis.GetVector();
        if (vector_parameter.size() != 3) {
            KRATOS_ERROR << "'rotation_axis' must be of size 3, but got" << vector_parameter.size();
        }

        array_1d<double,3> rotation_axis;
        rotation_axis[0] = vector_parameter[0];
        rotation_axis[1] = vector_parameter[1];
        rotation_axis[2] = vector_parameter[2];

        const double rotation_angle = rRotationAngle.GetDouble();

        p_transform = std::make_shared<LinearTransform>(
            rotation_axis,
            rotation_angle,
            reference_point,
            translation_vector);
    } // rotation_definition == "rotation_axis"

    else if (rRotationDefinition == "euler_angles") {
        vector_parameter = rEulerAngles.GetVector();
        if (vector_parameter.size() != 3) {
            KRATOS_ERROR << "'euler_angles' must be of size 3, but got" << vector_parameter.size();
        }

        array_1d<double,3> euler_angles;
        euler_angles[0] = vector_parameter[0];
        euler_angles[1] = vector_parameter[1];
        euler_angles[2] = vector_parameter[2];

        p_transform = std::make_shared<LinearTransform>(
            euler_angles,
            reference_point,
            translation_vector);
    } // rotation_definition == "euler_angles"

    else {
        KRATOS_ERROR << "Invalid 'rotation_definition': " << rRotationDefinition;
    } // unsupported rotation_definition

    // Set transform function
    mTransformFunctor = [p_transform, this](const Node<3>& rNode)
    {
        KRATOS_TRY
        return p_transform->Apply(rNode);
        KRATOS_CATCH("");
    };

    KRATOS_CATCH("");
}


void ImposeMeshMotionProcess::ParseAndSetParametricTransform(const std::string& rRotationDefinition,
                                                             const Parameters& rEulerAngles,
                                                             const Parameters& rRotationAxis,
                                                             const Parameters& rRotationAngle,
                                                             const Parameters& rReferencePoint,
                                                             const Parameters& rTranslationVector)
{
    KRATOS_TRY

    ParametricLinearTransform::Pointer p_transform;

        if (rRotationDefinition == "rotation_axis") {
            p_transform = std::make_shared<ParametricLinearTransform>(
                rRotationAxis,
                rRotationAngle,
                rReferencePoint,
                rTranslationVector);
        } // rotation_definition == "rotation_axis"

        else if (rRotationDefinition == "euler_angles") {
            p_transform = std::make_shared<ParametricLinearTransform>(
                rEulerAngles,
                rReferencePoint,
                rTranslationVector);
        } // rotation_definition == "euler_angles"

        else {
            KRATOS_ERROR << "Invalid 'rotation_definition': " << rRotationDefinition;
        } // unsupported rotation_definition

        mTransformFunctor = [p_transform, this](const Node<3>& rNode)
        {
            KRATOS_TRY

            const double time = mrModelPart.GetProcessInfo().GetValue(TIME);

            return p_transform->Apply(
                rNode,
                time,
                rNode.X0(),
                rNode.Y0(),
                rNode.Z0());

            KRATOS_CATCH("");
        };

    KRATOS_CATCH("");
}


void ImposeMeshMotionProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY;

    const double time = mrModelPart.GetProcessInfo().GetValue(TIME);

    if (mIntervalUtility.IsInInterval(time)) {
        block_for_each(mrModelPart.Nodes(),
            [this](Node<3>& rNode) {
                array_1d<double,3> transformed_point = mTransformFunctor(rNode);
                rNode.GetSolutionStepValue(MESH_DISPLACEMENT) = transformed_point - rNode;
            }
        ); // block_for_each
    } // if time in interval

    KRATOS_CATCH("");
}


const Parameters ImposeMeshMotionProcess::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "model_part_name"       : "",
        "interval"              : [0.0, "End"],
        "rotation_definition"   : "rotation_axis",
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