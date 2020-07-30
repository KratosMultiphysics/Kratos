//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes
#include <limits>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_utilities/rans_check_utilities.h"

// Include base h
#include "rans_vector_align_process.h"

namespace Kratos
{
RansVectorAlignProcess::RansVectorAlignProcess(Model& rModel, Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"         : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "input_variable_name"     : "PLEASE_SPECIFY_INPUT_VARIABLE_NAME",
            "output_variable_name"    : "PLEASE_SPECIFY_OUTPUT_VARIABLE_NAME",
            "alignment_variable_name" : "NORMAL",
            "is_tangential_alignment" : true,
            "echo_level"              : 0
        })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mModelPartName = mrParameters["model_part_name"].GetString();
    mInputVariableName = mrParameters["input_variable_name"].GetString();
    mOutputVariableName = mrParameters["output_variable_name"].GetString();
    mAlignmentVariableName = mrParameters["alignment_variable_name"].GetString();
    mIsTangentialAlignment = mrParameters["is_tangential_alignment"].GetBool();
    mEchoLevel = mrParameters["echo_level"].GetInt();

    KRATOS_CATCH("");
}

int RansVectorAlignProcess::Check()
{
    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    const Variable<array_1d<double, 3>>& input_variable =
        KratosComponents<Variable<array_1d<double, 3>>>::Get(mInputVariableName);
    const Variable<array_1d<double, 3>>& output_variable =
        KratosComponents<Variable<array_1d<double, 3>>>::Get(mOutputVariableName);
    const Variable<array_1d<double, 3>>& alignment_variable =
        KratosComponents<Variable<array_1d<double, 3>>>::Get(mAlignmentVariableName);

    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, input_variable);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, output_variable);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, alignment_variable);

    return 0.0;
}

void RansVectorAlignProcess::Execute()
{
    CalculateAlignmentVectors();
}

std::string RansVectorAlignProcess::Info() const
{
    return std::string("RansVectorAlignProcess");
}

void RansVectorAlignProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansVectorAlignProcess::PrintData(std::ostream& rOStream) const
{
}

void RansVectorAlignProcess::CalculateAlignmentVectors()
{
    KRATOS_TRY

    const Variable<array_1d<double, 3>>& r_input_variable =
        KratosComponents<Variable<array_1d<double, 3>>>::Get(mInputVariableName);
    const Variable<array_1d<double, 3>>& r_output_variable =
        KratosComponents<Variable<array_1d<double, 3>>>::Get(mOutputVariableName);
    const Variable<array_1d<double, 3>>& r_alignment_variable =
        KratosComponents<Variable<array_1d<double, 3>>>::Get(mAlignmentVariableName);

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    const int number_of_nodes = r_model_part.NumberOfNodes();

    if (mIsTangentialAlignment)
    {
#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_model_part.NodesBegin() + i_node);

            const array_1d<double, 3>& r_normal =
                r_node.FastGetSolutionStepValue(r_alignment_variable);
            const double normal_magnitude = norm_2(r_normal);

            KRATOS_ERROR_IF(normal_magnitude < std::numeric_limits<double>::epsilon())
                << "NORMAL magnitude is zero in node " << r_node.Id() << " at "
                << r_node.Coordinates() << " in " << mModelPartName
                << " [ NORMAL = " << r_normal << " ].\n";

            const array_1d<double, 3> r_unit_normal = r_normal / norm_2(r_normal);

            const array_1d<double, 3>& r_velocity =
                r_node.FastGetSolutionStepValue(r_input_variable);
            r_node.FastGetSolutionStepValue(r_output_variable) =
                r_velocity - r_unit_normal * inner_prod(r_velocity, r_unit_normal);
        }
    }
    else
    {
#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_model_part.NodesBegin() + i_node);

            const array_1d<double, 3>& r_normal =
                r_node.FastGetSolutionStepValue(r_alignment_variable);
            const double normal_magnitude = norm_2(r_normal);

            KRATOS_ERROR_IF(normal_magnitude < std::numeric_limits<double>::epsilon())
                << "NORMAL magnitude is zero in node " << r_node.Id() << " at "
                << r_node.Coordinates() << " in " << mModelPartName
                << " [ NORMAL = " << r_normal << " ].\n";

            const array_1d<double, 3> r_unit_normal = r_normal / norm_2(r_normal);

            const array_1d<double, 3>& r_velocity =
                r_node.FastGetSolutionStepValue(r_input_variable);
            r_node.FastGetSolutionStepValue(r_output_variable) =
                r_unit_normal * inner_prod(r_velocity, r_unit_normal);
        }
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << r_input_variable.Name() << " is aligned using "
        << r_alignment_variable.Name() << " and stored at "
        << r_output_variable.Name() << " for nodes in " << mModelPartName << ".\n";

    KRATOS_CATCH("");
}

} // namespace Kratos.
