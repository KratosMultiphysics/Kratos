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
#include <fstream>

// External includes

// Project includes
#include "includes/data_communicator.h"
#include "includes/define.h"
#include "input_output/logger_output.h"
#include "utilities/brute_force_point_locator.h"

// Application includes
#include "custom_utilities/rans_check_utilities.h"

// Include base h
#include "rans_line_output_process.h"

namespace Kratos
{
RansLineOutputProcess::RansLineOutputProcess(Model& rModel, Parameters rParameters)
    : Process(), mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"                   : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variable_names_list"               : [],
            "historical_value"                  : true,
            "start_point"                       : [0.0, 0.0, 0.0],
            "end_point"                         : [0.0, 0.0, 0.0],
            "number_of_sampling_points"         : 0,
            "output_file_name"                  : "",
            "output_step_control_variable_name" : "STEP",
            "output_step_interval"              : 1
        })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mModelPartName = mrParameters["model_part_name"].GetString();

    const Vector& r_start_point = mrParameters["start_point"].GetVector();
    KRATOS_ERROR_IF(r_start_point.size() != 3)
        << "\"start_point\" should be given as 3-valued vector. [ "
           "start_point = "
        << r_start_point << " ]\n";
    mStartPoint[0] = r_start_point[0];
    mStartPoint[1] = r_start_point[1];
    mStartPoint[2] = r_start_point[2];

    const Vector& r_end_point = mrParameters["end_point"].GetVector();
    KRATOS_ERROR_IF(r_end_point.size() != 3)
        << "\"end_point\" should be given as 3-valued vector. [ "
           "end_point = "
        << r_end_point << " ]\n";
    mEndPoint[0] = r_end_point[0];
    mEndPoint[1] = r_end_point[1];
    mEndPoint[2] = r_end_point[2];

    mNumberOfSamplingPoints = mrParameters["number_of_sampling_points"].GetInt();
    KRATOS_ERROR_IF(mNumberOfSamplingPoints <= 2)
        << "Sampling points should be greater than 2. [ "
           "number_of_sampling_points = "
        << mNumberOfSamplingPoints << " ]\n";

    mVariableNames = mrParameters["variable_names_list"].GetStringArray();
    mIsHistoricalValue = mrParameters["historical_value"].GetBool();

    // output controls
    mOutputFileName = mrParameters["output_file_name"].GetString();
    mOutputStepControlVariableName =
        mrParameters["output_step_control_variable_name"].GetString();
    mOutputStepInterval = mrParameters["output_step_interval"].GetInt();

    mPreviousStepValue = 0.0;

    KRATOS_CATCH("");
}

int RansLineOutputProcess::Check()
{
    KRATOS_TRY

        RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    if (mIsHistoricalValue)
    {
        for (std::string variable_name : mVariableNames)
        {
            if (KratosComponents<Variable<double>>::Has(variable_name))
            {
                const Variable<double>& r_variable =
                    KratosComponents<Variable<double>>::Get(variable_name);
                RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, r_variable);
            }
            else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name))
            {
                const Variable<array_1d<double, 3>>& r_variable =
                    KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_name);
                RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, r_variable);
            }
            else
            {
                KRATOS_ERROR
                    << "Variable name = " << variable_name
                    << " not found in double, 3d-double variable lists.\n";
            }
        }
    }

    const ProcessInfo& r_process_info = r_model_part.GetProcessInfo();

    if (KratosComponents<Variable<int>>::Has(mOutputStepControlVariableName))
    {
        const Variable<int>& r_variable =
            KratosComponents<Variable<int>>::Get(mOutputStepControlVariableName);
        KRATOS_ERROR_IF(!r_process_info.Has(r_variable))
            << "\"output_step_control_variable_name\" ( " << mOutputStepControlVariableName
            << " ) not found in process info of " << r_model_part.Name() << ".\n";
    }
    else if (KratosComponents<Variable<double>>::Has(mOutputStepControlVariableName))
    {
        const Variable<double>& r_variable =
            KratosComponents<Variable<double>>::Get(mOutputStepControlVariableName);
        KRATOS_ERROR_IF(!r_process_info.Has(r_variable))
            << "\"output_step_control_variable_name\" ( " << mOutputStepControlVariableName
            << " ) not found in process info of " << r_model_part.Name() << ".\n";
    }
    else if (KratosComponents<Variable<std::string>>::Has(mOutputStepControlVariableName))
    {
        const Variable<std::string>& r_variable =
            KratosComponents<Variable<std::string>>::Get(mOutputStepControlVariableName);
        KRATOS_ERROR_IF(!r_process_info.Has(r_variable))
            << "\"output_step_control_variable_name\" ( " << mOutputStepControlVariableName
            << " ) not found in process info of " << r_model_part.Name() << ".\n";
    }

    return 0;

    KRATOS_CATCH("");
}

void RansLineOutputProcess::ExecuteInitialize()
{
    KRATOS_TRY

    const array_1d<double, 3>& r_direction_vector = mEndPoint - mStartPoint;

    mSamplingPoints.resize(mNumberOfSamplingPoints * 3);
    mSamplingPointElementIds.resize(mNumberOfSamplingPoints, 0);
    mSamplingPointElementShapeFunctions.resize(mNumberOfSamplingPoints);

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
    const DataCommunicator& r_communicator =
        r_model_part.GetCommunicator().GetDataCommunicator();

    // calculate the sampling points in the rank zero
    if (r_communicator.Rank() == 0)
    {
        for (int i = 0; i < mNumberOfSamplingPoints; ++i)
        {
            const array_1d<double, 3> current_point =
                mStartPoint + r_direction_vector *
                                  (static_cast<double>(i) /
                                   static_cast<double>(mNumberOfSamplingPoints - 1));
            mSamplingPoints[i * 3] = current_point[0];
            mSamplingPoints[i * 3 + 1] = current_point[1];
            mSamplingPoints[i * 3 + 2] = current_point[2];
        }
    }

    r_communicator.Broadcast(mSamplingPoints, 0);

    BruteForcePointLocator brute_force_point_locator(r_model_part);
    for (int i = 0; i < mNumberOfSamplingPoints; ++i)
    {
        Point current_point(mSamplingPoints[i * 3], mSamplingPoints[i * 3 + 1],
                            mSamplingPoints[i * 3 + 2]);
        mSamplingPointElementIds[i] = brute_force_point_locator.FindElement(
            current_point, mSamplingPointElementShapeFunctions[i]);
        if (mSamplingPointElementIds[i] > -1)
        {
            mSamplePointLocalIndexList.push_back(i);
        }
    }

    mSamplePointLocalIndexListMaster =
        r_communicator.Gatherv(mSamplePointLocalIndexList, 0);
    if (r_communicator.Rank() == 0)
    {
        mFoundGlobalPoints.resize(mNumberOfSamplingPoints, -1);
        for (std::vector<int> current_index_list : mSamplePointLocalIndexListMaster)
        {
            for (int sample_index : current_index_list)
            {
                KRATOS_ERROR_IF(sample_index > mNumberOfSamplingPoints)
                    << "Sampling index error.\n";
                KRATOS_ERROR_IF(mFoundGlobalPoints[sample_index] != -1)
                    << "Two or more partitions found enclosed element for "
                       "sample point at [ "
                    << mSamplingPoints[sample_index * 3] << ", "
                    << mSamplingPoints[sample_index * 3 + 1] << ", "
                    << mSamplingPoints[sample_index * 3 + 2] << " ].\n";
                mFoundGlobalPoints[sample_index] = 1;
            }
        }
        for (int sample_index = 0; sample_index < mNumberOfSamplingPoints; ++sample_index)
        {
            KRATOS_WARNING_IF(this->Info(), mFoundGlobalPoints[sample_index] != 1)
                << "Element not found for sample point at [ "
                << mSamplingPoints[sample_index * 3] << ", "
                << mSamplingPoints[sample_index * 3 + 1] << ", "
                << mSamplingPoints[sample_index * 3 + 2] << " ].\n";
        }
    }

    KRATOS_CATCH("");
}

void RansLineOutputProcess::Execute()
{
    if (this->IsOutputStep())
        this->WriteOutputFile();
}

void RansLineOutputProcess::ExecuteFinalizeSolutionStep()
{
    if (this->IsOutputStep())
        this->WriteOutputFile();
}

std::string RansLineOutputProcess::Info() const
{
    return std::string("RansLineOutputProcess");
}

void RansLineOutputProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansLineOutputProcess::PrintData(std::ostream& rOStream) const
{
}

bool RansLineOutputProcess::IsOutputStep()
{
    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
    const ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    double current_step_value = 0.0;

    if (KratosComponents<Variable<int>>::Has(mOutputStepControlVariableName))
    {
        const Variable<int>& r_variable =
            KratosComponents<Variable<int>>::Get(mOutputStepControlVariableName);
        current_step_value = static_cast<double>(r_process_info[r_variable]);
    }
    else if (KratosComponents<Variable<double>>::Has(mOutputStepControlVariableName))
    {
        const Variable<double>& r_variable =
            KratosComponents<Variable<double>>::Get(mOutputStepControlVariableName);
        current_step_value = r_process_info[r_variable];
    }
    else
    {
        return true;
    }

    mCurrentStepCount += (current_step_value - mPreviousStepValue);
    mPreviousStepValue = current_step_value;

    if (mCurrentStepCount >= mOutputStepInterval)
    {
        mCurrentStepCount = 0.0;
        return true;
    }
    else
    {
        return false;
    }
}

void RansLineOutputProcess::WriteOutputFile()
{
    KRATOS_TRY

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
    const DataCommunicator& r_communicator =
        r_model_part.GetCommunicator().GetDataCommunicator();

    const int my_pid = r_communicator.Rank();
    const int number_of_local_points = mSamplePointLocalIndexList.size();

    std::vector<std::vector<double>> global_sample_point_variables_value_list;
    std::vector<std::string> global_variable_names_list;

    for (std::string variable_name : mVariableNames)
    {
        std::vector<double> sample_point_values_list;
        int variable_length = 0;
        if (KratosComponents<Variable<double>>::Has(variable_name))
        {
            sample_point_values_list.resize(number_of_local_points);
            variable_length = 1;
            global_variable_names_list.push_back(variable_name);

            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(variable_name);

            for (int i = 0; i < number_of_local_points; ++i)
            {
                const double sample_point_value =
                    InterpolateVariable(r_variable, mSamplePointLocalIndexList[i]);
                sample_point_values_list[i] = sample_point_value;
            }
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name))
        {
            sample_point_values_list.resize(number_of_local_points * 3);
            variable_length = 3;
            global_variable_names_list.push_back(variable_name + "_X");
            global_variable_names_list.push_back(variable_name + "_Y");
            global_variable_names_list.push_back(variable_name + "_Z");
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_name);

            for (int i = 0; i < number_of_local_points; ++i)
            {
                const array_1d<double, 3>& sample_point_value =
                    InterpolateVariable(r_variable, mSamplePointLocalIndexList[i]);
                sample_point_values_list[i * 3] = sample_point_value[0];
                sample_point_values_list[i * 3 + 1] = sample_point_value[1];
                sample_point_values_list[i * 3 + 2] = sample_point_value[2];
            }
        }
        else
        {
            KRATOS_ERROR << "Variable name = " << variable_name
                         << " not found in double, 3d-double variable lists.\n";
        }

        auto global_sample_point_values_list =
            r_communicator.Gatherv(sample_point_values_list, 0);

        if (my_pid == 0)
        {
            for (int i_variable_index = 0; i_variable_index < variable_length; ++i_variable_index)
            {
                std::vector<double> current_variable_values(mNumberOfSamplingPoints);
                for (int i_process = 0;
                     i_process <
                     static_cast<int>(global_sample_point_values_list.size());
                     ++i_process)
                {
                    const std::vector<double>& local_sample_point_values_list =
                        global_sample_point_values_list[i_process];
                    const std::vector<int>& local_sample_point_indices =
                        mSamplePointLocalIndexListMaster[i_process];
                    const int number_of_local_sample_points =
                        local_sample_point_indices.size();
                    for (int i_local_index = 0;
                         i_local_index < number_of_local_sample_points; ++i_local_index)
                    {
                        current_variable_values[local_sample_point_indices[i_local_index]] =
                            local_sample_point_values_list[i_local_index * variable_length + i_variable_index];
                    }
                }
                global_sample_point_variables_value_list.push_back(current_variable_values);
            }

            // writing the output to a file
            std::ofstream output_file;
            output_file.open(GetOutputFileName());

            // writing the header
            WriteOutputFileHeader(output_file);
            output_file << "\n#POINT_ID,POINT_X,POINT_Y,POINT_Z";

            const int number_of_variables =
                static_cast<int>(global_sample_point_variables_value_list.size());
            for (int i_variable_index = 0; i_variable_index < number_of_variables; ++i_variable_index)
            {
                output_file << "," << global_variable_names_list[i_variable_index];
            }
            output_file << "\n";

            // writing the sample point data
            for (int i_point_index = 0; i_point_index < mNumberOfSamplingPoints; ++i_point_index)
            {
                if (mFoundGlobalPoints[i_point_index] == 1)
                {
                    output_file << (i_point_index + 1) << ","
                                << mSamplingPoints[i_point_index * 3] << ","
                                << mSamplingPoints[i_point_index * 3 + 1] << ","
                                << mSamplingPoints[i_point_index * 3 + 2];

                    for (int i_variable_index = 0;
                         i_variable_index < number_of_variables; ++i_variable_index)
                    {
                        output_file
                            << ","
                            << global_sample_point_variables_value_list[i_variable_index][i_point_index];
                    }
                    output_file << "\n";
                }
            }
            output_file << "# End of file";
            output_file.close();
        }
    }

    KRATOS_CATCH("");
}

void RansLineOutputProcess::WriteOutputFileHeader(std::ofstream& rOutputFileStream) const
{
    std::stringstream kratos_header;
    LoggerOutput current_logger_output(kratos_header);
    current_logger_output.WriteHeader();
    std::string kratos_header_str = kratos_header.str();
    std::string commented_kratos_header = "";
    for (char current_char : kratos_header_str)
    {
        if (current_char == '\n')
        {
            commented_kratos_header += "\n# ";
        }
        else
        {
            commented_kratos_header += current_char;
        }
    }

    rOutputFileStream << "# "
                         "-------------------------------------------------"
                         "-----------------\n# ";
    rOutputFileStream << commented_kratos_header;
    rOutputFileStream << "\n# ------------------ Summary of the line "
                         "settings ------------------\n";
    rOutputFileStream << "# Model part name                   : " << mModelPartName << "\n";
    rOutputFileStream << "# Line start location               : " << mStartPoint[0]
                      << ", " << mStartPoint[1] << ", " << mStartPoint[2] << "\n";
    rOutputFileStream << "# Line end location                 : " << mEndPoint[0]
                      << ", " << mEndPoint[1] << ", " << mEndPoint[2] << "\n";
    rOutputFileStream << "# Number of sampling points         : " << mNumberOfSamplingPoints
                      << "\n";
    rOutputFileStream << "# Output step control variable name : " << mOutputStepControlVariableName
                      << "\n";

    const ProcessInfo& r_process_info =
        mrModel.GetModelPart(mModelPartName).GetProcessInfo();

    std::stringstream output_step_control_variable_value;
    if (KratosComponents<Variable<int>>::Has(mOutputStepControlVariableName))
    {
        const Variable<int>& r_variable =
            KratosComponents<Variable<int>>::Get(mOutputStepControlVariableName);

        output_step_control_variable_value << r_process_info[r_variable];
    }
    else if (KratosComponents<Variable<double>>::Has(mOutputStepControlVariableName))
    {
        const Variable<double>& r_variable =
            KratosComponents<Variable<double>>::Get(mOutputStepControlVariableName);
        output_step_control_variable_value << r_process_info[r_variable];
    }
    else if (KratosComponents<Variable<std::string>>::Has(mOutputStepControlVariableName))
    {
        const Variable<std::string>& r_variable =
            KratosComponents<Variable<std::string>>::Get(mOutputStepControlVariableName);
        output_step_control_variable_value << r_process_info[r_variable];
    }
    else
    {
        KRATOS_ERROR << "Output step control variable name not found in "
                        "int, double or string variables list. [ "
                        "output_step_control_variable_name = "
                     << mOutputStepControlVariableName << " ].\n";
    }

    rOutputFileStream << "# Output step control variable value: "
                      << output_step_control_variable_value.str() << "\n";
    rOutputFileStream << "# Output step frequency             : " << mOutputStepInterval
                      << "\n";
    rOutputFileStream << "# output historical values          : "
                      << (mIsHistoricalValue ? "true" : "false") << "\n";

    rOutputFileStream << "# -------------------- End of line output "
                         "settings ----------------\n";
}

double RansLineOutputProcess::InterpolateVariable(const Variable<double>& rVariable,
                                                  const int SamplingIndex) const
{
    double value = 0.0;

    if (mSamplingPointElementIds[SamplingIndex] > -1)
    {
        const ModelPart& r_modelpart = mrModel.GetModelPart(mModelPartName);
        const ModelPart::ElementType::GeometryType& r_geometry =
            r_modelpart.GetElement(mSamplingPointElementIds[SamplingIndex]).GetGeometry();
        const Vector& r_shape_functions = mSamplingPointElementShapeFunctions[SamplingIndex];

        for (int i = 0; i < static_cast<int>(r_geometry.PointsNumber()); ++i)
        {
            if (mIsHistoricalValue)
            {
                value += r_geometry[i].FastGetSolutionStepValue(rVariable) *
                         r_shape_functions[i];
            }
            else
            {
                value += r_geometry[i].GetValue(rVariable) * r_shape_functions[i];
            }
        }
    }

    return value;
}

array_1d<double, 3> RansLineOutputProcess::InterpolateVariable(
    const Variable<array_1d<double, 3>>& rVariable, const int SamplingIndex) const
{
    array_1d<double, 3> value;
    value.clear();

    if (mSamplingPointElementIds[SamplingIndex] > -1)
    {
        const ModelPart& r_modelpart = mrModel.GetModelPart(mModelPartName);
        const ModelPart::ElementType::GeometryType& r_geometry =
            r_modelpart.GetElement(mSamplingPointElementIds[SamplingIndex]).GetGeometry();
        const Vector& r_shape_functions = mSamplingPointElementShapeFunctions[SamplingIndex];

        for (int i = 0; i < static_cast<int>(r_geometry.PointsNumber()); ++i)
        {
            if (mIsHistoricalValue)
            {
                value += r_geometry[i].FastGetSolutionStepValue(rVariable) *
                         r_shape_functions[i];
            }
            else
            {
                value += r_geometry[i].GetValue(rVariable) * r_shape_functions[i];
            }
        }
    }

    return value;
}

std::string RansLineOutputProcess::GetOutputFileName() const
{
    const ProcessInfo& r_process_info =
        mrModel.GetModelPart(mModelPartName).GetProcessInfo();

    if (KratosComponents<Variable<int>>::Has(mOutputStepControlVariableName))
    {
        const Variable<int>& r_variable =
            KratosComponents<Variable<int>>::Get(mOutputStepControlVariableName);
        std::stringstream output_name;
        output_name << mOutputFileName << "_" << r_process_info[r_variable] << ".csv";
        return output_name.str();
    }
    else if (KratosComponents<Variable<double>>::Has(mOutputStepControlVariableName))
    {
        const Variable<double>& r_variable =
            KratosComponents<Variable<double>>::Get(mOutputStepControlVariableName);
        std::stringstream output_name;
        output_name << mOutputFileName << "_" << r_process_info[r_variable] << ".csv";
        return output_name.str();
    }
    else if (KratosComponents<Variable<std::string>>::Has(mOutputStepControlVariableName))
    {
        const Variable<std::string>& r_variable =
            KratosComponents<Variable<std::string>>::Get(mOutputStepControlVariableName);
        std::stringstream output_name;
        output_name << mOutputFileName << "_" << r_process_info[r_variable] << ".csv";
        return output_name.str();
    }
    else
    {
        KRATOS_ERROR << "Output step control variable name not found in "
                        "int, double or string variables list. [ "
                        "output_step_control_variable_name = "
                     << mOutputStepControlVariableName << " ].\n";
    }

    return "";
}

} // namespace Kratos.
