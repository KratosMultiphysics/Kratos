//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <fstream>
#include <functional>

// External includes

// Project includes
#include "includes/data_communicator.h"
#include "includes/define.h"
#include "input_output/logger_output.h"
#include "utilities/brute_force_point_locator.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"

// Include base h
#include "rans_line_output_process.h"

namespace Kratos
{
// VariableDataCollector type methods

// *************** double *********************
template <>
void LineOutputProcessUtilities::VariableDataCollector<double>::AddToValuesVector(
    std::vector<double>& rValuesVector,
    const double& rValue,
    const SizeType StartIndex)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(StartIndex >= rValuesVector.size())
        << "StartIndex is greater than or equal to rValuesVector size. [ "
        << StartIndex << " >= " << rValuesVector.size() << " ].\n";

    rValuesVector[StartIndex] += rValue;

    KRATOS_CATCH("");
}

template <>
void LineOutputProcessUtilities::VariableDataCollector<double>::AddNamesToVector(
    std::vector<std::string>& rNamesVector,
    const Variable<double>& rVariable,
    const double&,
    const SizeType StartIndex)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(StartIndex >= rNamesVector.size())
        << "StartIndex is greater than or equal to rNamesVector size. [ "
        << StartIndex << " >= " << rNamesVector.size() << " ].\n";

    rNamesVector[StartIndex] = rVariable.Name();

    KRATOS_CATCH("");
}

template <>
std::size_t LineOutputProcessUtilities::VariableDataCollector<double>::GetVariableDataLength(
    const double&)
{
    return 1;
}

// *************** array_1d<double, 3> *********************
template <>
void LineOutputProcessUtilities::VariableDataCollector<array_1d<double, 3>>::AddToValuesVector(
    std::vector<double>& rValuesVector,
    const array_1d<double, 3>& rValue,
    const SizeType StartIndex)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(StartIndex >= rValuesVector.size() + 3)
        << "rValuesVector size is not enough to allocate values. rValuesVector "
           "vector should have atleast 3 elements from(including) StartIndex value. [ "
           "StartIndex = "
        << StartIndex << ", rValuesVector.size() = " << rValuesVector.size() << " ].\n";

    rValuesVector[StartIndex] += rValue[0];
    rValuesVector[StartIndex + 1] += rValue[1];
    rValuesVector[StartIndex + 2] += rValue[2];

    KRATOS_CATCH("");
}

template <>
void LineOutputProcessUtilities::VariableDataCollector<array_1d<double, 3>>::AddNamesToVector(
    std::vector<std::string>& rNamesVector,
    const Variable<array_1d<double, 3>>& rVariable,
    const array_1d<double, 3>&,
    const SizeType StartIndex)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(StartIndex >= rNamesVector.size() + 3)
        << "rNamesVector size is not enough to allocate values. rNamesVector "
           "vector should have atleast 3 elements from(including) StartIndex "
           "value. [ StartIndex = "
        << StartIndex << ", rNamesVector.size() = " << rNamesVector.size() << " ].\n";

    rNamesVector[StartIndex] = rVariable.Name() + "_X";
    rNamesVector[StartIndex + 1] = rVariable.Name() + "_Y";
    rNamesVector[StartIndex + 2] = rVariable.Name() + "_Z";

    KRATOS_CATCH("");
}

template <>
std::size_t LineOutputProcessUtilities::VariableDataCollector<array_1d<double, 3>>::GetVariableDataLength(
    const array_1d<double, 3>&)
{
    return 3;
}

// *************** Matrix *********************
template <>
void LineOutputProcessUtilities::VariableDataCollector<Matrix>::AddToValuesVector(
    std::vector<double>& rValuesVector,
    const Matrix& rValue,
    const SizeType StartIndex)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(StartIndex >=
                          rValuesVector.size() + rValue.size1() * rValue.size2())
        << "rValuesVector size is not enough to allocate values. "
           "rValuesVector vector should have atleast "
        << rValue.size1() * rValue.size2()
        << " elements from(including) StartIndex value. [ StartIndex = " << StartIndex
        << ", rValuesVector.size() = " << rValuesVector.size() << " ].\n";

    SizeType local_index = 0;
    for (std::size_t i = 0; i < rValue.size1(); ++i) {
        for (std::size_t j = 0; j < rValue.size2(); ++j) {
            rValuesVector[StartIndex + local_index++] += rValue(i, j);
        }
    }

    KRATOS_CATCH("");
}

template <>
void LineOutputProcessUtilities::VariableDataCollector<Matrix>::AddNamesToVector(
    std::vector<std::string>& rNamesVector,
    const Variable<Matrix>& rVariable,
    const Matrix& rValue,
    const SizeType StartIndex)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(StartIndex >=
                          rNamesVector.size() + rValue.size1() * rValue.size2())
        << "rNamesVector size is not enough to allocate values. "
           "rNamesVector vector should have atleast "
        << rValue.size1() * rValue.size2()
        << " elements from(including) StartIndex value. [ StartIndex = " << StartIndex
        << ", rNamesVector.size() = " << rNamesVector.size() << " ].\n";

    SizeType local_index = 0;
    for (SizeType i = 0; i < rValue.size1(); ++i) {
        for (SizeType j = 0; j < rValue.size2(); ++j) {
            rNamesVector[StartIndex + local_index++] =
                rVariable.Name() + "(" + std::to_string(i + 1) + "|" +
                std::to_string(j + 1) + ")";
        }
    }

    KRATOS_CATCH("");
}

template <>
std::size_t LineOutputProcessUtilities::VariableDataCollector<Matrix>::GetVariableDataLength(
    const Matrix& rValue)
{
    return rValue.size1() * rValue.size2();
}

RansLineOutputProcess::RansLineOutputProcess(
    Model& rModel,
    Parameters rParameters)
    : Process(), mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mModelPartName = rParameters["model_part_name"].GetString();
    mWriteHeader = rParameters["write_header_information"].GetBool();
    mUpdatePointsEachStep = rParameters["update_output_points_each_step"].GetBool();

    const Vector& r_start_point = rParameters["start_point"].GetVector();
    KRATOS_ERROR_IF(r_start_point.size() != 3)
        << "\"start_point\" should be given as 3-valued vector. [ "
           "start_point = "
        << r_start_point << " ]\n";
    const Vector& r_end_point = rParameters["end_point"].GetVector();
    KRATOS_ERROR_IF(r_end_point.size() != 3)
        << "\"end_point\" should be given as 3-valued vector. [ "
           "end_point = "
        << r_end_point << " ]\n";
    noalias(mStartPoint) = r_start_point;
    noalias(mEndPoint) = r_end_point;

    mNumberOfSamplingPoints = rParameters["number_of_sampling_points"].GetInt();
    KRATOS_ERROR_IF(mNumberOfSamplingPoints <= 2)
        << "Sampling points should be greater than 2. [ "
           "number_of_sampling_points = "
        << mNumberOfSamplingPoints << " ]\n";

    mVariableNames = rParameters["variable_names_list"].GetStringArray();
    mIsHistoricalValue = rParameters["historical_value"].GetBool();

    // output controls
    mOutputFileName = rParameters["output_file_name"].GetString();
    mOutputStepControlVariableName =
        rParameters["output_step_control_variable_name"].GetString();
    mOutputStepInterval = rParameters["output_step_interval"].GetInt();

    mPreviousStepValue = 0.0;

    KRATOS_CATCH("");
}

int RansLineOutputProcess::Check()
{
    KRATOS_TRY

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    for (const auto& variable_name : mVariableNames) {
        KRATOS_ERROR_IF(!(
            CheckAndAddVariableToList(mDoubleVariablesList, r_model_part, variable_name) ||
            CheckAndAddVariableToList(mArray3VariablesList, r_model_part, variable_name) ||
            CheckAndAddVariableToList(mArray4VariablesList, r_model_part, variable_name) ||
            CheckAndAddVariableToList(mArray6VariablesList, r_model_part, variable_name) ||
            CheckAndAddVariableToList(mArray9VariablesList, r_model_part, variable_name) ||
            CheckAndAddVariableToList(mVectorVariablesList, r_model_part, variable_name) ||
            CheckAndAddVariableToList(mMatrixVariablesList, r_model_part, variable_name)))
            << variable_name << " not found in following list of variable types: "
            << "\n\t double"
            << "\n\t Array3"
            << "\n\t Array4"
            << "\n\t Array6"
            << "\n\t Array9"
            << "\n\t Vector"
            << "\n\t Matrix.";
    }

    return 0;

    KRATOS_CATCH("");
}

void RansLineOutputProcess::ExecuteInitialize()
{
    KRATOS_TRY

    CheckAndGetOutputVariableValue<int, double>(mOutputStepControlVariableName);

    UpdateSamplePoints();

    KRATOS_CATCH("");
}

void RansLineOutputProcess::UpdateSamplePoints()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    const double coeff = 1.0 / (mNumberOfSamplingPoints - 1.0);
    const array_1d<double, 3>& r_direction_vector = (mEndPoint - mStartPoint) * coeff;

    mSamplingPoints.resize(mNumberOfSamplingPoints * 3);
    mSamplingPointElementIds.resize(mNumberOfSamplingPoints, 0);
    mSamplingPointElementShapeFunctions.resize(mNumberOfSamplingPoints);

    const auto& r_communicator = r_model_part.GetCommunicator().GetDataCommunicator();

    // calculate the sampling points in the rank zero
    if (r_communicator.Rank() == 0) {
        SizeType local_index = 0;
        for (SizeType i = 0; i < mNumberOfSamplingPoints; ++i) {
            const auto& current_point = mStartPoint + (r_direction_vector * i);
            mSamplingPoints[local_index++] = current_point[0];
            mSamplingPoints[local_index++] = current_point[1];
            mSamplingPoints[local_index++] = current_point[2];
        }
    }

    r_communicator.Broadcast(mSamplingPoints, 0);

    // MPI call to FindElement fails with OMP_NUM_THREADS > 1 in clang full
    // debug. Therefore this is done without OMP.
    BruteForcePointLocator brute_force_point_locator(r_model_part);
    for (SizeType i = 0; i < mNumberOfSamplingPoints; ++i) {
        Point current_point(mSamplingPoints[i * 3], mSamplingPoints[i * 3 + 1],
                            mSamplingPoints[i * 3 + 2]);
        mSamplingPointElementIds[i] = brute_force_point_locator.FindElement(
            current_point, mSamplingPointElementShapeFunctions[i]);
        if (mSamplingPointElementIds[i] > -1) {
            mSamplePointLocalIndexList.push_back(i);
        }
    }

    mSamplePointLocalIndexListMaster =
        r_communicator.Gatherv(mSamplePointLocalIndexList, 0);
    if (r_communicator.Rank() == 0) {
        std::vector<int> found_global_points(mNumberOfSamplingPoints, -1);
        for (const auto& current_index_list : mSamplePointLocalIndexListMaster) {
            for (const auto sample_index : current_index_list) {
                KRATOS_ERROR_IF(static_cast<SizeType>(sample_index) > mNumberOfSamplingPoints)
                    << "Sampling index error.\n";
                KRATOS_ERROR_IF(found_global_points[sample_index] != -1)
                    << "Two or more partitions found enclosed element for "
                       "sample point at [ "
                    << mSamplingPoints[sample_index * 3] << ", "
                    << mSamplingPoints[sample_index * 3 + 1] << ", "
                    << mSamplingPoints[sample_index * 3 + 2] << " ].\n";
                found_global_points[sample_index] = 1;
            }
        }
        for (SizeType sample_index = 0; sample_index < mNumberOfSamplingPoints; ++sample_index) {
            KRATOS_WARNING_IF(this->Info(), found_global_points[sample_index] != 1)
                << "Element not found for sample point at [ "
                << mSamplingPoints[sample_index * 3] << ", "
                << mSamplingPoints[sample_index * 3 + 1] << ", "
                << mSamplingPoints[sample_index * 3 + 2] << " ].\n";
        }
    }

    KRATOS_CATCH("");
}

void RansLineOutputProcess::ExecuteInitializeSolutionStep()
{
    if (mUpdatePointsEachStep)
        this->UpdateSamplePoints();
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
    const double& current_step_value =
        CheckAndGetOutputVariableValue<int, double>(mOutputStepControlVariableName);

    mCurrentStepCount += (current_step_value - mPreviousStepValue);
    mPreviousStepValue = current_step_value;

    if (mCurrentStepCount >= mOutputStepInterval) {
        mCurrentStepCount = 0.0;
        return true;
    } else {
        return false;
    }
}

void RansLineOutputProcess::WriteOutputFile() const
{
    KRATOS_TRY

    using namespace LineOutputProcessUtilities;

    const auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    const auto& r_communicator = r_model_part.GetCommunicator().GetDataCommunicator();

    const int my_pid = r_communicator.Rank();
    const SizeType number_of_local_points = mSamplePointLocalIndexList.size();

    // create getter methods
    const auto get_double = (mIsHistoricalValue) ? &GetHistoricalValue<double> : &GetNonHistoricalValue<double>;
    const auto get_array3 = (mIsHistoricalValue) ? &GetHistoricalValue<array_1d<double, 3>> : &GetNonHistoricalValue<array_1d<double, 3>>;
    const auto get_array4 = (mIsHistoricalValue) ? &GetHistoricalValue<array_1d<double, 4>> : &GetNonHistoricalValue<array_1d<double, 4>>;
    const auto get_array6 = (mIsHistoricalValue) ? &GetHistoricalValue<array_1d<double, 6>> : &GetNonHistoricalValue<array_1d<double, 6>>;
    const auto get_array9 = (mIsHistoricalValue) ? &GetHistoricalValue<array_1d<double, 9>> : &GetNonHistoricalValue<array_1d<double, 9>>;
    const auto get_vector = (mIsHistoricalValue) ? &GetHistoricalValue<Vector> : &GetNonHistoricalValue<Vector>;
    const auto get_matrix = (mIsHistoricalValue) ? &GetHistoricalValue<Matrix> : &GetNonHistoricalValue<Matrix>;

    // create and initialize variable start indices vectors
    const auto& r_node_begin = *(r_model_part.NodesBegin());
    SizeType total_doubles_length = 0;
    const IndicesVector& double_start_indices = GetVariableDataStartIndices<double>(r_node_begin, mDoubleVariablesList, get_double, total_doubles_length);
    const IndicesVector& array3_start_indices = GetVariableDataStartIndices<array_1d<double, 3>>(r_node_begin, mArray3VariablesList, get_array3, total_doubles_length);
    const IndicesVector& array4_start_indices = GetVariableDataStartIndices<array_1d<double, 4>>(r_node_begin, mArray4VariablesList, get_array4, total_doubles_length);
    const IndicesVector& array6_start_indices = GetVariableDataStartIndices<array_1d<double, 6>>(r_node_begin, mArray6VariablesList, get_array6, total_doubles_length);
    const IndicesVector& array9_start_indices = GetVariableDataStartIndices<array_1d<double, 9>>(r_node_begin, mArray9VariablesList, get_array9, total_doubles_length);
    const IndicesVector& vector_start_indices = GetVariableDataStartIndices<Vector>(r_node_begin, mVectorVariablesList, get_vector, total_doubles_length);
    const IndicesVector& matrix_start_indices = GetVariableDataStartIndices<Matrix>(r_node_begin, mMatrixVariablesList, get_matrix, total_doubles_length);

    // check sizes of dynamic variables in nodes
    const int number_of_vector_variables = mVectorVariablesList.size();
    const int number_of_matrix_variables = mMatrixVariablesList.size();
    IndexPartition<int>(r_model_part.NumberOfNodes()).for_each([&](const int iNode)
    {
        // check for vector sizes
        const auto& r_node = *(r_model_part.NodesBegin() + iNode);
        int local_index = vector_start_indices[0];
        for (int i = 0; i < number_of_vector_variables; ++i) {
            const auto& r_variable = *(mVectorVariablesList[i]);
            const auto& r_value = (*get_vector)(r_node, r_variable);
            local_index += VariableDataCollector<Vector>::GetVariableDataLength(r_value);
            KRATOS_ERROR_IF(local_index != vector_start_indices[i + 1])
                << "Size mismatch in " << r_variable.Name() << " in Node at "
                << r_node.Coordinates() << " [ value at node = " << r_value << " ].\n";
        }

        // check for matrix sizes
        local_index = matrix_start_indices[0];
        for (int i = 0; i < number_of_matrix_variables; ++i) {
            const auto& r_variable = *(mMatrixVariablesList[i]);
            const auto& r_value = (*get_matrix)(r_node, r_variable);
            local_index += VariableDataCollector<Matrix>::GetVariableDataLength(r_value);
            KRATOS_ERROR_IF(local_index != matrix_start_indices[i + 1])
                << "Size mismatch in " << r_variable.Name() << " in Node at "
                << r_node.Coordinates() << " [ value at node = " << r_value << " ].\n";
        }
    });

    // create and initialize local vectors to hold flatten variable data
    std::vector<double> local_sampled_values(number_of_local_points * total_doubles_length, 0.0);

    // calculate local sampling point values
    IndexPartition<int>(number_of_local_points).for_each([&](const int SamplingIndex) {
        const SizeType global_index = mSamplePointLocalIndexList[SamplingIndex];
        const auto& r_geometry =
            r_model_part.GetElement(mSamplingPointElementIds[global_index]).GetGeometry();
        const Vector& r_shape_functions = mSamplingPointElementShapeFunctions[global_index];
        const SizeType local_sample_point_offset = static_cast<SizeType>(SamplingIndex) * total_doubles_length;

        InterpolateVariables(
            local_sampled_values, r_geometry, r_shape_functions, local_sample_point_offset,
            std::tie(double_start_indices, mDoubleVariablesList, get_double),
            std::tie(array3_start_indices, mArray3VariablesList, get_array3),
            std::tie(array4_start_indices, mArray4VariablesList, get_array4),
            std::tie(array6_start_indices, mArray6VariablesList, get_array6),
            std::tie(array9_start_indices, mArray9VariablesList, get_array9),
            std::tie(vector_start_indices, mVectorVariablesList, get_vector),
            std::tie(matrix_start_indices, mMatrixVariablesList, get_matrix));

        KRATOS_DEBUG_ERROR_IF(mSamplingPointElementIds[global_index] == -1) << "Sampling point not found.";
    });

    // gather everything to rank zero
    auto global_sample_point_double_values = r_communicator.Gatherv(local_sampled_values, 0);

    if (my_pid == 0) {
        // put all doubles in continuos double vector
        const int number_of_processes = global_sample_point_double_values.size();
        std::vector<double> global_values(mNumberOfSamplingPoints * total_doubles_length);
        for (int rank = 0; rank < number_of_processes; ++rank) {
            const auto& local_indices = mSamplePointLocalIndexListMaster[rank];
            IndexPartition<int>(local_indices.size()).for_each([&](const int i) {
                const SizeType global_offset = local_indices[i] * total_doubles_length;
                const SizeType local_offset = static_cast<SizeType>(i) * total_doubles_length;
                for (SizeType j = 0; j < total_doubles_length; ++j) {
                    global_values[global_offset + j] =
                        global_sample_point_double_values[rank][local_offset + j];
                }
            });
        }

        // writing the output to a file
        std::ofstream output_file;
        output_file.open(GetOutputFileName());


        // writing the header
        if (mWriteHeader) {
            WriteOutputFileHeader(output_file);
        }

        output_file << "#\n#POINT_ID,POINT_X,POINT_Y,POINT_Z";

        std::vector<std::string> variable_names(total_doubles_length);
        AddVariablesListNamesToVector(variable_names, r_node_begin, mDoubleVariablesList, double_start_indices, get_double);
        AddVariablesListNamesToVector(variable_names, r_node_begin, mArray3VariablesList, array3_start_indices, get_array3);
        AddVariablesListNamesToVector(variable_names, r_node_begin, mArray4VariablesList, array4_start_indices, get_array4);
        AddVariablesListNamesToVector(variable_names, r_node_begin, mArray6VariablesList, array6_start_indices, get_array6);
        AddVariablesListNamesToVector(variable_names, r_node_begin, mArray9VariablesList, array9_start_indices, get_array9);
        AddVariablesListNamesToVector(variable_names, r_node_begin, mVectorVariablesList, vector_start_indices, get_vector);
        AddVariablesListNamesToVector(variable_names, r_node_begin, mMatrixVariablesList, matrix_start_indices, get_matrix);

        const int number_of_variables = variable_names.size();
        for (int i = 0; i < number_of_variables; ++i) {
            output_file << "," << variable_names[i];
        }
        output_file << "\n";

        // writing the sample point data
        SizeType local_index = 0;
        for (SizeType i = 0; i < mNumberOfSamplingPoints; ++i) {
            output_file << (i + 1) << "," << mSamplingPoints[i * 3] << ","
                        << mSamplingPoints[i * 3 + 1] << ","
                        << mSamplingPoints[i * 3 + 2];

            for (SizeType j = 0; j < total_doubles_length; ++j) {
                output_file << "," << global_values[local_index++];
            }
            output_file << "\n";
        }
        output_file << "# End of file";
        output_file.close();
    }

    KRATOS_CATCH("");
}

void RansLineOutputProcess::WriteOutputFileHeader(
    std::ofstream& rOutputFileStream) const
{
    KRATOS_TRY

    std::stringstream kratos_header;
    LoggerOutput current_logger_output(kratos_header);
    current_logger_output.WriteHeader();
    std::string kratos_header_str = kratos_header.str();
    std::string commented_kratos_header = "";
    for (char current_char : kratos_header_str) {
        if (current_char == '\n') {
            commented_kratos_header += "\n# ";
        } else {
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

    const std::string& output_value = std::to_string(
        CheckAndGetOutputVariableValue<int, double>(mOutputStepControlVariableName));

    rOutputFileStream << "# Output step control variable value: "
                      << output_value << "\n";
    rOutputFileStream << "# Output step frequency             : " << mOutputStepInterval
                      << "\n";
    rOutputFileStream << "# output historical values          : "
                      << (mIsHistoricalValue ? "true" : "false") << "\n";

    rOutputFileStream << "# -------------------- End of line output "
                         "settings ----------------\n";

    KRATOS_CATCH("");
}

std::string RansLineOutputProcess::GetOutputFileName() const
{
    KRATOS_TRY

    const std::string& output_value = std::to_string(
        CheckAndGetOutputVariableValue<int, double>(mOutputStepControlVariableName));

    std::stringstream output_name;
    output_name << mOutputFileName << "_" << output_value << ".csv";

    return output_name.str();

    KRATOS_CATCH("");
}

const Parameters RansLineOutputProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name"                   : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variable_names_list"               : [],
            "historical_value"                  : true,
            "start_point"                       : [0.0, 0.0, 0.0],
            "end_point"                         : [0.0, 0.0, 0.0],
            "number_of_sampling_points"         : 0,
            "output_file_name"                  : "",
            "output_step_control_variable_name" : "STEP",
            "output_step_interval"              : 1,
            "write_header_information"          : true,
            "update_output_points_each_step"    : false
        })");

    return default_parameters;
}

} // namespace Kratos.
