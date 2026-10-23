//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "processes/json_output_process.h"
#include "includes/variables.h"
#include "containers/model.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

namespace Kratos
{
JsonOutputProcess::JsonOutputProcess(
    Model& rModel,
    Parameters rSettings
    ) : Process(),
        mSettings(rSettings)
{
    // Validate and assign defaults
    Parameters default_parameters = GetDefaultParameters();
    mSettings.ValidateAndAssignDefaults(default_parameters);

    // We get the submodelpart
    const std::string& r_model_part_name = mSettings["model_part_name"].GetString();
    const std::string& r_sub_model_part_name = mSettings["sub_model_part_name"].GetString();
    if (r_sub_model_part_name != "") {
        mpSubModelPart = &(rModel.GetModelPart(r_model_part_name).GetSubModelPart(r_sub_model_part_name));
    } else {
        mpSubModelPart = &(rModel.GetModelPart(r_model_part_name));
    }
}

/***********************************************************************************/
/***********************************************************************************/

JsonOutputProcess::JsonOutputProcess(
    ModelPart& rModelPart,
    Parameters rSettings
    ) : Process(),
        mpSubModelPart(&rModelPart),
        mSettings(rSettings)
{
    // Validate and assign defaults
    Parameters default_parameters = GetDefaultParameters();
    mSettings.ValidateAndAssignDefaults(default_parameters);
}

/***********************************************************************************/
/***********************************************************************************/

void JsonOutputProcess::ExecuteInitialize()
{
    // Check if we are in MPI
    KRATOS_ERROR_IF(mpSubModelPart->GetCommunicator().GetDataCommunicator().Size() > 1) << "This process cannot be used for writing output in MPI!" << std::endl;

    // If we consider any flag
    const std::string& r_flag_name = mSettings["check_for_flag"].GetString();
    if (r_flag_name != "") {
        mpFlag = &KratosComponents<Flags>::Get(r_flag_name);
    }

    // Assign some member variables
    mOutputFileName = mSettings["output_file_name"].GetString();
    mFrequency = mSettings["time_frequency"].GetDouble();
    mResultantSolution = mSettings["resultant_solution"].GetBool();
    mHistoricalValue = mSettings["historical_value"].GetBool();
    mUseNodeCoordinates = mSettings["use_node_coordinates"].GetBool();

    // The variables to output
    ParseVariables(mSettings["output_variables"], mOutputVariables, mOutputArrayVariables, mOutputVectorVariables);
    ParseVariables(mSettings["gauss_points_output_variables"], mGaussPointsOutputVariables, mGaussPointsOutputArrayVariables, mGaussPointsOutputVectorVariables);
}

/***********************************************************************************/
/***********************************************************************************/

void JsonOutputProcess::ExecuteBeforeSolutionLoop()
{
    InitializeJson();
}

/***********************************************************************************/
/***********************************************************************************/

void JsonOutputProcess::ExecuteFinalizeSolutionStep()
{
    // We write the json file if the time counter is greater than the frequency
    const double dt = mpSubModelPart->GetProcessInfo().GetValue(DELTA_TIME);
    mTimeCounter += dt;
    if (mTimeCounter > mFrequency) {
        mTimeCounter = 0.0;
        WriteJson();
    }
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters JsonOutputProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "help"                          : "This process generates a json file containing the solution of a list of variables from a given submodelpart",
        "output_variables"              : [],
        "gauss_points_output_variables" : [],
        "output_file_name"              : "",
        "model_part_name"               : "",
        "sub_model_part_name"           : "",
        "check_for_flag"                : "",
        "time_frequency"                : 1.00,
        "historical_value"              : true,
        "resultant_solution"            : false,
        "use_node_coordinates"          : false
    })" );
    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

bool JsonOutputProcess::CheckFlag(const Kratos::shared_ptr<const Geometry<Node>>& pGeometry)
{
    if (mpFlag != nullptr) {
        for (auto& r_node : *pGeometry) {
            if (r_node.IsNot(*mpFlag)) {
                return false;
            }
        }
    }
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

bool JsonOutputProcess::CheckFlag(const Node& rNode)
{
    if (mpFlag != nullptr) {
        if (rNode.IsNot(*mpFlag)) {
            return false;
        }
    }
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

std::string JsonOutputProcess::GetNodeIdentifier(const Node& rNode)
{
    if (mUseNodeCoordinates) {
        return "X_" + std::to_string(rNode.X0()) + "_Y_" + std::to_string(rNode.Y0()) + "_Z_" + std::to_string(rNode.Z0());
    } else {
        return std::to_string(rNode.Id());
    }
}

/***********************************************************************************/
/***********************************************************************************/

void JsonOutputProcess::ParseVariables(
    const Parameters rOutputVariables,
    std::vector<const Variable<double>*>& rDoubleVariables,
    std::vector<const Variable<array_1d<double, 3>>*>& rArray1dVariables,
    std::vector<const Variable<Vector>*>& rVectorVariables
    )
{
    for (unsigned int i = 0; i < rOutputVariables.size(); ++i) {
        const std::string& r_variable_name = rOutputVariables[i].GetString();
        if (KratosComponents<Variable<double>>::Has(r_variable_name)) {
            const auto* p_variable = &KratosComponents<Variable<double>>::Get(r_variable_name);
            rDoubleVariables.push_back(p_variable);
        } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable_name)) {
            const auto* p_variable = &KratosComponents<Variable<array_1d<double, 3>>>::Get(r_variable_name);
            rArray1dVariables.push_back(p_variable);
        } else if (KratosComponents<Variable<Vector>>::Has(r_variable_name)) {
            const auto* p_variable = &KratosComponents<Variable<Vector>>::Get(r_variable_name);
            rVectorVariables.push_back(p_variable);
        } else {
            KRATOS_WARNING("JsonOutputProcess") << "Variable " << r_variable_name << " is not a double, array_1d or vector variable, skipping it" << std::endl;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void JsonOutputProcess::InitializeJson()
{
    // Initialize the json file and add the time array
    mJsonFile = Parameters(R"({})");
    mJsonFile.AddEmptyArray("TIME");

    // Nodal values
    if (mOutputVariables.size() + mOutputArrayVariables.size() + mOutputVectorVariables.size() > 0) {
        int count = 0;
        for (auto& r_node : mpSubModelPart->Nodes()) {
            if (CheckFlag(r_node)) {
                const std::string node_identifier = "NODE_" + GetNodeIdentifier(r_node);
                if (!mResultantSolution) {
                    auto json_node = mJsonFile.AddEmptyValue(node_identifier);
                    for (const auto* p_variable : mOutputVariables) {
                        json_node.AddEmptyArray(p_variable->Name());
                    }
                    for (const auto* p_variable : mOutputArrayVariables) {
                        const auto& r_variable_name = p_variable->Name();
                        json_node.AddEmptyArray(r_variable_name + "_X");
                        json_node.AddEmptyArray(r_variable_name + "_Y");
                        json_node.AddEmptyArray(r_variable_name + "_Z");
                    }
                    for (const auto* p_variable : mOutputVectorVariables) {
                        json_node.AddEmptyArray(p_variable->Name());
                    }
                } else {
                    if (count == 0) {
                        auto json_node = mJsonFile.AddEmptyValue("RESULTANT");
                        for (const auto* p_variable : mOutputVariables) {
                            json_node.AddEmptyArray(p_variable->Name());
                        }
                        for (const auto* p_variable : mOutputArrayVariables) {
                            const auto& r_variable_name = p_variable->Name();
                            json_node.AddEmptyArray(r_variable_name + "_X");
                            json_node.AddEmptyArray(r_variable_name + "_Y");
                            json_node.AddEmptyArray(r_variable_name + "_Z");
                        }
                        for (const auto* p_variable : mOutputVectorVariables) {
                            json_node.AddEmptyArray(p_variable->Name());
                        }
                    }
                }
                ++count;
            }
        }
    }

    // Gauss points values
    if (mGaussPointsOutputVariables.size() + mGaussPointsOutputArrayVariables.size() + mGaussPointsOutputVectorVariables.size() > 0) {
        const auto& r_process_info = mpSubModelPart->GetProcessInfo();
        if (!mResultantSolution) {
            for (auto& r_elem : mpSubModelPart->Elements()) {
                if (CheckFlag(r_elem.pGetGeometry())) {
                    const std::string element_identifier = "ELEMENT_" + std::to_string(r_elem.Id());
                    mJsonFile.AddEmptyValue(element_identifier);
                    auto json_element = mJsonFile[element_identifier];
                    for (const auto* p_variable : mGaussPointsOutputVariables) {
                        std::vector<double> values;
                        r_elem.CalculateOnIntegrationPoints(*p_variable, values, r_process_info);
                        const std::string& r_variable_name = p_variable->Name();
                        json_element.AddEmptyValue(r_variable_name);
                        auto json_variable = json_element[r_variable_name];
                        for (unsigned int i = 0; i < values.size(); ++i) {
                            json_variable.AddEmptyArray(std::to_string(i));
                        }
                    }
                    for (const auto* p_variable : mGaussPointsOutputArrayVariables) {
                        std::vector<array_1d<double, 3>> values;
                        r_elem.CalculateOnIntegrationPoints(*p_variable, values, r_process_info);
                        const std::string& r_variable_name = p_variable->Name();
                        json_element.AddEmptyValue(r_variable_name + "_X");
                        json_element.AddEmptyValue(r_variable_name + "_Y");
                        json_element.AddEmptyValue(r_variable_name + "_Z");
                        auto json_variable_x = json_element[r_variable_name + "_X"];
                        auto json_variable_y = json_element[r_variable_name + "_Y"];
                        auto json_variable_z = json_element[r_variable_name + "_Z"];
                        for (unsigned int i = 0; i < values.size(); ++i) {
                            const std::string i_string = std::to_string(i);
                            json_variable_x.AddEmptyArray(i_string);
                            json_variable_y.AddEmptyArray(i_string);
                            json_variable_z.AddEmptyArray(i_string);
                        }
                    }
                    for (const auto* p_variable : mGaussPointsOutputVectorVariables) {
                        std::vector<Vector> values;
                        r_elem.CalculateOnIntegrationPoints(*p_variable, values, r_process_info);
                        const std::string& r_variable_name = p_variable->Name();
                        json_element.AddEmptyValue(r_variable_name);
                        auto json_variable = json_element[r_variable_name];
                        for (unsigned int i = 0; i < values.size(); ++i) {
                            json_variable.AddEmptyArray(std::to_string(i));
                        }
                    }
                }
            }
        } else {
            mJsonFile.AddEmptyValue("RESULTANT");
            auto json_element = mJsonFile["RESULTANT"];
            auto it_elem_begin = mpSubModelPart->ElementsBegin();
            for (const auto* p_variable : mGaussPointsOutputVariables) {
                std::vector<double> values;
                it_elem_begin->CalculateOnIntegrationPoints(*p_variable, values, r_process_info);
                const std::string& r_variable_name = p_variable->Name();
                json_element.AddEmptyValue(r_variable_name);
                auto json_variable = json_element[r_variable_name];
                for (unsigned int i = 0; i < values.size(); ++i) {
                    json_variable.AddEmptyArray(std::to_string(i));
                }
            }
            for (const auto* p_variable : mGaussPointsOutputArrayVariables) {
                std::vector<array_1d<double, 3>> values;
                it_elem_begin->CalculateOnIntegrationPoints(*p_variable, values, r_process_info);
                const std::string& r_variable_name = p_variable->Name();
                json_element.AddEmptyValue(r_variable_name + "_X");
                json_element.AddEmptyValue(r_variable_name + "_Y");
                json_element.AddEmptyValue(r_variable_name + "_Z");
                auto json_variable_x = json_element[r_variable_name + "_X"];
                auto json_variable_y = json_element[r_variable_name + "_Y"];
                auto json_variable_z = json_element[r_variable_name + "_Z"];
                for (unsigned int i = 0; i < values.size(); ++i) {
                    const std::string i_string = std::to_string(i);
                    json_variable_x.AddEmptyArray(i_string);
                    json_variable_y.AddEmptyArray(i_string);
                    json_variable_z.AddEmptyArray(i_string);
                }
            }
            for (const auto* p_variable : mGaussPointsOutputVectorVariables) {
                std::vector<Vector> values;
                it_elem_begin->CalculateOnIntegrationPoints(*p_variable, values, r_process_info);
                const std::string& r_variable_name = p_variable->Name();
                json_element.AddEmptyValue(r_variable_name);
                auto json_variable = json_element[r_variable_name];
                for (unsigned int i = 0; i < values.size(); ++i) {
                    json_variable.AddEmptyArray(std::to_string(i));
                }
            }
        }
    }

    // Write the file
    std::ofstream output_file(mOutputFileName);
    output_file << mJsonFile.PrettyPrintJsonString();
}

/***********************************************************************************/
/***********************************************************************************/

void JsonOutputProcess::WriteJson()
{
    // We write the time
    const auto& r_process_info = mpSubModelPart->GetProcessInfo();
    const double time = r_process_info.GetValue(TIME);
    mJsonFile["TIME"].Append(time);

    // Nodal values
    if (mOutputVariables.size() + mOutputArrayVariables.size() + mOutputVectorVariables.size() > 0) {
        if (!mResultantSolution) {
            for (auto& r_node : mpSubModelPart->Nodes()) {
                if (CheckFlag(r_node)) {
                    const std::string node_identifier = "NODE_" + GetNodeIdentifier(r_node);
                    auto json_node = mJsonFile[node_identifier];
                    for (const auto* p_variable : mOutputVariables) {
                        const double value = mHistoricalValue ? r_node.GetSolutionStepValue(*p_variable, 0) : r_node.GetValue(*p_variable);
                        json_node[p_variable->Name()].Append(value);
                    }
                    for (const auto* p_variable : mOutputArrayVariables) {
                        const array_1d<double, 3>& r_value = mHistoricalValue ? r_node.GetSolutionStepValue(*p_variable, 0) : r_node.GetValue(*p_variable);
                        const auto& r_variable_name = p_variable->Name();
                        json_node[r_variable_name + "_X"].Append(r_value[0]);
                        json_node[r_variable_name + "_Y"].Append(r_value[1]);
                        json_node[r_variable_name + "_Z"].Append(r_value[2]);
                    }
                    for (const auto* p_variable : mOutputVectorVariables) {
                        const Vector& r_value = mHistoricalValue ? r_node.GetSolutionStepValue(*p_variable, 0) : r_node.GetValue(*p_variable);
                        json_node[p_variable->Name()].Append(r_value);
                    }
                }
            }
        } else {
            auto json_node = mJsonFile["RESULTANT"];
            auto& r_nodes_array = mpSubModelPart->Nodes();
            for (const auto* p_variable : mOutputVariables) {
                const double resultant_value = block_for_each<SumReduction<double>>(r_nodes_array,[&](Node& r_node){
                    if (CheckFlag(r_node)) {
                        return mHistoricalValue ? r_node.GetSolutionStepValue(*p_variable, 0) : r_node.GetValue(*p_variable);
                    } else {
                        return 0.0;
                    }
                });
                json_node[p_variable->Name()].Append(resultant_value);
            }
            for (const auto* p_variable : mOutputArrayVariables) {
                const array_1d<double, 3> zero_vector = ZeroVector(3);
                const array_1d<double, 3> resultant_value = block_for_each<SumReduction<array_1d<double, 3>>>(r_nodes_array,[&](Node& r_node){
                    if (CheckFlag(r_node)) {
                        return mHistoricalValue ? r_node.GetSolutionStepValue(*p_variable, 0) : r_node.GetValue(*p_variable);
                    } else {
                        return zero_vector;
                    }
                });
                const auto& r_variable_name = p_variable->Name();
                json_node[r_variable_name + "_X"].Append(resultant_value[0]);
                json_node[r_variable_name + "_Y"].Append(resultant_value[1]);
                json_node[r_variable_name + "_Z"].Append(resultant_value[2]);
            }
            for (const auto* p_variable : mOutputVectorVariables) {
                auto it_node_begin = mpSubModelPart->NodesBegin();
                const std::size_t vector_size = mHistoricalValue ? it_node_begin->GetSolutionStepValue(*p_variable, 0).size() : it_node_begin->GetValue(*p_variable).size();
                Vector resultant_value = ZeroVector(vector_size);
                for (auto& r_node : mpSubModelPart->Nodes()) {
                    if (CheckFlag(r_node)) {
                        resultant_value += mHistoricalValue ? r_node.GetSolutionStepValue(*p_variable, 0) : r_node.GetValue(*p_variable);
                    }
                }
                json_node[p_variable->Name()].Append(resultant_value);
            }
        }
    }

    // Gauss points values
    if (mGaussPointsOutputVariables.size() + mGaussPointsOutputArrayVariables.size() + mGaussPointsOutputVectorVariables.size() > 0) {
        if (!mResultantSolution) {
            for (auto& r_elem : mpSubModelPart->Elements()) {
                if (CheckFlag(r_elem.pGetGeometry())) {
                    const std::string element_identifier = "ELEMENT_" + std::to_string(r_elem.Id());
                    auto json_element = mJsonFile[element_identifier];
                    for (const auto* p_variable : mGaussPointsOutputVariables) {
                        std::vector<double> values;
                        r_elem.CalculateOnIntegrationPoints(*p_variable, values, r_process_info);
                        const std::string& r_variable_name = p_variable->Name();
                        auto json_variable = json_element[r_variable_name];
                        for (unsigned int i = 0; i < values.size(); ++i) {
                            json_variable[std::to_string(i)].Append(values[i]);
                        }
                    }
                    for (const auto* p_variable : mGaussPointsOutputArrayVariables) {
                        std::vector<array_1d<double, 3>> values;
                        r_elem.CalculateOnIntegrationPoints(*p_variable, values, r_process_info);
                        const std::string& r_variable_name = p_variable->Name();
                        auto json_variable_x = json_element[r_variable_name + "_X"];
                        auto json_variable_y = json_element[r_variable_name + "_Y"];
                        auto json_variable_z = json_element[r_variable_name + "_Z"];
                        for (unsigned int i = 0; i < values.size(); ++i) {
                            json_variable_x[std::to_string(i)].Append(values[i][0]);
                            json_variable_y[std::to_string(i)].Append(values[i][1]);
                            json_variable_z[std::to_string(i)].Append(values[i][2]);
                        }
                    }
                    for (const auto* p_variable : mGaussPointsOutputVectorVariables) {
                        std::vector<Vector> values;
                        r_elem.CalculateOnIntegrationPoints(*p_variable, values, r_process_info);
                        const std::string& r_variable_name = p_variable->Name();
                        auto json_variable = json_element[r_variable_name];
                        for (unsigned int i = 0; i < values.size(); ++i) {
                            json_variable[std::to_string(i)].Append(values[i]);
                        }
                    }
                }
            }
        } else {
            auto json_element = mJsonFile["RESULTANT"];
            for (const auto* p_variable : mGaussPointsOutputVariables) {
                auto it_elem_begin = mpSubModelPart->ElementsBegin();
                std::vector<double> values;
                it_elem_begin->CalculateOnIntegrationPoints(*p_variable, values, r_process_info);
                std::vector<double> resultant_values(values.size(), 0.0);
                for (auto& r_elem : mpSubModelPart->Elements()) {
                    if (CheckFlag(r_elem.pGetGeometry())) {
                        r_elem.CalculateOnIntegrationPoints(*p_variable, values, r_process_info);
                        for (unsigned int i = 0; i < values.size(); ++i) {
                            resultant_values[i] += values[i];
                        }
                    }
                }
                // Add to the json
                auto json_variable = json_element[p_variable->Name()];
                for (unsigned int i = 0; i < resultant_values.size(); ++i) {
                    json_variable[std::to_string(i)].Append(resultant_values[i]);
                }
            }
            for (const auto* p_variable : mGaussPointsOutputArrayVariables) {
                auto it_elem_begin = mpSubModelPart->ElementsBegin();
                std::vector<array_1d<double, 3>> values;
                it_elem_begin->CalculateOnIntegrationPoints(*p_variable, values, r_process_info);
                std::vector<array_1d<double, 3>> resultant_values(values.size(), ZeroVector(3));
                for (auto& r_elem : mpSubModelPart->Elements()) {
                    if (CheckFlag(r_elem.pGetGeometry())) {
                        r_elem.CalculateOnIntegrationPoints(*p_variable, values, r_process_info);
                        for (unsigned int i = 0; i < values.size(); ++i) {
                            resultant_values[i] += values[i];
                        }
                    }
                }
                // Add to the json
                const std::string& r_variable_name = p_variable->Name();
                auto json_variable_x = json_element[r_variable_name + "_X"];
                auto json_variable_y = json_element[r_variable_name + "_Y"];
                auto json_variable_z = json_element[r_variable_name + "_Z"];
                for (unsigned int i = 0; i < resultant_values.size(); ++i) {
                    const std::string i_string = std::to_string(i);
                    json_variable_x[i_string].Append(resultant_values[i][0]);
                    json_variable_y[i_string].Append(resultant_values[i][1]);
                    json_variable_z[i_string].Append(resultant_values[i][2]);
                }
            }
            for (const auto* p_variable : mGaussPointsOutputVectorVariables) {
                auto it_elem_begin = mpSubModelPart->ElementsBegin();
                std::vector<Vector> values;
                it_elem_begin->CalculateOnIntegrationPoints(*p_variable, values, r_process_info);
                std::vector<Vector> resultant_values(values.size(), ZeroVector(values[0].size()));
                for (auto& r_elem : mpSubModelPart->Elements()) {
                    if (CheckFlag(r_elem.pGetGeometry())) {
                        r_elem.CalculateOnIntegrationPoints(*p_variable, values, mpSubModelPart->GetProcessInfo());
                        for (unsigned int i = 0; i < values.size(); ++i) {
                            resultant_values[i] += values[i];
                        }
                    }
                }
                // Add to the json
                auto json_variable = json_element[p_variable->Name()];
                for (unsigned int i = 0; i < resultant_values.size(); ++i) {
                    json_variable[std::to_string(i)].Append(resultant_values[i]);
                }
            }
        }
    }

    // Write the file
    std::ofstream output_file(mOutputFileName);
    output_file << mJsonFile.PrettyPrintJsonString();
}

}  // namespace Kratos.