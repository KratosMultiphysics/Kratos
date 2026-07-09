//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//
//

// System includes

// External includes

// Project includes
#include "input_output/unv_output.h"
#include "includes/kratos_filesystem.h"

namespace Kratos {

void UnvOutput::VariablesLists::Initialize(
    const std::vector<std::string>& rVariableNames,
    std::unordered_map<std::size_t, int>& rUnvVariableKeys
    )
{
    // Iterate through the provided variable names and populate the corresponding variable lists
    int double_counter = 0, array_counter = 0;
    for (auto& r_variable_name : rVariableNames) {
        if (KratosComponents<Variable<bool>>::Has(r_variable_name)) {
            const Variable<bool>& r_variable = KratosComponents<Variable<bool>>::Get(r_variable_name);
            mBoolVariables.push_back(&r_variable);
            
            // Assign UNV variable keys based on known variable names
            const std::size_t variable_key = r_variable.Key();
            const auto& r_pair = unv_scalar_variables.begin() + double_counter;
            KRATOS_WARNING("UnvOutput") << "Unknown variable: " << r_variable_name << ". Using UNV scalar " << r_pair->second << " variable name id: " << r_pair->first << std::endl;
            rUnvVariableKeys[variable_key] = r_pair->first;
            ++double_counter;
        } else if (KratosComponents<Variable<int>>::Has(r_variable_name)) {
            const Variable<int>& r_variable = KratosComponents<Variable<int>>::Get(r_variable_name);
            mIntVariables.push_back(&r_variable);
            
            // Assign UNV variable keys based on known variable names
            const std::size_t variable_key = r_variable.Key();
            const auto& r_pair = unv_scalar_variables.begin() + double_counter;
            KRATOS_WARNING("UnvOutput") << "Unknown variable: " << r_variable_name << ". Using UNV scalar " << r_pair->second << " variable name id: " << r_pair->first << std::endl;
            rUnvVariableKeys[variable_key] = r_pair->first;
            ++double_counter;
        } else if (KratosComponents<Variable<double>>::Has(r_variable_name)) {
            const Variable<double>& r_variable = KratosComponents<Variable<double>>::Get(r_variable_name);
            mDoubleVariables.push_back(&r_variable);
            
            // Assign UNV variable keys based on known variable names
            const std::size_t variable_key = r_variable.Key();
            if (r_variable_name == "TEMPERATURE") {
                rUnvVariableKeys[variable_key] = 5;
            } else if (StringUtilities::ContainsPartialString(r_variable_name, "PRESSURE")) {
                rUnvVariableKeys[variable_key] = 117;
            } else {
                const auto& r_pair = unv_scalar_variables.begin() + double_counter;
                KRATOS_WARNING("UnvOutput") << "Unknown variable: " << r_variable_name << ". Using UNV scalar " << r_pair->second << " variable name id: " << r_pair->first << std::endl;
                rUnvVariableKeys[variable_key] = r_pair->first;
                ++double_counter;
            }

        } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable_name)) {
            const Variable<array_1d<double, 3>>& r_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(r_variable_name);
            mArray1DVariables.push_back(&r_variable);

            // Assign UNV variable keys based on known variable names
            const std::size_t variable_key = r_variable.Key();
            if (r_variable_name == "VELOCITY"){
                rUnvVariableKeys[variable_key] = 11;
            } else if (r_variable_name == "DISPLACEMENT") {
                rUnvVariableKeys[variable_key] = 8;
            } else if (r_variable_name == "ACCELERATION") {
                rUnvVariableKeys[variable_key] = 12;
            } else if (r_variable_name == "ACCELERATION") {
                rUnvVariableKeys[variable_key] = 12;
            } else if (r_variable_name == "REACTION") {
                rUnvVariableKeys[variable_key] = 8;
            } else {
                const auto& r_pair = unv_vector3_variables.begin() + array_counter;
                KRATOS_WARNING("UnvOutput") << "Unknown variable: " << r_variable_name << ". Using UNV 3 components vector " << r_pair->second << " variable name id: " << r_pair->first << std::endl;
                rUnvVariableKeys[variable_key] = r_pair->first;
                ++array_counter;
            }
        }
        // else if (const auto* p_vector_var = dynamic_cast<const Variable<Vector>*>(&r_variable))
        //     mVectorVariables.push_back(p_vector_var); // NOTE: Current unsupported
        // else if (const auto* p_matrix_var = dynamic_cast<const Variable<Matrix>*>(&r_variable))
        //     mMatrixVariables.push_back(p_matrix_var); // NOTE: Current unsupported
    }
}

UnvOutput::UnvOutput(
    Model& rModel,
    Parameters Settings
    ) :  mrOutputModelPart(rModel.GetModelPart(Settings["model_part_name"].GetString()))
{

    // Validate and assign default parameters
    const Parameters default_parameters = GetDefaultParameters();
    Settings.ValidateAndAssignDefaults(default_parameters);

    // The output file name is constructed using the output path, custom name prefix, model part name, and custom name postfix.
    const std::string output_path = Settings["output_path"].GetString();
    const std::string custom_name_prefix = Settings["custom_name_prefix"].GetString();
    const std::string custom_name_postfix = Settings["custom_name_postfix"].GetString();
    const bool save_output_files_in_folder = Settings["save_output_files_in_folder"].GetBool();
    if (save_output_files_in_folder) {
        mOutputFileName = (std::filesystem::path(output_path) / (custom_name_prefix + mrOutputModelPart.Name() + custom_name_postfix + ".unv")).string();
    } else {
        mOutputFileName = custom_name_prefix + mrOutputModelPart.Name() + custom_name_postfix + ".unv";
    }

    // Set the entity type based on the provided settings. The entity type determines whether to print elements, conditions, or automatically determine based on the model part.
    const std::string entity_type_str = Settings["entity_type"].GetString();
    if (entity_type_str == "automatic") {
        mEntityType = EntityType::AUTOMATIC;
    } else if (entity_type_str == "elements") {
        mEntityType = EntityType::ELEMENTS;
    } else if (entity_type_str == "conditions") {
        mEntityType = EntityType::CONDITIONS;
    } else {
        KRATOS_ERROR << "Invalid entity_type: " << entity_type_str << ". Valid options are 'automatic', 'elements', or 'conditions'." << std::endl;
    }

    // Initialize the lists of historical and non-historical variables to be printed based on the provided settings.
    mHistoricalVariables.Initialize(Settings["nodal_solution_step_data_variables"].GetStringArray(), mUnvVariableKeys);
    mNonHistoricalVariables.Initialize(Settings["nodal_data_value_variables"].GetStringArray(), mUnvVariableKeys);
}

UnvOutput::UnvOutput(
    Kratos::ModelPart& rModelPart, 
    const std::string& rOutFileWithoutExtension
    ) : mrOutputModelPart(rModelPart),
        mOutputFileName(rOutFileWithoutExtension + ".unv") {
}

void UnvOutput::WriteMesh() {
    // Check if the output file has been initialized. If not, initialize it.
    if (!mInitializedOutputFile) {
        KRATOS_WARNING("UnvOutput") << "Output file has not been initialized yet. Initializing output file before writing mesh." << std::endl;
        InitializeOutputFile();
    }

    // Write the nodes
    WriteNodes();

    // Write the geometry (elements or conditions) based on availability
    if (mEntityType == EntityType::AUTOMATIC) {
        if (mrOutputModelPart.Elements().size() > 0) {
            // Write the elements if they exist in the model part
            WriteElements();
        } else if (mrOutputModelPart.Conditions().size() > 0) {
            KRATOS_WARNING("UnvOutput") << "No elements found in the model part. Writing conditions instead." << std::endl;
            
            // Write the conditions if they exist in the model part
            WriteConditions();
        } else {
            KRATOS_WARNING("UnvOutput") << "No elements or conditions found in the model part. No mesh will be written." << std::endl;
        }
    } else if (mEntityType == EntityType::ELEMENTS) {
        if (mrOutputModelPart.Elements().size() > 0) {
            // Write the elements if they exist in the model part
            WriteElements();
        } else {
            KRATOS_WARNING("UnvOutput") << "No elements found in the model part. No mesh will be written." << std::endl;
        }
    } else if (mEntityType == EntityType::CONDITIONS) {
        if (mrOutputModelPart.Conditions().size() > 0) {
            // Write the conditions if they exist in the model part
            WriteConditions();
        } else {
            KRATOS_WARNING("UnvOutput") << "No conditions found in the model part. No mesh will be written." << std::endl;
        }
    }

    // Set the flag to indicate that the mesh has been written
    mMeshWritten = true;
}

void UnvOutput::PrintOutput() {
    // Check if the mesh has been written before printing output. If not, write the mesh first.
    if (!mMeshWritten) {
        KRATOS_WARNING("UnvOutput") << "Mesh has not been written yet. Writing mesh before printing output." << std::endl;
        WriteMesh();
    }

    // Get the current time step from the model part's process info
    const auto& r_process_info = mrOutputModelPart.GetProcessInfo();
    const double time_step = r_process_info.GetValue(TIME);

    // Print output for nodal historical variables
    for (const auto& p_variable : mHistoricalVariables.mBoolVariables) {
        WriteNodalResults(*p_variable, time_step);
    }
    for (const auto& p_variable : mHistoricalVariables.mIntVariables) {
        WriteNodalResults(*p_variable, time_step);
    }
    for (const auto& p_variable : mHistoricalVariables.mDoubleVariables) {
        WriteNodalResults(*p_variable, time_step);
    }
    for (const auto& p_variable : mHistoricalVariables.mArray1DVariables) {
        WriteNodalResults(*p_variable, time_step);
    }

    // Print output for nodal non-historical variables
    for (const auto& p_variable : mNonHistoricalVariables.mBoolVariables) {
        WriteNodalNonHistoricalResults(*p_variable, time_step);
    }
    for (const auto& p_variable : mNonHistoricalVariables.mIntVariables) {
        WriteNodalNonHistoricalResults(*p_variable, time_step);
    }
    for (const auto& p_variable : mNonHistoricalVariables.mDoubleVariables) {
        WriteNodalNonHistoricalResults(*p_variable, time_step);
    }
    for (const auto& p_variable : mNonHistoricalVariables.mArray1DVariables) {
        WriteNodalNonHistoricalResults(*p_variable, time_step);
    }
}

void UnvOutput::InitializeOutputFile() {
    // Open the output file in truncation mode to clear any existing content and prepare it for writing.
    std::ofstream rOutputFile;
    rOutputFile.open(mOutputFileName, std::ios::out | std::ios::trunc);
    rOutputFile.close();

    // Set the flag to indicate that the output file has been initialized
    mInitializedOutputFile = true;
}

void UnvOutput::WriteNodes() {
    std::ofstream rOutputFile;
    rOutputFile.open(mOutputFileName, std::ios::out | std::ios::app);

    rOutputFile << std::scientific;
    rOutputFile << std::setprecision(15);

    const int export_coordinate_system = 0;
    const int displacement_coordinate_system_number = 0;
    const int color = 0;

    rOutputFile << std::setw(6) << "-1" << "\n";
    rOutputFile << std::setw(6) << as_integer(DatasetID::NODES_DATASET) << "\n";

    int node_label;
    double x_coordinate, y_coordinate, z_coordinate;
    array_1d<double, 3> node_displacement = ZeroVector(3);
    for (auto& r_node : mrOutputModelPart.Nodes()) {
        node_label = r_node.Id();
        // Get the node displacement if the deformed configuration is to be written
        if (mWriteDeformedConfiguration) {
            noalias(node_displacement) = r_node.FastGetSolutionStepValue(DISPLACEMENT);
        }
        x_coordinate = r_node.X() + node_displacement[0];
        y_coordinate = r_node.Y() + node_displacement[1];
        z_coordinate = r_node.Z() + node_displacement[2];
        rOutputFile << std::setw(10) << node_label << std::setw(10) << export_coordinate_system
                    << std::setw(10)
                    << displacement_coordinate_system_number << std::setw(10) << color << "\n";
        rOutputFile << std::setw(25) << x_coordinate << std::setw(25) << y_coordinate << std::setw(25)
                    << z_coordinate << "\n";
    }
    rOutputFile << std::setw(6) << "-1" << "\n";

    rOutputFile.close();
}

void UnvOutput::WriteElements() {
    std::ofstream rOutputFile;
    rOutputFile.open(mOutputFileName, std::ios::out | std::ios::app);

    const int physical_property_table_number = 1;
    const int material_property_table_number = 1;
    const int color = 0;

    rOutputFile << std::setw(6) << "-1" << "\n";
    rOutputFile << std::setw(6) << as_integer(DatasetID::ELEMENTS_DATASET) << "\n";

    for (auto& r_element : mrOutputModelPart.Elements()) {
        const int element_label = r_element.Id();
        auto& r_geometry = r_element.GetGeometry();
        // Write triangles
        if (r_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3) {
            const int fe_descriptor_id = 41; // Plane Stress Linear Triangle
            const int number_of_nodes = 3;
            rOutputFile << std::setw(10) << element_label;
            rOutputFile << std::setw(10) << fe_descriptor_id;
            rOutputFile << std::setw(10) << physical_property_table_number;
            rOutputFile << std::setw(10) << material_property_table_number;
            rOutputFile << std::setw(10) << color;
            rOutputFile << std::setw(10) << number_of_nodes << "\n";
            rOutputFile << std::setw(10) << r_geometry[0].Id();
            rOutputFile << std::setw(10) << r_geometry[1].Id();
            rOutputFile << std::setw(10) << r_geometry[2].Id() << "\n";
        } else if (r_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) {
            const int fe_descriptor_id = 111; // Solid linear tetrahedron
            const int number_of_nodes = 4;
            rOutputFile << std::setw(10) << element_label;
            rOutputFile << std::setw(10) << fe_descriptor_id;
            rOutputFile << std::setw(10) << physical_property_table_number;
            rOutputFile << std::setw(10) << material_property_table_number;
            rOutputFile << std::setw(10) << color;
            rOutputFile << std::setw(10) << number_of_nodes << "\n";
            rOutputFile << std::setw(10) << r_geometry[0].Id();
            rOutputFile << std::setw(10) << r_geometry[1].Id();
            rOutputFile << std::setw(10) << r_geometry[2].Id();
            rOutputFile << std::setw(10) << r_geometry[3].Id() << "\n";
        } else {
            KRATOS_ERROR << "Element with ID " << element_label << " has unsupported geometry for UNV output." << std::endl;
        }
    }
    rOutputFile << std::setw(6) << "-1" << "\n";
    rOutputFile.close();
}

void UnvOutput::WriteConditions() {
    std::ofstream rOutputFile;
    rOutputFile.open(mOutputFileName, std::ios::out | std::ios::app);

    const int physical_property_table_number = 1;
    const int material_property_table_number = 1;
    const int color = 0;

    rOutputFile << std::setw(6) << "-1" << "\n";
    rOutputFile << std::setw(6) << as_integer(DatasetID::ELEMENTS_DATASET) << "\n";

    for (auto& r_condition : mrOutputModelPart.Conditions()) {
        const int element_label = r_condition.Id();
        auto& r_geometry = r_condition.GetGeometry();
        // Write lines and triangles
        if (r_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) {
            const int fe_descriptor_id = 21; // Linear beam
            const int number_of_nodes = 2;
            rOutputFile << std::setw(10) << element_label;
            rOutputFile << std::setw(10) << fe_descriptor_id;
            rOutputFile << std::setw(10) << physical_property_table_number;
            rOutputFile << std::setw(10) << material_property_table_number;
            rOutputFile << std::setw(10) << color;
            rOutputFile << std::setw(10) << number_of_nodes << "\n";
            rOutputFile << std::setw(10) << r_geometry[0].Id();
            rOutputFile << std::setw(10) << r_geometry[1].Id();
        } else if (r_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) {
            const int fe_descriptor_id = 91; // Thin Shell Linear Triangle
            const int number_of_nodes = 3;
            rOutputFile << std::setw(10) << element_label;
            rOutputFile << std::setw(10) << fe_descriptor_id;
            rOutputFile << std::setw(10) << physical_property_table_number;
            rOutputFile << std::setw(10) << material_property_table_number;
            rOutputFile << std::setw(10) << color;
            rOutputFile << std::setw(10) << number_of_nodes << "\n";
            rOutputFile << std::setw(10) << r_geometry[0].Id();
            rOutputFile << std::setw(10) << r_geometry[1].Id();
            rOutputFile << std::setw(10) << r_geometry[2].Id() << "\n";
        } else {
            KRATOS_ERROR << "Element with ID " << element_label << " has unsupported geometry for UNV output." << std::endl;
        }
    }
    rOutputFile << std::setw(6) << "-1" << "\n";
    rOutputFile.close();
}

void UnvOutput::WriteNodalResults(const Variable<bool>& rVariable, const double TimeStep) {
    WriteNodalResultRecords<Variable<bool>, WriteType::HISTORICAL>(rVariable, 1, TimeStep);
}

void UnvOutput::WriteNodalResults(const Variable<int>& rVariable, const double TimeStep) {
    WriteNodalResultRecords<Variable<int>, WriteType::HISTORICAL>(rVariable, 1, TimeStep);
}   

void UnvOutput::WriteNodalResults(const Variable<double>& rVariable, const double TimeStep) {
    WriteNodalResultRecords<Variable<double>, WriteType::HISTORICAL>(rVariable, 1, TimeStep);
} 

void UnvOutput::WriteNodalResults(const Variable<array_1d<double,3>>& rVariable, const double TimeStep) {
    WriteNodalResultRecords<Variable<array_1d<double,3>>, WriteType::HISTORICAL>(rVariable, 3, TimeStep);
}

void UnvOutput::WriteNodalResults(const Variable<Vector>& rVariable, const double TimeStep) {
    KRATOS_ERROR << "Dynamic Vector results are not yet supported in UNV" << std::endl;
}

void UnvOutput::WriteNodalResults(const Variable<Matrix>& rVariable, const double TimeStep) {
    KRATOS_ERROR << "Matrix results are not yet supported in UNV" << std::endl;
}

void UnvOutput::WriteNodalNonHistoricalResults(const Variable<bool>& rVariable, const double TimeStep) {
    WriteNodalResultRecords<Variable<bool>, WriteType::NON_HISTORICAL>(rVariable, 1, TimeStep);
}

void UnvOutput::WriteNodalNonHistoricalResults(const Variable<int>& rVariable, const double TimeStep) {
    WriteNodalResultRecords<Variable<int>, WriteType::NON_HISTORICAL>(rVariable, 1, TimeStep);
}   

void UnvOutput::WriteNodalNonHistoricalResults(const Variable<double>& rVariable, const double TimeStep) {
    WriteNodalResultRecords<Variable<double>, WriteType::NON_HISTORICAL>(rVariable, 1, TimeStep);
} 

void UnvOutput::WriteNodalNonHistoricalResults(const Variable<array_1d<double,3>>& rVariable, const double TimeStep) {
    WriteNodalResultRecords<Variable<array_1d<double,3>>, WriteType::NON_HISTORICAL>(rVariable, 3, TimeStep);
}

void UnvOutput::WriteNodalNonHistoricalResults(const Variable<Vector>& rVariable, const double TimeStep) {
    KRATOS_ERROR << "Dynamic Vector results are not yet supported in UNV" << std::endl;
}

void UnvOutput::WriteNodalNonHistoricalResults(const Variable<Matrix>& rVariable, const double TimeStep) {
    KRATOS_ERROR << "Matrix results are not yet supported in UNV" << std::endl;
}

UnvOutput::DataCharacteristics UnvOutput::GetDataType(const Variable<bool>& rVariable) {
    return UnvOutput::DataCharacteristics::SCALAR;
}

UnvOutput::DataCharacteristics UnvOutput::GetDataType(const Variable<int>& rVariable) {
    return UnvOutput::DataCharacteristics::SCALAR;
}   

UnvOutput::DataCharacteristics UnvOutput::GetDataType(const Variable<double>& rVariable) {
    return UnvOutput::DataCharacteristics::SCALAR;
} 

UnvOutput::DataCharacteristics UnvOutput::GetDataType(const Variable<array_1d<double,3>>& rVariable) {
    return UnvOutput::DataCharacteristics::D3_DOF_GLOBAL_TRANSLATION_VECTOR;
}

UnvOutput::DataCharacteristics UnvOutput::GetDataType(const Variable<Vector>& rVariable) {
    KRATOS_ERROR << "Dynamic Vector results are not yet supported in UNV" << std::endl;
} 

UnvOutput::DataCharacteristics UnvOutput::GetDataType(const Variable<Matrix>& rVariable) {
    KRATOS_ERROR << "Matrix results are not yet supported in UNV" << std::endl;
}

void UnvOutput::WriteNodalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<bool>& rVariable) {
    rOutputFile << std::setw(13) << rNode.FastGetSolutionStepValue(rVariable) << "\n";
}

void UnvOutput::WriteNodalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<int>& rVariable) {
    rOutputFile << std::setw(13) << rNode.FastGetSolutionStepValue(rVariable) << "\n";
}

void UnvOutput::WriteNodalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<double>& rVariable) {
    rOutputFile << std::setw(13) << rNode.FastGetSolutionStepValue(rVariable) << "\n";
}

void UnvOutput::WriteNodalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<array_1d<double,3>>& rVariable) {
    const auto& r_temp = rNode.FastGetSolutionStepValue(rVariable);

    rOutputFile << std::setw(13) << r_temp[0];
    rOutputFile << std::setw(13) << r_temp[1];
    rOutputFile << std::setw(13) << r_temp[2];
    rOutputFile << "\n";
}

void UnvOutput::WriteNodalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<Vector>& rVariable) {
    KRATOS_ERROR << "Dynamic Vector results are not yet supported by in UNV" << std::endl;
}

void UnvOutput::WriteNodalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<Matrix>& rVariable) {
    KRATOS_ERROR << "Matrix results are not yet supported by in UNV" << std::endl;
}

void UnvOutput::WriteNodalNonHistoricalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<bool>& rVariable) {
    rOutputFile << std::setw(13) << rNode.GetValue(rVariable) << "\n";
}

void UnvOutput::WriteNodalNonHistoricalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<int>& rVariable) {
    rOutputFile << std::setw(13) << rNode.GetValue(rVariable) << "\n";
}

void UnvOutput::WriteNodalNonHistoricalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<double>& rVariable) {
    rOutputFile << std::setw(13) << rNode.GetValue(rVariable) << "\n";
}

void UnvOutput::WriteNodalNonHistoricalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<array_1d<double,3>>& rVariable) {
    const auto& r_temp = rNode.GetValue(rVariable);

    rOutputFile << std::setw(13) << r_temp[0];
    rOutputFile << std::setw(13) << r_temp[1];
    rOutputFile << std::setw(13) << r_temp[2];
    rOutputFile << "\n";
}

void UnvOutput::WriteNodalNonHistoricalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<Vector>& rVariable) {
    KRATOS_ERROR << "Dynamic Vector results are not yet supported by in UNV" << std::endl;
}

void UnvOutput::WriteNodalNonHistoricalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<Matrix>& rVariable) {
    KRATOS_ERROR << "Matrix results are not yet supported by in UNV" << std::endl;
}

int UnvOutput::GetUnvVariableName(const VariableData& rVariable) {
    const std::size_t variable_key = rVariable.Key();
    if (mUnvVariableKeys.find(variable_key) != mUnvVariableKeys.end()) {
        return mUnvVariableKeys[variable_key];
    }
        
    return variable_key + 1000; // If not found, return a value greater than or equal to 1000
}

const Parameters UnvOutput::GetDefaultParameters() const
{
    Parameters default_parameters(R"(
    {
        "model_part_name"                             : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "output_control_type"                         : "step",
        "output_interval"                             : 1.0,
        "output_path"                                 : "UNV_Output",
        "custom_name_prefix"                          : "",
        "custom_name_postfix"                         : "",                                     
        "save_output_files_in_folder"                 : true,
        "entity_type"                                 : "automatic",
        "write_deformed_configuration"                : false,
        "nodal_solution_step_data_variables"          : [],
        "nodal_data_value_variables"                  : []
    })");

    return default_parameters;
}

} // namespace Kratos
