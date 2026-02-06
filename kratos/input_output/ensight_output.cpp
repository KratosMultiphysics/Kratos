//
//   |   /         |
//   ' /    __| _` | __|  _ \   __|
//   . \   |   (   | |    (   |\__ `
//  _|\_\_|  \__,_|\__|\___/ ____/
//         Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <iomanip>

// External includes

// Project includes
#include "input_output/ensight_output.h"
#include "containers/model.h"
#include "includes/kratos_filesystem.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/input_output_utilities.h"

namespace Kratos
{

Parameters EnSightOutput::GetDefaultParameters()
{
    KRATOS_TRY

    return Parameters(R"(
    {
        "model_part_name"                             : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "ensight_file_format"                         : "gold", // Options: "6", "gold"
        "file_format"                                 : "ascii", // Options: "ascii", "binary" // TODO: Check if binary is supported
        "output_precision"                            : 6,
        "step_label_precision"                        : 4,
        "output_control_type"                         : "step",
        "output_interval"                             : 1.0,
        "output_sub_model_parts"                      : false,
        "output_path"                                 : "EnSight_Output",
        "custom_name_prefix"                          : "",
        "custom_name_postfix"                         : "",
        "entity_type"                                 : "automatic",
        "save_output_files_in_folder"                 : true,
        "evolving_geometry"                           : true, // For meshes with evolving geometry
        "nodal_solution_step_data_variables"          : [],
        "nodal_data_value_variables"                  : [],
        "nodal_flags"                                 : [],
        "element_data_value_variables"                : [],
        "element_flags"                               : [],
        "condition_data_value_variables"              : [],
        "condition_flags"                             : [],
        "gauss_point_variables_extrapolated_to_nodes" : [],
        "gauss_point_variables_in_elements"           : []
    })" );

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

EnSightOutput::EnSightOutput(ModelPart& rModelPart, Parameters ThisParameters)
    : mrModelPart(rModelPart),
      mOutputSettings(ThisParameters)
{
    KRATOS_TRY

    // Validate and assign default parameters
    mOutputSettings.ValidateAndAssignDefaults(GetDefaultParameters());

    // Initialize other variables
    const std::string ensight_file_format = mOutputSettings["ensight_file_format"].GetString();
    if (ensight_file_format == "6") {
        mEnSightFileFormat = EnSightFileFormat::EnSight6;
    } else if (ensight_file_format == "gold") {
        mEnSightFileFormat = EnSightFileFormat::EnSightGold;
    } else {
        KRATOS_ERROR << "Unknown ensight_file_format \"" << ensight_file_format << "\". Available options are \"6\", \"gold\"" << std::endl;
    }

    const std::string file_format = mOutputSettings["file_format"].GetString();
    if (file_format == "ascii") {
        mFileFormat = FileFormat::ASCII;
    } else if (file_format == "binary") {
        mFileFormat = FileFormat::BINARY;
    } else {
        KRATOS_ERROR << "Option for \"file_format\": " << file_format << " not recognised! Possible options are: \"ascii\", \"binary\"" << std::endl;
    }

    // Set the default precision
    mDefaultPrecision = mOutputSettings["output_precision"].GetInt();
    KRATOS_ERROR_IF(mDefaultPrecision < 0) << "\"output_precision\" must be a non-negative integer." << std::endl;
    KRATOS_ERROR_IF(mDefaultPrecision > 6) << "\"output_precision\" cannot be higher than 6 for EnSight format due to floating point precision limitations." << std::endl;

    // Set the step label precision
    mStepLabelPrecision = mOutputSettings["step_label_precision"].GetInt();

    // Adding GP variables to nodal data variables list
    if(mOutputSettings["gauss_point_variables_extrapolated_to_nodes"].size() > 0) {
        Parameters gauss_intergration_param_non_hist = Parameters(R"(
        {
            "echo_level"                 : 0,
            "area_average"               : true,
            "average_variable"           : "NODAL_AREA",
            "list_of_variables"          : [],
            "extrapolate_non_historical" : true
        })");

        gauss_intergration_param_non_hist.SetValue("list_of_variables", mOutputSettings["gauss_point_variables_extrapolated_to_nodes"]);

        for(auto const& gauss_var : mOutputSettings["gauss_point_variables_extrapolated_to_nodes"]) {
            mOutputSettings["nodal_data_value_variables"].Append(gauss_var);
        }

        // Making the gauss point to nodes process if any gauss point result is requested for
        mpGaussToNodesProcess = Kratos::make_unique<IntegrationValuesExtrapolationToNodesProcess>(rModelPart, gauss_intergration_param_non_hist);
    }

    // Get the entity type
    const std::string entity_type = mOutputSettings["entity_type"].GetString();
    if (entity_type == "element") {
        mEntityType = EntityType::ELEMENT;
    } else if (entity_type == "condition") {
        mEntityType = EntityType::CONDITION;
    } else if (entity_type == "automatic") {
        mEntityType = EntityType::AUTOMATIC;

        const std::size_t num_elements = rModelPart.GetCommunicator().GlobalNumberOfElements();
        const std::size_t num_conditions = rModelPart.GetCommunicator().GlobalNumberOfConditions();

        KRATOS_WARNING_IF("EnSightOutput", num_elements > 0 && num_conditions > 0) << "Modelpart \"" << rModelPart.Name() << "\" has both elements and conditions.\nGiving precedence to elements and writing only elements!" << std::endl;
    } else {
        KRATOS_ERROR << "Unknown entity_type \"" << entity_type << "\". Available options are \"element\", \"condition\", \"automatic\"" << std::endl;
    }

    // Create part data
    if (!mOutputSettings["evolving_geometry"].GetBool()) {
        UpdatePartData();

        // If no output filename is provided, use the model part name
        const std::string& r_custom_name_prefix = mOutputSettings["custom_name_prefix"].GetString();
        const std::string& r_custom_name_postfix = mOutputSettings["custom_name_postfix"].GetString();
        const std::string geometry_filename = r_custom_name_prefix + mrModelPart.Name() + r_custom_name_postfix;

        // Write the geometry file
        WriteGeometryFile(geometry_filename);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void EnSightOutput::PrintOutput(const std::string& rOutputFilename)
{
    KRATOS_TRY

    // Update part data
    if (mOutputSettings["evolving_geometry"].GetBool()) {
        UpdatePartData();
    }

    // Prepare the variable types
    PrepareVariableTypeMap();

    // Extrapolate GP results to nodes if requested
    PrepareGaussPointResults();

    // If no output filename is provided, use the model part name
    const std::string& r_custom_name_prefix = mOutputSettings["custom_name_prefix"].GetString();
    const std::string& r_custom_name_postfix = mOutputSettings["custom_name_postfix"].GetString();
    std::string base_filename = rOutputFilename == "" ? r_custom_name_prefix + mrModelPart.Name() + r_custom_name_postfix : r_custom_name_prefix + rOutputFilename + r_custom_name_postfix;

    // Write the geometry file for the current step
    const std::string step_label = InputOutputUtilities::GenerateStepLabel(mTimeValues.size(), mStepLabelPrecision);
    std::string geometry_filename = base_filename + "." + step_label;
    WriteGeometryFile(geometry_filename);

    // Write each variable file for the current step

    // Nodal historical variables
    for (const auto& r_variable_name_param : mOutputSettings["nodal_solution_step_data_variables"]) {
        const std::string r_variable_name = r_variable_name_param.GetString();
        const std::string var_filename = base_filename + "." + r_variable_name + "." + step_label + ".node";
        WriteNodalVariableToFile(r_variable_name, true, var_filename);
    }

    // Nodal non-historical variables
    for (const auto& r_variable_name_param : mOutputSettings["nodal_data_value_variables"]) {
        const std::string r_variable_name = r_variable_name_param.GetString();
        const std::string var_filename = base_filename + "." + r_variable_name + "." + step_label + ".node";
        WriteNodalVariableToFile(r_variable_name, false, var_filename);
    }

    // Nodal flags
    for (const auto& r_variable_name_param : mOutputSettings["nodal_flags"]) {
        const std::string r_flag_name = r_variable_name_param.GetString();
        const std::string var_filename = base_filename + "." + r_flag_name + "." + step_label + ".node";
        WriteNodalFlagToFile(r_flag_name, var_filename);
    }

    // Elemental variables
    for (const auto& r_variable_name_param : mOutputSettings["element_data_value_variables"]) {
        const std::string r_variable_name = r_variable_name_param.GetString();
        const std::string var_filename = base_filename + "." + r_variable_name + "." + step_label + ".elem";
        WriteGeometricalVariableToFile(r_variable_name, var_filename);
    }

    // Elemental flags
    for (const auto& r_variable_name_param : mOutputSettings["element_flags"]) {
        const std::string r_flag_name = r_variable_name_param.GetString();
        const std::string var_filename = base_filename + "." + r_flag_name + "." + step_label + ".elem";
        WriteGeometricalFlagToFile(r_flag_name, var_filename);
    }

    // Gauss points in elements
    for (const auto& r_variable_name_param : mOutputSettings["gauss_point_variables_in_elements"]) {
        const std::string r_variable_name = r_variable_name_param.GetString();
        const std::string var_filename = base_filename + "." + r_variable_name + "." + step_label + ".gp_elem";
        WriteGeometricalGaussVariableToFile(r_variable_name, var_filename);
    }

    // Conditional variables
    for (const auto& r_variable_name_param : mOutputSettings["condition_data_value_variables"]) {
        const std::string r_variable_name = r_variable_name_param.GetString();
        const std::string var_filename = base_filename + "." + r_variable_name + "." + step_label + ".cond";
        WriteGeometricalVariableToFile(r_variable_name, var_filename, false);
    }

    // Conditional flags
    for (const auto& r_variable_name_param : mOutputSettings["condition_flags"]) {
        const std::string r_flag_name = r_variable_name_param.GetString();
        const std::string var_filename = base_filename + "." + r_flag_name + "." + step_label + ".cond";
        WriteGeometricalFlagToFile(r_flag_name, var_filename, false);
    }

    // Determine the base filename for this time step
    const double current_time = mrModelPart.GetProcessInfo()[TIME];
    mTimeValues.push_back(current_time);

    // (Re)Write the main .case file to link everything
    WriteCaseFile();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

EnSightOutput::EntityType EnSightOutput::GetEntityType(const ModelPart& rModelPart) const
{
    KRATOS_TRY

    const std::size_t num_elements = rModelPart.GetCommunicator().GlobalNumberOfElements();
    const std::size_t num_conditions = rModelPart.GetCommunicator().GlobalNumberOfConditions();

    if (mEntityType == EntityType::ELEMENT) {
        return (num_elements > 0) ? EntityType::ELEMENT : EntityType::NONE;
    } else if (mEntityType == EntityType::CONDITION) {
        return (num_conditions > 0) ? EntityType::CONDITION : EntityType::NONE;
    }

    // Automatic: elements take precedence over conditions
    if (num_elements > 0) {
        return EntityType::ELEMENT;
    } else if(num_conditions > 0) {
        return EntityType::CONDITION;
    } else {
        return EntityType::NONE;
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void EnSightOutput::PrepareVariableTypeMap()
{
    KRATOS_TRY

    // Clear the existing map
    mVariableTypeMap.clear();

    // Populate the map with variable names and their types

    // Populate nodal solution step data variables
    for (const auto& r_variable_name_param : mOutputSettings["nodal_solution_step_data_variables"]) {
        const std::string& r_variable_name = r_variable_name_param.GetString();
        if (mVariableTypeMap.find(r_variable_name) == mVariableTypeMap.end()) {
            mVariableTypeMap[r_variable_name] = GetVariableType(r_variable_name, &mrModelPart, EntityType::NODE, true);
        }
    }

    // Populate nodal data value variables
    for (const auto& r_variable_name_param : mOutputSettings["nodal_data_value_variables"]) {
        const std::string& r_variable_name = r_variable_name_param.GetString();
        if (mVariableTypeMap.find(r_variable_name) == mVariableTypeMap.end()) {
            mVariableTypeMap[r_variable_name] = GetVariableType(r_variable_name, &mrModelPart, EntityType::NODE, false);
        }
    }

    // Populate nodal flags
    for (const auto& r_variable_name_param : mOutputSettings["nodal_flags"]) {
        const std::string& r_variable_name = r_variable_name_param.GetString();
        if (mVariableTypeMap.find(r_variable_name) == mVariableTypeMap.end()) {
            mVariableTypeMap[r_variable_name] = VariableType::SCALAR;
        }
    }

    // Populate element data value variables
    for (const auto& r_variable_name_param : mOutputSettings["element_data_value_variables"]) {
        const std::string& r_variable_name = r_variable_name_param.GetString();
        if (mVariableTypeMap.find(r_variable_name) == mVariableTypeMap.end()) {
            mVariableTypeMap[r_variable_name] = GetVariableType(r_variable_name, &mrModelPart, EntityType::ELEMENT, false);
        }
    }

    // Populate elemental flags
    for (const auto& r_variable_name_param : mOutputSettings["element_flags"]) {
        const std::string& r_variable_name = r_variable_name_param.GetString();
        if (mVariableTypeMap.find(r_variable_name) == mVariableTypeMap.end()) {
            mVariableTypeMap[r_variable_name] = VariableType::SCALAR;
        }
    }

    // Populate element Gauss point variables
    for (const auto& r_variable_name_param : mOutputSettings["gauss_point_variables_in_elements"]) {
        const std::string& r_variable_name = r_variable_name_param.GetString();
        if (mVariableTypeMap.find(r_variable_name) == mVariableTypeMap.end()) {
            mVariableTypeMap[r_variable_name] = GetVariableType(r_variable_name, &mrModelPart, EntityType::ELEMENT, true);
        }
    }

    // Populate condition data value variables
    for (const auto& r_variable_name_param : mOutputSettings["condition_data_value_variables"]) {
        const std::string& r_variable_name = r_variable_name_param.GetString();
        if (mVariableTypeMap.find(r_variable_name) == mVariableTypeMap.end()) {
            mVariableTypeMap[r_variable_name] = GetVariableType(r_variable_name, &mrModelPart, EntityType::CONDITION, false);
        }
    }

    // Populate conditional flags
    for (const auto& r_variable_name_param : mOutputSettings["condition_flags"]) {
        const std::string& r_variable_name = r_variable_name_param.GetString();
        if (mVariableTypeMap.find(r_variable_name) == mVariableTypeMap.end()) {
            mVariableTypeMap[r_variable_name] = VariableType::SCALAR;
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void EnSightOutput::PrepareGaussPointResults()
{
    KRATOS_TRY

    // Check if any Gauss point variables are requested for extrapolation
    if (mOutputSettings["gauss_point_variables_extrapolated_to_nodes"].size() > 0) {
        // Ensure the Gauss point to nodes process is initialized
        KRATOS_ERROR_IF_NOT(mpGaussToNodesProcess) << "Gauss point extrapolation process is not initialized!" << std::endl;
        // Execute the process to extrapolate Gauss point results to nodes
        mpGaussToNodesProcess->Execute();
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void EnSightOutput::WriteCaseFile()
{
    KRATOS_TRY

    // Get the model part name from the settings
    const std::string base_filename = mrModelPart.Name();
    const std::string case_filename = GetOutputFileName(base_filename, ".case");

    // Open the case file for writing
    std::ofstream case_file(case_filename);

    // Check if the file was opened successfully
    KRATOS_ERROR_IF_NOT(case_file.is_open()) << "File \"" << case_filename << "\" could not be opened!" << std::endl;

    // Example with evolving geometry:
    /* ENSIGHT 6 */
    // FORMAT
    // Type: ensight
    //
    // GEOMETRY
    // model:            1                 example2.geo**
    //
    // VARIABLE
    // scalar per node:  1   Stress        example2.scl**
    // vector per node:  1   Displacement  example2.dis**
    //
    // TIME
    // time set:               1
    // number of steps:        3
    // filename start number:  0
    // filename increment:     1
    // time values:            1.0 2.0 3.0

    /* ENSIGHT GOLD */
    // FORMAT
    // type:       ensight gold
    //
    // GEOMETRY
    // model:             1                   exgold2.geo**
    //
    // VARIABLE
    // scalar per node:   1   Stress          exgold2.scl**
    // vector per node:   1   Displacement    exgold2.vct**
    //
    // TIME
    // time set:               1
    // number of steps:        3
    // filename start number:  0
    // filename increment:     1
    // time values:            1.0    2.0    3.0

    // FORMAT Section
    case_file << "FORMAT\n";
    if (mEnSightFileFormat == EnSightFileFormat::EnSightGold) {
        case_file << "type: ensight gold\n\n";
    } else {
        case_file << "type: ensight\n\n";
    }

    // GEOMETRY Section
    case_file << "GEOMETRY\n";
    // Using wildcards for transient geometry
    std::string asterisk_label = "";
    for (std::size_t i = 0; i < mStepLabelPrecision; ++i) {
        asterisk_label += "*";
    }

    // Write the geometry filename line
    if (!mOutputSettings["evolving_geometry"].GetBool()) {
        // If no output filename is provided, use the model part name
        const std::string& r_custom_name_prefix = mOutputSettings["custom_name_prefix"].GetString();
        const std::string& r_custom_name_postfix = mOutputSettings["custom_name_postfix"].GetString();
        case_file << "model: 1 " << r_custom_name_prefix << mrModelPart.Name() << r_custom_name_postfix << ".geo\n\n";
    } else {
        case_file << "model: 1 " << base_filename << "." << asterisk_label << ".geo\n\n";
    }

    // VARIABLE Section
    case_file << "VARIABLE\n";

    // Nodal solution step variables
    for (const auto& r_variable_name_param : mOutputSettings["nodal_solution_step_data_variables"]) {
        const std::string& r_variable_name = r_variable_name_param.GetString();
        const auto variable_type = mVariableTypeMap[r_variable_name];
        const std::string file_extension = GetExtensionFile(variable_type);
        if (variable_type == VariableType::SCALAR) {
            case_file << "scalar per node: 1 " << r_variable_name << " " << base_filename << "." << r_variable_name << "." << asterisk_label << ".node" << file_extension << "\n";
        } else if (variable_type == VariableType::VECTOR) {
            case_file << "vector per node: 1 " << r_variable_name << " " << base_filename << "." << r_variable_name << "." << asterisk_label << ".node" << file_extension << "\n";
        } else if (variable_type == VariableType::TENSOR_SYMMETRIC) {
            case_file << "tensor symm per node: 1 " << r_variable_name << " " << base_filename << "." << r_variable_name << "." << asterisk_label << ".node" << file_extension << "\n";
        } else {
            KRATOS_ERROR << "Unknown variable type for element data variable: " << r_variable_name << std::endl;
        }
    }

    // Nodal solution data variables
    for (const auto& r_variable_name_param : mOutputSettings["nodal_data_value_variables"]) {
        const std::string& r_variable_name = r_variable_name_param.GetString();
        const auto variable_type = mVariableTypeMap[r_variable_name];
        const std::string file_extension = GetExtensionFile(variable_type);
        if (variable_type == VariableType::SCALAR) {
            case_file << "scalar per node: 1 " << r_variable_name << " " << base_filename << "." << r_variable_name << "." << asterisk_label << ".node" << file_extension << "\n";
        } else if (variable_type == VariableType::VECTOR) {
            case_file << "vector per node: 1 " << r_variable_name << " " << base_filename << "." << r_variable_name << "." << asterisk_label << ".node" << file_extension << "\n";
        } else if (variable_type == VariableType::TENSOR_SYMMETRIC) {
            case_file << "tensor symm per node: 1 " << r_variable_name << " " << base_filename << "." << r_variable_name << "." << asterisk_label << ".node" << file_extension << "\n";
        } else {
            KRATOS_ERROR << "Unknown variable type for element data variable: " << r_variable_name << std::endl;
        }
    }

    // Nodal flags (always VariableType::SCALAR)
    for (const auto& r_flag_name_param : mOutputSettings["nodal_flags"]) {
        const std::string& r_flag_name = r_flag_name_param.GetString();
        case_file << "scalar per node: 1 " << r_flag_name << " " << base_filename << "." << r_flag_name << "." << asterisk_label << ".node" << ".scl\n";
    }

    // Element data variables
    for (const auto& r_variable_name_param : mOutputSettings["element_data_value_variables"]) {
        const std::string& r_variable_name = r_variable_name_param.GetString();
        const auto variable_type = mVariableTypeMap[r_variable_name];
        const std::string file_extension = GetExtensionFile(variable_type);
        if (variable_type == VariableType::SCALAR) {
            case_file << "scalar per element: 1 " << r_variable_name << " " << base_filename << "." << r_variable_name << "." << asterisk_label << ".elem" << file_extension << "\n";
        } else if (variable_type == VariableType::VECTOR) {
            case_file << "vector per element: 1 " << r_variable_name << " " << base_filename << "." << r_variable_name << "." << asterisk_label << ".elem" << file_extension << "\n";
        } else if (variable_type == VariableType::TENSOR_SYMMETRIC) {
            case_file << "tensor symm per element: 1 " << r_variable_name << " " << base_filename << "." << r_variable_name << "." << asterisk_label << ".elem" << file_extension << "\n";
        } else {
            KRATOS_ERROR << "Unknown variable type for element data variable: " << r_variable_name << std::endl;
        }
    }

    // Element flags (always VariableType::SCALAR)
    for (const auto& r_flag_name_param : mOutputSettings["element_flags"]) {
        const std::string& r_flag_name = r_flag_name_param.GetString();
        case_file << "scalar per element: 1 " << r_flag_name << " " << base_filename << "." << r_flag_name << "." << asterisk_label << ".elem" << ".scl\n";
    }

    // Gauss point variables for elements
    for (const auto& r_variable_name_param : mOutputSettings["gauss_point_variables_in_elements"]) {
        const std::string& r_variable_name = r_variable_name_param.GetString();
        const auto variable_type = mVariableTypeMap[r_variable_name];
        const std::string file_extension = GetExtensionFile(variable_type);
        if (variable_type == VariableType::SCALAR) {
            case_file << "scalar per element: 1 " << r_variable_name << " " << base_filename << "." << r_variable_name << "." << asterisk_label << ".gp_elem" << file_extension << "\n";
        } else if (variable_type == VariableType::VECTOR) {
            case_file << "vector per element: 1 " << r_variable_name << " " << base_filename << "." << r_variable_name << "." << asterisk_label << ".gp_elem" << file_extension << "\n";
        } else if (variable_type == VariableType::TENSOR_SYMMETRIC) {
            case_file << "tensor symm per element: 1 " << r_variable_name << " " << base_filename << "." << r_variable_name << "." << asterisk_label << ".gp_elem" << file_extension << "\n";
        } else {
            KRATOS_ERROR << "Unknown variable type for element data variable: " << r_variable_name << std::endl;
        }
    }

    // Condition data variables
    for (const auto& r_variable_name_param : mOutputSettings["condition_data_value_variables"]) {
        const std::string& r_variable_name = r_variable_name_param.GetString();
        const auto variable_type = mVariableTypeMap[r_variable_name];
        const std::string file_extension = GetExtensionFile(variable_type);
        if (variable_type == VariableType::SCALAR) {
            case_file << "scalar per element: 1 " << r_variable_name << " " << base_filename << "." << r_variable_name << "." << asterisk_label << ".cond" << file_extension << "\n";
        } else if (variable_type == VariableType::VECTOR) {
            case_file << "vector per element: 1 " << r_variable_name << " " << base_filename << "." << r_variable_name << "." << asterisk_label << ".cond" << file_extension << "\n";
        } else if (variable_type == VariableType::TENSOR_SYMMETRIC) {
            case_file << "tensor symm per element: 1 " << r_variable_name << " " << base_filename << "." << r_variable_name << "." << asterisk_label << ".cond" << file_extension << "\n";
        } else {
            KRATOS_ERROR << "Unknown variable type for condition data variable: " << r_variable_name << std::endl;
        }
    }

    // Condition flags (always VariableType::SCALAR)
    for (const auto& r_flag_name_param : mOutputSettings["condition_flags"]) {
        const std::string& r_flag_name = r_flag_name_param.GetString();
        case_file << "scalar per element: 1 " << r_flag_name << " " << base_filename << "." << r_flag_name << "." << asterisk_label << ".cond" << ".scl\n";
    }

    case_file << "\n";

    // TIME Section
    case_file << "TIME\n";
    case_file << "time set: 1\n";
    case_file << "number of steps: " << mTimeValues.size() << "\n";
    case_file << "filename start number: 0\n"; // NOTE: Assuming steps start at 0
    case_file << "filename increment: 1\n";    // NOTE: Assuming steps increment by 1
    case_file << "time values: ";
    for(const double time : mTimeValues) {
        case_file << time << " ";
    }
    case_file << std::endl;

    case_file.close();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void EnSightOutput::WriteGeometryFile(const std::string& rFileName)
{
    KRATOS_TRY

    // Open the geometry file for writing
    const std::string step_label = InputOutputUtilities::GenerateStepLabel(mTimeValues.size(), mStepLabelPrecision);
    const std::string full_path = GetOutputFileName(rFileName, ".geo");
    std::ofstream geo_file(full_path, mFileFormat == FileFormat::BINARY ? std::ios::binary : std::ios::out);

    KRATOS_ERROR_IF_NOT(geo_file.is_open()) << "File \"" << full_path << "\" could not be opened!" << std::endl;
    if (mFileFormat == FileFormat::ASCII) {
        geo_file << std::scientific << std::setprecision(mDefaultPrecision); /// Set precision for ASCII output
    } else {
        geo_file << std::ios::binary;
    }

    // Write Header
    const bool is_ensight_6 = mEnSightFileFormat == EnSightFileFormat::EnSight6;
    if (is_ensight_6) {
        WriteString(geo_file, "EnSight 6 Geometry File");
        // Example EnSight6 geometry file writing:
        // This is the 1st description line of the EnSight6 geometry
        // This is the 2nd description line of the EnSight6 geometry
        // node id given
        // element id given
        // coordinates
        //         11
        //         15 4.00000e+00 0.00000e+00 0.00000e+00
        //         31 3.00000e+00 0.00000e+00 0.00000e+00
        //         20 5.00000e+00 0.00000e+00 0.00000e+00
        //         40 6.00000e+00 0.00000e+00 0.00000e+00
        //         22 5.00000e+00 1.00000e+00 0.00000e+00
        //         44 6.00000e+00 1.00000e+00 0.00000e+00
        //         55 6.00000e+00 3.00000e+00 0.00000e+00
        //         60 5.00000e+00 0.00000e+00 2.00000e+00
        //         61 6.00000e+00 0.00000e+00 2.00000e+00
        //         62 6.00000e+00 1.00000e+00 2.00000e+00
        //         63 5.00000e+00 1.00000e+00 2.00000e+00
        // part      1
        // 2D uns-elements  (description line for part 1)
        // tria3
        //         2
        //         102         15         20         22
        //         103         22         44         55
        // hexa8
        //         1
        //         104         20         40         44         22         60         61         62         63
        // part      2
        // 1D uns-elements  (description line for part 2)
        // bar2
        //         1
        //         101         31         15
    } else { // Otherwise we consider EnSight Gold format
        WriteString(geo_file, "EnSight Gold Geometry File");
        // This is the 1st description line of the EnSight Gold geometry
        // This is the 2nd description line of the EnSight Gold geometry
        // node id given
        // element id given
        // extents
        // 0.00000e+00 2.00000e+00
        // 0.00000e+00 2.00000e+00
        // 0.00000e+00 2.00000e+00
        // part
        //         1
        // 2D uns-elements (description line for part 1)
        // coordinates
        //     10
        //     15
        //     20
        //     40
        //     22
        //     44
        //     55
        //     60
        //     61
        //     62
        //     63
        // 4.00000e+00
        // 5.00000e+00
        // 6.00000e+00
        // 5.00000e+00
        // 6.00000e+00
        // 6.00000e+00
        // 5.00000e+00
        // 6.00000e+00
        // 6.00000e+00
        // 5.00000e+00
        // 0.00000e+00
        // 0.00000e+00
        // 0.00000e+00
        // 1.00000e+00
        // 1.00000e+00
        // 3.00000e+00
        // 0.00000e+00
        // 0.00000e+00
        // 1.00000e+00
        // 1.00000e+00
        // 0.00000e+00
        // 0.00000e+00
        // 0.00000e+00
        // 0.00000e+00
        // 0.00000e+00
        // 2.00000e+00
        // 2.00000e+00
        // 2.00000e+00
        // 2.00000e+00
        // 2.00000e+00
        // tria3
        //         2
        //     102
        //     103
        //         1         2         4
        //         4         5         6
        // hexa8
        //         1
        //     104
        //         2         3         5         4         7         8         9        10
        // part
        //         2
        // 1D uns-elements (description line for part 2)
        // coordinates
        //         2
        //         15
        //         31
        // 4.00000e+00
        // 3.00000e+00
        // 0.00000e+00
        // 0.00000e+00
        // 0.00000e+00
        // 0.00000e+00
        // bar2
        //         1
        //     101
        //         2         1
    }
    WriteString(geo_file, "Written by Kratos Multi-Physics");
    // Options:
    // node id <off/given/assign/ignore>
    // element id <off/given/assign/ignore>
    WriteString(geo_file, "node id given");    // TODO: Add option to use local IDs
    WriteString(geo_file, "element id given"); // TODO: Add option to use local IDs

    // // Define extents
    // extents
    // xmin xmax
    // ymin ymax
    // zmin zmax
    if (!is_ensight_6) {
        const auto bb = ComputeBoundingBox(mrModelPart, 1.0 + 1.0e-12);
        const auto& r_point_min = bb.GetMinPoint();
        const auto& r_point_max = bb.GetMaxPoint();

        WriteScalarData(geo_file, r_point_min.X(), false, true);
        WriteScalarData(geo_file, r_point_max.X());
        WriteScalarData(geo_file, r_point_min.Y(), false, true);
        WriteScalarData(geo_file, r_point_max.Y());
        WriteScalarData(geo_file, r_point_min.Z(), false, true);
        WriteScalarData(geo_file, r_point_max.Z());
    }

    // Coordinates for EnSight 6 are written before the parts
    if (is_ensight_6) {
        // Write coordinates
        WriteString(geo_file, "coordinates");
        WriteScalarData(geo_file, mrModelPart.NumberOfNodes());

        // Write IDs, then X, then Y, then Z blocks
        for (const auto& r_node : mrModelPart.Nodes()) {
            WriteScalarData(geo_file, r_node.Id(), false);
            WriteString(geo_file, " ", false);
            WriteScalarData(geo_file, r_node.X(),  false);
            WriteString(geo_file, " ", false);
            WriteScalarData(geo_file, r_node.Y(),  false);
            WriteString(geo_file, " ", false);
            WriteScalarData(geo_file, r_node.Z());
        }

        // // Write node coordinates // TODO: Add option to write with local IDs
        // for (const auto& r_node : mrModelPart.Nodes()) {
        //     WriteVectorData(geo_file, r_node.Coordinates());
        // }
    }

    // Iterate over the parts
    for (const auto& r_part_data : mPartDatas) {
        // Write part header
        if (is_ensight_6) {
            WriteString(geo_file, "part\t" + std::to_string(r_part_data.PartId));
        } else {
            WriteString(geo_file, "part");
            WriteString(geo_file, "\t\t", false);
            WriteScalarData(geo_file, r_part_data.PartId);
        }
        WriteString(geo_file, "Submodelpart\t" + r_part_data.PartName);

        // Write coordinates for EnSight Gold in each part
        if (mEnSightFileFormat == EnSightFileFormat::EnSightGold) {
            // Write coordinates
            WriteString(geo_file, "coordinates");
            WriteScalarData(geo_file, r_part_data.PartNodes.size(), false, true);

            // Write IDs, then X, then Y, then Z blocks
            for (const auto* p_node : r_part_data.PartNodes) {
                WriteScalarData(geo_file, p_node->Id(), false, true);
            }
            for (const auto* p_node : r_part_data.PartNodes) {
                WriteScalarData(geo_file, p_node->X());
            }
            for (const auto* p_node : r_part_data.PartNodes) {
                WriteScalarData(geo_file, p_node->Y());
            }
            for (const auto* p_node : r_part_data.PartNodes) {
                WriteScalarData(geo_file, p_node->Z());
            }
        }

        // Elements by Type
        std::vector<std::size_t> connectivity;
        for (const auto& it : r_part_data.PartGeometricalObjects) {
            const std::string& r_ensight_element_type = it.first;
            const auto& r_geometrical_objects = it.second;

            WriteString(geo_file, r_ensight_element_type);
            WriteString(geo_file, "\t", false);
            WriteScalarData(geo_file, r_geometrical_objects.size());

            for (const auto* p_geometrical_object : r_geometrical_objects) {
                WriteScalarData(geo_file, p_geometrical_object->Id(), !is_ensight_6, is_ensight_6);
                GetGeometryConnectivity(*p_geometrical_object, connectivity);
                for (const auto index : connectivity) {
                    // WriteScalarData(geo_file, r_part_data.KratosIdToLocalId.at(r_node.Id()));
                    WriteString(geo_file, "\t\t", false);
                    WriteScalarData(geo_file, index, false);
                }
                WriteString(geo_file, "");
            }
        }
    }
    geo_file.close();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void EnSightOutput::WriteNodalVariableToFile(
    const std::string& rVariableName,
    const bool IsHistorical,
    const std::string& rFileName
    )
{
    KRATOS_TRY

    // TODO: Use IndexPartition in tensors

    // Get the variable type and its extension
    const auto variable_type = mVariableTypeMap[rVariableName];
    const std::string file_extension = GetExtensionFile(variable_type);

    // Get the current step
    const std::string full_path = GetOutputFileName(rFileName, file_extension);
    std::ofstream var_file(full_path, mFileFormat == FileFormat::BINARY ? std::ios::binary : std::ios::out);

    KRATOS_ERROR_IF_NOT(var_file.is_open()) << "File \"" << full_path << "\" could not be opened!" << std::endl;
    if (mFileFormat == FileFormat::ASCII) {
        var_file << std::scientific << std::setprecision(mDefaultPrecision); // Set precision for ASCII output
    } else {
        var_file << std::ios::binary;
    }

    // Write variable information
    const bool is_ensight_gold = mEnSightFileFormat == EnSightFileFormat::EnSightGold;
    const std::string label = is_ensight_gold ? "EnSightGold" : "EnSight6";
    const std::string type_label = GetTypeLabel(variable_type);
    WriteString(var_file, "Per_node " + type_label + " values for the " + label + " geometry. Nodal variable: " + rVariableName);

    // For EnSightGold
    if (is_ensight_gold) {
        // Example for EnSight Gold scalar variable:
        // Per_node scalar values for the EnSight Gold geometry
        // part
        // 1
        // coordinates
        // 1.00000E+00
        // 3.00000E+00
        // 4.00000E+00
        // 5.00000E+00
        // 6.00000E+00
        // 7.00000E+00
        // 8.00000E+00
        // 9.00000E+00
        // 1.00000E+01
        // 1.10000E+01
        // part
        // 2
        // coordinates
        // 1.00000E+00
        // 2.00000E+00

        // Example for EnSight Gold vector variable:
        // ‘vx_m1 vx_m2 ... vx_mm’
        // ‘vy_m1 vy_m2 ... vy_mm’
        // ‘vz_m1 vz_m2 ... vz_mm’

        // Per_node vector values for the EnSight Gold geometry
        // part
        // 1
        // coordinates
        // 1.10000E+00
        // 3.10000E+00
        // 4.10000E+00
        // 5.10000E+00
        // 6.10000E+00
        // 7.10000E+00
        // 8.10000E+00
        // 9.10000E+00
        // 1.01000E+01
        // 1.11000E+01
        // 1.20000E+00
        // 3.20000E+00
        // 4.20000E+00
        // 5.20000E+00
        // 6.20000E+00
        // 7.20000E+00
        // 8.20000E+00
        // 9.20000E+00
        // 1.02000E+01
        // 1.12000E+01
        // 1.30000E+00
        // 3.30000E+00
        // 4.30000E+00
        // 5.30000E+00
        // 6.30000E+00
        // 7.30000E+00
        // 8.30000E+00
        // 9.30000E+00
        // 1.03000E+01
        // 1.13000E+01
        // part
        // 2
        // coordinates
        // 1.10000E+00
        // 2.10000E+00
        // 1.20000E+00
        // 2.20000E+00
        // 1.30000E+00
        // 2.30000E+00

        // Example for EnSight Gold symmetric tensor variable:
        // ‘v11_m1 v11_m2 ... v11_mm’
        // ‘v22_m1 v22_m2 ... v22_mm’
        // ‘v33_m1 v33_m2 ... v33_mm’
        // ‘v12_m1 v12_m2 ... v12_mm’
        // ‘v13_m1 v13_m2 ... v13_mm’
        // ‘v23_m1 v23_m2 ... v23_mm’

        // Per_node symmetric tensor values for the EnSight Gold geometry
        // part
        // 1
        // coordinates
        // 1.10000E+00
        // 3.10000E+00
        // 4.10000E+00
        // 5.10000E+00
        // 6.10000E+00
        // 7.10000E+00
        // 8.10000E+00
        // 9.10000E+00
        // 1.01000E+01
        // 1.11000E+01
        // 1.20000E+00
        // 3.20000E+00
        // 4.20000E+00
        // 5.20000E+00
        // 6.20000E+00
        // 7.20000E+00
        // 8.20000E+00
        // 9.20000E+00
        // 1.02000E+01
        // 1.12000E+01
        // 1.30000E+00
        // 3.30000E+00
        // 4.30000E+00
        // 5.30000E+00
        // 6.30000E+00
        // 7.30000E+00
        // 8.30000E+00
        // 9.30000E+00
        // 1.03000E+01
        // 1.13000E+01
        // 1.40000E+00
        // 3.40000E+00
        // 4.40000E+00
        // 5.40000E+00
        // 6.40000E+00
        // 7.40000E+00
        // 8.40000E+00
        // 9.40000E+00
        // 1.04000E+01
        // 1.14000E+01
        // 1.50000E+00
        // 3.50000E+00
        // 4.50000E+00
        // 5.50000E+00
        // 6.50000E+00
        // 7.50000E+00
        // 8.50000E+00
        // 9.50000E+00
        // 1.05000E+01
        // 1.15000E+01
        // 1.60000E+00
        // 3.60000E+00
        // 4.60000E+00
        // 5.60000E+00
        // 6.60000E+00
        // 7.60000E+00
        // 8.60000E+00
        // 9.60000E+00
        // 1.06000E+01
        // 1.16000E+01
        // part
        // 2
        // coordinates
        // 1.10000E+00
        // 2.10000E+00
        // 1.20000E+00
        // 2.20000E+00
        // 1.30000E+00
        // 2.30000E+00
        // 1.40000E+00
        // 2.40000E+00
        // 1.50000E+00
        // 2.50000E+00
        // 1.60000E+00
        // 2.60000E+00

        // Write historical data
        if (IsHistorical) {
            for (const auto& r_part_data : mPartDatas) {
                // Write part data
                WriteString(var_file, "part");
                WriteScalarData(var_file, r_part_data.PartId);
                WriteString(var_file, "coordinates");
                if (KratosComponents<Variable<bool>>::Has(rVariableName)) {
                    const auto& r_variable = KratosComponents<Variable<bool>>::Get(rVariableName);
                    for (const auto* p_node : r_part_data.PartNodes) {
                        WriteScalarData(var_file, p_node->FastGetSolutionStepValue(r_variable));
                    }
                } else if (KratosComponents<Variable<int>>::Has(rVariableName)) {
                    const auto& r_variable = KratosComponents<Variable<int>>::Get(rVariableName);
                    for (const auto* p_node : r_part_data.PartNodes) {
                        WriteScalarData(var_file, p_node->FastGetSolutionStepValue(r_variable));
                    }
                } else if (KratosComponents<Variable<double>>::Has(rVariableName)) {
                    const auto& r_variable = KratosComponents<Variable<double>>::Get(rVariableName);
                    for (const auto* p_node : r_part_data.PartNodes) {
                        WriteScalarData(var_file, p_node->FastGetSolutionStepValue(r_variable));
                    }
                } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)) {
                    const auto& r_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
                    std::array<std::vector<double>, 3> components;
                    const std::size_t number_of_nodes = r_part_data.PartNodes.size();
                    components[0].reserve(number_of_nodes);
                    components[1].reserve(number_of_nodes);
                    components[2].reserve(number_of_nodes);
                    for (const auto* p_node : r_part_data.PartNodes) {
                        const auto& r_value = p_node->FastGetSolutionStepValue(r_variable);
                        components[0].push_back(r_value[0]);
                        components[1].push_back(r_value[1]);
                        components[2].push_back(r_value[2]);
                    }
                    for (unsigned int i = 0; i < 3; ++i) {
                        for (std::size_t j = 0; j < number_of_nodes; ++j) {
                            WriteScalarData(var_file, components[i][j]);
                        }
                    }
                } else if (KratosComponents<Variable<Vector>>::Has(rVariableName)) {
                    const auto& r_variable = KratosComponents<Variable<Vector>>::Get(rVariableName);
                    if (variable_type == VariableType::VECTOR) {
                        std::array<std::vector<double>, 3> components;
                        const std::size_t number_of_nodes = r_part_data.PartNodes.size();
                        components[0].reserve(number_of_nodes);
                        components[1].reserve(number_of_nodes);
                        components[2].reserve(number_of_nodes);
                        for (const auto* p_node : r_part_data.PartNodes) {
                            const auto& r_value = p_node->FastGetSolutionStepValue(r_variable);
                            components[0].push_back(r_value[0]);
                            components[1].push_back(r_value[1]);
                            components[2].push_back(r_value[2]);
                        }
                        for (unsigned int i = 0; i < 3; ++i) {
                            for (std::size_t j = 0; j < number_of_nodes; ++j) {
                                WriteScalarData(var_file, components[i][j]);
                            }
                        }
                    } else if (variable_type == VariableType::TENSOR_SYMMETRIC) {
                        std::array<std::vector<double>, 6> components;
                        const std::size_t number_of_nodes = r_part_data.PartNodes.size();
                        components[0].reserve(number_of_nodes);
                        components[1].reserve(number_of_nodes);
                        components[2].reserve(number_of_nodes);
                        components[3].reserve(number_of_nodes);
                        components[4].reserve(number_of_nodes);
                        components[5].reserve(number_of_nodes);
                        for (const auto* p_node : r_part_data.PartNodes) {
                            const auto& r_value = p_node->FastGetSolutionStepValue(r_variable);
                            components[0].push_back(r_value[0]);
                            components[1].push_back(r_value[1]);
                            components[2].push_back(r_value[2]);
                            components[3].push_back(r_value[3]);
                            components[4].push_back(r_value[5]);
                            components[5].push_back(r_value[4]);
                        }
                        for (unsigned int i = 0; i < 6; ++i) {
                            for (std::size_t j = 0; j < number_of_nodes; ++j) {
                                WriteScalarData(var_file, components[i][j]);
                            }
                        }
                    } else {
                        KRATOS_ERROR << "Unknown variable type for vector: " << rVariableName << std::endl;
                    }
                } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)) {
                    const auto& r_variable = KratosComponents<Variable<array_1d<double, 6>>>::Get(rVariableName);
                    std::array<std::vector<double>, 6> components;
                    const std::size_t number_of_nodes = r_part_data.PartNodes.size();
                    components[0].reserve(number_of_nodes);
                    components[1].reserve(number_of_nodes);
                    components[2].reserve(number_of_nodes);
                    components[3].reserve(number_of_nodes);
                    components[4].reserve(number_of_nodes);
                    components[5].reserve(number_of_nodes);
                    for (const auto* p_node : r_part_data.PartNodes) {
                        const auto& r_value = p_node->FastGetSolutionStepValue(r_variable);
                        components[0].push_back(r_value[0]);
                        components[1].push_back(r_value[1]);
                        components[2].push_back(r_value[2]);
                        components[3].push_back(r_value[3]);
                        components[4].push_back(r_value[5]);
                        components[5].push_back(r_value[4]);
                    }
                    for (unsigned int i = 0; i < 6; ++i) {
                        for (std::size_t j = 0; j < number_of_nodes; ++j) {
                            WriteScalarData(var_file, components[i][j]);
                        }
                    }
                } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(rVariableName)) {
                    // TODO: Implement asymmetric tensor.
                    // const auto& r_variable = KratosComponents<Variable<array_1d<double, 9>>>::Get(rVariableName);
                    KRATOS_ERROR << "Asymmetric tensor output is not implemented yet." << std::endl;
                } else if (KratosComponents<Variable<Matrix>>::Has(rVariableName)) {
                    const auto& r_variable = KratosComponents<Variable<Matrix>>::Get(rVariableName);
                    if (variable_type == VariableType::TENSOR_SYMMETRIC) {
                        std::array<std::vector<double>, 6> components;
                        const std::size_t number_of_nodes = r_part_data.PartNodes.size();
                        components[0].reserve(number_of_nodes);
                        components[1].reserve(number_of_nodes);
                        components[2].reserve(number_of_nodes);
                        components[3].reserve(number_of_nodes);
                        components[4].reserve(number_of_nodes);
                        components[5].reserve(number_of_nodes);
                        for (const auto* p_node : r_part_data.PartNodes) {
                            const auto& r_value = p_node->FastGetSolutionStepValue(r_variable);
                            components[0].push_back(r_value(0, 0));
                            components[1].push_back(r_value(1, 1));
                            components[2].push_back(r_value(2, 2));
                            components[3].push_back(r_value(0, 1));
                            components[4].push_back(r_value(0, 2));
                            components[5].push_back(r_value(1, 2));
                        }
                        for (unsigned int i = 0; i < 6; ++i) {
                            for (std::size_t j = 0; j < number_of_nodes; ++j) {
                                WriteScalarData(var_file, components[i][j]);
                            }
                        }
                    } else if (variable_type == VariableType::TENSOR_ASYMMETRIC) {
                        // TODO: Implement asymmetric tensor.
                        KRATOS_ERROR << "Asymmetric tensor output is not implemented yet." << std::endl;
                    } else {
                        KRATOS_ERROR << "Unknown variable type for vector: " << rVariableName << std::endl;
                    }
                }
                WriteString(var_file, "");
            }
        } else { // Non-historical data
            for (const auto& r_part_data : mPartDatas) {
                // Write part data
                WriteString(var_file, "part");
                WriteScalarData(var_file, r_part_data.PartId);
                WriteString(var_file, "coordinates");
                if (KratosComponents<Variable<bool>>::Has(rVariableName)) {
                    const auto& r_variable = KratosComponents<Variable<bool>>::Get(rVariableName);
                    for (const auto* p_node : r_part_data.PartNodes) {
                        WriteScalarData(var_file, p_node->GetValue(r_variable));
                    }
                } else if (KratosComponents<Variable<int>>::Has(rVariableName)) {
                    const auto& r_variable = KratosComponents<Variable<int>>::Get(rVariableName);
                    for (const auto* p_node : r_part_data.PartNodes) {
                        WriteScalarData(var_file, p_node->GetValue(r_variable));
                    }
                } else if (KratosComponents<Variable<double>>::Has(rVariableName)) {
                    const auto& r_variable = KratosComponents<Variable<double>>::Get(rVariableName);
                    for (const auto* p_node : r_part_data.PartNodes) {
                        WriteScalarData(var_file, p_node->GetValue(r_variable));
                    }
                } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)) {
                    const auto& r_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
                    std::array<std::vector<double>, 3> components;
                    const std::size_t number_of_nodes = r_part_data.PartNodes.size();
                    components[0].reserve(number_of_nodes);
                    components[1].reserve(number_of_nodes);
                    components[2].reserve(number_of_nodes);
                    for (const auto* p_node : r_part_data.PartNodes) {
                        const auto& r_value = p_node->GetValue(r_variable);
                        components[0].push_back(r_value[0]);
                        components[1].push_back(r_value[1]);
                        components[2].push_back(r_value[2]);
                    }
                    for (unsigned int i = 0; i < 3; ++i) {
                        for (std::size_t j = 0; j < number_of_nodes; ++j) {
                            WriteScalarData(var_file, components[i][j]);
                        }
                    }
                } else if (KratosComponents<Variable<Vector>>::Has(rVariableName)) {
                    const auto& r_variable = KratosComponents<Variable<Vector>>::Get(rVariableName);
                    if (variable_type == VariableType::VECTOR) {
                        std::array<std::vector<double>, 3> components;
                        const std::size_t number_of_nodes = r_part_data.PartNodes.size();
                        components[0].reserve(number_of_nodes);
                        components[1].reserve(number_of_nodes);
                        components[2].reserve(number_of_nodes);
                        for (const auto* p_node : r_part_data.PartNodes) {
                            const auto& r_value = p_node->GetValue(r_variable);
                            components[0].push_back(r_value[0]);
                            components[1].push_back(r_value[1]);
                            components[2].push_back(r_value[2]);
                        }
                        for (unsigned int i = 0; i < 3; ++i) {
                            for (std::size_t j = 0; j < number_of_nodes; ++j) {
                                WriteScalarData(var_file, components[i][j]);
                            }
                        }
                    } else if (variable_type == VariableType::TENSOR_SYMMETRIC) {
                        std::array<std::vector<double>, 6> components;
                        const std::size_t number_of_nodes = r_part_data.PartNodes.size();
                        components[0].reserve(number_of_nodes);
                        components[1].reserve(number_of_nodes);
                        components[2].reserve(number_of_nodes);
                        components[3].reserve(number_of_nodes);
                        components[4].reserve(number_of_nodes);
                        components[5].reserve(number_of_nodes);
                        for (const auto* p_node : r_part_data.PartNodes) {
                            const auto& r_value = p_node->GetValue(r_variable);
                            components[0].push_back(r_value[0]);
                            components[1].push_back(r_value[1]);
                            components[2].push_back(r_value[2]);
                            components[3].push_back(r_value[3]);
                            components[4].push_back(r_value[5]);
                            components[5].push_back(r_value[4]);
                        }
                        for (unsigned int i = 0; i < 6; ++i) {
                            for (std::size_t j = 0; j < number_of_nodes; ++j) {
                                WriteScalarData(var_file, components[i][j]);
                            }
                        }
                    } else {
                        KRATOS_ERROR << "Unknown variable type for vector: " << rVariableName << std::endl;
                    }
                } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)) {
                    const auto& r_variable = KratosComponents<Variable<array_1d<double, 6>>>::Get(rVariableName);
                    std::array<std::vector<double>, 6> components;
                    const std::size_t number_of_nodes = r_part_data.PartNodes.size();
                    components[0].reserve(number_of_nodes);
                    components[1].reserve(number_of_nodes);
                    components[2].reserve(number_of_nodes);
                    components[3].reserve(number_of_nodes);
                    components[4].reserve(number_of_nodes);
                    components[5].reserve(number_of_nodes);
                    for (const auto* p_node : r_part_data.PartNodes) {
                        const auto& r_value = p_node->GetValue(r_variable);
                        components[0].push_back(r_value[0]);
                        components[1].push_back(r_value[1]);
                        components[2].push_back(r_value[2]);
                        components[3].push_back(r_value[3]);
                        components[4].push_back(r_value[5]);
                        components[5].push_back(r_value[4]);
                    }
                    for (unsigned int i = 0; i < 6; ++i) {
                        for (std::size_t j = 0; j < number_of_nodes; ++j) {
                            WriteScalarData(var_file, components[i][j]);
                        }
                    }
                } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(rVariableName)) {
                    // TODO: Implement asymmetric tensor.
                    // const auto& r_variable = KratosComponents<Variable<array_1d<double, 9>>>::Get(rVariableName);
                    KRATOS_ERROR << "Asymmetric tensor output is not implemented yet." << std::endl;
                } else if (KratosComponents<Variable<Matrix>>::Has(rVariableName)) {
                    const auto& r_variable = KratosComponents<Variable<Matrix>>::Get(rVariableName);
                    if (variable_type == VariableType::TENSOR_SYMMETRIC) {
                        std::array<std::vector<double>, 6> components;
                        const std::size_t number_of_nodes = r_part_data.PartNodes.size();
                        components[0].reserve(number_of_nodes);
                        components[1].reserve(number_of_nodes);
                        components[2].reserve(number_of_nodes);
                        components[3].reserve(number_of_nodes);
                        components[4].reserve(number_of_nodes);
                        components[5].reserve(number_of_nodes);
                        for (const auto* p_node : r_part_data.PartNodes) {
                            const auto& r_value = p_node->GetValue(r_variable);
                            components[0].push_back(r_value(0, 0));
                            components[1].push_back(r_value(1, 1));
                            components[2].push_back(r_value(2, 2));
                            components[3].push_back(r_value(0, 1));
                            components[4].push_back(r_value(0, 2));
                            components[5].push_back(r_value(1, 2));
                        }
                        for (unsigned int i = 0; i < 6; ++i) {
                            for (std::size_t j = 0; j < number_of_nodes; ++j) {
                                WriteScalarData(var_file, components[i][j]);
                            }
                        }
                    } else if (variable_type == VariableType::TENSOR_ASYMMETRIC) {
                        // TODO: Implement asymmetric tensor.
                        KRATOS_ERROR << "Asymmetric tensor output is not implemented yet." << std::endl;
                    } else {
                        KRATOS_ERROR << "Unknown variable type for vector: " << rVariableName << std::endl;
                    }
                }
                WriteString(var_file, "");
            }
        }
    } else { // For Ensight 6
        // Example for EnSight6 scalar variable:
        // Per_node scalar values for the EnSight6 geometry
        // 1.00000E+00 2.00000E+00 3.00000E+00 4.00000E+00 5.00000E+00 6.00000E+00
        // 7.00000E+00 8.00000E+00 9.00000E+00 1.00000E+01 1.10000E+01

        // Example for EnSight6 vector variable:
        // Per_node vector values for the EnSight6 geometry
        // 1.10000E+00 1.20000E+00 1.30000E+00 2.10000E+00 2.20000E+00 2.30000E+00
        // 3.10000E+00 3.20000E+00 3.30000E+00 4.10000E+00 4.20000E+00 4.30000E+00
        // 5.10000E+00 5.20000E+00 5.30000E+00 6.10000E+00 6.20000E+00 6.30000E+00
        // 7.10000E+00 7.20000E+00 7.30000E+00 8.10000E+00 8.20000E+00 8.30000E+00
        // 9.l0000E+00 9.20000E+00 9.30000E+00 1.01000E+01 1.02000E+01 1.03000E+01
        // 1.11000E+01 1.12000E+01 1.13000E+01

        // Example for EnSight6 symmetric tensor variable:
        // Per_node symmetric tensor values for the EnSight6 geometry
        // 1.10000E+00 1.20000E+00 1.30000E+00 1.40000E+00 1.50000E+00 1.60000E+00
        // 2.10000E+00 2.20000E+00 2.30000E+00 2.40000E+00 2.50000E+00 2.60000E+00
        // 3.10000E+00 3.20000E+00 3.30000E+00 3.40000E+00 3.50000E+00 3.60000E+00
        // 4.10000E+00 4.20000E+00 4.30000E+00 4.40000E+00 4.50000E+00 4.60000E+00
        // 5.10000E+00 5.20000E+00 5.30000E+00 5.40000E+00 5.50000E+00 5.60000E+00
        // 6.10000E+00 6.20000E+00 6.30000E+00 6.40000E+00 6.50000E+00 6.60000E+00
        // 7.10000E+00 7.20000E+00 7.30000E+00 7.40000E+00 7.50000E+00 7.60000E+00
        // 8.10000E+00 8.20000E+00 8.30000E+00 8.40000E+00 8.50000E+00 8.60000E+00
        // 9.10000E+00 9.20000E+00 9.30000E+00 9.40000E+00 9.50000E+00 9.60000E+00
        // 1.01000E+01 1.02000E+01 1.03000E+01 1.04000E+01 1.05000E+01 1.06000E+01
        // 1.11000E+01 1.12000E+01 1.13000E+01 1.14000E+01 1.15000E+01 1.16000E+01

        // Write lambda that check that the counter is 6 for ensight 6
        auto check_counter = [](unsigned int& rCounter) -> bool {
            const bool check = rCounter == 6;
            if (check) {
                rCounter = 0;
            }
            return check;
        };

        // Write node data
        unsigned int counter = 0;
        bool new_line = true;
        if (IsHistorical) {
            // Write historical data
            if (KratosComponents<Variable<bool>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<bool>>::Get(rVariableName);
                for (const auto& r_node : mrModelPart.Nodes()) {
                    ++counter;
                    new_line = check_counter(counter);
                    WriteScalarData(var_file, r_node.FastGetSolutionStepValue(r_variable), new_line, false, !new_line);
                }
            } else if (KratosComponents<Variable<int>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<int>>::Get(rVariableName);
                for (const auto& r_node : mrModelPart.Nodes()) {
                    ++counter;
                    new_line = check_counter(counter);
                    WriteScalarData(var_file, r_node.FastGetSolutionStepValue(r_variable), new_line, false, !new_line);
                }
            } else if (KratosComponents<Variable<double>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<double>>::Get(rVariableName);
                for (const auto& r_node : mrModelPart.Nodes()) {
                    ++counter;
                    new_line = check_counter(counter);
                    WriteScalarData(var_file, r_node.FastGetSolutionStepValue(r_variable), new_line, false, !new_line);
                }
            } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
                for (const auto& r_node : mrModelPart.Nodes()) {
                    counter += 3;
                    new_line = check_counter(counter);
                    WriteVectorData(var_file, r_node.FastGetSolutionStepValue(r_variable), new_line, false, !new_line);
                }
            } else if (KratosComponents<Variable<Vector>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<Vector>>::Get(rVariableName);
                if (variable_type == VariableType::VECTOR) {
                    for (const auto& r_node : mrModelPart.Nodes()) {
                        counter += 3;
                        new_line = check_counter(counter);
                        WriteVectorData(var_file, r_node.FastGetSolutionStepValue(r_variable), new_line, false, !new_line);
                    }
                } else if (variable_type == VariableType::TENSOR_SYMMETRIC) {
                    for (const auto& r_node : mrModelPart.Nodes()) {
                        WriteSymmetricTensorData(var_file, r_node.FastGetSolutionStepValue(r_variable));
                    }
                } else {
                    KRATOS_ERROR << "Unknown variable type for vector: " << rVariableName << std::endl;
                }
            } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<array_1d<double, 6>>>::Get(rVariableName);
                for (const auto& r_node : mrModelPart.Nodes()) {
                    WriteSymmetricTensorData(var_file, r_node.FastGetSolutionStepValue(r_variable));
                }
            } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(rVariableName)) {
                // TODO: Implement asymmetric tensor.
                // const auto& r_variable = KratosComponents<Variable<array_1d<double, 9>>>::Get(rVariableName);
                KRATOS_ERROR << "Asymmetric tensor output is not implemented yet." << std::endl;
            } else if (KratosComponents<Variable<Matrix>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<Matrix>>::Get(rVariableName);
                if (variable_type == VariableType::TENSOR_SYMMETRIC) {
                    for (const auto& r_node : mrModelPart.Nodes()) {
                        WriteSymmetricTensorData(var_file, r_node.FastGetSolutionStepValue(r_variable));
                    }
                } else if (variable_type == VariableType::TENSOR_ASYMMETRIC) {
                    // TODO: Implement asymmetric tensor.
                    KRATOS_ERROR << "Asymmetric tensor output is not implemented yet." << std::endl;
                } else {
                    KRATOS_ERROR << "Unknown variable type for vector: " << rVariableName << std::endl;
                }
            }
        } else { // Non-historical data
            if (KratosComponents<Variable<bool>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<bool>>::Get(rVariableName);
                for (const auto& r_node : mrModelPart.Nodes()) {
                    ++counter;
                    new_line = check_counter(counter);
                    WriteScalarData(var_file, r_node.GetValue(r_variable), new_line, false, !new_line);
                }
            } else if (KratosComponents<Variable<int>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<int>>::Get(rVariableName);
                for (const auto& r_node : mrModelPart.Nodes()) {
                    ++counter;
                    new_line = check_counter(counter);
                    WriteScalarData(var_file, r_node.GetValue(r_variable), new_line, false, !new_line);
                }
            } else if (KratosComponents<Variable<double>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<double>>::Get(rVariableName);
                for (const auto& r_node : mrModelPart.Nodes()) {
                    ++counter;
                    new_line = check_counter(counter);
                    WriteScalarData(var_file, r_node.GetValue(r_variable), new_line, false, !new_line);
                }
            } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
                for (const auto& r_node : mrModelPart.Nodes()) {
                    counter += 3;
                    new_line = check_counter(counter);
                    WriteVectorData(var_file, r_node.GetValue(r_variable), new_line, false, !new_line);
                }
            } else if (KratosComponents<Variable<Vector>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<Vector>>::Get(rVariableName);
                if (variable_type == VariableType::VECTOR) {
                    for (const auto& r_node : mrModelPart.Nodes()) {
                        counter += 3;
                        new_line = check_counter(counter);
                        WriteVectorData(var_file, r_node.GetValue(r_variable), new_line, false, !new_line);
                    }
                } else if (variable_type == VariableType::TENSOR_SYMMETRIC) {
                    for (const auto& r_node : mrModelPart.Nodes()) {
                        WriteSymmetricTensorData(var_file, r_node.GetValue(r_variable));
                    }
                } else {
                    KRATOS_ERROR << "Unknown variable type for vector: " << rVariableName << std::endl;
                }
            } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<array_1d<double, 6>>>::Get(rVariableName);
                for (const auto& r_node : mrModelPart.Nodes()) {
                    WriteSymmetricTensorData(var_file, r_node.GetValue(r_variable));
                }
            } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(rVariableName)) {
                // TODO: Implement asymmetric tensor.
                // const auto& r_variable = KratosComponents<Variable<array_1d<double, 9>>>::Get(rVariableName);
                KRATOS_ERROR << "Asymmetric tensor output is not implemented yet." << std::endl;
            } else if (KratosComponents<Variable<Matrix>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<Matrix>>::Get(rVariableName);
                if (variable_type == VariableType::TENSOR_SYMMETRIC) {
                    for (const auto& r_node : mrModelPart.Nodes()) {
                        WriteSymmetricTensorData(var_file, r_node.GetValue(r_variable));
                    }
                } else if (variable_type == VariableType::TENSOR_ASYMMETRIC) {
                    // TODO: Implement asymmetric tensor.
                    KRATOS_ERROR << "Asymmetric tensor output is not implemented yet." << std::endl;
                } else {
                    KRATOS_ERROR << "Unknown variable type for vector: " << rVariableName << std::endl;
                }
            }
        }
    }

    // Close the variable file
    var_file.close();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void EnSightOutput::WriteNodalFlagToFile(
    const std::string& rFlagName,
    const std::string& rFileName
    )
{
    KRATOS_TRY

    // Get the variable type and its extension
    const auto variable_type = VariableType::SCALAR;
    const std::string file_extension = GetExtensionFile(variable_type);

    // Get the current step
    const std::string full_path = GetOutputFileName(rFileName, file_extension);
    std::ofstream var_file(full_path, mFileFormat == FileFormat::BINARY ? std::ios::binary : std::ios::out);

    KRATOS_ERROR_IF_NOT(var_file.is_open()) << "File \"" << full_path << "\" could not be opened!" << std::endl;
    if (mFileFormat == FileFormat::ASCII) {
        var_file << std::scientific << std::setprecision(mDefaultPrecision); // Set precision for ASCII output
    } else {
        var_file << std::ios::binary;
    }

    // Write variable information
    const bool is_ensight_gold = mEnSightFileFormat == EnSightFileFormat::EnSightGold;
    const std::string label = is_ensight_gold ? "EnSightGold" : "EnSight6";
    const std::string type_label = GetTypeLabel(variable_type);
    WriteString(var_file, "Per_node " + type_label + " values for the " + label + " geometry. Nodal variable: " + rFlagName);

    // Get flag
    const auto& r_flag = KratosComponents<Flags>::Get(rFlagName);

    // For EnSightGold
    if (is_ensight_gold) {
        for (const auto& r_part_data : mPartDatas) {
            // Write part data
            WriteString(var_file, "part");
            WriteScalarData(var_file, r_part_data.PartId);
            WriteString(var_file, "coordinates");
            for (const auto* p_node : r_part_data.PartNodes) {
                const int flag = p_node->IsDefined(r_flag) ? static_cast<int>(p_node->Is(r_flag)) : -1; // Default to -1 if not defined
                WriteScalarData(var_file, flag);
            }
            WriteString(var_file, "");
        }
    } else {
        // Write lambda that check that the counter is 6 for ensight 6
        auto check_counter = [](unsigned int& rCounter) -> bool {
            const bool check = rCounter == 6;
            if (check) {
                rCounter = 0;
            }
            return check;
        };

        // Write node data
        unsigned int counter = 0;
        bool new_line = true;
        for (const auto& r_node : mrModelPart.Nodes()) {
            ++counter;
            new_line = check_counter(counter);
            const int flag = r_node.IsDefined(r_flag) ? static_cast<int>(r_node.Is(r_flag)) : -1; // Default to -1 if not defined
            WriteScalarData(var_file, flag, new_line, false, !new_line);
        }
    }

    // Close the variable file
    var_file.close();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void EnSightOutput::WriteGeometricalVariableToFile(
    const std::string& rVariableName,
    const std::string& rFileName,
    const bool IsElement
    )
{
    KRATOS_TRY

    // TODO: Use IndexPartition in tensors

    // Get the variable type and its extension
    const auto variable_type = mVariableTypeMap[rVariableName];
    const std::string file_extension = GetExtensionFile(variable_type);

    // Open the file for writing
    const std::string full_path = GetOutputFileName(rFileName, file_extension);
    std::ofstream var_file(full_path, mFileFormat == FileFormat::BINARY ? std::ios::binary : std::ios::out);

    KRATOS_ERROR_IF_NOT(var_file.is_open()) << "File \"" << full_path << "\" could not be opened!" << std::endl;
    if (mFileFormat == FileFormat::ASCII) {
        var_file << std::scientific << std::setprecision(mDefaultPrecision); // Set precision for ASCII output
    } else {
        var_file << std::ios::binary;
    }

    // Write variable information
    const bool is_ensight_gold = mEnSightFileFormat == EnSightFileFormat::EnSightGold;
    const bool is_ensight_6 = mEnSightFileFormat == EnSightFileFormat::EnSight6;
    const std::string label = is_ensight_gold ? "EnSightGold" : "EnSight6";
    const std::string type_label = GetTypeLabel(variable_type);
    const std::string entity_label = IsElement ? "Elemental" : "Conditional";
    WriteString(var_file, "Per_elem " + type_label + " values for the " + label + " geometry. " + entity_label + " variable: " + rVariableName);

    /* EnSight6 */

    // Example for EnSight6 scalar variable:
    // Per_elem scalar values for the EnSight6 geometry
    // part 1
    // tria3
    //     2.00000E+00 3.00000E+00
    // hexa8
    //     4.00000E+00
    // part 2
    // bar2
    //     1.00000E+00

    // Example for EnSight6 vector variable:
    // Per_elem vector values for the EnSight6 geometry
    // part 1
    // tria3
    //     2.10000E+00 2.20000E+00 2.30000E+00 3.10000E+00 3.20000E+00 3.30000E+00
    // hexa8
    //     4.10000E+00 4.20000E+00 4.30000E+00
    // part 2
    // bar2
    //     1.10000E+00 1.20000E+00 1.30000E+00

    // Example for EnSight6 symmetric tensor variable:
    // Per_elem symmetric tensor values for the EnSight6 geometry
    // part 1
    // tria3
    //     2.10000E+00 2.20000E+00 2.30000E+00 2.40000E+00 2.50000E+00 2.60000E+00
    //     3.10000E+00 3.20000E+00 3.30000E+00 3.40000E+00 3.50000E+00 3.60000E+00
    // hexa8
    //     4.10000E+00 4.20000E+00 4.30000E+00 4.40000E+00 4.50000E+00 4.60000E+00
    // part 2
    // bar2
    //     1.10000E+00 1.20000E+00 1.30000E+00 1.40000E+00 1.50000E+00 1.60000E+00

    /* EnSight Gold */

    // Example for EnSight Goldscalar variable:
    // Per_elem scalar values for the EnSight Gold geometry
    // part
    // 1
    // tria3
    // 2.00000E+00
    // 3.00000E+00
    // hexa8
    // 4.00000E+00
    // part
    // 2
    // bar2
    // 1.00000E+00

    // Example for EnSight Gold vector variable:
    // Per_elem vector values for the EnSight Gold geometry
    // part
    // 1
    // tria3
    // 2.10000E+00
    // 3.10000E+00
    // 2.20000E+00
    // 3.20000E+00
    // 2.30000E+00
    // 3.30000E+00
    // hexa8
    // 4.10000E+00
    // 4.20000E+00
    // 4.30000E+00
    // part
    // 2
    // bar2
    // 1.10000E+00
    // 1.20000E+00
    // 1.30000E+00

    // Example for EnSight Gold symmetric tensor variable:
    // Per_elem symmetric tensor values for the EnSight Gold geometry
    // part
    // 1
    // tria3
    // 2.10000E+00
    // 3.10000E+00
    // 2.20000E+00
    // 3.20000E+00
    // 2.30000E+00
    // 3.30000E+00
    // 2.40000E+00
    // 3.40000E+00
    // 2.50000E+00
    // 3.50000E+00
    // 2.60000E+00
    // 3.60000E+00
    // hexa8
    // 4.10000E+00
    // 4.20000E+00
    // 4.30000E+00
    // 4.40000E+00
    // 4.50000E+00
    // 4.60000E+00
    // part
    // 2
    // bar2
    // 1.10000E+00
    // 1.20000E+00
    // 1.30000E+00
    // 1.40000E+00
    // 1.50000E+00
    // 1.60000E+00

    // Write lambda that check that the counter is 6 for ensight 6 and gold
    auto check_counter_ensight_6 = [](unsigned int& rCounter) -> bool {
        const bool check = rCounter == 6;
        if (check) {
            rCounter = 0;
        }
        return check;
    };
    auto check_counter_ensight_gold = [](unsigned int& rCounter) -> bool {
        (void)rCounter; // suppress unused warning
        return true;
    };
    std::function<bool(unsigned int&)> check_counter = is_ensight_gold
        ? std::function<bool(unsigned int&)>(check_counter_ensight_gold)
        : std::function<bool(unsigned int&)>(check_counter_ensight_6);

    // Write the geometrical data for the parts
    for (const auto& r_part_data : mPartDatas) {
        // Check if the part is of the correct type (element or condition)
        if (r_part_data.PartElements != IsElement) continue;

        // Write part data
        bool new_line = true;
        for (const auto& it : r_part_data.PartGeometricalObjects) {
            unsigned int counter = 0;
            const std::string& r_ensight_element_type = it.first;
            const auto& r_geometrical_objects = it.second;

            // Write part header
            if (is_ensight_6) {
                WriteString(var_file, "part\t" + std::to_string(r_part_data.PartId));
            } else { // Ensight Gold
                WriteString(var_file, "part");
                WriteScalarData(var_file, r_part_data.PartId);
            }
            WriteString(var_file, r_ensight_element_type);

            if (KratosComponents<Variable<bool>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<bool>>::Get(rVariableName);
                for (const auto* p_geometrical_object : r_geometrical_objects) {
                    ++counter;
                    new_line = check_counter(counter);
                    WriteScalarData(var_file, p_geometrical_object->GetValue(r_variable), new_line, !is_ensight_gold && counter == 1, !new_line);
                }
            } else if (KratosComponents<Variable<int>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<int>>::Get(rVariableName);
                for (const auto* p_geometrical_object : r_geometrical_objects) {
                    ++counter;
                    new_line = check_counter(counter);
                    WriteScalarData(var_file, p_geometrical_object->GetValue(r_variable), new_line, !is_ensight_gold && counter == 1, !new_line);
                }
            } else if (KratosComponents<Variable<double>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<double>>::Get(rVariableName);
                for (const auto* p_geometrical_object : r_geometrical_objects) {
                    ++counter;
                    new_line = check_counter(counter);
                    WriteScalarData(var_file, p_geometrical_object->GetValue(r_variable), new_line, !is_ensight_gold && counter == 1, !new_line);
                }
            } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
                if (is_ensight_gold) {
                    std::array<std::vector<double>, 3> components;
                    const std::size_t number_of_geometries = r_geometrical_objects.size();
                    components[0].reserve(number_of_geometries);
                    components[1].reserve(number_of_geometries);
                    components[2].reserve(number_of_geometries);
                    for (const auto* p_geometrical_object : r_geometrical_objects) {
                        const auto& r_value = p_geometrical_object->GetValue(r_variable);
                        components[0].push_back(r_value[0]);
                        components[1].push_back(r_value[1]);
                        components[2].push_back(r_value[2]);
                    }
                    for (unsigned int i = 0; i < 3; ++i) {
                        for (std::size_t j = 0; j < number_of_geometries; ++j) {
                            WriteScalarData(var_file, components[i][j]);
                        }
                    }
                } else {
                    for (const auto* p_geometrical_object : r_geometrical_objects) {
                        counter += 3;
                        new_line = check_counter(counter);
                        WriteVectorData(var_file, p_geometrical_object->GetValue(r_variable), new_line, counter == 1, !new_line);
                    }
                }
            } else if (KratosComponents<Variable<Vector>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<Vector>>::Get(rVariableName);
                if (is_ensight_gold) {
                    if (variable_type == VariableType::VECTOR) {
                        std::array<std::vector<double>, 3> components;
                        const std::size_t number_of_geometries = r_geometrical_objects.size();
                        components[0].reserve(number_of_geometries);
                        components[1].reserve(number_of_geometries);
                        components[2].reserve(number_of_geometries);
                        for (const auto* p_geometrical_object : r_geometrical_objects) {
                            const auto& r_value = p_geometrical_object->GetValue(r_variable);
                            components[0].push_back(r_value[0]);
                            components[1].push_back(r_value[1]);
                            components[2].push_back(r_value[2]);
                        }
                        for (unsigned int i = 0; i < 3; ++i) {
                            for (std::size_t j = 0; j < number_of_geometries; ++j) {
                                WriteScalarData(var_file, components[i][j]);
                            }
                        }
                    } else if (variable_type == VariableType::TENSOR_SYMMETRIC) {
                        std::array<std::vector<double>, 6> components;
                        const std::size_t number_of_geometries = r_geometrical_objects.size();
                        components[0].reserve(number_of_geometries);
                        components[1].reserve(number_of_geometries);
                        components[2].reserve(number_of_geometries);
                        components[3].reserve(number_of_geometries);
                        components[4].reserve(number_of_geometries);
                        components[5].reserve(number_of_geometries);
                        for (const auto* p_geometrical_object : r_geometrical_objects) {
                            const auto& r_value = p_geometrical_object->GetValue(r_variable);
                            components[0].push_back(r_value[0]);
                            components[1].push_back(r_value[1]);
                            components[2].push_back(r_value[2]);
                            components[3].push_back(r_value[3]);
                            components[4].push_back(r_value[5]);
                            components[5].push_back(r_value[4]);
                        }
                        for (unsigned int i = 0; i < 6; ++i) {
                            for (std::size_t j = 0; j < number_of_geometries; ++j) {
                                WriteScalarData(var_file, components[i][j]);
                            }
                        }
                    } else {
                        KRATOS_ERROR << "Unknown variable type for vector: " << rVariableName << std::endl;
                    }
                } else {
                    if (variable_type == VariableType::VECTOR) {
                        for (const auto* p_geometrical_object : r_geometrical_objects) {
                            counter += 3;
                            new_line = check_counter(counter);
                            WriteVectorData(var_file, p_geometrical_object->GetValue(r_variable), new_line, counter == 1, !new_line);
                        }
                    } else if (variable_type == VariableType::TENSOR_SYMMETRIC) {
                        for (const auto* p_geometrical_object : r_geometrical_objects) {
                            WriteSymmetricTensorData(var_file, p_geometrical_object->GetValue(r_variable), true, true);
                        }
                    } else {
                        KRATOS_ERROR << "Unknown variable type for vector: " << rVariableName << std::endl;
                    }
                }
            } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<array_1d<double, 6>>>::Get(rVariableName);
                if (is_ensight_gold) {
                    std::array<std::vector<double>, 6> components;
                    const std::size_t number_of_geometries = r_geometrical_objects.size();
                    components[0].reserve(number_of_geometries);
                    components[1].reserve(number_of_geometries);
                    components[2].reserve(number_of_geometries);
                    components[3].reserve(number_of_geometries);
                    components[4].reserve(number_of_geometries);
                    components[5].reserve(number_of_geometries);
                    for (const auto* p_geometrical_object : r_geometrical_objects) {
                        const auto& r_value = p_geometrical_object->GetValue(r_variable);
                        components[0].push_back(r_value[0]);
                        components[1].push_back(r_value[1]);
                        components[2].push_back(r_value[2]);
                        components[3].push_back(r_value[3]);
                        components[4].push_back(r_value[5]);
                        components[5].push_back(r_value[4]);
                    }
                    for (unsigned int i = 0; i < 6; ++i) {
                        for (std::size_t j = 0; j < number_of_geometries; ++j) {
                            WriteScalarData(var_file, components[i][j]);
                        }
                    }
                } else {
                    for (const auto* p_geometrical_object : r_geometrical_objects) {
                        WriteSymmetricTensorData(var_file, p_geometrical_object->GetValue(r_variable), true, true);
                    }
                }
            } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(rVariableName)) {
                // TODO: Implement asymmetric tensor.
                // const auto& r_variable = KratosComponents<Variable<array_1d<double, 9>>>::Get(rVariableName);
                KRATOS_ERROR << "Asymmetric tensor output is not implemented yet." << std::endl;
            } else if (KratosComponents<Variable<Matrix>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<Matrix>>::Get(rVariableName);
                if (is_ensight_gold) {
                    if (variable_type == VariableType::TENSOR_SYMMETRIC) {
                        std::array<std::vector<double>, 6> components;
                        const std::size_t number_of_geometries = r_geometrical_objects.size();
                        components[0].reserve(number_of_geometries);
                        components[1].reserve(number_of_geometries);
                        components[2].reserve(number_of_geometries);
                        components[3].reserve(number_of_geometries);
                        components[4].reserve(number_of_geometries);
                        components[5].reserve(number_of_geometries);
                        for (const auto* p_geometrical_object : r_geometrical_objects) {
                            const auto& r_value = p_geometrical_object->GetValue(r_variable);
                            components[0].push_back(r_value(0, 0));
                            components[1].push_back(r_value(1, 1));
                            components[2].push_back(r_value(2, 2));
                            components[3].push_back(r_value(0, 1));
                            components[4].push_back(r_value(0, 2));
                            components[5].push_back(r_value(1, 2));
                        }
                    for (unsigned int i = 0; i < 6; ++i) {
                        for (std::size_t j = 0; j < number_of_geometries; ++j) {
                            WriteScalarData(var_file, components[i][j]);
                        }
                    }
                    } else if (variable_type == VariableType::TENSOR_ASYMMETRIC) {
                        // TODO: Implement asymmetric tensor
                        KRATOS_ERROR << "Asymmetric tensor output is not implemented yet." << std::endl;
                    } else {
                        KRATOS_ERROR << "Unknown variable type for matrix: " << rVariableName << std::endl;
                    }
                } else {
                    if (variable_type == VariableType::TENSOR_SYMMETRIC) {
                        for (const auto* p_geometrical_object : r_geometrical_objects) {
                            WriteSymmetricTensorData(var_file, p_geometrical_object->GetValue(r_variable), true, true);
                        }
                    } else if (variable_type == VariableType::TENSOR_ASYMMETRIC) {
                        // TODO: Implement asymmetric tensor
                        KRATOS_ERROR << "Asymmetric tensor output is not implemented yet." << std::endl;
                    } else {
                        KRATOS_ERROR << "Unknown variable type for matrix: " << rVariableName << std::endl;
                    }
                }
            }

            // Adding new line if needed
            if (!new_line) {
                WriteString(var_file, "");
            }
        }
    }

    // Close the variable file
    var_file.close();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void EnSightOutput::WriteGeometricalFlagToFile(
    const std::string& rFlagName,
    const std::string& rFileName,
    const bool IsElement
    )
{
    KRATOS_TRY

    // Get the variable type and its extension
    const auto variable_type = VariableType::SCALAR;
    const std::string file_extension = GetExtensionFile(variable_type);

    // Open the file for writing
    const std::string full_path = GetOutputFileName(rFileName, file_extension);
    std::ofstream var_file(full_path, mFileFormat == FileFormat::BINARY ? std::ios::binary : std::ios::out);

    KRATOS_ERROR_IF_NOT(var_file.is_open()) << "File \"" << full_path << "\" could not be opened!" << std::endl;
    if (mFileFormat == FileFormat::ASCII) {
        var_file << std::scientific << std::setprecision(mDefaultPrecision); // Set precision for ASCII output
    } else {
        var_file << std::ios::binary;
    }

    // Write variable information
    const bool is_ensight_gold = mEnSightFileFormat == EnSightFileFormat::EnSightGold;
    const bool is_ensight_6 = mEnSightFileFormat == EnSightFileFormat::EnSight6;
    const std::string label = is_ensight_gold ? "EnSightGold" : "EnSight6";
    const std::string type_label = GetTypeLabel(variable_type);
    const std::string entity_label = IsElement ? "Elemental" : "Conditional";
    WriteString(var_file, "Per_elem " + type_label + " values for the " + label + " geometry. " + entity_label + " variable: " + rFlagName);

    // Write lambda that check that the counter is 6 for ensight 6 and gold
    auto check_counter_ensight_6 = [](unsigned int& rCounter) -> bool {
        const bool check = rCounter == 6;
        if (check) {
            rCounter = 0;
        }
        return check;
    };
    auto check_counter_ensight_gold = [](unsigned int& rCounter) -> bool {
        (void)rCounter; // suppress unused warning
        return true;
    };
    std::function<bool(unsigned int&)> check_counter = is_ensight_gold
        ? std::function<bool(unsigned int&)>(check_counter_ensight_gold)
        : std::function<bool(unsigned int&)>(check_counter_ensight_6);

    // Write the geometrical data for the parts
    for (const auto& r_part_data : mPartDatas) {
        // Check if the part is of the correct type (element or condition)
        if (r_part_data.PartElements != IsElement) continue;

        // Write part data
        bool new_line = true;
        for (const auto& it : r_part_data.PartGeometricalObjects) {
            unsigned int counter = 0;
            const std::string& r_ensight_element_type = it.first;
            const auto& r_geometrical_objects = it.second;

            // Write part header
            if (is_ensight_6) {
                WriteString(var_file, "part\t" + std::to_string(r_part_data.PartId));
            } else { // Ensight Gold
                WriteString(var_file, "part");
                WriteScalarData(var_file, r_part_data.PartId);
            }
            WriteString(var_file, r_ensight_element_type);

            const auto& r_flag = KratosComponents<Flags>::Get(rFlagName);
            for (const auto* p_geometrical_object : r_geometrical_objects) {
                ++counter;
                new_line = check_counter(counter);
                const int flag = p_geometrical_object->IsDefined(r_flag) ? static_cast<int>(p_geometrical_object->Is(r_flag)) : -1; // Default to -1 if not defined
                WriteScalarData(var_file, flag, new_line, is_ensight_6 && counter == 1, !new_line);
            }

            // Adding new line if needed
            if (!new_line) {
                WriteString(var_file, "");
            }
        }
    }

    // Close the variable file
    var_file.close();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void EnSightOutput::WriteGeometricalGaussVariableToFile(
    const std::string& rVariableName,
    const std::string& rFileName,
    const bool IsElement
    )
{
    KRATOS_TRY

    // Check if the variable is defined on the element
    if (!IsElement) return;

    // Get the variable type and its extension
    const auto variable_type = mVariableTypeMap[rVariableName];
    const std::string file_extension = GetExtensionFile(variable_type);

    // Open the file for writing
    const std::string full_path = GetOutputFileName(rFileName, file_extension);
    std::ofstream var_file(full_path, mFileFormat == FileFormat::BINARY ? std::ios::binary : std::ios::out);

    KRATOS_ERROR_IF_NOT(var_file.is_open()) << "File \"" << full_path << "\" could not be opened!" << std::endl;
    if (mFileFormat == FileFormat::ASCII) {
        var_file << std::scientific << std::setprecision(mDefaultPrecision); // Set precision for ASCII output
    } else {
        var_file << std::ios::binary;
    }

    // Write variable information
    const bool is_ensight_gold = mEnSightFileFormat == EnSightFileFormat::EnSightGold;
    const bool is_ensight_6 = mEnSightFileFormat == EnSightFileFormat::EnSight6;
    const std::string label = is_ensight_gold ? "EnSightGold" : "EnSight6";
    const std::string type_label = GetTypeLabel(variable_type);
    const std::string entity_label = IsElement ? "Elemental" : "Conditional";
    WriteString(var_file, "Per_elem " + type_label + " values for the " + label + " geometry. " + entity_label + " variable: " + rVariableName);

    // Write lambda that check that the counter is 6 for ensight 6 and gold
    auto check_counter_ensight_6 = [](unsigned int& rCounter) -> bool {
        const bool check = rCounter == 6;
        if (check) {
            rCounter = 0;
        }
        return check;
    };
    auto check_counter_ensight_gold = [](unsigned int& rCounter) -> bool {
        (void)rCounter; // suppress unused warning
        return true;
    };
    std::function<bool(unsigned int&)> check_counter = is_ensight_gold
        ? std::function<bool(unsigned int&)>(check_counter_ensight_gold)
        : std::function<bool(unsigned int&)>(check_counter_ensight_6);

    // Write the geometrical data for the parts
    for (const auto& r_part_data : mPartDatas) {
        // Check if the part is of the correct type (element or condition)
        if (r_part_data.PartElements != IsElement) continue;

        // Write part data
        bool new_line = true;
        for (const auto& it : r_part_data.PartGeometricalObjects) {
            unsigned int counter = 0;
            const std::string& r_ensight_element_type = it.first;
            const auto& r_geometrical_objects = it.second;

            // Write part header
            if (is_ensight_6) {
                WriteString(var_file, "part\t" + std::to_string(r_part_data.PartId));
            } else { // Ensight Gold
                WriteString(var_file, "part");
                WriteScalarData(var_file, r_part_data.PartId);
            }
            WriteString(var_file, r_ensight_element_type);

            if (KratosComponents<Variable<bool>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<bool>>::Get(rVariableName);
                for (const auto* p_geometrical_object : r_geometrical_objects) {
                    ++counter;
                    new_line = check_counter(counter);
                    const Element& r_element = dynamic_cast<const Element&>(*p_geometrical_object);

                    // Non-const reference to the element
                    Element& r_non_const_element = const_cast<Element&>(r_element);

                    // Auxiliary values
                    auto& r_this_geometry_begin = r_element.GetGeometry();
                    const GeometryData::IntegrationMethod this_integration_method = r_element.GetIntegrationMethod();
                    const auto& r_integration_points = r_this_geometry_begin.IntegrationPoints(this_integration_method);
                    const std::size_t integration_points_number = r_integration_points.size();

                    // Just if number of GP is greater than 0
                    bool value = false;
                    if (integration_points_number > 0) {
                        // Auxiliary values
                        const auto& r_process_info = mrModelPart.GetProcessInfo();

                        std::vector<bool> aux_result(integration_points_number);
                        r_non_const_element.CalculateOnIntegrationPoints(r_variable, aux_result, r_process_info);
                        for (const bool aux_value : aux_result) {
                            if (aux_value) {
                                value = true;
                                break;
                            }
                        }
                    }
                    WriteScalarData(var_file, value, new_line, is_ensight_6 && counter == 1, !new_line);
                }
            } else if (KratosComponents<Variable<int>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<int>>::Get(rVariableName);
                for (const auto* p_geometrical_object : r_geometrical_objects) {
                    ++counter;
                    new_line = check_counter(counter);
                    const Element& r_element = dynamic_cast<const Element&>(*p_geometrical_object);
                    const int value = GetAverageIntegrationValue<int>(r_element, r_variable);
                    WriteScalarData(var_file, value, new_line, is_ensight_6 && counter == 1, !new_line);
                }
            } else if (KratosComponents<Variable<double>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<double>>::Get(rVariableName);
                for (const auto* p_geometrical_object : r_geometrical_objects) {
                    ++counter;
                    new_line = check_counter(counter);
                    const Element& r_element = dynamic_cast<const Element&>(*p_geometrical_object);
                    const double value = GetAverageIntegrationValue<double>(r_element, r_variable);
                    WriteScalarData(var_file, value, new_line, is_ensight_6 && counter == 1, !new_line);
                }
            } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
                if (is_ensight_gold) {
                    std::array<std::vector<double>, 3> components;
                    const std::size_t number_of_geometries = r_geometrical_objects.size();
                    components[0].reserve(number_of_geometries);
                    components[1].reserve(number_of_geometries);
                    components[2].reserve(number_of_geometries);
                    for (const auto* p_geometrical_object : r_geometrical_objects) {
                        const Element& r_element = dynamic_cast<const Element&>(*p_geometrical_object);
                        const array_1d<double, 3> value = GetAverageIntegrationValue<array_1d<double, 3>>(r_element, r_variable);
                        components[0].push_back(value[0]);
                        components[1].push_back(value[1]);
                        components[2].push_back(value[2]);
                    }
                    for (unsigned int i = 0; i < 3; ++i) {
                        for (std::size_t j = 0; j < number_of_geometries; ++j) {
                            WriteScalarData(var_file, components[i][j]);
                        }
                    }
                } else {
                    for (const auto* p_geometrical_object : r_geometrical_objects) {
                        counter += 3;
                        new_line = check_counter(counter);
                        const Element& r_element = dynamic_cast<const Element&>(*p_geometrical_object);
                        const array_1d<double, 3> value = GetAverageIntegrationValue<array_1d<double, 3>>(r_element, r_variable);
                        WriteVectorData(var_file, value, new_line, counter == 1, !new_line);
                    }
                }
            } else if (KratosComponents<Variable<Vector>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<Vector>>::Get(rVariableName);
                if (is_ensight_gold) {
                    if (variable_type == VariableType::VECTOR) {
                        std::array<std::vector<double>, 3> components;
                        const std::size_t number_of_geometries = r_geometrical_objects.size();
                        components[0].reserve(number_of_geometries);
                        components[1].reserve(number_of_geometries);
                        components[2].reserve(number_of_geometries);
                        for (const auto* p_geometrical_object : r_geometrical_objects) {
                            const Element& r_element = dynamic_cast<const Element&>(*p_geometrical_object);
                            const Vector value = GetAverageIntegrationValue<Vector>(r_element, r_variable);
                            components[0].push_back(value[0]);
                            components[1].push_back(value[1]);
                            components[2].push_back(value[2]);
                        }
                        for (unsigned int i = 0; i < 3; ++i) {
                            for (std::size_t j = 0; j < number_of_geometries; ++j) {
                                WriteScalarData(var_file, components[i][j]);
                            }
                        }
                    } else if (variable_type == VariableType::TENSOR_SYMMETRIC) {
                        std::array<std::vector<double>, 6> components;
                        const std::size_t number_of_geometries = r_geometrical_objects.size();
                        components[0].reserve(number_of_geometries);
                        components[1].reserve(number_of_geometries);
                        components[2].reserve(number_of_geometries);
                        components[3].reserve(number_of_geometries);
                        components[4].reserve(number_of_geometries);
                        components[5].reserve(number_of_geometries);
                        for (const auto* p_geometrical_object : r_geometrical_objects) {
                            const Element& r_element = dynamic_cast<const Element&>(*p_geometrical_object);
                            const Vector value = GetAverageIntegrationValue<Vector>(r_element, r_variable);
                            components[0].push_back(value[0]);
                            components[1].push_back(value[1]);
                            components[2].push_back(value[2]);
                            components[3].push_back(value[3]);
                            components[4].push_back(value[5]);
                            components[5].push_back(value[4]);
                        }
                        for (unsigned int i = 0; i < 6; ++i) {
                            for (std::size_t j = 0; j < number_of_geometries; ++j) {
                                WriteScalarData(var_file, components[i][j]);
                            }
                        }
                    } else {
                        KRATOS_ERROR << "Unknown variable type for vector: " << rVariableName << std::endl;
                    }
                } else {
                    if (variable_type == VariableType::VECTOR) {
                        for (const auto* p_geometrical_object : r_geometrical_objects) {
                            counter += 3;
                            new_line = check_counter(counter);
                            const Element& r_element = dynamic_cast<const Element&>(*p_geometrical_object);
                            const Vector value = GetAverageIntegrationValue<Vector>(r_element, r_variable);
                            WriteVectorData(var_file, value, new_line, counter == 1, !new_line);
                        }
                    } else if (variable_type == VariableType::TENSOR_SYMMETRIC) {
                        for (const auto* p_geometrical_object : r_geometrical_objects) {
                            const Element& r_element = dynamic_cast<const Element&>(*p_geometrical_object);
                            const Vector value = GetAverageIntegrationValue<Vector>(r_element, r_variable);
                            WriteSymmetricTensorData(var_file, value, true, true);
                        }
                    } else {
                        KRATOS_ERROR << "Unknown variable type for vector: " << rVariableName << std::endl;
                    }
                }
            } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<array_1d<double, 6>>>::Get(rVariableName);
                if (is_ensight_gold) {
                    std::array<std::vector<double>, 6> components;
                    const std::size_t number_of_geometries = r_geometrical_objects.size();
                    components[0].reserve(number_of_geometries);
                    components[1].reserve(number_of_geometries);
                    components[2].reserve(number_of_geometries);
                    components[3].reserve(number_of_geometries);
                    components[4].reserve(number_of_geometries);
                    components[5].reserve(number_of_geometries);
                    for (const auto* p_geometrical_object : r_geometrical_objects) {
                        const Element& r_element = dynamic_cast<const Element&>(*p_geometrical_object);
                        const array_1d<double, 6> value = GetAverageIntegrationValue<array_1d<double, 6>>(r_element, r_variable);
                        components[0].push_back(value[0]);
                        components[1].push_back(value[1]);
                        components[2].push_back(value[2]);
                        components[3].push_back(value[3]);
                        components[4].push_back(value[5]);
                        components[5].push_back(value[4]);
                    }
                    for (unsigned int i = 0; i < 6; ++i) {
                        for (std::size_t j = 0; j < number_of_geometries; ++j) {
                            WriteScalarData(var_file, components[i][j]);
                        }
                    }
                } else {
                    for (const auto* p_geometrical_object : r_geometrical_objects) {
                        const Element& r_element = dynamic_cast<const Element&>(*p_geometrical_object);
                        const array_1d<double, 6> value = GetAverageIntegrationValue<array_1d<double, 6>>(r_element, r_variable);
                        WriteSymmetricTensorData(var_file, value, true, true);
                    }
                }
            } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(rVariableName)) {
                // TODO: Implement asymmetric tensor.
                // const auto& r_variable = KratosComponents<Variable<array_1d<double, 9>>>::Get(rVariableName);
                KRATOS_ERROR << "Asymmetric tensor output is not implemented yet." << std::endl;
            } else if (KratosComponents<Variable<Matrix>>::Has(rVariableName)) {
                const auto& r_variable = KratosComponents<Variable<Matrix>>::Get(rVariableName);
                if (is_ensight_gold) {
                    if (variable_type == VariableType::TENSOR_SYMMETRIC) {
                        std::array<std::vector<double>, 6> components;
                        const std::size_t number_of_geometries = r_geometrical_objects.size();
                        components[0].reserve(number_of_geometries);
                        components[1].reserve(number_of_geometries);
                        components[2].reserve(number_of_geometries);
                        components[3].reserve(number_of_geometries);
                        components[4].reserve(number_of_geometries);
                        components[5].reserve(number_of_geometries);
                        for (const auto* p_geometrical_object : r_geometrical_objects) {
                            const Element& r_element = dynamic_cast<const Element&>(*p_geometrical_object);
                            const Matrix value = GetAverageIntegrationValue<Matrix>(r_element, r_variable);
                            components[0].push_back(value(0, 0));
                            components[1].push_back(value(1, 1));
                            components[2].push_back(value(2, 2));
                            components[3].push_back(value(0, 1));
                            components[4].push_back(value(0, 2));
                            components[5].push_back(value(1, 2));
                        }
                    for (unsigned int i = 0; i < 6; ++i) {
                        for (std::size_t j = 0; j < number_of_geometries; ++j) {
                            WriteScalarData(var_file, components[i][j]);
                        }
                    }
                    } else if (variable_type == VariableType::TENSOR_ASYMMETRIC) {
                        // TODO: Implement asymmetric tensor
                        KRATOS_ERROR << "Asymmetric tensor output is not implemented yet." << std::endl;
                    } else {
                        KRATOS_ERROR << "Unknown variable type for matrix: " << rVariableName << std::endl;
                    }
                } else {
                    if (variable_type == VariableType::TENSOR_SYMMETRIC) {
                        for (const auto* p_geometrical_object : r_geometrical_objects) {
                            const Element& r_element = dynamic_cast<const Element&>(*p_geometrical_object);
                            const Matrix value = GetAverageIntegrationValue<Matrix>(r_element, r_variable);
                            WriteSymmetricTensorData(var_file, value, true, true);
                        }
                    } else if (variable_type == VariableType::TENSOR_ASYMMETRIC) {
                        // TODO: Implement asymmetric tensor
                        KRATOS_ERROR << "Asymmetric tensor output is not implemented yet." << std::endl;
                    } else {
                        KRATOS_ERROR << "Unknown variable type for matrix: " << rVariableName << std::endl;
                    }
                }
            }

            // Adding new line if needed
            if (!new_line) {
                WriteString(var_file, "");
            }
        }
    }

    // Close the variable file
    var_file.close();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::string EnSightOutput::GetOutputFileName(
    const std::string& rFileLabel,
    const std::string& rFileExtension
    ) const
{
    KRATOS_TRY

    std::string output_file_name = "";

    if (rFileLabel != "") { // user specified file name externally
        output_file_name = rFileLabel;
    } else {
        const int rank = mrModelPart.GetCommunicator().MyPID();
        const std::string model_part_name = mrModelPart.Name();

        std::string label;
        std::stringstream ss;
        const std::string output_control = mOutputSettings["output_control_type"].GetString();
        if (output_control == "step") {
            ss << std::fixed << std::setprecision(mDefaultPrecision) << std::setfill('0')
            << mrModelPart.GetProcessInfo()[STEP];
            label = ss.str();
        } else if(output_control == "time") {
            ss << std::fixed << std::setprecision(mDefaultPrecision) << std::setfill('0')
            << mrModelPart.GetProcessInfo()[TIME];
            label = ss.str();
        } else {
            KRATOS_ERROR << "Option for \"output_control_type\": " << output_control
                <<" not recognised!\nPossible output_control_type options "
                << "are: \"step\", \"time\"" << std::endl;
        }

        const std::string& r_custom_name_prefix = mOutputSettings["custom_name_prefix"].GetString();
        const std::string& r_custom_name_postfix = mOutputSettings["custom_name_postfix"].GetString();
        output_file_name += r_custom_name_prefix + model_part_name + r_custom_name_postfix + "_" + std::to_string(rank) + "_" + label;
    }

    output_file_name += rFileExtension;

    if (mOutputSettings["save_output_files_in_folder"].GetBool()) {
        const std::filesystem::path output_path = mOutputSettings["output_path"].GetString();

        // Create folder if it doesn't exist before
        FilesystemExtensions::MPISafeCreateDirectories(output_path);

        output_file_name = (output_path / output_file_name).string();
    }

    return output_file_name;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void EnSightOutput::CollectPartData(
    const ModelPart& rSubModelPart,
    PartData& rPartData
    ) const
{
    KRATOS_TRY

    // Collect unique node IDs from the sub-model part
    std::unordered_set<std::size_t> unique_node_ids;

    // Get the entity type for the current model part
    const auto entity_type = GetEntityType(rSubModelPart);

    // Collect unique node IDs based on the entity type
    if (entity_type == EntityType::ELEMENT) {
        for (const auto& r_element : rSubModelPart.Elements()) {
            if (InputOutputUtilities::SkippableEntity(r_element, "EnSightOutput")) continue;
            for (const auto& r_node : r_element.GetGeometry()) {
                unique_node_ids.insert(r_node.Id());
            }
        }
    } else if (entity_type == EntityType::CONDITION) {
        for (const auto& r_condition : rSubModelPart.Conditions()) {
            if (InputOutputUtilities::SkippableEntity(r_condition, "EnSightOutput")) continue;
            for (const auto& r_node : r_condition.GetGeometry()) {
                unique_node_ids.insert(r_node.Id());
            }
        }
    } else {
        KRATOS_ERROR << "Entity type not supported for EnSight Gold output: " << static_cast<int>(entity_type) << std::endl;
    }

    // Store the nodes in the part data structure
    rPartData.PartNodes.reserve(unique_node_ids.size());
    for (const std::size_t id : unique_node_ids) {
        rPartData.PartNodes.push_back(&mrModelPart.GetNode(id));
    }

    // Map Kratos IDs to local IDs
    for (size_t i = 0; i < rPartData.PartNodes.size(); ++i) {
        rPartData.KratosIdToLocalId[rPartData.PartNodes[i]->Id()] = i + 1; // 1-based index
    }

    // Fill PartGeometricalObjects
    if (entity_type == EntityType::ELEMENT) {
        for (const auto& r_element : rSubModelPart.Elements()) {
            if (InputOutputUtilities::SkippableEntity(r_element, "EnSightOutput")) continue;
            const std::string name = GetEnSightName(r_element.GetGeometry().GetGeometryType());
            if (name != "unsupported_type") {
                rPartData.PartGeometricalObjects[name].push_back(&r_element);
            }
        }
    } else if (entity_type == EntityType::CONDITION) {
        for (const auto& r_condition : rSubModelPart.Conditions()) {
            if (InputOutputUtilities::SkippableEntity(r_condition, "EnSightOutput")) continue;
            const std::string name = GetEnSightName(r_condition.GetGeometry().GetGeometryType());
            if (name != "unsupported_type") {
                rPartData.PartGeometricalObjects[name].push_back(&r_condition);
            }
        }
    } else {
        KRATOS_ERROR << "Entity type not supported for EnSight Gold output: " << static_cast<int>(entity_type) << std::endl;
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template <typename TData>
void EnSightOutput::WriteScalarData(
    std::ofstream& rFileStream,
    const TData& rData,
    const bool EndOfLine,
    const bool AddInitialTabulation,
    const bool AddEndTabulation
    ) const
{
    KRATOS_TRY

    if (mFileFormat == FileFormat::ASCII) {
        if (AddInitialTabulation) rFileStream << "\t";
        if constexpr (std::is_integral_v<TData>) {
            rFileStream << std::setw(10) << rData; // Integer types are formatted with fixed width
        } else {
            rFileStream << std::scientific << std::setprecision(mDefaultPrecision) << rData << std::fixed; // Float/Double handled by stream settings
        }
        if (EndOfLine) rFileStream << "\n";
        if (AddEndTabulation) rFileStream << "\t";
    } else { // Binary
        rFileStream.write(reinterpret_cast<const char*>(&rData), sizeof(TData));
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template <typename TData>
void EnSightOutput::WriteVectorData(
    std::ofstream& rFileStream,
    const TData& rData,
    const bool EndOfLine,
    const bool AddInitialTabulation,
    const bool AddEndTabulation
    ) const
{
    KRATOS_TRY

    if (mFileFormat == FileFormat::ASCII) {
        if (AddInitialTabulation) rFileStream << "\t";
        rFileStream << std::scientific << std::setprecision(mDefaultPrecision) << rData[0] << "\t" << rData[1] << "\t" << rData[2] << std::fixed;
        if (EndOfLine) rFileStream << "\n";
        if (AddEndTabulation) rFileStream << "\t";
    } else { // Binary
        double buffer[3] = {static_cast<double>(rData[0]), static_cast<double>(rData[1]), static_cast<double>(rData[2])};
        rFileStream.write(reinterpret_cast<const char*>(buffer), sizeof(buffer));
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template <typename TData>
void EnSightOutput::WriteSymmetricTensorData(
    std::ofstream& rFileStream,
    const TData& rData,
    const bool EndOfLine,
    const bool AddInitialTabulation,
    const bool AddEndTabulation
    ) const
{
    KRATOS_TRY

    // Ordering of Voigt notation 11 22 33 12 23 13 (Kratos) -> 11 22 33 12 13 23 (EnSight)
    if (mFileFormat == FileFormat::ASCII) {
        if (AddInitialTabulation) rFileStream << "\t";
        rFileStream << std::scientific << std::setprecision(mDefaultPrecision);
        if constexpr (std::is_same_v<TData, Matrix>) { // Matrix type
            rFileStream << rData(0, 0) << "\t" << rData(1, 1) << "\t" << rData(2, 2) << "\t" << rData(0, 1) << "\t" << rData(0, 2) << "\t" << rData(1, 2) << std::fixed;
        } else { // Vector or array, last 2 indexes swapped
            rFileStream << rData[0] << "\t" << rData[1] << "\t" << rData[2] << "\t" << rData[3] << "\t" << rData[5] << "\t" << rData[4] << std::fixed;
        }
        if (EndOfLine) rFileStream << "\n";
        if (AddEndTabulation) rFileStream << "\t";
    } else { // Binary
        double buffer[6];
        if constexpr (std::is_same_v<TData, Matrix>) { // Matrix type
            buffer[0] = static_cast<double>(rData(0, 0));
            buffer[1] = static_cast<double>(rData(1, 1));
            buffer[2] = static_cast<double>(rData(2, 2));
            buffer[3] = static_cast<double>(rData(0, 1));
            buffer[4] = static_cast<double>(rData(0, 2));
            buffer[5] = static_cast<double>(rData(1, 2));
        } else { // Vector or array, last 2 indexes swapped
            buffer[0] = static_cast<double>(rData[0]);
            buffer[1] = static_cast<double>(rData[1]);
            buffer[2] = static_cast<double>(rData[2]);
            buffer[3] = static_cast<double>(rData[3]);
            buffer[4] = static_cast<double>(rData[5]);
            buffer[5] = static_cast<double>(rData[4]);
        }
        rFileStream.write(reinterpret_cast<const char*>(buffer), sizeof(buffer));
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void EnSightOutput::WriteString(
    std::ofstream& rFileStream,
    const std::string& rString,
    const bool EndOfLine,
    const bool AddInitialTabulation,
    const bool AddEndTabulation
    ) const
{
    KRATOS_TRY

    if (mFileFormat == FileFormat::ASCII) {
        if (AddInitialTabulation) rFileStream << "\t";
        rFileStream << rString;
        if (EndOfLine) rFileStream << "\n";
        if (AddEndTabulation) rFileStream << "\t";
    } else { // Binary
        char buffer[80] = {}; // Zero-initialize
        rString.copy(buffer, std::min(rString.length(), static_cast<size_t>(79)));
        rFileStream.write(buffer, 80);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::string EnSightOutput::GetEnSightName(const GeometryData::KratosGeometryType& rGeometryType) const
{
    KRATOS_TRY

    switch (rGeometryType) {
        // 0D Geometries
        case GeometryData::KratosGeometryType::Kratos_Point2D:               return "point";
        case GeometryData::KratosGeometryType::Kratos_Point3D:               return "point";

        // 1D Geometries
        case GeometryData::KratosGeometryType::Kratos_Line2D2:               return "bar2";
        case GeometryData::KratosGeometryType::Kratos_Line3D2:               return "bar2";
        case GeometryData::KratosGeometryType::Kratos_Line2D3:               return "bar3";
        case GeometryData::KratosGeometryType::Kratos_Line3D3:               return "bar3";

        // 2D Geometries
        case GeometryData::KratosGeometryType::Kratos_Triangle2D3:           return "tria3";
        case GeometryData::KratosGeometryType::Kratos_Triangle3D3:           return "tria3";
        case GeometryData::KratosGeometryType::Kratos_Triangle2D6:           return "tria6";
        case GeometryData::KratosGeometryType::Kratos_Triangle3D6:           return "tria6";
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4:      return "quad4";
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4:      return "quad4";
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8:      return "quad8";
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8:      return "quad8";

        // 3D Geometries
        case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:         return "tetra4";
        case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10:        return "tetra10";
        case GeometryData::KratosGeometryType::Kratos_Pyramid3D5:            return "pyramid5";
        case GeometryData::KratosGeometryType::Kratos_Pyramid3D13:           return "pyramid13";
        case GeometryData::KratosGeometryType::Kratos_Prism3D6:              return "penta6";
        case GeometryData::KratosGeometryType::Kratos_Prism3D15:             return "penta15";
        case GeometryData::KratosGeometryType::Kratos_Hexahedra3D8:          return "hexa8";
        case GeometryData::KratosGeometryType::Kratos_Hexahedra3D20:         return "hexa20";

        default:
            KRATOS_ERROR << "Geometry type not supported for EnSight Gold output: " << static_cast<int>(rGeometryType) << std::endl;
            return "unsupported_type"; // Fallback for unsupported types
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void EnSightOutput::UpdatePartData()
{
    KRATOS_TRY

    // Define model part pointer
    std::vector<const ModelPart*> p_model_parts;

    // Check if we need to output sub-model parts
    const bool output_sub_model_parts = mOutputSettings["output_sub_model_parts"].GetBool();

    // Check if the model part has sub-model parts
    const auto& r_sub_model_parts = mrModelPart.SubModelParts();
    if (!r_sub_model_parts.empty() && output_sub_model_parts) {
        // First we push back the main model part
        p_model_parts.reserve(r_sub_model_parts.size() + 1);
        p_model_parts.push_back(&mrModelPart);

        // Now add all submodel parts
        for (const auto& r_sub_model_part : r_sub_model_parts) {
            p_model_parts.push_back(&r_sub_model_part);
        }

        // Sort all entries except the first one (main model part)
        std::sort(p_model_parts.begin() + 1, p_model_parts.end(), 
            [](const ModelPart* a, const ModelPart* b) {
                return a->Name() < b->Name();
            }
        );
    } else {
        // If no sub-model parts, use the main model part
        p_model_parts.push_back(&mrModelPart);
    }

    // Prepare parts data
    mPartDatas.clear();
    mPartDatas.resize(p_model_parts.size());
    std::size_t part_index = 1; // EnSight part index starts at 1 (at least Hypervgiew does not like 0, Paraview is OK with 0)
    for (const auto& p_model_part : p_model_parts) {
        const ModelPart& r_sub_model_part = *p_model_part;
        // Skip empty parts
        const std::size_t number_of_elements = r_sub_model_part.NumberOfElements();
        const std::size_t number_of_conditions = r_sub_model_part.NumberOfConditions();
        if (number_of_elements == 0 && number_of_conditions == 0) continue;

        // Collect data for the part
        PartData& r_part_data = mPartDatas[part_index - 1];

        // If elements are zero, we set PartElements to false to know we have only conditions
        if (number_of_elements == 0) {
            r_part_data.PartElements = false;
        }

        // Collect data for the part
        CollectPartData(r_sub_model_part, r_part_data);

        // Set part name and ID
        r_part_data.PartName = r_sub_model_part.Name();
        r_part_data.PartId = part_index++;
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

EnSightOutput::VariableType EnSightOutput::GetVariableType(
    const std::string& rVariableName,
    ModelPart* pModelPart,
    const EntityType Entity,
    const bool IsHistorical
    ) const
{
    KRATOS_TRY

    if (KratosComponents<Variable<bool>>::Has(rVariableName)) {
        return VariableType::SCALAR;
    } else if (KratosComponents<Variable<int>>::Has(rVariableName)) {
        return VariableType::SCALAR;
    }  else if (KratosComponents<Variable<double>>::Has(rVariableName)) {
        return VariableType::SCALAR;
    } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)) {
        return VariableType::VECTOR;
    } else if (KratosComponents<Variable<Vector>>::Has(rVariableName)) {
        // return VariableType::VECTOR_UNDEFINED;
        // We need to check the variable size manually
        KRATOS_ERROR_IF_NOT(pModelPart) << "Model part pointer is null for variable: " << rVariableName << std::endl;
        KRATOS_ERROR_IF(Entity == EntityType::UNDEFINED) << "Undefined entity type for variable: " << rVariableName << std::endl;
        // Check if the variable is a vector with size 3 or 6
        const auto& r_variable = KratosComponents<Variable<Vector>>::Get(rVariableName);
        Vector value;
        // Get the vector value
        if (Entity == EntityType::NODE) {
            auto it_node_begin = pModelPart->Nodes().begin();
            if (IsHistorical) {
                noalias(value) = it_node_begin->FastGetSolutionStepValue(r_variable);
            } else {
                noalias(value) = it_node_begin->GetValue(r_variable);
            }
        } else if (Entity == EntityType::ELEMENT) {
            auto it_elem_begin = pModelPart->Elements().begin();
            if (IsHistorical) { // Gauss point
                // Auxiliary values
                const auto& r_process_info = mrModelPart.GetProcessInfo();
                auto& r_this_geometry_begin = it_elem_begin->GetGeometry();
                const GeometryData::IntegrationMethod this_integration_method = it_elem_begin->GetIntegrationMethod();
                const auto& r_integration_points = r_this_geometry_begin.IntegrationPoints(this_integration_method);
                const SizeType integration_points_number = r_integration_points.size();
                std::vector<Vector> aux_result(integration_points_number);
                it_elem_begin->CalculateOnIntegrationPoints(r_variable, aux_result, r_process_info);
                noalias(value) = aux_result[0];
            } else {
                noalias(value) = it_elem_begin->GetValue(r_variable);
            }
        } else if (Entity == EntityType::CONDITION) {
            auto it_cond_begin = pModelPart->Conditions().begin();
            noalias(value) = it_cond_begin->GetValue(r_variable);
        } else {
            KRATOS_ERROR << "Unknown entity type for variable: " << rVariableName << std::endl;
        }
        // Check size
        if (value.size() == 3) {
            return VariableType::VECTOR;
        } else if (value.size() == 6) {
            return VariableType::TENSOR_SYMMETRIC;
        } else if (value.size() == 9) {
            return VariableType::TENSOR_ASYMMETRIC;
        } else {
            KRATOS_ERROR << "Unknown vector size for variable: " << rVariableName << std::endl;
        }
    } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)) {
        return VariableType::TENSOR_SYMMETRIC;
    } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(rVariableName)) {
        return VariableType::TENSOR_ASYMMETRIC;
    } else if (KratosComponents<Variable<Matrix>>::Has(rVariableName)) {
        // return VariableType::TENSOR_UNDEFINED;
        // We need to check the variable size manually
        KRATOS_ERROR_IF_NOT(pModelPart) << "Model part pointer is null for variable: " << rVariableName << std::endl;
        KRATOS_ERROR_IF(Entity == EntityType::UNDEFINED) << "Undefined entity type for variable: " << rVariableName << std::endl;
        // Check if the variable is a vector with size 6 or 9
        const auto& r_variable = KratosComponents<Variable<Matrix>>::Get(rVariableName);
        Matrix value;
        // Get the matrix value
        if (Entity == EntityType::NODE) {
            auto it_node_begin = pModelPart->Nodes().begin();
            if (IsHistorical) {
                noalias(value) = it_node_begin->FastGetSolutionStepValue(r_variable);
            } else {
                noalias(value) = it_node_begin->GetValue(r_variable);
            }
        } else if (Entity == EntityType::ELEMENT) {
            auto it_elem_begin = pModelPart->Elements().begin();
            if (IsHistorical) { // Gauss point
                // Auxiliary values
                const auto& r_process_info = mrModelPart.GetProcessInfo();
                auto& r_this_geometry_begin = it_elem_begin->GetGeometry();
                const GeometryData::IntegrationMethod this_integration_method = it_elem_begin->GetIntegrationMethod();
                const auto& r_integration_points = r_this_geometry_begin.IntegrationPoints(this_integration_method);
                const SizeType integration_points_number = r_integration_points.size();
                std::vector<Matrix> aux_result(integration_points_number);
                it_elem_begin->CalculateOnIntegrationPoints(r_variable, aux_result, r_process_info);
                noalias(value) = aux_result[0];
            } else {
                noalias(value) = it_elem_begin->GetValue(r_variable);
            }
        } else if (Entity == EntityType::CONDITION) {
            auto it_cond_begin = pModelPart->Conditions().begin();
            noalias(value) = it_cond_begin->GetValue(r_variable);
        } else {
            KRATOS_ERROR << "Unknown entity type for variable: " << rVariableName << std::endl;
        }
        // Lambda to check if a matrix is symmetric
        auto is_symmetric = [](const Matrix& rMatrix) {
            KRATOS_ERROR_IF(rMatrix.size1() != 3 || rMatrix.size2() != 3) << "Matrix is not 3x3: " << rMatrix.size1() << "x" << rMatrix.size2() << std::endl;
            for (std::size_t i = 0; i < 3; ++i) {
                for (std::size_t j = i + 1; j < 3; ++j) {
                    if (std::abs(rMatrix(i, j) - rMatrix(j, i)) > 1e-12) return false;
                }
            }
            return true;
        };
        // Check the symmetry of the matrix
        if (is_symmetric(value)) {
            return VariableType::TENSOR_SYMMETRIC;
        } else {
            return VariableType::TENSOR_ASYMMETRIC;
        }
    } else {
        KRATOS_ERROR << "Unknown variable type for EnSight output: " << rVariableName << std::endl;
        return VariableType::UNKNOWN;
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::string EnSightOutput::GetExtensionFile(const VariableType Type) const
{
    KRATOS_TRY

    switch (Type) {
        case VariableType::SCALAR: return ".scl";
        case VariableType::VECTOR: return ".vec";
        case VariableType::VECTOR_UNDEFINED: return ".vec";  // NOTE: Should not happen
        case VariableType::TENSOR_SYMMETRIC: return ".ten";
        case VariableType::TENSOR_ASYMMETRIC: return ".ten";
        case VariableType::TENSOR_UNDEFINED: return ".ten";  // NOTE: Should not happen
        case VariableType::COMPLEX_SCALAR: return ".cplx";   // TODO: Currently unsupported
        case VariableType::COMPLEX_VECTOR: return ".cplx";   // TODO: Currently unsupported
        default:
            KRATOS_ERROR << "Unknown variable type for EnSight Gold output: " << static_cast<int>(Type) << std::endl;
    }
    return ".unknown"; // Fallback for unknown types

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::string EnSightOutput::GetTypeLabel(const VariableType Type) const
{
    KRATOS_TRY

    switch (Type) {
        case VariableType::SCALAR: return "scalar";
        case VariableType::VECTOR: return "vector";
        case VariableType::VECTOR_UNDEFINED: return "vector undefined";   // NOTE: Should not happen
        case VariableType::TENSOR_SYMMETRIC: return "tensor symmetric";
        case VariableType::TENSOR_ASYMMETRIC: return "tensor asymmetric"; // TODO: Check if this is correct
        case VariableType::TENSOR_UNDEFINED: return "tensor undefined";   // NOTE: Should not happen
        case VariableType::COMPLEX_SCALAR: return "complex scalar";       // TODO: Currently unsupported
        case VariableType::COMPLEX_VECTOR: return "complex vector";       // TODO: Currently unsupported
        default:
            KRATOS_ERROR << "Unknown variable type for EnSight Gold output: " << static_cast<int>(Type) << std::endl;
    }
    return ".unknown"; // Fallback for unknown types

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

BoundingBox<Point> EnSightOutput::ComputeBoundingBox(
    const ModelPart& rModelPart,
    const double Coefficient
    )
{
    KRATOS_TRY

    // The bounding box is the same as the bins bounding box
    using MultipleReduction = CombinedReduction<MinReduction<double>, MinReduction  <double>, MinReduction<double>, MaxReduction<double>, MaxReduction<double>, MaxReduction<double>>;
    BoundingBox<Point> bb;
    auto& r_min = bb.GetMinPoint();
    auto& r_max = bb.GetMaxPoint();
    r_min[0] = std::numeric_limits<double>::max();
    r_min[1] = std::numeric_limits<double>::max();
    r_min[2] = std::numeric_limits<double>::max();
    r_max[0] = std::numeric_limits<double>::lowest();
    r_max[1] = std::numeric_limits<double>::lowest();
    r_max[2] = std::numeric_limits<double>::lowest();
    std::tie(r_min[0], r_min[1], r_min[2], r_max[0], r_max[1], r_max[2]) = block_for_each<MultipleReduction>(rModelPart.Nodes(), [](Node& rNode) {
        return std::make_tuple(rNode.X(), rNode.Y(), rNode.Z(), rNode.X(), rNode.Y(), rNode.Z());
    });
    r_max[0] = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(r_max[0]);
    r_max[1] = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(r_max[1]);
    r_max[2] = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(r_max[2]);

    r_min[0] = rModelPart.GetCommunicator().GetDataCommunicator().MinAll(r_min[0]);
    r_min[1] = rModelPart.GetCommunicator().GetDataCommunicator().MinAll(r_min[1]);
    r_min[2] = rModelPart.GetCommunicator().GetDataCommunicator().MinAll(r_min[2]);

    // Modify the bounding box size in case of computing the complementary distance
    if (Coefficient != 1.0) {
        // The length of the sides of the bounding box
        const double dx = r_max[0] - r_min[0];
        const double dy = r_max[1] - r_min[1];
        const double dz = r_max[2] - r_min[2];

        // Increase the bounding box
        r_max[0] += (Coefficient - 1.0) * dx;
        r_max[1] += (Coefficient - 1.0) * dy;
        r_max[2] += (Coefficient - 1.0) * dz;
        r_min[0] -= (Coefficient - 1.0) * dx;
        r_min[1] -= (Coefficient - 1.0) * dy;
        r_min[2] -= (Coefficient - 1.0) * dz;
    }

    return bb;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void EnSightOutput::GetGeometryConnectivity(
    const GeometricalObject& rGeometricalObject,
    std::vector<std::size_t>& rConnectivity
    )
{
    KRATOS_TRY

    const auto& r_geometry = rGeometricalObject.GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    if (rConnectivity.size() != number_of_nodes) {
        rConnectivity.resize(number_of_nodes);
    }
    const auto geometry_type = r_geometry.GetGeometryType();

    switch (geometry_type) {
        case GeometryData::KratosGeometryType::Kratos_Prism3D15:
            // Handle quadratic prism geometry
            // First 9 nodes are the same
            for (std::size_t i = 0; i < 9; ++i) {
                rConnectivity[i] = r_geometry[i].Id();
            }
            // Last 6 nodes are swapped
            rConnectivity[9] = r_geometry[12].Id();
            rConnectivity[10] = r_geometry[13].Id();
            rConnectivity[11] = r_geometry[14].Id();
            rConnectivity[12] = r_geometry[9].Id();
            rConnectivity[13] = r_geometry[10].Id();
            rConnectivity[14] = r_geometry[11].Id();
            break;
        case GeometryData::KratosGeometryType::Kratos_Hexahedra3D20:
            // Handle quadratic hexahedra geometry
            // First 12 nodes are the same
            for (std::size_t i = 0; i < 12; ++i) {
                rConnectivity[i] = r_geometry[i].Id();
            }
            // Last 8 nodes are swapped
            rConnectivity[12] = r_geometry[16].Id();
            rConnectivity[13] = r_geometry[17].Id();
            rConnectivity[14] = r_geometry[18].Id();
            rConnectivity[15] = r_geometry[19].Id();
            rConnectivity[16] = r_geometry[12].Id();
            rConnectivity[17] = r_geometry[13].Id();
            rConnectivity[18] = r_geometry[14].Id();
            rConnectivity[19] = r_geometry[15].Id();
            break;
        default:
            for (std::size_t i = 0; i < number_of_nodes; ++i) {
                rConnectivity[i] = r_geometry[i].Id();
            }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TData>
TData EnSightOutput::GetAverageIntegrationValue(
    const Element& rElement,
    const Variable<TData>& rVariable
    )
{
    KRATOS_TRY

    // Average value
    TData value = TData();

    // Non-const reference to the element
    Element& r_non_const_element = const_cast<Element&>(rElement);

    // Auxiliary values
    auto& r_this_geometry_begin = rElement.GetGeometry();
    const GeometryData::IntegrationMethod this_integration_method = rElement.GetIntegrationMethod();
    const auto& r_integration_points = r_this_geometry_begin.IntegrationPoints(this_integration_method);
    const std::size_t integration_points_number = r_integration_points.size();

    // Just if number of GP is greater than 0
    if (integration_points_number > 0) {
        // Auxiliary values
        const auto& r_process_info = mrModelPart.GetProcessInfo();

        std::vector<TData> aux_result(integration_points_number);
        r_non_const_element.CalculateOnIntegrationPoints(rVariable, aux_result, r_process_info);
        value = aux_result[0];
        for (unsigned int i = 1; i < integration_points_number; ++i) {
            value += aux_result[i];
        }
        value /= static_cast<double>(integration_points_number);
    }
    return value;

    KRATOS_CATCH("")
}

} // namespace Kratos