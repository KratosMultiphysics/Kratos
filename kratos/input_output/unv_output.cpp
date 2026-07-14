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
//  Contributors:    Vicente Mataix Ferrandiz
//
//

// System includes
#include <set>
#include <functional>

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
            } else if (r_variable_name == "REACTION") {
                rUnvVariableKeys[variable_key] = 8;
            } else {
                const auto& r_pair = unv_vector3_variables.begin() + array_counter;
                KRATOS_WARNING("UnvOutput") << "Unknown variable: " << r_variable_name << ". Using UNV 3 components vector " << r_pair->second << " variable name id: " << r_pair->first << std::endl;
                rUnvVariableKeys[variable_key] = r_pair->first;
                ++array_counter;
            }
        }
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

    // The file format: UNV is a text (ASCII) format only, so anything else is rejected.
    const std::string file_format = Settings["file_format"].GetString();
    KRATOS_ERROR_IF(file_format != "ascii")
        << "Invalid file_format: '" << file_format << "'. UNV is an ASCII-only format; use 'ascii'." << std::endl;

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

    // General output settings
    mDefaultPrecision = Settings["output_precision"].GetInt();
    mDeformationFactor = Settings["deformation_factor"].GetDouble();
    mDecomposeQuadraticIntoLinear = Settings["decompose_quadratic_into_linear"].GetBool();
    mWriteDeformedConfiguration = Settings["write_deformed_configuration"].GetBool();

    // Unit system declared in the UNV file (validated against the known systems)
    mUnitSystem = Settings["unit_system"].GetString();
    UnitSystemData unit_system_data;
    KRATOS_ERROR_IF_NOT(TryGetUnitSystem(mUnitSystem, unit_system_data))
        << "Unknown unit_system: '" << mUnitSystem << "'. Currently supported: 'SI'." << std::endl;
    mWriteIds = Settings["write_ids"].GetBool();
    mOutputSubModelParts = Settings["output_sub_model_parts"].GetBool();

    // Initialize the lists of historical and non-historical nodal variables to be printed based on the provided settings.
    mHistoricalVariables.Initialize(Settings["nodal_solution_step_data_variables"].GetStringArray(), mUnvVariableKeys);
    mNonHistoricalVariables.Initialize(Settings["nodal_data_value_variables"].GetStringArray(), mUnvVariableKeys);

    // Initialize the lists of element and condition (non-historical) variables to be printed.
    mElementVariables.Initialize(Settings["element_data_value_variables"].GetStringArray(), mUnvVariableKeys);
    mConditionVariables.Initialize(Settings["condition_data_value_variables"].GetStringArray(), mUnvVariableKeys);

    // Initialize the list of gauss point variables to be averaged in elements.
    mGaussPointVariables.Initialize(Settings["gauss_point_variables_in_elements"].GetStringArray(), mUnvVariableKeys);

    // Setup the gauss point extrapolation to nodes process, if requested. The extrapolated values are
    // written as non-historical nodal data (mirroring VtkOutput).
    const std::vector<std::string> gauss_to_nodes = Settings["gauss_point_variables_extrapolated_to_nodes"].GetStringArray();
    if (gauss_to_nodes.size() > 0) {
        Parameters extrapolation_parameters(R"({
            "echo_level"                 : 0,
            "area_average"               : true,
            "average_variable"           : "NODAL_AREA",
            "list_of_variables"          : [],
            "extrapolate_non_historical" : true
        })");
        for (const auto& r_name : gauss_to_nodes) {
            extrapolation_parameters["list_of_variables"].Append(r_name);
        }
        mpGaussToNodesProcess = Kratos::make_unique<IntegrationValuesExtrapolationToNodesProcess>(mrOutputModelPart, extrapolation_parameters);
        mNonHistoricalVariables.Initialize(gauss_to_nodes, mUnvVariableKeys);
    }

    // Initialize the flag lists.
    InitializeFlags(Settings["nodal_flags"].GetStringArray(), mNodalFlags);
    InitializeFlags(Settings["element_flags"].GetStringArray(), mElementFlags);
    InitializeFlags(Settings["condition_flags"].GetStringArray(), mConditionFlags);
}

UnvOutput::UnvOutput(
    Kratos::ModelPart& rModelPart,
    const std::string& rOutFileWithoutExtension
    ) : mrOutputModelPart(rModelPart),
        mOutputFileName(rOutFileWithoutExtension + ".unv") {
}

void UnvOutput::InitializeFlags(
    const std::vector<std::string>& rFlagNames,
    std::vector<std::pair<const Flags*, std::string>>& rFlagList
    )
{
    for (const auto& r_flag_name : rFlagNames) {
        KRATOS_ERROR_IF_NOT(KratosComponents<Flags>::Has(r_flag_name))
            << "The flag '" << r_flag_name << "' is not registered in KratosComponents." << std::endl;
        rFlagList.emplace_back(&KratosComponents<Flags>::Get(r_flag_name), r_flag_name);
    }
}

void UnvOutput::WriteMesh() {
    // Check if the output file has been initialized. If not, initialize it.
    if (!mInitializedOutputFile) {
        KRATOS_WARNING("UnvOutput") << "Output file has not been initialized yet. Initializing output file before writing mesh." << std::endl;
        InitializeOutputFile();
    }

    // Write the units dataset (declares the model unit system, e.g. SI)
    WriteUnitsDataset();

    // Write the nodes
    WriteNodes();

    // Write the geometry (elements or conditions) based on availability
    if (mEntityType == EntityType::AUTOMATIC) {
        if (mrOutputModelPart.Elements().size() > 0) {
            // Write the elements if they exist in the model part
            WriteEntities(mrOutputModelPart.Elements());
        } else if (mrOutputModelPart.Conditions().size() > 0) {
            KRATOS_WARNING("UnvOutput") << "No elements found in the model part. Writing conditions instead." << std::endl;

            // Write the conditions if they exist in the model part
            WriteEntities(mrOutputModelPart.Conditions());
        } else {
            KRATOS_WARNING("UnvOutput") << "No elements or conditions found in the model part. No mesh will be written." << std::endl;
        }
    } else if (mEntityType == EntityType::ELEMENTS) {
        if (mrOutputModelPart.Elements().size() > 0) {
            // Write the elements if they exist in the model part
            WriteEntities(mrOutputModelPart.Elements());
        } else {
            KRATOS_WARNING("UnvOutput") << "No elements found in the model part. No mesh will be written." << std::endl;
        }
    } else if (mEntityType == EntityType::CONDITIONS) {
        if (mrOutputModelPart.Conditions().size() > 0) {
            // Write the conditions if they exist in the model part
            WriteEntities(mrOutputModelPart.Conditions());
        } else {
            KRATOS_WARNING("UnvOutput") << "No conditions found in the model part. No mesh will be written." << std::endl;
        }
    }

    // Write the sub model parts as UNV groups, if requested
    if (mOutputSubModelParts) {
        WriteGroups();
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

    // Run the gauss point extrapolation to nodes process, if configured.
    PrepareGaussPointResults();

    // Get the current time step from the model part's process info
    const auto& r_process_info = mrOutputModelPart.GetProcessInfo();
    const double time_step = r_process_info.GetValue(TIME);

    // Print output for nodal historical variables
    for (const auto& p_variable : mHistoricalVariables.mBoolVariables) {
        WriteResultRecords<Variable<bool>, WriteType::HISTORICAL, ResultLocation::NODES>(*p_variable, 1, time_step);
    }
    for (const auto& p_variable : mHistoricalVariables.mIntVariables) {
        WriteResultRecords<Variable<int>, WriteType::HISTORICAL, ResultLocation::NODES>(*p_variable, 1, time_step);
    }
    for (const auto& p_variable : mHistoricalVariables.mDoubleVariables) {
        WriteResultRecords<Variable<double>, WriteType::HISTORICAL, ResultLocation::NODES>(*p_variable, 1, time_step);
    }
    for (const auto& p_variable : mHistoricalVariables.mArray1DVariables) {
        WriteResultRecords<Variable<array_1d<double,3>>, WriteType::HISTORICAL, ResultLocation::NODES>(*p_variable, 3, time_step);
    }

    // Print output for nodal non-historical variables
    for (const auto& p_variable : mNonHistoricalVariables.mBoolVariables) {
        WriteResultRecords<Variable<bool>, WriteType::NON_HISTORICAL, ResultLocation::NODES>(*p_variable, 1, time_step);
    }
    for (const auto& p_variable : mNonHistoricalVariables.mIntVariables) {
        WriteResultRecords<Variable<int>, WriteType::NON_HISTORICAL, ResultLocation::NODES>(*p_variable, 1, time_step);
    }
    for (const auto& p_variable : mNonHistoricalVariables.mDoubleVariables) {
        WriteResultRecords<Variable<double>, WriteType::NON_HISTORICAL, ResultLocation::NODES>(*p_variable, 1, time_step);
    }
    for (const auto& p_variable : mNonHistoricalVariables.mArray1DVariables) {
        WriteResultRecords<Variable<array_1d<double,3>>, WriteType::NON_HISTORICAL, ResultLocation::NODES>(*p_variable, 3, time_step);
    }

    // Print output for element (non-historical) variables
    for (const auto& p_variable : mElementVariables.mBoolVariables) {
        WriteResultRecords<Variable<bool>, WriteType::NON_HISTORICAL, ResultLocation::ELEMENTS>(*p_variable, 1, time_step);
    }
    for (const auto& p_variable : mElementVariables.mIntVariables) {
        WriteResultRecords<Variable<int>, WriteType::NON_HISTORICAL, ResultLocation::ELEMENTS>(*p_variable, 1, time_step);
    }
    for (const auto& p_variable : mElementVariables.mDoubleVariables) {
        WriteResultRecords<Variable<double>, WriteType::NON_HISTORICAL, ResultLocation::ELEMENTS>(*p_variable, 1, time_step);
    }
    for (const auto& p_variable : mElementVariables.mArray1DVariables) {
        WriteResultRecords<Variable<array_1d<double,3>>, WriteType::NON_HISTORICAL, ResultLocation::ELEMENTS>(*p_variable, 3, time_step);
    }

    // Print output for condition (non-historical) variables
    for (const auto& p_variable : mConditionVariables.mBoolVariables) {
        WriteResultRecords<Variable<bool>, WriteType::NON_HISTORICAL, ResultLocation::CONDITIONS>(*p_variable, 1, time_step);
    }
    for (const auto& p_variable : mConditionVariables.mIntVariables) {
        WriteResultRecords<Variable<int>, WriteType::NON_HISTORICAL, ResultLocation::CONDITIONS>(*p_variable, 1, time_step);
    }
    for (const auto& p_variable : mConditionVariables.mDoubleVariables) {
        WriteResultRecords<Variable<double>, WriteType::NON_HISTORICAL, ResultLocation::CONDITIONS>(*p_variable, 1, time_step);
    }
    for (const auto& p_variable : mConditionVariables.mArray1DVariables) {
        WriteResultRecords<Variable<array_1d<double,3>>, WriteType::NON_HISTORICAL, ResultLocation::CONDITIONS>(*p_variable, 3, time_step);
    }

    // Print output for gauss point variables averaged in elements
    for (const auto& p_variable : mGaussPointVariables.mDoubleVariables) {
        WriteGaussPointElementResults<double>(*p_variable, 1, time_step);
    }
    for (const auto& p_variable : mGaussPointVariables.mArray1DVariables) {
        WriteGaussPointElementResults<array_1d<double,3>>(*p_variable, 3, time_step);
    }

    // Print the flags
    for (const auto& r_flag : mNodalFlags) {
        WriteFlagRecords<ResultLocation::NODES>(*r_flag.first, r_flag.second, mrOutputModelPart.Nodes(), time_step);
    }
    for (const auto& r_flag : mElementFlags) {
        WriteFlagRecords<ResultLocation::ELEMENTS>(*r_flag.first, r_flag.second, mrOutputModelPart.Elements(), time_step);
    }
    for (const auto& r_flag : mConditionFlags) {
        WriteFlagRecords<ResultLocation::CONDITIONS>(*r_flag.first, r_flag.second, mrOutputModelPart.Conditions(), time_step);
    }

    // Print the entity ids, if requested
    if (mWriteIds) {
        WriteIdRecords<ResultLocation::NODES>("NODE_ID", mrOutputModelPart.Nodes(), time_step);
        WriteIdRecords<ResultLocation::ELEMENTS>("ELEMENT_ID", mrOutputModelPart.Elements(), time_step);
        WriteIdRecords<ResultLocation::CONDITIONS>("CONDITION_ID", mrOutputModelPart.Conditions(), time_step);
    }
}

void UnvOutput::PrepareGaussPointResults() {
    if (mpGaussToNodesProcess != nullptr) {
        mpGaussToNodesProcess->Execute();
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

bool UnvOutput::TryGetUnitSystem(
    const std::string& rName,
    UnitSystemData& rOut
    )
{
    if (rName == "SI") {
        // SI: Meter (newton). All conversion factors are 1.0 (values are already in SI).
        rOut = {1, "SI: Meter (newton)", 1.0, 1.0, 1.0, 273.15, 2};
        return true;
    }
    // Additional systems (e.g. mm) can be added here.
    return false;
}

void UnvOutput::WriteUnitsDataset() {
    UnitSystemData unit_system_data;
    KRATOS_ERROR_IF_NOT(TryGetUnitSystem(mUnitSystem, unit_system_data))
        << "Unknown unit_system: '" << mUnitSystem << "'." << std::endl;

    std::ofstream rOutputFile;
    rOutputFile.open(mOutputFileName, std::ios::out | std::ios::app);
    rOutputFile << std::scientific << std::setprecision(15);

    rOutputFile << std::setw(6) << "-1" << "\n";
    rOutputFile << std::setw(6) << as_integer(DatasetID::UNITS_DATASET) << "\n";

    // Record 1: units code (I10), units description (20A1, left-justified), temperature mode (I10)
    rOutputFile << std::setw(10) << unit_system_data.code << "  " << std::left << std::setw(20)
                << unit_system_data.description << std::right << std::setw(10) << unit_system_data.temperature_mode << "\n";

    // Record 2: length, force and temperature conversion factors to SI
    rOutputFile << std::setw(25) << unit_system_data.length_factor;
    rOutputFile << std::setw(25) << unit_system_data.force_factor;
    rOutputFile << std::setw(25) << unit_system_data.temperature_factor << "\n";

    // Record 3: temperature offset
    rOutputFile << std::setw(25) << unit_system_data.temperature_offset << "\n";

    rOutputFile << std::setw(6) << "-1" << "\n";
    rOutputFile.close();
}

void UnvOutput::WriteNodes() {
    std::ofstream rOutputFile;
    rOutputFile.open(mOutputFileName, std::ios::out | std::ios::app);

    rOutputFile << std::scientific;
    rOutputFile << std::setprecision(15);

    const int export_coordinate_system = 0;
    const int displacement_coordinate_system_number = 0;
    const int color = 0;

    // Only read DISPLACEMENT if the deformed configuration is requested and the variable is available.
    const bool write_deformed = mWriteDeformedConfiguration && mrOutputModelPart.HasNodalSolutionStepVariable(DISPLACEMENT);
    KRATOS_WARNING_IF("UnvOutput", mWriteDeformedConfiguration && !write_deformed)
        << "write_deformed_configuration is set but DISPLACEMENT is not part of the nodal solution step variables. Writing the reference configuration." << std::endl;

    rOutputFile << std::setw(6) << "-1" << "\n";
    rOutputFile << std::setw(6) << as_integer(DatasetID::NODES_DATASET) << "\n";

    int node_label;
    double x_coordinate, y_coordinate, z_coordinate;
    array_1d<double, 3> node_displacement = ZeroVector(3);
    for (auto& r_node : mrOutputModelPart.Nodes()) {
        node_label = r_node.Id();
        // Get the node displacement if the deformed configuration is to be written (scaled by the deformation factor)
        if (write_deformed) {
            noalias(node_displacement) = mDeformationFactor * r_node.FastGetSolutionStepValue(DISPLACEMENT);
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

bool UnvOutput::TryGetUnvDescriptor(
    const GeometryData::KratosGeometryType Type,
    UnvElementDescriptor& rOut
    )
{
    using GeometryType = GeometryData::KratosGeometryType;
    switch (Type) {
        // Lines / beams (need the extra orientation record)
        case GeometryType::Kratos_Line2D2:
        case GeometryType::Kratos_Line3D2:
            rOut = {21, 21, 2, true, {0, 1}, false};
            return true;
        case GeometryType::Kratos_Line2D3:
        case GeometryType::Kratos_Line3D3:
            rOut = {24, 24, 3, true, {0, 2, 1}, false};
            return true;
        // Triangles
        case GeometryType::Kratos_Triangle2D3:
            rOut = {41, 91, 3, false, {0, 1, 2}, false};
            return true;
        case GeometryType::Kratos_Triangle3D3:
            rOut = {41, 91, 3, false, {0, 1, 2}, false};
            return true;
        case GeometryType::Kratos_Triangle2D6:
        case GeometryType::Kratos_Triangle3D6:
            rOut = {42, 92, 6, false, {0, 3, 1, 4, 2, 5}, false};
            return true;
        // Quadrilaterals
        case GeometryType::Kratos_Quadrilateral2D4:
        case GeometryType::Kratos_Quadrilateral3D4:
            rOut = {44, 94, 4, false, {0, 1, 2, 3}, false};
            return true;
        case GeometryType::Kratos_Quadrilateral2D8:
        case GeometryType::Kratos_Quadrilateral3D8:
            rOut = {45, 95, 8, false, {0, 4, 1, 5, 2, 6, 3, 7}, false};
            return true;
        case GeometryType::Kratos_Quadrilateral2D9:
        case GeometryType::Kratos_Quadrilateral3D9:
            // UNV has no 9-node quad: degrade to the 8-node parabolic quad (drop the center node 8).
            rOut = {45, 95, 8, false, {0, 4, 1, 5, 2, 6, 3, 7}, true};
            return true;
        // Tetrahedra
        case GeometryType::Kratos_Tetrahedra3D4:
            rOut = {-1, 111, 4, false, {0, 1, 2, 3}, false};
            return true;
        case GeometryType::Kratos_Tetrahedra3D10:
            rOut = {-1, 118, 10, false, {0, 4, 1, 5, 2, 6, 7, 8, 9, 3}, false};
            return true;
        // Hexahedra
        case GeometryType::Kratos_Hexahedra3D8:
            rOut = {-1, 115, 8, false, {0, 1, 2, 3, 4, 5, 6, 7}, false};
            return true;
        case GeometryType::Kratos_Hexahedra3D20:
            rOut = {-1, 116, 20, false, {0, 8, 1, 9, 2, 10, 3, 11, 12, 13, 14, 15, 4, 16, 5, 17, 6, 18, 7, 19}, false};
            return true;
        case GeometryType::Kratos_Hexahedra3D27:
            // UNV has no 27-node brick: degrade to the 20-node parabolic brick (drop nodes 20-26).
            rOut = {-1, 116, 20, false, {0, 8, 1, 9, 2, 10, 3, 11, 12, 13, 14, 15, 4, 16, 5, 17, 6, 18, 7, 19}, true};
            return true;
        // Prisms / wedges
        case GeometryType::Kratos_Prism3D6:
            rOut = {-1, 112, 6, false, {0, 1, 2, 3, 4, 5}, false};
            return true;
        case GeometryType::Kratos_Prism3D15:
            // Parabolic wedge (FE 113); best-effort ordering (not in the FEconv reference).
            rOut = {-1, 113, 15, false, {0, 6, 1, 7, 2, 8, 9, 10, 11, 3, 12, 4, 13, 5, 14}, false};
            return true;
        default:
            // Pyramids and any other geometry have no UNV representation.
            return false;
    }
}

bool UnvOutput::TryGetLinearDecomposition(
    const GeometryData::KratosGeometryType Type,
    UnvLinearDecomposition& rOut
    )
{
    using GeometryType = GeometryData::KratosGeometryType;
    switch (Type) {
        // Lines / beams
        case GeometryType::Kratos_Line2D2:
        case GeometryType::Kratos_Line3D2:
            rOut = {21, 21, 2, true, false, {{0, 1}}};
            return true;
        case GeometryType::Kratos_Line2D3:
        case GeometryType::Kratos_Line3D3:
            // Split the 3-node line into two linear segments.
            rOut = {21, 21, 2, true, false, {{0, 2}, {2, 1}}};
            return true;
        // Triangles
        case GeometryType::Kratos_Triangle2D3:
        case GeometryType::Kratos_Triangle3D3:
            rOut = {41, 91, 3, false, false, {{0, 1, 2}}};
            return true;
        case GeometryType::Kratos_Triangle2D6:
        case GeometryType::Kratos_Triangle3D6:
            // Split the 6-node triangle into four linear triangles.
            rOut = {41, 91, 3, false, false, {{0, 3, 5}, {3, 1, 4}, {5, 4, 2}, {3, 4, 5}}};
            return true;
        // Quadrilaterals
        case GeometryType::Kratos_Quadrilateral2D4:
        case GeometryType::Kratos_Quadrilateral3D4:
            rOut = {44, 94, 4, false, false, {{0, 1, 2, 3}}};
            return true;
        case GeometryType::Kratos_Quadrilateral2D8:
        case GeometryType::Kratos_Quadrilateral3D8:
            // Serendipity quad (no center node): split into six linear triangles.
            rOut = {41, 91, 3, false, false, {{0, 4, 7}, {4, 1, 5}, {5, 2, 6}, {6, 3, 7}, {4, 5, 6}, {4, 6, 7}}};
            return true;
        case GeometryType::Kratos_Quadrilateral2D9:
        case GeometryType::Kratos_Quadrilateral3D9:
            // Has a center node: split into four linear quadrilaterals.
            rOut = {44, 94, 4, false, false, {{0, 4, 8, 7}, {4, 1, 5, 8}, {8, 5, 2, 6}, {7, 8, 6, 3}}};
            return true;
        // Tetrahedra
        case GeometryType::Kratos_Tetrahedra3D4:
            rOut = {-1, 111, 4, false, false, {{0, 1, 2, 3}}};
            return true;
        case GeometryType::Kratos_Tetrahedra3D10:
            // Standard subdivision into eight linear tetrahedra (four corner tets plus the inner
            // octahedron split along the 4-9 edge).
            rOut = {-1, 111, 4, false, false, {
                {0, 4, 6, 7}, {4, 1, 5, 8}, {6, 5, 2, 9}, {7, 8, 9, 3},
                {4, 9, 6, 7}, {4, 9, 7, 8}, {4, 9, 8, 5}, {4, 9, 5, 6}}};
            return true;
        // Hexahedra
        case GeometryType::Kratos_Hexahedra3D8:
            rOut = {-1, 115, 8, false, false, {{0, 1, 2, 3, 4, 5, 6, 7}}};
            return true;
        case GeometryType::Kratos_Hexahedra3D20:
        case GeometryType::Kratos_Hexahedra3D27:
            // No conforming linear refinement without interior nodes: reduce to the linear corner brick.
            rOut = {-1, 115, 8, false, true, {{0, 1, 2, 3, 4, 5, 6, 7}}};
            return true;
        // Prisms / wedges
        case GeometryType::Kratos_Prism3D6:
            rOut = {-1, 112, 6, false, false, {{0, 1, 2, 3, 4, 5}}};
            return true;
        case GeometryType::Kratos_Prism3D15:
            // Reduce to the linear corner wedge (mid-side nodes dropped).
            rOut = {-1, 112, 6, false, true, {{0, 1, 2, 3, 4, 5}}};
            return true;
        default:
            // Pyramids and any other geometry have no UNV representation.
            return false;
    }
}

int UnvOutput::GetSubElementCount(const GeometryData::KratosGeometryType Type) const
{
    if (!mDecomposeQuadraticIntoLinear) {
        return 1;
    }
    UnvLinearDecomposition decomposition;
    if (TryGetLinearDecomposition(Type, decomposition)) {
        return static_cast<int>(decomposition.sub_elements.size());
    }
    return 1;
}

template<typename TContainerType>
void UnvOutput::WriteEntities(const TContainerType& rContainer) {
    std::ofstream rOutputFile;
    rOutputFile.open(mOutputFileName, std::ios::out | std::ios::app);

    const int physical_property_table_number = 1;
    const int material_property_table_number = 1;
    const int color = 0;

    rOutputFile << std::setw(6) << "-1" << "\n";
    rOutputFile << std::setw(6) << as_integer(DatasetID::ELEMENTS_DATASET) << "\n";

    // Helper writing a single 2412 element record (header, optional beam record and connectivity).
    auto write_record = [&](const long long Label, const int FeDescriptorId, const int NumberOfNodes,
                            const bool IsBeam, const std::vector<int>& rConnectivity, const auto& rGeometry) {
        rOutputFile << std::setw(10) << Label;
        rOutputFile << std::setw(10) << FeDescriptorId;
        rOutputFile << std::setw(10) << physical_property_table_number;
        rOutputFile << std::setw(10) << material_property_table_number;
        rOutputFile << std::setw(10) << color;
        rOutputFile << std::setw(10) << NumberOfNodes << "\n";

        // Beam entities require the extra orientation record between the header and the connectivity.
        if (IsBeam) {
            rOutputFile << std::setw(10) << 0 << std::setw(10) << 0 << std::setw(10) << 0 << "\n";
        }

        // Connectivity record: node ids wrapped at 8 per line.
        int count = 0;
        for (const int local_index : rConnectivity) {
            rOutputFile << std::setw(10) << static_cast<int>(rGeometry[local_index].Id());
            if (++count % 8 == 0) {
                rOutputFile << "\n";
            }
        }
        if (count % 8 != 0) {
            rOutputFile << "\n";
        }
    };

    std::set<int> warned_geometries;
    for (auto& r_entity : rContainer) {
        const int entity_label = static_cast<int>(r_entity.Id());
        const auto& r_geometry = r_entity.GetGeometry();
        const auto geometry_type = r_geometry.GetGeometryType();
        const bool is_3d = r_geometry.WorkingSpaceDimension() == 3;

        if (mDecomposeQuadraticIntoLinear) {
            // Decompose each geometry into one or more linear sub-elements (stricter reader support).
            UnvLinearDecomposition decomposition;
            KRATOS_ERROR_IF_NOT(TryGetLinearDecomposition(geometry_type, decomposition))
                << "Entity with ID " << entity_label << " has a geometry that cannot be represented in UNV "
                << "(e.g. pyramids, which have no UNV element type)." << std::endl;

            const int fe_descriptor_id = (is_3d && decomposition.fe_descriptor_3d != -1) ? decomposition.fe_descriptor_3d
                : (decomposition.fe_descriptor_2d != -1 ? decomposition.fe_descriptor_2d : decomposition.fe_descriptor_3d);

            if (decomposition.drops_nodes && warned_geometries.insert(static_cast<int>(geometry_type)).second) {
                KRATOS_WARNING("UnvOutput") << "Geometry of entity " << entity_label
                    << " has no linear refinement; writing the linear corner element (mid-side nodes are dropped)." << std::endl;
            }

            const auto labels = GetEntityLabels(r_entity);
            for (std::size_t k = 0; k < decomposition.sub_elements.size(); ++k) {
                write_record(labels[k], fe_descriptor_id, decomposition.num_nodes,
                    decomposition.is_beam, decomposition.sub_elements[k], r_geometry);
            }
        } else {
            // Write the geometry as-is (quadratic descriptors preserved).
            UnvElementDescriptor descriptor;
            KRATOS_ERROR_IF_NOT(TryGetUnvDescriptor(geometry_type, descriptor))
                << "Entity with ID " << entity_label << " has a geometry that cannot be represented in UNV "
                << "(e.g. pyramids, which have no UNV element type)." << std::endl;

            const int fe_descriptor_id = (is_3d && descriptor.fe_descriptor_3d != -1) ? descriptor.fe_descriptor_3d
                : (descriptor.fe_descriptor_2d != -1 ? descriptor.fe_descriptor_2d : descriptor.fe_descriptor_3d);

            if (descriptor.degrade && warned_geometries.insert(static_cast<int>(geometry_type)).second) {
                KRATOS_WARNING("UnvOutput") << "Geometry of entity " << entity_label
                    << " has no exact UNV representation; writing a degraded element with "
                    << descriptor.num_nodes << " nodes (interior nodes are dropped)." << std::endl;
            }

            write_record(static_cast<long long>(entity_label), fe_descriptor_id, descriptor.num_nodes,
                descriptor.is_beam, descriptor.reorder, r_geometry);
        }
    }
    rOutputFile << std::setw(6) << "-1" << "\n";
    rOutputFile.close();
}

// Explicit instantiations for the element and condition containers.
template void UnvOutput::WriteEntities<ModelPart::ElementsContainerType>(const ModelPart::ElementsContainerType&);
template void UnvOutput::WriteEntities<ModelPart::ConditionsContainerType>(const ModelPart::ConditionsContainerType&);

void UnvOutput::WriteGroups() {
    // Collect the sub model parts recursively.
    std::vector<const ModelPart*> sub_model_parts;
    std::function<void(const ModelPart&)> collect = [&](const ModelPart& rModelPart) {
        for (const auto& r_sub_model_part : rModelPart.SubModelParts()) {
            sub_model_parts.push_back(&r_sub_model_part);
            collect(r_sub_model_part);
        }
    };
    collect(mrOutputModelPart);

    if (sub_model_parts.empty()) {
        return;
    }

    std::ofstream rOutputFile;
    rOutputFile.open(mOutputFileName, std::ios::out | std::ios::app);

    rOutputFile << std::setw(6) << "-1" << "\n";
    rOutputFile << std::setw(6) << 2467 << "\n";

    int group_number = 1;
    for (const ModelPart* p_sub_model_part : sub_model_parts) {
        const auto& r_sub_model_part = *p_sub_model_part;

        // Count the total number of UNV entities (accounting for the decomposition into sub-elements).
        std::size_t number_of_entities = 0;
        for (const auto& r_element : r_sub_model_part.Elements()) {
            number_of_entities += GetEntityLabels(r_element).size();
        }
        for (const auto& r_condition : r_sub_model_part.Conditions()) {
            number_of_entities += GetEntityLabels(r_condition).size();
        }

        // Record 1: group header (group number, 6 active-set ids = 0, number of entities)
        rOutputFile << std::setw(10) << group_number;
        rOutputFile << std::setw(10) << 0 << std::setw(10) << 0 << std::setw(10) << 0;
        rOutputFile << std::setw(10) << 0 << std::setw(10) << 0 << std::setw(10) << 0;
        rOutputFile << std::setw(10) << static_cast<int>(number_of_entities) << "\n";

        // Record 2: group name
        rOutputFile << r_sub_model_part.Name() << "\n";

        // Records 3..N: entity tuples (entity type code, tag, 0, 0), two per line.
        // NOTE: elements and conditions are both written into the 2412 dataset with entity type code 8.
        // Kratos elements and conditions have independent id spaces, so a group referencing both may be
        // ambiguous if an element and a condition share the same id.
        int tuple_count = 0;
        auto write_tuple = [&](const int TypeCode, const long long Tag) {
            rOutputFile << std::setw(10) << TypeCode << std::setw(10) << Tag << std::setw(10) << 0 << std::setw(10) << 0;
            if (++tuple_count % 2 == 0) {
                rOutputFile << "\n";
            }
        };
        for (const auto& r_element : r_sub_model_part.Elements()) {
            for (const long long label : GetEntityLabels(r_element)) {
                write_tuple(8, label);
            }
        }
        for (const auto& r_condition : r_sub_model_part.Conditions()) {
            for (const long long label : GetEntityLabels(r_condition)) {
                write_tuple(8, label);
            }
        }
        if (tuple_count % 2 != 0) {
            rOutputFile << "\n";
        }
        ++group_number;
    }
    rOutputFile << std::setw(6) << "-1" << "\n";
    rOutputFile.close();
}

void UnvOutput::WriteResultDatasetHeader(
    std::ofstream& rOutputFile,
    const std::string& rLabel,
    const std::string& rName,
    const int Location,
    const int DataCharacteristic,
    const int UnvVariableId,
    const int NumberOfComponents,
    const double TimeStep
    )
{
    rOutputFile << std::setw(6) << "-1" << "\n";                                    // Begin block
    rOutputFile << std::setw(6) << as_integer(DatasetID::RESULTS_DATASET) << "\n";  // DatasetID

    rOutputFile << std::setw(10) << rLabel << "\n";                                 // Record 1 - Label
    rOutputFile << std::setw(6)  << rName << "\n";                                  // Record 2 - Name
    rOutputFile << std::setw(10) << Location << "\n";                               // Record 3 - Location

    // Records 4-8: free text
    rOutputFile << "\n" << "\n" << "\n" << "\n" << "\n";

    // Record 9: ModelType, AnalysisType, DataCharacteristic, ResultType, DataType, NumberOfDataValues
    rOutputFile << std::setw(10) << as_integer(ModelType::STRUCTURAL);
    rOutputFile << std::setw(10) << as_integer(AnalysisType::TRANSIENT);
    rOutputFile << std::setw(10) << DataCharacteristic;
    rOutputFile << std::setw(10) << UnvVariableId;
    rOutputFile << std::setw(10) << as_integer(DataType::SINGLE_PRECISION_FLOATING_POINT);
    rOutputFile << std::setw(10) << NumberOfComponents;
    rOutputFile << "\n";

    // Record 10: DesignSetId, IterationNumber, SolutionSetId, BoundaryCondition, LoadSet, ModeNumber, TimeStampNumber, FrequencyNumber
    rOutputFile << std::setw(10) << 0;
    rOutputFile << std::setw(10) << 0;
    rOutputFile << std::setw(10) << 0;
    rOutputFile << std::setw(10) << 0;
    rOutputFile << std::setw(10) << 1;
    rOutputFile << std::setw(10) << 1;
    rOutputFile << std::setw(10) << 1;
    rOutputFile << std::setw(10) << 0;
    rOutputFile << "\n";

    // Record 11: CreationOption, (Unknown)
    rOutputFile << std::setw(10) << 0;
    rOutputFile << std::setw(10) << 0;
    rOutputFile << "\n";

    // Record 12: Time, Frequency, Eigenvalue, NodalMass, ViscousDampingRatio, HystereticDampingRatio
    rOutputFile << std::scientific << std::setprecision(mDefaultPrecision);
    rOutputFile << std::setw(static_cast<int>(mDefaultPrecision) + 9) << TimeStep;
    rOutputFile << std::setw(13) << "0.00000E+00";
    rOutputFile << std::setw(13) << "0.00000E+00";
    rOutputFile << std::setw(13) << "0.00000E+00";
    rOutputFile << std::setw(13) << "0.00000E+00";
    rOutputFile << std::setw(13) << "0.00000E+00";
    rOutputFile << "\n";

    // Record 13: Eigenvalue_re, Eigenvalue_im, ModalA_re, ModalA_im, ModalB_re, ModalB_im
    rOutputFile << std::setw(13) << "0.00000E+00";
    rOutputFile << std::setw(13) << "0.00000E+00";
    rOutputFile << std::setw(13) << "0.00000E+00";
    rOutputFile << std::setw(13) << "0.00000E+00";
    rOutputFile << std::setw(13) << "0.00000E+00";
    rOutputFile << std::setw(13) << "0.00000E+00";
    rOutputFile << "\n";
}

void UnvOutput::WriteNodalResults(const Variable<bool>& rVariable, const double TimeStep) {
    WriteResultRecords<Variable<bool>, WriteType::HISTORICAL, ResultLocation::NODES>(rVariable, 1, TimeStep);
}

void UnvOutput::WriteNodalResults(const Variable<int>& rVariable, const double TimeStep) {
    WriteResultRecords<Variable<int>, WriteType::HISTORICAL, ResultLocation::NODES>(rVariable, 1, TimeStep);
}

void UnvOutput::WriteNodalResults(const Variable<double>& rVariable, const double TimeStep) {
    WriteResultRecords<Variable<double>, WriteType::HISTORICAL, ResultLocation::NODES>(rVariable, 1, TimeStep);
}

void UnvOutput::WriteNodalResults(const Variable<array_1d<double,3>>& rVariable, const double TimeStep) {
    WriteResultRecords<Variable<array_1d<double,3>>, WriteType::HISTORICAL, ResultLocation::NODES>(rVariable, 3, TimeStep);
}

void UnvOutput::WriteNodalResults(const Variable<Vector>& rVariable, const double TimeStep) {
    KRATOS_ERROR << "Dynamic Vector results are not yet supported in UNV" << std::endl;
}

void UnvOutput::WriteNodalResults(const Variable<Matrix>& rVariable, const double TimeStep) {
    KRATOS_ERROR << "Matrix results are not yet supported in UNV" << std::endl;
}

void UnvOutput::WriteNodalNonHistoricalResults(const Variable<bool>& rVariable, const double TimeStep) {
    WriteResultRecords<Variable<bool>, WriteType::NON_HISTORICAL, ResultLocation::NODES>(rVariable, 1, TimeStep);
}

void UnvOutput::WriteNodalNonHistoricalResults(const Variable<int>& rVariable, const double TimeStep) {
    WriteResultRecords<Variable<int>, WriteType::NON_HISTORICAL, ResultLocation::NODES>(rVariable, 1, TimeStep);
}

void UnvOutput::WriteNodalNonHistoricalResults(const Variable<double>& rVariable, const double TimeStep) {
    WriteResultRecords<Variable<double>, WriteType::NON_HISTORICAL, ResultLocation::NODES>(rVariable, 1, TimeStep);
}

void UnvOutput::WriteNodalNonHistoricalResults(const Variable<array_1d<double,3>>& rVariable, const double TimeStep) {
    WriteResultRecords<Variable<array_1d<double,3>>, WriteType::NON_HISTORICAL, ResultLocation::NODES>(rVariable, 3, TimeStep);
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
        "file_format"                                 : "ascii",
        "output_precision"                            : 7,
        "output_control_type"                         : "step",
        "output_interval"                             : 1.0,
        "output_sub_model_parts"                      : false,
        "output_path"                                 : "UNV_Output",
        "custom_name_prefix"                          : "",
        "custom_name_postfix"                         : "",
        "save_output_files_in_folder"                 : true,
        "entity_type"                                 : "automatic",
        "unit_system"                                 : "SI",
        "decompose_quadratic_into_linear"             : false,
        "write_deformed_configuration"                : false,
        "deformation_factor"                          : 1.0,
        "write_ids"                                   : false,
        "nodal_solution_step_data_variables"          : [],
        "nodal_data_value_variables"                  : [],
        "nodal_flags"                                 : [],
        "element_data_value_variables"                : [],
        "element_flags"                               : [],
        "condition_data_value_variables"              : [],
        "condition_flags"                             : [],
        "gauss_point_variables_extrapolated_to_nodes" : [],
        "gauss_point_variables_in_elements"           : []
    })");

    return default_parameters;
}

} // namespace Kratos
