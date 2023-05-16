//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Crescenzio
//

// System includes

// External includes

// Project includes
#include "particle_mechanics_application_variables.h"
#include "particle_vtk_output.h"
#include "includes/kratos_filesystem.h"

namespace Kratos
{

Parameters ParticleVtkOutput::GetDefaultParameters()
{
    // IMPORTANT: when "output_control_type" is "time", then paraview will not be able to group them
    Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"                             : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "output_entity"                               : "elements",
        "file_format"                                 : "binary",
        "output_precision"                            : 7,
        "output_control_type"                         : "step",
        "output_interval"                             : 1.0,
        "output_sub_model_parts"                      : false,
        "output_path"                                 : "MPM_VTK_Output",
        "custom_name_prefix"                          : "",
        "custom_name_postfix"                         : "",
        "save_output_files_in_folder"                 : true,
        "write_ids"                                   : false,
        "particle_data_value_variables"               : [],
        "particle_flags"                              : []
    })" );

    return default_parameters;
}

ParticleVtkOutput::ParticleVtkOutput(
    ModelPart& rModelPart,
    Parameters ThisParameters
    ) : VtkOutput(rModelPart)
{
    // Validate and assign default parameters for current class
    Parameters default_parameters = GetDefaultParameters();
    ThisParameters.ValidateAndAssignDefaults(default_parameters);
    mOutputSettings = ThisParameters;

    // Initialize other variables
    mDefaultPrecision = mOutputSettings["output_precision"].GetInt();
    const std::string file_format = mOutputSettings["file_format"].GetString();
    if (file_format == "ascii") {
        mFileFormat = VtkOutput::FileFormat::VTK_ASCII;
    } else if (file_format == "binary") {
        mFileFormat = VtkOutput::FileFormat::VTK_BINARY;
        // test for endian-format
        int num = 1;
        if (*(char *)&num == 1) {
            mShouldSwap = true;
        }
    } else {
        KRATOS_ERROR << "Option for \"file_format\": " << file_format << " not recognised!\n Possible output formats options are: \"ascii\", \"binary\"" << std::endl;
    }

    const std::string output_entities = mOutputSettings["output_entity"].GetString();
    if (output_entities == "elements") {
        mOutputEntities = ParticleVtkOutput::OutputEntities::ELEMENTS;
    } else if (output_entities == "conditions") {
        mOutputEntities = ParticleVtkOutput::OutputEntities::CONDITIONS;
    } else {
        KRATOS_ERROR << "ParticleVtkOutput. Option for \"output_entity\": " << output_entities <<
            " not recognised!\n Admissible options are: \"elements\", \"conditions\"" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ParticleVtkOutput::PrintOutput(const std::string& rOutputFilename)
{
    // For whole model part
    WriteModelPartToFile(mrModelPart, false, rOutputFilename);
    // For sub model parts
    const bool print_sub_model_parts = mOutputSettings["output_sub_model_parts"].GetBool();
    if (print_sub_model_parts) {
        for (auto& r_sub_model_part : mrModelPart.SubModelParts()) {
            WriteModelPartToFile(r_sub_model_part, true, rOutputFilename);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ParticleVtkOutput::WriteModelPartToFile(
    const ModelPart& rModelPart,
    const bool IsSubModelPart,
    const std::string& rOutputFilename
    )
{
    // Make the file stream object
    const std::string output_file_name = GetOutputFileName(rModelPart, IsSubModelPart, rOutputFilename);
    std::ofstream output_file;
    if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) {
        output_file.open(output_file_name, std::ios::out | std::ios::trunc);
        output_file << std::scientific;
        output_file << std::setprecision(mDefaultPrecision);
    } else if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_BINARY) {
        output_file.open(output_file_name, std::ios::out | std::ios::binary | std::ios::trunc);
    }

    KRATOS_ERROR_IF_NOT(output_file.is_open()) << "File \"" << output_file_name << "\" could not be opened!" << std::endl;

    WriteHeaderToFile(rModelPart, output_file);
    WriteMeshToFile(rModelPart, output_file);
    if (mOutputEntities == ParticleVtkOutput::OutputEntities::ELEMENTS) {
        WriteElementResultsToFile(rModelPart, output_file);
    } else if (mOutputEntities == ParticleVtkOutput::OutputEntities::CONDITIONS) {
        WriteConditionResultsToFile(rModelPart, output_file);
    }

    output_file.close();
}

/***********************************************************************************/
/***********************************************************************************/

void ParticleVtkOutput::WriteMeshToFile(
    const ModelPart& rModelPart,
    std::ofstream& rFileStream
    ) const
{
    WriteConditionsAndElementsToFile(rModelPart, rFileStream);
}

/***********************************************************************************/
/***********************************************************************************/

void ParticleVtkOutput::WriteConditionsAndElementsToFile(
    const ModelPart& rModelPart,
    std::ofstream& rFileStream
    ) const
{
    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();

    if (mOutputEntities == ParticleVtkOutput::OutputEntities::ELEMENTS) {
        const SizeType num_elements = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(static_cast<int>(r_local_mesh.NumberOfElements()));
        // Write points data
        rFileStream << "POINTS " << num_elements << " float\n";
        for (auto itr_element = r_local_mesh.ElementsBegin(); itr_element != r_local_mesh.ElementsEnd(); itr_element++) {
            std::vector<array_1d<double, 3>> mp_coord = { ZeroVector(3) };
            itr_element->CalculateOnIntegrationPoints(MP_COORD, mp_coord, rModelPart.GetProcessInfo());
            WriteVectorDataToFile(mp_coord[0], rFileStream);
            if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
        }
        // Write cells data
        rFileStream << "\nCELLS " << num_elements << " " << 2*num_elements << "\n";
        for (IndexType itr_element = 0; itr_element < num_elements; ++itr_element) {
            WriteScalarDataToFile((unsigned int)1, rFileStream);
            if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) rFileStream << " ";
            WriteScalarDataToFile((int)itr_element, rFileStream);
            if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
        }
        // Write cell types data
        rFileStream << "\nCELL_TYPES " << num_elements << "\n";
        for (IndexType itr_element = 0; itr_element < num_elements; ++itr_element) {
            WriteScalarDataToFile((unsigned int)1, rFileStream);
            if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
        }
    } else if (mOutputEntities == ParticleVtkOutput::OutputEntities::CONDITIONS) {
        const SizeType num_conditions = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(static_cast<int>(r_local_mesh.NumberOfConditions()));
        // Write points data
        rFileStream << "POINTS " << num_conditions << " float\n";
        for (auto itr_condition = r_local_mesh.ConditionsBegin(); itr_condition != r_local_mesh.ConditionsEnd(); itr_condition++) {
            std::vector<array_1d<double, 3>> mpc_coord = { ZeroVector(3) };
            itr_condition->CalculateOnIntegrationPoints(MPC_COORD, mpc_coord, rModelPart.GetProcessInfo());
            WriteVectorDataToFile(mpc_coord[0], rFileStream);
            if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
        }
        // Write cells data
        rFileStream << "\nCELLS " << num_conditions << " " << 2*num_conditions << "\n";
        for (IndexType itr_condition = 0; itr_condition < num_conditions; ++itr_condition) {
            WriteScalarDataToFile((unsigned int)1, rFileStream);
            if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) rFileStream << " ";
            WriteScalarDataToFile((int)itr_condition, rFileStream);
            if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
        }
        // Write cell types data
        rFileStream << "\nCELL_TYPES " << num_conditions << "\n";
        for (IndexType itr_condition = 0; itr_condition < num_conditions; ++itr_condition) {
            WriteScalarDataToFile((unsigned int)1, rFileStream);
            if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ParticleVtkOutput::WriteElementResultsToFile(
    const ModelPart& rModelPart,
    std::ofstream& rFileStream
    )
{
    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();
    const int num_elements = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(static_cast<int>(r_local_mesh.NumberOfElements()));

    Parameters element_data_value_variables = mOutputSettings["particle_data_value_variables"];
    Parameters element_flags = mOutputSettings["particle_flags"];

    SizeType counter_element_data_value_variables = 0;
    for (IndexType entry = 0; entry < element_data_value_variables.size(); ++entry) {
        const std::string& r_element_result_name = element_data_value_variables[entry].GetString();
        if (IsCompatibleVariable(r_element_result_name)) {
            ++counter_element_data_value_variables;
        }
    }

    if (num_elements > 0) {
        rFileStream << "POINT_DATA " << r_local_mesh.NumberOfElements() << "\n";
        const bool write_ids = mOutputSettings["write_ids"].GetBool();
        rFileStream << "FIELD FieldData " << counter_element_data_value_variables + element_flags.size() + (write_ids ? 2 : 0) << "\n";

        for (IndexType entry = 0; entry < element_data_value_variables.size(); ++entry) {
            const std::string& r_element_result_name = element_data_value_variables[entry].GetString();
            WriteGeometricalContainerResults(r_element_result_name,r_local_mesh.Elements(),rFileStream);
        }

        if (element_flags.size() > 0) {
            mrModelPart.GetCommunicator().SynchronizeElementalFlags();
        }

        for (IndexType entry = 0; entry < element_flags.size(); ++entry) {
            const std::string& r_element_result_name = element_flags[entry].GetString();
            const Flags flag = KratosComponents<Flags>::Get(r_element_result_name);
            WriteFlagContainerVariable(r_local_mesh.Elements(), flag, r_element_result_name, rFileStream);
        }

        if (write_ids) {
            WritePropertiesIdsToFile(r_local_mesh.Elements(), rFileStream);
            WriteIdsToFile(r_local_mesh.Elements(), "KRATOS_ELEMENT_ID", rFileStream);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ParticleVtkOutput::WriteConditionResultsToFile(
    const ModelPart& rModelPart,
    std::ofstream& rFileStream
    )
{
    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();
    const int num_conditions = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(static_cast<int>(r_local_mesh.NumberOfConditions()));

    Parameters condition_data_value_variables = mOutputSettings["particle_data_value_variables"];
    Parameters condition_flags = mOutputSettings["particle_flags"];

    SizeType counter_condition_data_value_variables = 0;
    for (IndexType entry = 0; entry < condition_data_value_variables.size(); ++entry) {
        const std::string& r_condition_result_name = condition_data_value_variables[entry].GetString();
        if (IsCompatibleVariable(r_condition_result_name)) {
            ++counter_condition_data_value_variables;
        }
    }

    if (num_conditions > 0) {
        rFileStream << "POINT_DATA " << r_local_mesh.NumberOfConditions() << "\n";
        const bool write_ids = mOutputSettings["write_ids"].GetBool();
        rFileStream << "FIELD FieldData " << counter_condition_data_value_variables + condition_flags.size() + (write_ids ? 2 : 0) << "\n";

        for (IndexType entry = 0; entry < condition_data_value_variables.size(); ++entry) {
            const std::string& r_condition_result_name = condition_data_value_variables[entry].GetString();
            WriteGeometricalContainerResults(r_condition_result_name,r_local_mesh.Conditions(),rFileStream);
        }

        for (IndexType entry = 0; entry < condition_flags.size(); ++entry) {
            const std::string& r_condition_result_name = condition_flags[entry].GetString();
            const Flags flag = KratosComponents<Flags>::Get(r_condition_result_name);
            WriteFlagContainerVariable(r_local_mesh.Conditions(), flag, r_condition_result_name, rFileStream);
        }

        if (write_ids) {
            WritePropertiesIdsToFile(r_local_mesh.Conditions(), rFileStream);
            WriteIdsToFile(r_local_mesh.Conditions(), "KRATOS_CONDITION_ID", rFileStream);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType>
void ParticleVtkOutput::WriteGeometricalContainerResults(
    const std::string& rVariableName,
    const TContainerType& rContainer,
    std::ofstream& rFileStream
    ) const
{
    if (KratosComponents<Variable<double>>::Has(rVariableName)){
        const Variable<double>& var_to_write = KratosComponents<Variable<double>>::Get(rVariableName);
        WriteScalarContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<bool>>::Has(rVariableName)){
        const Variable<bool>& var_to_write = KratosComponents<Variable<bool>>::Get(rVariableName);
        WriteScalarContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<int>>::Has(rVariableName)){
        const Variable<int>& var_to_write = KratosComponents<Variable<int>>::Get(rVariableName);
        WriteScalarContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)){
        const Variable<array_1d<double, 3>>& var_to_write = KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
        WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<Vector>>::Has(rVariableName)){
        const Variable<Vector>& var_to_write = KratosComponents<Variable<Vector>>::Get(rVariableName);
        WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
    } else {
        KRATOS_WARNING_ONCE(rVariableName) << mrModelPart.GetCommunicator().GetDataCommunicator()
            << "Variable \"" << rVariableName << "\" is "
            << "not suitable for ParticleVtkOutput, skipping it" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType>
void ParticleVtkOutput::WritePropertiesIdsToFile(
    const TContainerType& rContainer,
    std::ofstream& rFileStream
    ) const
{
    rFileStream << "PROPERTIES_ID" << " 1 " << rContainer.size() << "  int\n";
    for (const auto& r_entity : rContainer) {
        WriteScalarDataToFile((int)r_entity.GetProperties().Id(), rFileStream);
        if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) {
            rFileStream <<"\n";
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType>
void ParticleVtkOutput::WriteIdsToFile(
    const TContainerType& rContainer,
    const std::string& DataName,
    std::ofstream& rFileStream
    ) const
{
    rFileStream << DataName << " 1 " << rContainer.size() << "  int\n";
    for (const auto& r_entity : rContainer) {
        WriteScalarDataToFile((int)r_entity.Id(), rFileStream);
        if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) {
            rFileStream <<"\n";
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType, class TVarType>
void ParticleVtkOutput::WriteScalarContainerVariable(
    const TContainerType& rContainer,
    const Variable<TVarType>& rVariable,
    std::ofstream& rFileStream
    ) const
{
    if (rContainer.size() == 0) {
        return;
    }

    rFileStream << rVariable.Name() << " 1 " << rContainer.size() << "  float\n";

    std::vector<TVarType> result;
    const ProcessInfo dummy_process_info = ProcessInfo();

    for (auto& r_entity : rContainer) {
        r_entity.CalculateOnIntegrationPoints(rVariable, result, dummy_process_info);
        WriteScalarDataToFile((float)result[0], rFileStream);
        if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) {
            rFileStream <<"\n";
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType, class TVarType>
void ParticleVtkOutput::WriteVectorContainerVariable(
    const TContainerType& rContainer,
    const Variable<TVarType>& rVariable,
    std::ofstream& rFileStream
    ) const
{
    if (rContainer.size() == 0) {
        return;
    }

    const auto& r_process_info = mrModelPart.GetProcessInfo();
    std::vector<TVarType> result;
    rContainer.begin()->CalculateOnIntegrationPoints(rVariable, result, r_process_info);
    const int res_size = result[0].size();

    rFileStream << rVariable.Name() << " " << res_size << " " << rContainer.size() << "  float\n";

    for (auto& r_entity : rContainer) {
        r_entity.CalculateOnIntegrationPoints(rVariable, result, r_process_info);
        WriteVectorDataToFile(result[0], rFileStream);
        if (mFileFormat == ParticleVtkOutput::FileFormat::VTK_ASCII) {
            rFileStream <<"\n";
        }
    }
}

} // namespace Kratos
