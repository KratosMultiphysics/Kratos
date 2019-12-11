//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala, Philipp Bucher
//  Collaborator:    Vicente Mataix Ferrandiz
//
//

// System includes
#include <map>

// External includes

// Project includes
#include "vtk_output.h"
#include "containers/model.h"
#include "processes/fast_transfer_between_model_parts_process.h"

namespace Kratos
{
VtkOutput::VtkOutput(
    ModelPart& rModelPart,
    Parameters ThisParameters
    ) : mrModelPart(rModelPart),
        mOutputSettings(ThisParameters)
{
    // The default parameters
    Parameters default_parameters = GetDefaultParameters();
    mOutputSettings.ValidateAndAssignDefaults(default_parameters);

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

        for(auto const& gauss_var : mOutputSettings["gauss_point_variables_extrapolated_to_nodes"])
            mOutputSettings["nodal_data_value_variables"].Append(gauss_var);

        // Making the gauss point to nodes process if any gauss point result is requested for
        mpGaussToNodesProcess = Kratos::make_unique<IntegrationValuesExtrapolationToNodesProcess>(rModelPart, gauss_intergration_param_non_hist);
    }

    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();
    const auto& r_data_comm = rModelPart.GetCommunicator().GetDataCommunicator();

    const int num_elements = r_data_comm.SumAll(static_cast<int>(r_local_mesh.NumberOfElements()));
    const int num_conditions = r_data_comm.SumAll(static_cast<int>(r_local_mesh.NumberOfConditions()));

    KRATOS_WARNING_IF("VtkOutput", num_elements > 0 && num_conditions > 0) << r_data_comm << "Modelpart \"" << rModelPart.Name() << "\" has both elements and conditions.\nGiving precedence to elements and writing only elements!" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void VtkOutput::PrepareGaussPointResults()
{
    if(mOutputSettings["gauss_point_variables_extrapolated_to_nodes"].size() > 0){
        mpGaussToNodesProcess->Execute();
    }
}

/***********************************************************************************/
/***********************************************************************************/

void VtkOutput::PrintOutput()
{
    // For Gauss point results
    PrepareGaussPointResults();

    // For whole model part
    WriteModelPartToFile(mrModelPart, false);

    // For sub model parts
    const bool print_sub_model_parts = mOutputSettings["output_sub_model_parts"].GetBool();
    if (print_sub_model_parts) {
        for (auto& r_sub_model_part : mrModelPart.SubModelParts()) {

            const auto& r_local_mesh = mrModelPart.GetCommunicator().LocalMesh();
            const auto& r_data_comm = mrModelPart.GetCommunicator().GetDataCommunicator();

            const int num_nodes = r_data_comm.SumAll(static_cast<int>(r_local_mesh.NumberOfNodes()));
            const int num_elements = r_data_comm.SumAll(static_cast<int>(r_local_mesh.NumberOfElements()));
            const int num_conditions = r_data_comm.SumAll(static_cast<int>(r_local_mesh.NumberOfConditions()));

            if (num_nodes == 0 && (num_elements != 0 || num_conditions != 0)) {
                WriteModelPartWithoutNodesToFile(r_sub_model_part);
            } else if (num_nodes != 0) {
                WriteModelPartToFile(r_sub_model_part, true);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void VtkOutput::WriteModelPartToFile(const ModelPart& rModelPart, const bool IsSubModelPart)
{
    Initialize(rModelPart);

    // Make the file stream object
    const std::string output_file_name = GetOutputFileName(rModelPart, IsSubModelPart);
    std::ofstream output_file;
    if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) {
        output_file.open(output_file_name, std::ios::out | std::ios::trunc);
        output_file << std::scientific;
        output_file << std::setprecision(mDefaultPrecision);
    } else if (mFileFormat == VtkOutput::FileFormat::VTK_BINARY) {
        output_file.open(output_file_name, std::ios::out | std::ios::binary | std::ios::trunc);
    }

    WriteHeaderToFile(rModelPart, output_file);
    WriteMeshToFile(rModelPart, output_file);
    WriteNodalResultsToFile(rModelPart, output_file);
    WriteElementResultsToFile(rModelPart, output_file);
    WriteConditionResultsToFile(rModelPart, output_file);

    output_file.close();
}

/***********************************************************************************/
/***********************************************************************************/

std::string VtkOutput::GetOutputFileName(const ModelPart& rModelPart, const bool IsSubModelPart)
{
    const int rank = rModelPart.GetCommunicator().MyPID();
    std::string model_part_name;

    if (IsSubModelPart) {
        model_part_name = rModelPart.GetParentModelPart()->Name() + "_" + rModelPart.Name();
    } else {
        model_part_name = rModelPart.Name();
    }

    std::string label;
    std::stringstream ss;
    const std::string output_control = mOutputSettings["output_control_type"].GetString();
    if (output_control == "step") {
        ss << std::fixed << std::setprecision(mDefaultPrecision)<< std::setfill('0')
           << rModelPart.GetProcessInfo()[STEP];
        label = ss.str();
    } else if(output_control == "time") {
        ss << std::fixed << std::setprecision(mDefaultPrecision) << std::setfill('0')
           << rModelPart.GetProcessInfo()[TIME];
        label = ss.str();
    } else {
        KRATOS_ERROR << "Option for \"output_control_type\": " << output_control
            <<" not recognised!\nPossible output_control_type options "
            << "are: \"step\", \"time\"" << std::endl;
    }

    // Putting everything together
    std::string output_file_name;
    if (mOutputSettings["save_output_files_in_folder"].GetBool()) {
        output_file_name += mOutputSettings["folder_name"].GetString() + "/";
    }
    const std::string& r_custom_name_prefix = mOutputSettings["custom_name_prefix"].GetString();
    const std::string& r_custom_name_postfix = mOutputSettings["custom_name_postfix"].GetString();
    output_file_name += r_custom_name_prefix + model_part_name + r_custom_name_postfix + "_" + std::to_string(rank) + "_" + label + ".vtk";

    return output_file_name;
}

/***********************************************************************************/
/***********************************************************************************/

void VtkOutput::Initialize(const ModelPart& rModelPart)
{
    CreateMapFromKratosIdToVTKId(rModelPart);
}

/***********************************************************************************/
/***********************************************************************************/

void VtkOutput::CreateMapFromKratosIdToVTKId(const ModelPart& rModelPart)
{
    int vtk_id = 0;
    for(const auto& r_node : rModelPart.Nodes()) {
        mKratosIdToVtkId[r_node.Id()] = vtk_id++;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void VtkOutput::WriteHeaderToFile(const ModelPart& rModelPart, std::ofstream& rFileStream) const
{
    rFileStream << "# vtk DataFile Version 4.0"
                << "\n"
                << "vtk output"
                << "\n";
    if(mFileFormat == VtkOutput::FileFormat::VTK_ASCII) {
        rFileStream << "ASCII" << "\n";
    }
    else if (mFileFormat == VtkOutput::FileFormat::VTK_BINARY) {
        rFileStream << "BINARY" << "\n";
    }
    rFileStream << "DATASET UNSTRUCTURED_GRID"
                << "\n";
}

/***********************************************************************************/
/***********************************************************************************/

void VtkOutput::WriteMeshToFile(const ModelPart& rModelPart, std::ofstream& rFileStream) const
{
    WriteNodesToFile(rModelPart, rFileStream);
    WriteConditionsAndElementsToFile(rModelPart, rFileStream);
}

/***********************************************************************************/
/***********************************************************************************/

void VtkOutput::WriteNodesToFile(const ModelPart& rModelPart, std::ofstream& rFileStream) const
{
    // NOTE: also in MPI all nodes (local and ghost) have to be written, because
    // they might be needed by the elements/conditions due to the connectivity

    // Write nodes header
    rFileStream << "POINTS " << rModelPart.NumberOfNodes() << " float\n";

    // Write nodes
    if (mOutputSettings["write_deformed_configuration"].GetBool()) {
        for (const auto& r_node : rModelPart.Nodes()) {
            WriteVectorDataToFile(r_node.Coordinates(), rFileStream);
            if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
        }
    } else {
        for (const auto& r_node : rModelPart.Nodes()) {
            WriteVectorDataToFile(r_node.GetInitialPosition(), rFileStream);
            if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void VtkOutput::WriteConditionsAndElementsToFile(const ModelPart& rModelPart, std::ofstream& rFileStream) const
{
    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();

    const int num_elements = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(static_cast<int>(r_local_mesh.NumberOfElements()));
    const int num_conditions = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(static_cast<int>(r_local_mesh.NumberOfConditions()));

    if (num_elements > 0) {
        // write cells header
        rFileStream << "\nCELLS " << r_local_mesh.NumberOfElements() << " "
            << DetermineVtkCellListSize(r_local_mesh.Elements()) << "\n";
        WriteConnectivity(r_local_mesh.Elements(), rFileStream);
        // write cell types header
        rFileStream << "\nCELL_TYPES " << r_local_mesh.NumberOfElements() << "\n";
        WriteCellType(r_local_mesh.Elements(), rFileStream);
    } else if (num_conditions > 0) {
        // write cells header
        rFileStream << "\nCELLS " << r_local_mesh.NumberOfConditions() << " "
            << DetermineVtkCellListSize(r_local_mesh.Conditions()) << "\n";
        WriteConnectivity(r_local_mesh.Conditions(), rFileStream);
        // write cell types header
        rFileStream << "\nCELL_TYPES " << r_local_mesh.NumberOfConditions() << "\n";
        WriteCellType(r_local_mesh.Conditions(), rFileStream);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType>
unsigned int VtkOutput::DetermineVtkCellListSize(const TContainerType& rContainer) const
{
    unsigned int vtk_cell_list_size_container = 0;

    const auto container_begin = rContainer.begin();
    const int num_entities = static_cast<int>(rContainer.size());
    #pragma omp parallel for reduction(+:vtk_cell_list_size_container)
    for (int i=0; i<num_entities; ++i) {
        const auto entity_i = container_begin + i;
        vtk_cell_list_size_container += entity_i->GetGeometry().PointsNumber() + 1;
    }

    return vtk_cell_list_size_container;
}

/***********************************************************************************/
/***********************************************************************************/

template <typename TContainerType>
void VtkOutput::WriteConnectivity(const TContainerType& rContainer, std::ofstream& rFileStream) const
{
    // NOTE: also in MPI all nodes (local and ghost) have to be written, because
    // they might be needed by the elements/conditions due to the connectivity

    const auto& r_id_map = mKratosIdToVtkId; // const reference to not accidentially modify the map
    for (const auto& r_entity : rContainer) {
        const auto& r_geom = r_entity.GetGeometry();
        const unsigned int number_of_nodes = r_geom.size();

        WriteScalarDataToFile((unsigned int)number_of_nodes, rFileStream);
        for (const auto& r_node : r_geom) {
            if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << " ";
            int id = r_id_map.at(r_node.Id());
            WriteScalarDataToFile((int)id, rFileStream);
        }
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <typename TContainerType>
void VtkOutput::WriteCellType(const TContainerType& rContainer, std::ofstream& rFileStream) const
{
    // IMPORTANT: The map geo_type_vtk_cell_type_map is to be extended to support new geometries
    // NOTE: See https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    const std::map<GeometryData::KratosGeometryType, int> geo_type_vtk_cell_type_map = {
        { GeometryData::KratosGeometryType::Kratos_Point2D,          1 },
        { GeometryData::KratosGeometryType::Kratos_Point3D,          1 },
        { GeometryData::KratosGeometryType::Kratos_Line2D2,          3 },
        { GeometryData::KratosGeometryType::Kratos_Line3D2,          3 },
        { GeometryData::KratosGeometryType::Kratos_Triangle2D3,      5 },
        { GeometryData::KratosGeometryType::Kratos_Triangle3D3,      5 },
        { GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4, 9 },
        { GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4, 9 },
        { GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4,    10 },
        { GeometryData::KratosGeometryType::Kratos_Hexahedra3D8,     12 },
        { GeometryData::KratosGeometryType::Kratos_Prism3D6,         13 },
        { GeometryData::KratosGeometryType::Kratos_Line2D3,          21 },
        { GeometryData::KratosGeometryType::Kratos_Line3D3,          21 },
        { GeometryData::KratosGeometryType::Kratos_Triangle2D6,      22 },
        { GeometryData::KratosGeometryType::Kratos_Triangle3D6,      22 },
        { GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8, 23 },
        { GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8, 23 },
        { GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10,   24 }
//         { GeometryData::KratosGeometryType::Kratos_Hexahedra3D20,    25 } // NOTE: Quadratic hexahedra (20) requires a conversor, order does not coincide with VTK
    };
    // Write entity types
    for (const auto& r_entity : rContainer) {
        int cell_type = -1;
        const auto& r_kratos_cell = r_entity.GetGeometry().GetGeometryType();
        if (geo_type_vtk_cell_type_map.count(r_kratos_cell) > 0) {
            cell_type = geo_type_vtk_cell_type_map.at(r_kratos_cell);
        } else {
            const auto& r_kratos_cell = r_entity.GetGeometry().GetGeometryType();
            KRATOS_ERROR << "Modelpart contains elements or conditions with "
             << "geometries for which no VTK-output is implemented!" << std::endl
             << "Cell type: " << static_cast<int>(r_kratos_cell) << std::endl;
        }

        WriteScalarDataToFile( (int)cell_type, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool VtkOutput::IsCompatibleVariable(const std::string& rVariableName) const
{
    if (KratosComponents<Variable<double>>::Has(rVariableName)){
        return true;
    } else if (KratosComponents<Variable<bool>>::Has(rVariableName)){
        return true;
    } else if (KratosComponents<Variable<int>>::Has(rVariableName)){
        return true;
    } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)){
        return true;
    } else if (KratosComponents<Variable<Vector>>::Has(rVariableName)){
        return true;
    } else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(rVariableName)){
        return true;
    } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)){
        return true;
    } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(rVariableName)){
        return true;
    } else {
        return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void VtkOutput::WriteNodalResultsToFile(const ModelPart& rModelPart, std::ofstream& rFileStream)
{
    // NOTE: also in MPI all nodes (local and ghost) have to be written, because
    // they might be needed by the elements/conditions due to the connectivity
    // Paraview needs a result on every node, therefore all results are written
    // this is why the synchronization is necessary

    // TODO perform synchronization of nodal results at the same time to
    // improve performance in MPI

    // write nodal results header
    Parameters nodal_solution_step_results = mOutputSettings["nodal_solution_step_data_variables"];
    Parameters nodal_variable_data_results = mOutputSettings["nodal_data_value_variables"];
    Parameters nodal_flags = mOutputSettings["nodal_flags"];
    const bool write_ids = mOutputSettings["write_ids"].GetBool();

    // Checking nodal_solution_step_results
    std::size_t counter_nodal_solution_step_results = 0;
    for (IndexType entry = 0; entry < nodal_solution_step_results.size(); ++entry) {
        // write nodal results variable header
        const std::string& r_nodal_result_name = nodal_solution_step_results[entry].GetString();
        if (IsCompatibleVariable(r_nodal_result_name)) ++counter_nodal_solution_step_results;
    }

    // Checking nodal_variable_data_results
    std::size_t counter_nodal_variable_data_results = 0;
    for (IndexType entry = 0; entry < nodal_variable_data_results.size(); ++entry) {
        // write nodal results variable header
        const std::string& r_nodal_result_name = nodal_variable_data_results[entry].GetString();
        if (IsCompatibleVariable(r_nodal_result_name)) ++counter_nodal_variable_data_results;
    }

    rFileStream << "POINT_DATA " << rModelPart.NumberOfNodes() << "\n";
    rFileStream << "FIELD FieldData " << counter_nodal_solution_step_results + counter_nodal_variable_data_results + nodal_flags.size() + (write_ids ? 1 : 0)  << "\n";

    // Writing nodal_solution_step_results
    for (IndexType entry = 0; entry < nodal_solution_step_results.size(); ++entry) {
        // write nodal results variable header
        const std::string& r_nodal_result_name = nodal_solution_step_results[entry].GetString();
        WriteNodalContainerResults(r_nodal_result_name, rModelPart.Nodes(), true, rFileStream);
    }

    // Writing nodal_variable_data_results
    for (IndexType entry = 0; entry < nodal_variable_data_results.size(); ++entry) {
        // write nodal results variable header
        const std::string& r_nodal_result_name = nodal_variable_data_results[entry].GetString();
        WriteNodalContainerResults(r_nodal_result_name, rModelPart.Nodes(), false, rFileStream);
    }

    // Writing nodal_flags
    if (nodal_flags.size() > 0) {
        mrModelPart.GetCommunicator().SynchronizeNodalFlags();
    }
    for (IndexType entry = 0; entry < nodal_flags.size(); ++entry) {
        // write nodal results variable header
        const std::string& r_nodal_result_name = nodal_flags[entry].GetString();
        const Flags flag = KratosComponents<Flags>::Get(r_nodal_result_name);
        WriteFlagContainerVariable(rModelPart.Nodes(), flag, r_nodal_result_name, rFileStream);
    }

    // If we write ids
    if (write_ids) {
        WriteIdsToFile(rModelPart.Nodes(), "KRATOS_NODE_ID", rFileStream);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void VtkOutput::WriteElementResultsToFile(const ModelPart& rModelPart, std::ofstream& rFileStream)
{
    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();
    Parameters element_data_value_variables = mOutputSettings["element_data_value_variables"];
    Parameters element_flags = mOutputSettings["element_flags"];
    Parameters gauss_point_variables_in_elements = mOutputSettings["gauss_point_variables_in_elements"];

    // Checking element_data_value_variables
    std::size_t counter_element_data_value_variables = 0;
    for (IndexType entry = 0; entry < element_data_value_variables.size(); ++entry) {
        // write nodal results variable header
        const std::string& r_element_result_name = element_data_value_variables[entry].GetString();
        if (IsCompatibleVariable(r_element_result_name)) ++counter_element_data_value_variables;
    }

    // Checking gauss_point_variables_in_elements
    std::size_t counter_gauss_point_variables_in_elements = 0;
    for (IndexType entry = 0; entry < gauss_point_variables_in_elements.size(); ++entry) {
        // write nodal results variable header
        const std::string& r_element_result_name = gauss_point_variables_in_elements[entry].GetString();
        if (IsCompatibleVariable(r_element_result_name)) ++counter_gauss_point_variables_in_elements;
    }

    const int num_elements = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(static_cast<int>(r_local_mesh.NumberOfElements()));

    if (num_elements > 0) {
        // write cells header
        rFileStream << "CELL_DATA " << r_local_mesh.NumberOfElements() << "\n";
        const bool write_ids = mOutputSettings["write_ids"].GetBool();
        rFileStream << "FIELD FieldData " << counter_element_data_value_variables + element_flags.size() + (write_ids ? 2 : 0) + counter_gauss_point_variables_in_elements << "\n";
        for (IndexType entry = 0; entry < element_data_value_variables.size(); ++entry) {
            const std::string& r_element_result_name = element_data_value_variables[entry].GetString();
            WriteGeometricalContainerResults(r_element_result_name,r_local_mesh.Elements(),rFileStream);
        }

        // Writing element_flags
        if (element_flags.size() > 0) {
            mrModelPart.GetCommunicator().SynchronizeElementalFlags();
        }
        for (IndexType entry = 0; entry < element_flags.size(); ++entry) {
            // Write elemental flags results variable header
            const std::string& r_element_result_name = element_flags[entry].GetString();
            const Flags flag = KratosComponents<Flags>::Get(r_element_result_name);
            WriteFlagContainerVariable(r_local_mesh.Elements(), flag, r_element_result_name, rFileStream);
        }

        // If we write ids
        if (write_ids) {
            WritePropertiesIdsToFile(r_local_mesh.Elements(), rFileStream);
            WriteIdsToFile(r_local_mesh.Elements(), "KRATOS_ELEMENT_ID", rFileStream);
        }

        // Direct write GP values
        for (IndexType entry = 0; entry < gauss_point_variables_in_elements.size(); ++entry) {
            const std::string& r_element_result_name = gauss_point_variables_in_elements[entry].GetString();
            WriteGeometricalContainerIntegrationResults(r_element_result_name,r_local_mesh.Elements(),rFileStream);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void VtkOutput::WriteConditionResultsToFile(const ModelPart& rModelPart, std::ofstream& rFileStream)
{
    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();
    Parameters condition_results = mOutputSettings["condition_data_value_variables"];
    Parameters condition_flags = mOutputSettings["condition_flags"];
    Parameters gauss_point_variables_in_elements = mOutputSettings["gauss_point_variables_in_elements"];

    // Checking condition_results
    std::size_t counter_condition_results = 0;
    for (IndexType entry = 0; entry < condition_results.size(); ++entry) {
        // write nodal results variable header
        const std::string& r_condition_result_name = condition_results[entry].GetString();
        if (IsCompatibleVariable(r_condition_result_name)) ++counter_condition_results;
    }

    // Checking gauss_point_variables_in_elements
    std::size_t counter_gauss_point_variables_in_elements = 0;
    for (IndexType entry = 0; entry < gauss_point_variables_in_elements.size(); ++entry) {
        // write nodal results variable header
        const std::string& r_element_result_name = gauss_point_variables_in_elements[entry].GetString();
        if (IsCompatibleVariable(r_element_result_name)) ++counter_gauss_point_variables_in_elements;
    }

    const int num_elements = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(static_cast<int>(r_local_mesh.NumberOfElements()));
    const int num_conditions = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(static_cast<int>(static_cast<int>(r_local_mesh.NumberOfConditions())));

    if (num_elements == 0 && num_conditions > 0) { // TODO: Can we have conditions and elements at the same time?
        // Write cells header
        rFileStream << "CELL_DATA " << r_local_mesh.NumberOfConditions() << "\n";
        const bool write_ids = mOutputSettings["write_ids"].GetBool();
        rFileStream << "FIELD FieldData " << counter_condition_results + condition_flags.size() + (write_ids ? 2 : 0) + counter_gauss_point_variables_in_elements << "\n";
        for (IndexType entry = 0; entry < condition_results.size(); ++entry) {
            const std::string& r_condition_result_name = condition_results[entry].GetString();
            WriteGeometricalContainerResults(r_condition_result_name,r_local_mesh.Conditions(),rFileStream);
        }

        // Writing condition_flags
        if (condition_flags.size() > 0) {
            // mrModelPart.GetCommunicator().SynchronizeConditionFlags(); // TODO implement this if at some point ghost-conditions are used
        }
        for (IndexType entry = 0; entry < condition_flags.size(); ++entry) {
            // Write conditional flags results variable header
            const std::string& r_condition_result_name = condition_flags[entry].GetString();
            const Flags flag = KratosComponents<Flags>::Get(r_condition_result_name);
            WriteFlagContainerVariable(r_local_mesh.Conditions(), flag, r_condition_result_name, rFileStream);
        }

        // If we write properties_id
        if (write_ids) {
            WritePropertiesIdsToFile(r_local_mesh.Conditions(), rFileStream);
            WriteIdsToFile(r_local_mesh.Conditions(), "KRATOS_CONDITION_ID", rFileStream);
        }

        // Direct write GP values
        for (IndexType entry = 0; entry < gauss_point_variables_in_elements.size(); ++entry) {
            const std::string& r_condition_result_name = gauss_point_variables_in_elements[entry].GetString();
            WriteGeometricalContainerIntegrationResults(r_condition_result_name,r_local_mesh.Conditions(),rFileStream);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void VtkOutput::WriteNodalContainerResults(
    const std::string& rVariableName,
    const ModelPart::NodesContainerType& rNodes,
    const bool IsHistoricalValue,
    std::ofstream& rFileStream) const
{
    if (KratosComponents<Variable<double>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<double>>::Get(rVariableName);
        WriteNodalScalarValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    } else if (KratosComponents<Variable<bool>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<bool>>::Get(rVariableName);
        WriteNodalScalarValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    } else if (KratosComponents<Variable<int>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<int>>::Get(rVariableName);
        WriteNodalScalarValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
        WriteNodalVectorValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    } else if (KratosComponents<Variable<Vector>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<Vector>>::Get(rVariableName);
        WriteNodalVectorValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 4>>>::Get(rVariableName);
        WriteNodalVectorValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 6>>>::Get(rVariableName);
        WriteNodalVectorValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 9>>>::Get(rVariableName);
        WriteNodalVectorValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    } else {
        KRATOS_WARNING_ONCE(rVariableName) << mrModelPart.GetCommunicator().GetDataCommunicator() << "Variable \"" << rVariableName << "\" is "
            << "not suitable for VtkOutput, skipping it" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType>
void VtkOutput::WriteGeometricalContainerResults(
    const std::string& rVariableName,
    const TContainerType& rContainer,
    std::ofstream& rFileStream) const
{
    if (KratosComponents<Variable<double>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<double>>::Get(rVariableName);
        WriteScalarContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<bool>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<bool>>::Get(rVariableName);
        WriteScalarContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<int>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<int>>::Get(rVariableName);
        WriteScalarContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<Flags>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<Flags>>::Get(rVariableName);
        WriteScalarContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
        WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<Vector>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<Vector>>::Get(rVariableName);
        WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 4>>>::Get(rVariableName);
        WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 6>>>::Get(rVariableName);
        WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 9>>>::Get(rVariableName);
        WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
    } else {
        KRATOS_WARNING_ONCE(rVariableName) << mrModelPart.GetCommunicator().GetDataCommunicator() << "Variable \"" << rVariableName << "\" is "
            << "not suitable for VtkOutput, skipping it" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType>
void VtkOutput::WriteGeometricalContainerIntegrationResults(
    const std::string& rVariableName,
    const TContainerType& rContainer,
    std::ofstream& rFileStream) const
{
    if (KratosComponents<Variable<double>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<double>>::Get(rVariableName);
        WriteIntegrationScalarContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<bool>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<bool>>::Get(rVariableName);
        WriteIntegrationScalarContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<int>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<int>>::Get(rVariableName);
        WriteIntegrationScalarContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
        WriteIntegrationVectorContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<Vector>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<Vector>>::Get(rVariableName);
        WriteIntegrationVectorContainerVariable(rContainer, var_to_write, rFileStream);
    } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 6>>>::Get(rVariableName);
        WriteIntegrationVectorContainerVariable(rContainer, var_to_write, rFileStream);
    } else {
        KRATOS_WARNING_ONCE(rVariableName) << mrModelPart.GetCommunicator().GetDataCommunicator() << "Variable \"" << rVariableName << "\" is "
            << "not suitable for VtkOutput, skipping it" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TVarType>
void VtkOutput::WriteNodalScalarValues(
    const ModelPart::NodesContainerType& rNodes,
    const TVarType& rVariable,
    const bool IsHistoricalValue,
    std::ofstream& rFileStream) const
{
    if (IsHistoricalValue) {
        mrModelPart.GetCommunicator().SynchronizeVariable(rVariable);
        WriteScalarSolutionStepVariable(rNodes, rVariable, rFileStream);
    } else {
        mrModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rVariable);
        WriteScalarContainerVariable(rNodes, rVariable, rFileStream);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TVarType>
void VtkOutput::WriteNodalVectorValues(
    const ModelPart::NodesContainerType& rNodes,
    const TVarType& rVariable,
    const bool IsHistoricalValue,
    std::ofstream& rFileStream) const
{
    if (IsHistoricalValue) {
        mrModelPart.GetCommunicator().SynchronizeVariable(rVariable);
        WriteVectorSolutionStepVariable(rNodes, rVariable, rFileStream);
    } else {
        mrModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rVariable);
        WriteVectorContainerVariable(rNodes, rVariable, rFileStream);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType, class TVarType>
void VtkOutput::WriteScalarSolutionStepVariable(
    const TContainerType& rContainer,
    const TVarType& rVariable,
    std::ofstream& rFileStream) const
{
    rFileStream << rVariable.Name() << " 1 "
                << rContainer.size() << "  float\n";

    for (const auto& r_entity : rContainer) {
        const auto& r_result = r_entity.FastGetSolutionStepValue(rVariable);
        WriteScalarDataToFile((float)r_result, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType, class TVarType>
void VtkOutput::WriteVectorSolutionStepVariable(
    const TContainerType& rContainer,
    const TVarType& rVariable,
    std::ofstream& rFileStream) const
{
    if (rContainer.size() == 0) {
        KRATOS_WARNING("VtkOutput") << mrModelPart.GetCommunicator().GetDataCommunicator() << "Empty container!" << std::endl;
        return;
    }

    const int res_size = static_cast<int>((rContainer.begin()->FastGetSolutionStepValue(rVariable)).size());

    rFileStream << rVariable.Name() << " " << res_size
                << " " << rContainer.size() << "  float\n";

    for (const auto& r_entity : rContainer) {
        const auto& r_result = r_entity.FastGetSolutionStepValue(rVariable);
        WriteVectorDataToFile(r_result, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType>
void VtkOutput::WriteFlagContainerVariable(
    const TContainerType& rContainer,
    const Flags Flag,
    const std::string& rFlagName,
    std::ofstream& rFileStream) const
{
    rFileStream << rFlagName << " 1 "
                << rContainer.size() << "  float\n";

    for (const auto& r_entity : rContainer) {
        const float result = r_entity.IsDefined(Flag) ? float(r_entity.Is(Flag)) : -1.0;
        WriteScalarDataToFile(result, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType, class TVarType>
void VtkOutput::WriteScalarContainerVariable(
    const TContainerType& rContainer,
    const TVarType& rVariable,
    std::ofstream& rFileStream) const
{
    rFileStream << rVariable.Name() << " 1 "
                << rContainer.size() << "  float\n";

    for (const auto& r_entity : rContainer) {
        const double result = r_entity.GetValue(rVariable);
        WriteScalarDataToFile((float)result, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType, class TVarType>
void VtkOutput::WriteIntegrationScalarContainerVariable(
    const TContainerType& rContainer,
    const Variable<TVarType>& rVariable,
    std::ofstream& rFileStream) const
{
    rFileStream << rVariable.Name() << " 1 "
                << rContainer.size() << "  float\n";

    // Auxiliar values
    const auto& r_process_info = mrModelPart.GetProcessInfo();
    auto& r_this_geometry_begin = (rContainer.begin())->GetGeometry();
    const GeometryData::IntegrationMethod this_integration_method = (rContainer.begin())->GetIntegrationMethod();
    const auto& r_integration_points = r_this_geometry_begin.IntegrationPoints(this_integration_method);
    const SizeType integration_points_number = r_integration_points.size();

    double aux_value;
    for (auto& r_entity : rContainer) { // TODO: CalculateOnIntegrationPoints should be const methods
        aux_value = 0.0;
        std::vector<TVarType> aux_result(integration_points_number);
        r_entity.CalculateOnIntegrationPoints(rVariable, aux_result, r_process_info);
        for (const double value : aux_result) {
            aux_value += value;
        }
        aux_value /= static_cast<double>(integration_points_number);
        WriteScalarDataToFile((float)aux_value, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType, class TVarType>
void VtkOutput::WriteVectorContainerVariable(
    const TContainerType& rContainer,
    const TVarType& rVariable,
    std::ofstream& rFileStream) const
{
    if (rContainer.size() == 0) {
        KRATOS_WARNING("VtkOutput") << mrModelPart.GetCommunicator().GetDataCommunicator() << "Empty container!" << std::endl;
        return;
    }

    const int res_size = static_cast<int>((rContainer.begin()->GetValue(rVariable)).size());

    rFileStream << rVariable.Name() << " " << res_size << " " << rContainer.size() << "  float\n";

    for (const auto& r_entity : rContainer) {
        const auto& r_result = r_entity.GetValue(rVariable);
        WriteVectorDataToFile(r_result, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType, class TVarType>
void VtkOutput::WriteIntegrationVectorContainerVariable(
    const TContainerType& rContainer,
    const Variable<TVarType>& rVariable,
    std::ofstream& rFileStream) const
{
    if (rContainer.size() == 0) {
        KRATOS_WARNING("VtkOutput") << mrModelPart.GetCommunicator().GetDataCommunicator() << "Empty container!" << std::endl;
        return;
    }

    const int res_size = static_cast<int>((rContainer.begin()->GetValue(rVariable)).size());

    rFileStream << rVariable.Name() << " " << res_size << " " << rContainer.size() << "  float\n";

    // Auxiliar values
    const auto& r_process_info = mrModelPart.GetProcessInfo();
    auto& r_this_geometry_begin = (rContainer.begin())->GetGeometry();
    const GeometryData::IntegrationMethod this_integration_method = (rContainer.begin())->GetIntegrationMethod();
    const auto& r_integration_points = r_this_geometry_begin.IntegrationPoints(this_integration_method);
    const SizeType integration_points_number = r_integration_points.size();

    TVarType aux_value;
    for (auto& r_entity : rContainer) { // TODO: CalculateOnIntegrationPoints should be const methods
        aux_value = rVariable.Zero();
        std::vector<TVarType> aux_result(integration_points_number);
        r_entity.CalculateOnIntegrationPoints(rVariable, aux_result, r_process_info);
        for (const TVarType& r_value : aux_result) {
            noalias(aux_value) += r_value;
        }
        aux_value /= static_cast<double>(integration_points_number);
        WriteVectorDataToFile(aux_value, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

/***********************************************************************************/
/***********************************************************************************/

void VtkOutput::ForceBigEndian(unsigned char* pBytes) const
{
    if (mShouldSwap) {
        unsigned char tmp = pBytes[0];
        pBytes[0] = pBytes[3];
        pBytes[3] = tmp;
        tmp = pBytes[1];
        pBytes[1] = pBytes[2];
        pBytes[2] = tmp;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType>
void VtkOutput::WritePropertiesIdsToFile(
    const TContainerType& rContainer,
    std::ofstream& rFileStream) const
{
    rFileStream << "PROPERTIES_ID" << " 1 "
                << rContainer.size() << "  int\n";

    for (const auto& r_entity : rContainer) {
        WriteScalarDataToFile((int)r_entity.GetProperties().Id(), rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<typename TContainerType>
void VtkOutput::WriteIdsToFile(
    const TContainerType& rContainer,
    const std::string DataName,
    std::ofstream& rFileStream) const
{
    rFileStream << DataName << " 1 "
                << rContainer.size() << "  int\n";

    for (const auto& r_entity : rContainer) {
        WriteScalarDataToFile((int)r_entity.Id(), rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}


/***********************************************************************************/
/***********************************************************************************/

void VtkOutput::WriteModelPartWithoutNodesToFile(ModelPart& rModelPart)
{
    // Getting model and creating auxiliar model part
    auto& r_model = mrModelPart.GetModel();
    const std::string& r_name_model_part = rModelPart.Name();
    auto& r_auxiliar_model_part = r_model.CreateModelPart("AUXILIAR_" + r_name_model_part);

    // Tranfering entities of the submodelpart
    FastTransferBetweenModelPartsProcess(r_auxiliar_model_part, rModelPart).Execute();

    // Tranfering nodes from root model part
    FastTransferBetweenModelPartsProcess(r_auxiliar_model_part, mrModelPart, FastTransferBetweenModelPartsProcess::EntityTransfered::NODES).Execute();

    // Marking to remove the nodes
    for (auto& r_node : r_auxiliar_model_part.Nodes()) {
        r_node.Set(TO_ERASE, true);
    }

    // Checking nodes from conditions
    for (auto& r_cond : r_auxiliar_model_part.Conditions()) {
        auto& r_geometry = r_cond.GetGeometry();
        for (auto& r_node : r_geometry) {
            r_node.Set(TO_ERASE, false);
        }
    }

    // Checking nodes from elements
    for (auto& r_elem : r_auxiliar_model_part.Elements()) {
        auto& r_geometry = r_elem.GetGeometry();
        for (auto& r_node : r_geometry) {
            r_node.Set(TO_ERASE, false);
        }
    }

    // Removing unused nodes
    r_auxiliar_model_part.RemoveNodes(TO_ERASE);

    // Actually writing the
    WriteModelPartToFile(r_auxiliar_model_part, true);

    // Deletin auxiliar modek part
    r_model.DeleteModelPart("AUXILIAR_" + r_name_model_part);
}

/***********************************************************************************/
/***********************************************************************************/

Parameters VtkOutput::GetDefaultParameters()
{
    // IMPORTANT: when "output_control_type" is "time", then paraview will not be able to group them
    Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"                             : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "file_format"                                 : "ascii",
        "output_precision"                            : 7,
        "output_control_type"                         : "step",
        "output_frequency"                            : 1.0,
        "output_sub_model_parts"                      : false,
        "folder_name"                                 : "VTK_Output",
        "custom_name_prefix"                          : "",
        "custom_name_postfix"                         : "",
        "save_output_files_in_folder"                 : true,
        "write_deformed_configuration"                : false,
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
    })" );

    return default_parameters;
}

} // namespace Kratos
