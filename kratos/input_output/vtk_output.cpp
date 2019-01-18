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
//
//

// System includes
#include <map>

// External includes

// Project includes
#include "vtk_output.h"


namespace Kratos
{

VtkOutput::VtkOutput(ModelPart& rModelPart, Parameters Parameters)
    : mrModelPart(rModelPart), mOutputSettings(Parameters)
{
    mDefaultPrecision = mOutputSettings["output_precision"].GetInt();
    const std::string file_format = mOutputSettings["file_format"].GetString();
    if (file_format == "ascii") {
        mFileFormat = VtkOutput::FileFormat::VTK_ASCII;
    }
    else if (file_format == "binary") {
        mFileFormat = VtkOutput::FileFormat::VTK_BINARY;
        // test for endian-format
        int num = 1;
        if (*(char *)&num == 1) {
            mShouldSwap = true;
        }
    }
    else {
        KRATOS_ERROR << "Option for \"file_format\": " << file_format
            << " not recognised!\n Possible output formats "
            << "options are: \"ascii\", \"binary\"" << std::endl;
    }
}

void VtkOutput::PrintOutput()
{
    //For whole model part
    WriteModelPartToFile(mrModelPart, false);

    //For sub model parts
    const bool print_sub_model_parts = mOutputSettings["output_sub_model_parts"].GetBool();
    if(print_sub_model_parts) {
        for (const auto& r_sub_model_part : mrModelPart.SubModelParts()) {
            WriteModelPartToFile(r_sub_model_part, true);
        }
    }
}

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

std::string VtkOutput::GetOutputFileName(const ModelPart& rModelPart, const bool IsSubModelPart)
{
    const int rank = rModelPart.GetCommunicator().MyPID();
    std::string model_part_name;

    if (IsSubModelPart) {
        model_part_name = rModelPart.GetParentModelPart()->Name() + "_" + rModelPart.Name();
    }
    else {
        model_part_name = rModelPart.Name();
    }

    std::string label;
    std::stringstream ss;
    const std::string output_control = mOutputSettings["output_control_type"].GetString();
    if (output_control == "step") {
        ss << std::setprecision(mDefaultPrecision)<< std::setfill('0')
           << rModelPart.GetProcessInfo()[STEP];
        label = ss.str();
    } else if(output_control == "time") {
        ss << std::setprecision(mDefaultPrecision) << std::setfill('0')
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
    output_file_name += model_part_name + "_" + std::to_string(rank) + "_" + label + ".vtk";

    return output_file_name;
}

void VtkOutput::Initialize(const ModelPart& rModelPart)
{
    CreateMapFromKratosIdToVTKId(rModelPart);
}

void VtkOutput::CreateMapFromKratosIdToVTKId(const ModelPart& rModelPart)
{
    int vtk_id = 0;
    for(const auto& r_node : rModelPart.Nodes()) {
        mKratosIdToVtkId[r_node.Id()] = vtk_id++;
    }
}

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

void VtkOutput::WriteMeshToFile(const ModelPart& rModelPart, std::ofstream& rFileStream) const
{
    WriteNodesToFile(rModelPart, rFileStream);
    WriteConditionsAndElementsToFile(rModelPart, rFileStream);
}

void VtkOutput::WriteNodesToFile(const ModelPart& rModelPart, std::ofstream& rFileStream) const
{
    // NOTE: also in MPI all nodes (local and ghost) have to be written, because
    // they might be needed by the elements/conditions due to the connectivity

    // write nodes header
    rFileStream << "POINTS " << rModelPart.NumberOfNodes() << " float\n";

    // write nodes
    for (const auto& r_node : rModelPart.Nodes()) {
        WriteVectorDataToFile(r_node.Coordinates(), rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
    }
}

void VtkOutput::WriteConditionsAndElementsToFile(const ModelPart& rModelPart, std::ofstream& rFileStream) const
{
    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();

    int num_elements = r_local_mesh.NumberOfElements();
    rModelPart.GetCommunicator().SumAll(num_elements);

    int num_conditions = r_local_mesh.NumberOfConditions();
    rModelPart.GetCommunicator().SumAll(num_conditions);

    if (num_elements > 0) {
        if (num_conditions > 0 && rModelPart.GetCommunicator().MyPID() == 0) {
            KRATOS_WARNING/*_ONCE*/("VtkOutput")<<"Modelpart \"" << rModelPart.Name() // TODO
                << "\" has both elements and conditions.\nGiving precedence to "
                << "elements and writing only elements!" <<std::endl;
        }
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

template<typename TContainerType>
unsigned int VtkOutput::DetermineVtkCellListSize(const TContainerType& rContainer) const
{
    unsigned int vtk_cell_list_size_container = 0;

    const auto container_begin = rContainer.begin();
    const int num_entities = rContainer.size();
    #pragma omp parallel for reduction(+:vtk_cell_list_size_container)
    for (int i=0; i<num_entities; ++i) {
        const auto entity_i = container_begin + i;
        vtk_cell_list_size_container += entity_i->GetGeometry().PointsNumber() + 1;
    }

    return vtk_cell_list_size_container;
}

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

template <typename TContainerType>
void VtkOutput::WriteCellType(const TContainerType& rContainer, std::ofstream& rFileStream) const
{
    // IMPORTANT: The map geo_type_vtk_cell_type_map is to be extended to support new geometries
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
        { GeometryData::KratosGeometryType::Kratos_Hexahedra3D8,     12 }
    };
    // write entity types
    for (const auto& r_entity : rContainer) {
        int cell_type = -1;
        if (geo_type_vtk_cell_type_map.count(r_entity.GetGeometry().GetGeometryType()) > 0) {
            cell_type = geo_type_vtk_cell_type_map.at(r_entity.GetGeometry().GetGeometryType());
        }
        else {
            KRATOS_ERROR << "Modelpart contains elements or conditions with "
             << "geometries for which no VTK-output is implemented!" << std::endl;
        }

        WriteScalarDataToFile( (int)cell_type, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
    }
}

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
    rFileStream << "POINT_DATA " << rModelPart.NumberOfNodes() << "\n";
    rFileStream << "FIELD FieldData " << nodal_solution_step_results.size() + nodal_variable_data_results.size()<< "\n";

    // Writing nodal_solution_step_results
    for (unsigned int entry = 0; entry < nodal_solution_step_results.size(); ++entry) {
        // write nodal results variable header
        const std::string nodal_result_name = nodal_solution_step_results[entry].GetString();
        WriteNodalContainerResults(nodal_result_name, rModelPart.Nodes(), true, rFileStream);
    }

    // Writing nodal_variable_data_results
    for (unsigned int entry = 0; entry < nodal_variable_data_results.size(); ++entry) {
        // write nodal results variable header
        const std::string nodal_result_name = nodal_variable_data_results[entry].GetString();
        WriteNodalContainerResults(nodal_result_name, rModelPart.Nodes(), false, rFileStream);
    }
}

void VtkOutput::WriteElementResultsToFile(const ModelPart& rModelPart, std::ofstream& rFileStream)
{
    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();
    Parameters element_results = mOutputSettings["element_data_value_variables"];

    int num_elements = r_local_mesh.NumberOfElements();
    rModelPart.GetCommunicator().SumAll(num_elements);

    if (num_elements > 0) {
        // write cells header
        rFileStream << "CELL_DATA " << r_local_mesh.NumberOfElements() << "\n";
        rFileStream << "FIELD FieldData " << element_results.size() << "\n";
        for (unsigned int entry = 0; entry < element_results.size(); ++entry) {
            const std::string element_result_name = element_results[entry].GetString();
            WriteGeometricalContainerResults(element_result_name,r_local_mesh.Elements(),rFileStream);
        }
    }
}

void VtkOutput::WriteConditionResultsToFile(const ModelPart& rModelPart, std::ofstream& rFileStream)
{
    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();
    Parameters condition_results = mOutputSettings["condition_data_value_variables"];

    int num_elements = r_local_mesh.NumberOfElements();
    rModelPart.GetCommunicator().SumAll(num_elements);

    int num_conditions = r_local_mesh.NumberOfConditions();
    rModelPart.GetCommunicator().SumAll(num_conditions);

    if (num_elements == 0 && num_conditions > 0) {
        // write cells header
        rFileStream << "CELL_DATA " << r_local_mesh.NumberOfConditions() << "\n";
        rFileStream << "FIELD FieldData " << condition_results.size() << "\n";
        for (unsigned int entry = 0; entry < condition_results.size(); ++entry) {
            const std::string condition_result_name = condition_results[entry].GetString();
            WriteGeometricalContainerResults(condition_result_name,r_local_mesh.Conditions(),rFileStream);
        }
    }
}

void VtkOutput::WriteNodalContainerResults(
    const std::string& rVariableName,
    const ModelPart::NodesContainerType& rNodes,
    const bool IsHistoricalValue,
    std::ofstream& rFileStream) const
{
    if (KratosComponents<Variable<double>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<double>>::Get(rVariableName);
        WriteNodalScalarValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    }
    else if (KratosComponents<Variable<int>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<int>>::Get(rVariableName);
        WriteNodalScalarValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    }
    else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
        WriteNodalVectorValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    }
    else if (KratosComponents<Variable<Vector>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<Vector>>::Get(rVariableName);
        WriteNodalVectorValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    }
    else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 4>>>::Get(rVariableName);
        WriteNodalVectorValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    }
    else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 6>>>::Get(rVariableName);
        WriteNodalVectorValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    }
    else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 9>>>::Get(rVariableName);
        WriteNodalVectorValues(rNodes, var_to_write, IsHistoricalValue, rFileStream);
    }
    else {
        KRATOS_WARNING_ONCE(rVariableName) << "Variable \"" << rVariableName << "\" is "
            << "not suitable for VtkOutput, skipping it" << std::endl;
    }
}

template<typename TContainerType>
void VtkOutput::WriteGeometricalContainerResults(
    const std::string& rVariableName,
    const TContainerType& rContainer,
    std::ofstream& rFileStream) const
{
    if (KratosComponents<Variable<double>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<double>>::Get(rVariableName);
        WriteScalarContainerVariable(rContainer, var_to_write, rFileStream);
    }
    else if (KratosComponents<Variable<int>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<int>>::Get(rVariableName);
        WriteScalarContainerVariable(rContainer, var_to_write, rFileStream);
    }
    else if (KratosComponents<Variable<Flags>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<Flags>>::Get(rVariableName);
        WriteScalarContainerVariable(rContainer, var_to_write, rFileStream);
    }
    else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
        WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
    }
    else if (KratosComponents<Variable<Vector>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<Vector>>::Get(rVariableName);
        WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
    }
    else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 4>>>::Get(rVariableName);
        WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
    }
    else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 6>>>::Get(rVariableName);
        WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
    }
    else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(rVariableName)){
        const auto& var_to_write = KratosComponents<Variable<array_1d<double, 9>>>::Get(rVariableName);
        WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
    }
    else {
        KRATOS_WARNING_ONCE(rVariableName) << "Variable \"" << rVariableName << "\" is "
            << "not suitable for VtkOutput, skipping it" << std::endl;
    }
}

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
    }
    else {
        mrModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rVariable);
        WriteScalarContainerVariable(rNodes, rVariable, rFileStream);
    }
}

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
    }
    else {
        mrModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rVariable);
        WriteVectorContainerVariable(rNodes, rVariable, rFileStream);
    }
}

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

template<typename TContainerType, class TVarType>
void VtkOutput::WriteVectorSolutionStepVariable(
    const TContainerType& rContainer,
    const TVarType& rVariable,
    std::ofstream& rFileStream) const
{
    KRATOS_DEBUG_ERROR_IF(rContainer.size() == 0) << "Empty container!" << std::endl;

    const int res_size = (rContainer.begin()->FastGetSolutionStepValue(rVariable)).size();

    rFileStream << rVariable.Name() << " " << res_size
                << " " << rContainer.size() << "  float\n";

    for (const auto& r_entity : rContainer) {
        const auto& r_result = r_entity.FastGetSolutionStepValue(rVariable);
        WriteVectorDataToFile(r_result, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

template<typename TContainerType, class TVarType>
void VtkOutput::WriteScalarContainerVariable(
    const TContainerType& rContainer,
    const TVarType& rVariable,
    std::ofstream& rFileStream) const
{
    rFileStream << rVariable.Name() << " 1 "
                << rContainer.size() << "  float\n";

    for (const auto& r_entity : rContainer) {
        const auto& r_result = r_entity.GetValue(rVariable);
        WriteScalarDataToFile((float)r_result, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

template<typename TContainerType, class TVarType>
void VtkOutput::WriteVectorContainerVariable(
    const TContainerType& rContainer,
    const TVarType& rVariable,
    std::ofstream& rFileStream) const
{
    KRATOS_DEBUG_ERROR_IF(rContainer.size() == 0) << "Empty container!" << std::endl;

    const int res_size = (rContainer.begin()->GetValue(rVariable)).size();

    rFileStream << rVariable.Name() << " " << res_size
                << " " << rContainer.size() << "  float\n";

    for (const auto& r_entity : rContainer) {
        const auto& r_result = r_entity.GetValue(rVariable);
        WriteVectorDataToFile(r_result, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

template <class TData>
void VtkOutput::WriteScalarDataToFile(const TData& rData, std::ofstream& rFileStream) const
{
    if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) {
        rFileStream << rData;
    }
    else if (mFileFormat == VtkOutput::FileFormat::VTK_BINARY) {
        TData data = rData;
        ForceBigEndian(reinterpret_cast<unsigned char *>(&data));
        rFileStream.write(reinterpret_cast<char *>(&data), sizeof(TData));
    }
}

template <class TData>
void VtkOutput::WriteVectorDataToFile(const TData& rData, std::ofstream& rFileStream) const
{
    if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) {
        for (const auto& r_data_comp : rData) {
            rFileStream << r_data_comp << " ";
        }
    }
    else if (mFileFormat == VtkOutput::FileFormat::VTK_BINARY)
    {
        for (const auto& r_data_comp : rData ) {
            float data_comp_local = (float)r_data_comp; // should not be const or a reference for enforcing big endian
            ForceBigEndian(reinterpret_cast<unsigned char *>(&data_comp_local));
            rFileStream.write(reinterpret_cast<char *>(&data_comp_local), sizeof(float));
        }
    }
}

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

} // namespace Kratos
