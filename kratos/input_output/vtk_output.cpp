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

// project includes
#include "vtk_output.h"

namespace Kratos
{

VtkOutput::VtkOutput(ModelPart& rModelPart, Parameters rParameters)
    : mrModelPart(rModelPart), mOutputSettings(rParameters)
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
    WriteModelPart(mrModelPart, false);

    //For sub model parts
    const bool print_sub_model_parts = mOutputSettings["output_sub_model_parts"].GetBool();
    if(print_sub_model_parts) {
        for (const auto& r_sub_model_part : mrModelPart.SubModelParts()) {
            WriteModelPart(r_sub_model_part, true);
        }
    }
}

void VtkOutput::WriteModelPart(const ModelPart& rModelPart, const bool IsSubModelPart)
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

    WriteHeader(rModelPart, output_file);
    WriteMesh(rModelPart, output_file);
    WriteNodalResults(rModelPart, output_file);
    WriteElementResults(rModelPart, output_file);
    WriteConditionResults(rModelPart, output_file);

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

void VtkOutput::WriteHeader(const ModelPart& rModelPart, std::ofstream& rFileStream) const
{
    rFileStream << "# vtk DataFile Version 4.0\nvtk output\n";

    if(mFileFormat == VtkOutput::FileFormat::VTK_ASCII) {
        rFileStream << "ASCII";
    }
    else if (mFileFormat == VtkOutput::FileFormat::VTK_BINARY) {
        rFileStream << "BINARY";
    }

    rFileStream << "\nDATASET UNSTRUCTURED_GRID\n";
}

void VtkOutput::WriteMesh(const ModelPart& rModelPart, std::ofstream& rFileStream) const
{
    WriteNodes(rModelPart, rFileStream);
    WriteConditionsAndElements(rModelPart, rFileStream);
}

void VtkOutput::WriteNodes(const ModelPart& rModelPart, std::ofstream& rFileStream) const
{
    // write nodes header
    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();
    rFileStream << "POINTS " << r_local_mesh.NumberOfNodes() << " float\n";

    // write nodes
    for (const auto& r_node : r_local_mesh.Nodes()) {
        WriteVectorDataToFile(r_node.Coordinates(), rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
    }
}

void VtkOutput::WriteConditionsAndElements(const ModelPart& rModelPart, std::ofstream& rFileStream) const
{
    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();

    int num_elements = r_local_mesh.NumberOfElements();
    rModelPart.GetCommunicator().SumAll(num_elements);

    int num_conditions = r_local_mesh.NumberOfConditions();
    rModelPart.GetCommunicator().SumAll(num_conditions);

    if (num_elements > 0) {
        if (num_conditions > 0) {
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
    const std::unordered_map<GeometryData::KratosGeometryType, int> geo_type_vtk_cell_type_map = {
        { GeometryData::KratosGeometryType::Kratos_Triangle2D3, 5 },
        { GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4, 9 },
        { GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4, 10 },
        { GeometryData::KratosGeometryType::Kratos_Line2D2, 3 },
        { GeometryData::KratosGeometryType::Kratos_Line3D2, 3 },
        { GeometryData::KratosGeometryType::Kratos_Point2D, 1 },
        { GeometryData::KratosGeometryType::Kratos_Point3D, 1 }
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

void VtkOutput::WriteNodalResults(const ModelPart& rModelPart, std::ofstream& rFileStream)
{
    const auto& r_local_mesh = rModelPart.GetCommunicator().LocalMesh();
    // write nodal results header
    Parameters nodal_solution_step_results = mOutputSettings["nodal_solution_step_data_variables"];
    Parameters nodal_variable_data_results = mOutputSettings["nodal_data_value_variables"];
    rFileStream << "POINT_DATA " << r_local_mesh.NumberOfNodes() << "\n";
    rFileStream << "FIELD FieldData " << nodal_solution_step_results.size() + nodal_variable_data_results.size()<< "\n";

    // Writing nodal_solution_step_results
    for (unsigned int entry = 0; entry < nodal_solution_step_results.size(); ++entry) {
        // write nodal results variable header
        const std::string nodal_result_name = nodal_solution_step_results[entry].GetString();
        WriteContainerSolutionsStepResults(nodal_result_name,r_local_mesh.Nodes(),rFileStream);
    }

    // Writing nodal_variable_data_results
    for (unsigned int entry = 0; entry < nodal_variable_data_results.size(); ++entry) {
        // write nodal results variable header
        const std::string nodal_result_name = nodal_variable_data_results[entry].GetString();
        WriteContainerVariableResults(nodal_result_name,r_local_mesh.Nodes(),rFileStream);
    }
}

void VtkOutput::WriteElementResults(const ModelPart& rModelPart, std::ofstream& rFileStream)
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
            WriteContainerVariableResults(element_result_name,r_local_mesh.Elements(),rFileStream);
        }
    }
}

void VtkOutput::WriteConditionResults(const ModelPart& rModelPart, std::ofstream& rFileStream)
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
            WriteContainerVariableResults(condition_result_name,r_local_mesh.Conditions(),rFileStream);
        }
    }
}

// IMPORTANT :: ONLY The following two functions are to be externded to add a new type of result.
VtkOutput::WriteDataType VtkOutput::GetWriteDataType(const std::string& rVariableName) const
{
    if (KratosComponents<Variable<double>>::Has(rVariableName))
        return VtkOutput::WriteDataType::VTK_SCALAR;
    else if (KratosComponents<Flags>::Has(rVariableName))
        return VtkOutput::WriteDataType::VTK_SCALAR;
    else if (KratosComponents<Variable<int>>::Has(rVariableName))
        return VtkOutput::WriteDataType::VTK_SCALAR;
    else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName))
        return VtkOutput::WriteDataType::VTK_VECTOR_3;
    else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(rVariableName))
        return VtkOutput::WriteDataType::VTK_VECTOR_4;
    else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName))
        return VtkOutput::WriteDataType::VTK_VECTOR_6;
    else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(rVariableName))
        return VtkOutput::WriteDataType::VTK_VECTOR_9;
    else
        KRATOS_ERROR << "Variable: \"" << rVariableName << "\" not supported "
                     << "for VtkOutput" << std::endl;
}

template<typename TContainerType>
void VtkOutput::WriteContainerSolutionsStepResults(
    const std::string& rVariableName,
    const TContainerType &rContainer,
    std::ofstream& rFileStream) const
{
    VtkOutput::WriteDataType vtk_data_type = GetWriteDataType(rVariableName);
    rFileStream << rVariableName << " " << (int)vtk_data_type
                << " " << rContainer.size() << "  float\n";

    if (vtk_data_type == VtkOutput::WriteDataType::VTK_SCALAR) {
        if (KratosComponents<Variable<double>>::Has(rVariableName)){
            const auto& var_to_write = KratosComponents<Variable<double>>::Get(rVariableName);
            WriteScalarSolutionStepVariable(rContainer, var_to_write, rFileStream);
        }
        if (KratosComponents<Variable<int>>::Has(rVariableName)){
            const auto& var_to_write = KratosComponents<Variable<int>>::Get(rVariableName);
            WriteScalarSolutionStepVariable(rContainer, var_to_write, rFileStream);
        }
        if (KratosComponents<Variable<Flags>>::Has(rVariableName)){
            const auto& var_to_write = KratosComponents<Variable<Flags>>::Get(rVariableName);
            WriteScalarSolutionStepVariable(rContainer, var_to_write, rFileStream);
        }
    }
    else {
        if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)){
           const auto& var_to_write = KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
            WriteVectorSolutionStepVariable(rContainer, var_to_write, rFileStream);
        }
        if (KratosComponents<Variable<array_1d<double, 4>>>::Has(rVariableName)){
            const auto& var_to_write = KratosComponents<Variable<array_1d<double, 4>>>::Get(rVariableName);
            WriteVectorSolutionStepVariable(rContainer, var_to_write, rFileStream);
        }
        if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)){
            const auto& var_to_write = KratosComponents<Variable<array_1d<double, 6>>>::Get(rVariableName);
            WriteVectorSolutionStepVariable(rContainer, var_to_write, rFileStream);
        }
        if (KratosComponents<Variable<array_1d<double, 9>>>::Has(rVariableName)){
            const auto& var_to_write = KratosComponents<Variable<array_1d<double, 9>>>::Get(rVariableName);
            WriteVectorSolutionStepVariable(rContainer, var_to_write, rFileStream);
        }
    }
}

template<typename TContainerType>
void VtkOutput::WriteContainerVariableResults(
    const std::string& rVariableName,
    const TContainerType&rContainer,
    std::ofstream& rFileStream) const
{
    VtkOutput::WriteDataType vtk_data_type = GetWriteDataType(rVariableName);
    rFileStream << rVariableName << " " << (int)vtk_data_type
                << " " << rContainer.size() << "  float\n";

    if (vtk_data_type == VtkOutput::WriteDataType::VTK_SCALAR) {
        if (KratosComponents<Variable<double>>::Has(rVariableName)){
            const auto& var_to_write = KratosComponents<Variable<double>>::Get(rVariableName);
            WriteScalarContainerVariable(rContainer, var_to_write, rFileStream);
        }
        if (KratosComponents<Variable<int>>::Has(rVariableName)){
            const auto& var_to_write = KratosComponents<Variable<int>>::Get(rVariableName);
            WriteScalarContainerVariable(rContainer, var_to_write, rFileStream);
        }
        if (KratosComponents<Variable<Flags>>::Has(rVariableName)){
            const auto& var_to_write = KratosComponents<Variable<Flags>>::Get(rVariableName);
            WriteScalarContainerVariable(rContainer, var_to_write, rFileStream);
        }
    }
    else {
        if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName)){
           const auto& var_to_write = KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
            WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
        }
        if (KratosComponents<Variable<array_1d<double, 4>>>::Has(rVariableName)){
            const auto& var_to_write = KratosComponents<Variable<array_1d<double, 4>>>::Get(rVariableName);
            WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
        }
        if (KratosComponents<Variable<array_1d<double, 6>>>::Has(rVariableName)){
            const auto& var_to_write = KratosComponents<Variable<array_1d<double, 6>>>::Get(rVariableName);
            WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
        }
        if (KratosComponents<Variable<array_1d<double, 9>>>::Has(rVariableName)){
            const auto& var_to_write = KratosComponents<Variable<array_1d<double, 9>>>::Get(rVariableName);
            WriteVectorContainerVariable(rContainer, var_to_write, rFileStream);
        }
    }
}

template<typename TContainerType, class TVarType>
void VtkOutput::WriteScalarSolutionStepVariable(
    const TContainerType& rContainer,
    const TVarType& rVariable,
    std::ofstream& rFileStream) const
{
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
