//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//
//

// project includes
#include "vtk_output.h"
#include <typeinfo>
#include <sstream>
// #include <byteswap.h>

namespace Kratos
{

VtkOutput::VtkOutput(ModelPart &rModelPart, Parameters rParameters) : mrModelPart(rModelPart), mOutputSettings(rParameters)
{
    mDefaultPrecision = this->mOutputSettings["output_precision"].GetInt();;
    this->mDoneTest = false;
    this->mShouldSwap = false;
    std::string file_format = this->mOutputSettings["file_format"].GetString();
    if(file_format == "ascii")
        mFileFormat =  VtkOutput::FileFormat::VTK_ASCII;
    else if (file_format == "binary")
        mFileFormat = VtkOutput::FileFormat::VTK_BINARY;
    else
        KRATOS_ERROR<<"Option for file_format : "<<file_format<<" not recognised.!"<<std::endl
            <<"Possible output formats options for VTKOutput are :: ascii and binary "<<std::endl;
}

VtkOutput::~VtkOutput(){};


void VtkOutput::CreateMapFromKratosIdToVTKId(const ModelPart &rModelPart)
{
    int vtk_id = 0;
    for(const auto& node : rModelPart.Nodes())
    {
        mKratosIdToVtkId[node.Id()] = vtk_id++;
    }
}

template<typename TContainerType>
unsigned int VtkOutput::DetermineVtkCellListSize(const TContainerType &rContainer)
{
    unsigned int vtk_cell_list_size_container = 0;

    const auto container_begin = rContainer.begin();
    const int num_entities = rContainer.size();
    #pragma omp parallel for reduction(+:vtk_cell_list_size_container)
    for (int i=0; i<num_entities; ++i)
    {
        auto entity_i = container_begin + i;
        vtk_cell_list_size_container += entity_i->GetGeometry().PointsNumber() +1;
    }

    return vtk_cell_list_size_container;
}

void VtkOutput::Initialize(const ModelPart &rModelPart)
{
    CreateMapFromKratosIdToVTKId(rModelPart);
}

void VtkOutput::WriteHeader(const ModelPart &rModelPart, std::ofstream& rFileStream)
{
    rFileStream << "# vtk DataFile Version 4.0"
               << "\n";
    rFileStream << "vtk output"
               << "\n";
    if(mFileFormat == VtkOutput::FileFormat::VTK_ASCII)
        rFileStream << "ASCII"<< "\n";
    else if (mFileFormat == VtkOutput::FileFormat::VTK_BINARY)
        rFileStream << "BINARY"<< "\n";
    rFileStream << "DATASET UNSTRUCTURED_GRID"
               << "\n";
}

void VtkOutput::WriteMesh(const ModelPart &rModelPart, std::ofstream& rFileStream)
{
    WriteNodes(rModelPart, rFileStream);
    WriteConditionsAndElements(rModelPart, rFileStream);
}

void VtkOutput::WriteNodes(const ModelPart &rModelPart, std::ofstream& rFileStream)
{
    // write nodes header
    const auto& local_mesh = rModelPart.GetCommunicator().LocalMesh();
    rFileStream << "POINTS " << local_mesh.NumberOfNodes() << " float"
               << "\n";

    // write nodes
    for(const auto& node : local_mesh.Nodes())
    {
        const auto& coordinates = node.Coordinates();
        WriteVectorDataToFile(coordinates, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
    }
}

void VtkOutput::WriteConditionsAndElements(const ModelPart &rModelPart, std::ofstream& rFileStream)
{
    const auto& local_mesh = rModelPart.GetCommunicator().LocalMesh();
    if(rModelPart.NumberOfElements() > 0)
    {
        if(local_mesh.NumberOfConditions()>0)
            KRATOS_WARNING/*_ONCE*/("VtkOutput")<<"Modelpart \"" << rModelPart.Name() // TODO
                << "\" has both elements and conditions.\nGiving precedence to "
                << "elements and writing only elements!" <<std::endl;

        // write cells header
        rFileStream << "\nCELLS " <<local_mesh.NumberOfElements() << " " << DetermineVtkCellListSize(local_mesh.Elements()) << "\n";
        WriteConnectivity(local_mesh.Elements(), rFileStream);
        // write cell types header
        rFileStream << "\nCELL_TYPES " <<local_mesh.NumberOfElements() << "\n";
        WriteCellType(local_mesh.Elements(), rFileStream);
    } else if(local_mesh.NumberOfConditions()>0)
    {
        // write cells header
        rFileStream << "\nCELLS " << local_mesh.NumberOfConditions() << " " << DetermineVtkCellListSize(local_mesh.Conditions()) << "\n";
        WriteConnectivity(local_mesh.Conditions(), rFileStream);
        // write cell types header
        rFileStream << "\nCELL_TYPES " << local_mesh.NumberOfConditions() << "\n";
        WriteCellType(local_mesh.Conditions(), rFileStream);
    }
}

template <typename TContainerType>
void VtkOutput::WriteConnectivity(const TContainerType& rContainer, std::ofstream& rFileStream)
{
    // write Conditions
    for (const auto& entity : rContainer)
    {
        const auto& entity_geometry = entity.GetGeometry();
        const unsigned int number_of_nodes = entity_geometry.size();

        WriteScalarDataToFile((unsigned int)number_of_nodes, rFileStream);
        for (const auto& node : entity_geometry){
            if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << " ";
            int id = mKratosIdToVtkId[node.Id()];
            WriteScalarDataToFile((int)id, rFileStream);
        }
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
    }
}

// IMPORTANT :: The map geo_type_vtk_cell_type_map is to be extended to support new geometries
template <typename TContainerType>
void VtkOutput::WriteCellType(const TContainerType& rContainer, std::ofstream& rFileStream)
{
    std::unordered_map<GeometryData::KratosGeometryType, int> geo_type_vtk_cell_type_map = {
            {GeometryData::KratosGeometryType::Kratos_Triangle2D3, 5 },
            {GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4, 9 },
            {GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4, 10},
            {GeometryData::KratosGeometryType::Kratos_Line2D2, 3},
            {GeometryData::KratosGeometryType::Kratos_Line3D2, 3},
            {GeometryData::KratosGeometryType::Kratos_Point2D, 1},
            {GeometryData::KratosGeometryType::Kratos_Point3D, 1}
    };
    // write entity types
    for (const auto& entity : rContainer)
    {
        int cell_type = -1;
        if( geo_type_vtk_cell_type_map.count( entity.GetGeometry().GetGeometryType()) > 0  )
            cell_type = geo_type_vtk_cell_type_map[  entity.GetGeometry().GetGeometryType() ];
        else
            KRATOS_ERROR<<"Modelpart contains elements or conditions with geometries for which no VTK-output is implemented!"<<std::endl;

        WriteScalarDataToFile( (int)cell_type, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
    }
}


// IMPORTANT :: ONLY The following two functions are to be externded to add a new type of result.
VtkOutput::WriteDataType VtkOutput::GetDataCharacterstic(std::string VariableName)
{
    if (KratosComponents<Variable<double>>::Has(VariableName))
        return  VtkOutput::WriteDataType::VTK_SCALAR;
    else if (KratosComponents<Flags>::Has(VariableName))
        return  VtkOutput::WriteDataType::VTK_SCALAR;
    else if (KratosComponents<Variable<int>>::Has(VariableName))
        return  VtkOutput::WriteDataType::VTK_SCALAR;
    else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(VariableName))
        return  VtkOutput::WriteDataType::VTK_VECTOR_3;
    else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(VariableName))
        return  VtkOutput::WriteDataType::VTK_VECTOR_4;
    else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(VariableName))
        return  VtkOutput::WriteDataType::VTK_VECTOR_6;
    else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(VariableName))
        return  VtkOutput::WriteDataType::VTK_VECTOR_9;
    else
        KRATOS_ERROR<<"Variable : "<<VariableName<<" does not exist."<<std::endl;
}

template<typename TContainerType>
void VtkOutput::WriteContainerSolutionsStepResult(std::string NodalResultName, const TContainerType &rTContainer,  std::ofstream& rFileStream)
{
    VtkOutput::WriteDataType data_characteristic = GetDataCharacterstic(NodalResultName);
    rFileStream << NodalResultName <<" "<< (int)data_characteristic<<" "<< rTContainer.size() <<" "<<" float"<<"\n";

    if(data_characteristic == VtkOutput::WriteDataType::VTK_SCALAR)
    {
        // For double variable
        if (KratosComponents<Variable<double>>::Has(NodalResultName)){
            const auto& var_to_write = KratosComponents<Variable<double>>::Get(NodalResultName);
            for(const auto& entity : rTContainer)
            {
                const auto &nodal_result =  entity.FastGetSolutionStepValue(var_to_write);
                WriteScalarDataToFile((float)nodal_result, rFileStream);
                if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
            }
        }
        // For Int variable
        if (KratosComponents<Variable<int>>::Has(NodalResultName)){
            const auto& var_to_write = KratosComponents<Variable<int>>::Get(NodalResultName);
            for(const auto& entity : rTContainer)
            {
                const auto &nodal_result =  entity.FastGetSolutionStepValue(var_to_write);
                WriteScalarDataToFile((float)nodal_result, rFileStream);
                if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
            }
        }
        // For Flag variable
        if (KratosComponents<Variable<Flags>>::Has(NodalResultName)){
            const auto& var_to_write = KratosComponents<Variable<Flags>>::Get(NodalResultName);
            for(const auto& entity : rTContainer)
            {
                const auto &nodal_result =  entity.FastGetSolutionStepValue(var_to_write);
                WriteScalarDataToFile((float)nodal_result, rFileStream);
                if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
            }
        }
    }
    else
    {
        // For Vector variable 3
        if (KratosComponents<Variable<array_1d<double, 3>>>::Has(NodalResultName)){
           const auto& var_to_write = KratosComponents<Variable<array_1d<double, 3>>>::Get(NodalResultName);
            for(const auto& entity : rTContainer)
            {
                const auto &nodal_result =  entity.FastGetSolutionStepValue(var_to_write);
                WriteVectorDataToFile(nodal_result, rFileStream);
                if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
            }
        }
        // For Vector variable 4
        if (KratosComponents<Variable<array_1d<double, 4>>>::Has(NodalResultName)){
            const auto& var_to_write = KratosComponents<Variable<array_1d<double, 4>>>::Get(NodalResultName);
            for(const auto& entity : rTContainer)
            {
                const auto &nodal_result =  entity.FastGetSolutionStepValue(var_to_write);
                WriteVectorDataToFile(nodal_result, rFileStream);
                if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
            }
        }
        // For Vector variable 6
        if (KratosComponents<Variable<array_1d<double, 6>>>::Has(NodalResultName)){
            const auto& var_to_write = KratosComponents<Variable<array_1d<double, 6>>>::Get(NodalResultName);
            for(const auto& entity : rTContainer)
            {
                const auto &nodal_result =  entity.FastGetSolutionStepValue(var_to_write);
                WriteVectorDataToFile(nodal_result, rFileStream);
                if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
            }
        }
        // For Vector variable 9
        if (KratosComponents<Variable<array_1d<double, 9>>>::Has(NodalResultName)){
            const auto& var_to_write = KratosComponents<Variable<array_1d<double, 9>>>::Get(NodalResultName);
            for(const auto& entity : rTContainer)
            {
                const auto &nodal_result =  entity.FastGetSolutionStepValue(var_to_write);
                WriteVectorDataToFile(nodal_result, rFileStream);
                if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
            }
        }
    }
}

template<typename TContainerType>
void VtkOutput::WriteContainerVariableResults(std::string NodalResultName, const TContainerType &rTContainer,  std::ofstream& rFileStream)
{
    VtkOutput::WriteDataType data_characteristic = GetDataCharacterstic(NodalResultName);
    rFileStream << NodalResultName <<" "<< (int)data_characteristic<<" "<< rTContainer.size() <<" "<<" float"<<"\n";

    if(data_characteristic == VtkOutput::WriteDataType::VTK_SCALAR)
    {
        // For double variable
        if (KratosComponents<Variable<double>>::Has(NodalResultName)){
            const auto& var_to_write = KratosComponents<Variable<double>>::Get(NodalResultName);
            for(const auto& entity : rTContainer)
            {
                const auto &nodal_result =  entity.GetValue(var_to_write);
                WriteScalarDataToFile((float)nodal_result, rFileStream);
                if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
            }
        }
        // For Int variable
        if (KratosComponents<Variable<int>>::Has(NodalResultName)){
            const auto& var_to_write = KratosComponents<Variable<int>>::Get(NodalResultName);
            for(const auto& entity : rTContainer)
            {
                const auto &nodal_result =  entity.GetValue(var_to_write);
                WriteScalarDataToFile((float)nodal_result, rFileStream);
                if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
            }
        }
        // For Flag variable
        if (KratosComponents<Variable<Flags>>::Has(NodalResultName)){
            const auto& var_to_write = KratosComponents<Variable<Flags>>::Get(NodalResultName);
            for(const auto& entity : rTContainer)
            {
                const auto &nodal_result =  entity.GetValue(var_to_write);
                WriteScalarDataToFile((float)nodal_result, rFileStream);
                if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
            }
        }
    }
    else
    {
        // For Vector variable 3
        if (KratosComponents<Variable<array_1d<double, 3>>>::Has(NodalResultName)){
           const auto& var_to_write = KratosComponents<Variable<array_1d<double, 3>>>::Get(NodalResultName);
            for(const auto& entity : rTContainer)
            {
                const auto &nodal_result =  entity.GetValue(var_to_write);
                WriteVectorDataToFile(nodal_result, rFileStream);
                if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
            }
        }
        // For Vector variable 4
        if (KratosComponents<Variable<array_1d<double, 4>>>::Has(NodalResultName)){
            const auto& var_to_write = KratosComponents<Variable<array_1d<double, 4>>>::Get(NodalResultName);
            for(const auto& entity : rTContainer)
            {
                const auto &nodal_result =  entity.GetValue(var_to_write);
                WriteVectorDataToFile(nodal_result, rFileStream);
                if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
            }
        }
        // For Vector variable 6
        if (KratosComponents<Variable<array_1d<double, 6>>>::Has(NodalResultName)){
            const auto& var_to_write = KratosComponents<Variable<array_1d<double, 6>>>::Get(NodalResultName);
            for(const auto& entity : rTContainer)
            {
                const auto &nodal_result =  entity.GetValue(var_to_write);
                WriteVectorDataToFile(nodal_result, rFileStream);
                if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
            }
        }
        // For Vector variable 9
        if (KratosComponents<Variable<array_1d<double, 9>>>::Has(NodalResultName)){
            const auto& var_to_write = KratosComponents<Variable<array_1d<double, 9>>>::Get(NodalResultName);
            for(const auto& entity : rTContainer)
            {
                const auto &nodal_result =  entity.GetValue(var_to_write);
                WriteVectorDataToFile(nodal_result, rFileStream);
                if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
            }
        }
    }
}

void VtkOutput::WriteNodalResults(const ModelPart &rModelPart, std::ofstream& rFileStream)
{
    const auto& local_mesh = rModelPart.GetCommunicator().LocalMesh();
    // write nodal results header
    Parameters nodal_solution_step_results = this->mOutputSettings["nodal_solution_step_data_variables"];
    Parameters nodal_variable_data_results = this->mOutputSettings["nodal_data_value_variables"];
    rFileStream << "POINT_DATA " << local_mesh.NumberOfNodes() << "\n";
    rFileStream << "FIELD FieldData " << nodal_solution_step_results.size() + nodal_variable_data_results.size()<< "\n";

    // Writing nodal_solution_step_results
    for (unsigned int entry = 0; entry < nodal_solution_step_results.size(); entry++)
    {
        // write nodal results variable header
        std::string nodal_result_name = nodal_solution_step_results[entry].GetString();
        WriteContainerSolutionsStepResult(nodal_result_name,local_mesh.Nodes(),rFileStream);
    }

    // Writing nodal_variable_data_results
    for (unsigned int entry = 0; entry < nodal_variable_data_results.size(); entry++)
    {
        // write nodal results variable header
        std::string nodal_result_name = nodal_variable_data_results[entry].GetString();
        WriteContainerVariableResults(nodal_result_name,local_mesh.Nodes(),rFileStream);
    }
}

void VtkOutput::WriteElementResults(const ModelPart &rModelPart, std::ofstream& rFileStream)
{
    const auto& local_mesh = rModelPart.GetCommunicator().LocalMesh();
    Parameters element_results = this->mOutputSettings["element_data_value_variables"];// list of element results

    // write cells header
    if (local_mesh.NumberOfElements() > 0)
    {
        rFileStream << "CELL_DATA " << local_mesh.Elements().size() << "\n";
        rFileStream << "FIELD FieldData " << element_results.size() << "\n";
        for (unsigned int entry = 0; entry < element_results.size(); entry++)
        {
            std::string element_result_name = element_results[entry].GetString();
            WriteContainerVariableResults(element_result_name,local_mesh.Elements(),rFileStream);
        }
    }
}

void VtkOutput::WriteConditionResults(const ModelPart &rModelPart, std::ofstream& rFileStream)
{
    const auto& local_mesh = rModelPart.GetCommunicator().LocalMesh();
    Parameters condition_results = this->mOutputSettings["condition_data_value_variables"];// list of element results
    if(local_mesh.NumberOfElements() == 0 && rModelPart.NumberOfConditions()>0)
    {
        rFileStream << "CELL_DATA " << local_mesh.Conditions().size() << "\n";
        rFileStream << "FIELD FieldData " << condition_results.size() << "\n";
        for (unsigned int entry = 0; entry < condition_results.size(); entry++)
        {
            std::string condition_result_name = condition_results[entry].GetString();
            WriteContainerVariableResults(condition_result_name,local_mesh.Conditions(),rFileStream);
        }
    }
}

void VtkOutput::WriteModelPart(const ModelPart &rModelPart)
{
    Initialize(rModelPart);

    // Make the file stream object
    std::string output_file_name = GetOutputFileName(rModelPart);
    std::ofstream output_file;
    if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII)
    {
        output_file.open(output_file_name, std::ios::out | std::ios::trunc);
        output_file << std::scientific;
        output_file << std::setprecision(mDefaultPrecision);
    } else if (mFileFormat == VtkOutput::FileFormat::VTK_BINARY)
    {
        output_file.open(output_file_name, std::ios::out | std::ios::binary | std::ios::trunc);
    }

    WriteHeader(rModelPart, output_file);
    WriteMesh(rModelPart, output_file);
    WriteNodalResults(rModelPart, output_file);
    WriteElementResults(rModelPart, output_file);
    WriteConditionResults(rModelPart, output_file);

    output_file.close();
}

void VtkOutput::PrintOutput()
{
    //For whole model part
    WriteModelPart(mrModelPart);
    //For sub model parts
    bool print_sub_model_parts = this->mOutputSettings["output_sub_model_parts"].GetBool();

    if(print_sub_model_parts)
    {
        for (const auto& sub_model_part : mrModelPart.SubModelParts())
        {
            std::string output_file_name = GetOutputFileName(sub_model_part);
            WriteModelPart(sub_model_part);
        }
    }

}

void VtkOutput::ForceBigEndian(unsigned char *bytes)
{
    if (!mDoneTest)
    {
        int num = 1;
        if (*(char *)&num == 1)
            mShouldSwap = true;
        mDoneTest = true;
    }

    if (mShouldSwap)
    {
        unsigned char tmp = bytes[0];
        bytes[0] = bytes[3];
        bytes[3] = tmp;
        tmp = bytes[1];
        bytes[1] = bytes[2];
        bytes[2] = tmp;
    }
}


std::string VtkOutput::GetOutputFileName(const ModelPart &rModelPart)
{
    int rank = rModelPart.GetCommunicator().MyPID();
    std::string model_part_name;
    if(rModelPart.IsSubModelPart())
        model_part_name = rModelPart.GetParentModelPart()->Name() + "_" + rModelPart.Name();
    else
        model_part_name = rModelPart.Name();

    std::string label;
    std::stringstream ss;
    std::string output_control = mOutputSettings["output_control_type"].GetString();
    if(output_control == "step"){
        ss << std::setprecision(mDefaultPrecision)<< std::setfill('0') << rModelPart.GetProcessInfo()[STEP];
        label = ss.str();
    }else if(output_control == "time"){
        ss << std::setprecision(mDefaultPrecision) << std::setfill('0') << rModelPart.GetProcessInfo()[TIME];
        label = ss.str();
    }else
        KRATOS_ERROR<<"Option for output_control_type : "<<output_control<<" not recognised.!"<<std::endl
            <<"Possible output_control_type options for VTKOutput are :: step and time "<<std::endl;

    // Putting every thing together
    std::string output_file_name = mOutputSettings["folder_name"].GetString() +"/"+
                                    model_part_name + "_" + std::to_string(rank) + "_" + label + ".vtk";
    return output_file_name;
}


template <class TData>
void VtkOutput::WriteScalarDataToFile(const TData& rData, std::ofstream& rFileStream)
{
    if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII)
    {
        rFileStream << rData;
    }
    else if (mFileFormat == VtkOutput::FileFormat::VTK_BINARY)
    {
        TData data = rData;
        ForceBigEndian(reinterpret_cast<unsigned char *>(&data));
        rFileStream.write(reinterpret_cast<char *>(&data), sizeof(TData));
    }
}

template <class TData>
void VtkOutput::WriteVectorDataToFile(const TData& rData, std::ofstream& rFileStream)
{
    if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII)
    {
        for(const auto& data_comp : rData )
            rFileStream << data_comp << " ";
    }
    else if (mFileFormat == VtkOutput::FileFormat::VTK_BINARY)
    {
        for(const auto& data_comp : rData ){
            float data_comp_local = (float)data_comp; // should not be const or a reference for enforcing big endian
            ForceBigEndian(reinterpret_cast<unsigned char *>(&data_comp_local));
            rFileStream.write(reinterpret_cast<char *>(&data_comp_local), sizeof(float));
        }
    }
}

} // namespace Kratos
