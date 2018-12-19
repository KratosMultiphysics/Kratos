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

namespace Kratos
{

VtkOutput::VtkOutput(ModelPart &rModelPart, Parameters rParameters) : mrModelPart(rModelPart), mrOutputSettings(rParameters)
{
    mDefaultPrecision = this->mrOutputSettings["output_precision"].GetInt();;
    this->mDoneTest = false;
    this->mShouldSwap = false;
    std::string file_format = this->mrOutputSettings["file_format"].GetString();
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
        int KratosId = node.Id();
        mKratosIdToVtkId[KratosId] = vtk_id;
        vtk_id++;
    }
}

unsigned int VtkOutput::DetermineVtkCellListSize(const ModelPart &rModelPart)
{
    unsigned int vtk_cell_list_size = 0;
    const auto& local_mesh = rModelPart.GetCommunicator().LocalMesh();

    for (const auto& elem : local_mesh.Elements())
    {
        vtk_cell_list_size++;
        vtk_cell_list_size += elem.GetGeometry().size();
    }

    for (const auto& condition : local_mesh.Conditions())
    {
        vtk_cell_list_size++;
        vtk_cell_list_size += condition.GetGeometry().size();
    }

    return vtk_cell_list_size;
}

void VtkOutput::Initialize(const ModelPart &rModelPart)
{
    CreateMapFromKratosIdToVTKId(rModelPart);
    mVtkCellListSize = DetermineVtkCellListSize(rModelPart);
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
    WriteConditionAndElementTypes(rModelPart, rFileStream);
}

void VtkOutput::WriteNodes(const ModelPart &rModelPart, std::ofstream& rFileStream)
{
    // write nodes header
    const auto& local_mesh = rModelPart.GetCommunicator().LocalMesh();
    rFileStream << "POINTS " << local_mesh.NumberOfNodes() << " double"
               << "\n";

    // write nodes
    for(const auto& node : local_mesh.Nodes())
    {
        const auto& coordinates = node.Coordinates();
        WriteScalarDataToFile((double)coordinates[0], rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << " ";
        WriteScalarDataToFile((double)coordinates[1], rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << " ";
        WriteScalarDataToFile((double)coordinates[2], rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << " ";
    }
}

void VtkOutput::WriteConditionsAndElements(const ModelPart &rModelPart, std::ofstream& rFileStream)
{
    const auto& local_mesh = rModelPart.GetCommunicator().LocalMesh();
    // write cells header
    rFileStream << "\nCELLS " << local_mesh.NumberOfConditions() + local_mesh.NumberOfElements() << " " << mVtkCellListSize << "\n";

    // write elements
    for (const auto& elem : local_mesh.Elements())
    {
        const auto& elem_geometry = elem.GetGeometry();
        const unsigned int number_of_nodes = elem_geometry.size();

        WriteScalarDataToFile((unsigned int)number_of_nodes, rFileStream);
        for (const auto& node : elem_geometry){
            if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << " ";
            int id = mKratosIdToVtkId[node.Id()];
            WriteScalarDataToFile((int)id, rFileStream);
        }
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
    }

    // write Conditions
    for (const auto& condition : local_mesh.Conditions())
    {
        const auto& condition_geometry = condition.GetGeometry();
        const unsigned int number_of_nodes = condition_geometry.size();

        WriteScalarDataToFile((unsigned int)number_of_nodes, rFileStream);
        for (const auto& node : condition_geometry){
            if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << " ";
            int id = mKratosIdToVtkId[node.Id()];
            WriteScalarDataToFile((int)id, rFileStream);
        }
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
    }

}

void VtkOutput::WriteConditionAndElementTypes(const ModelPart &rModelPart, std::ofstream& rFileStream)
{
    const auto& local_mesh = rModelPart.GetCommunicator().LocalMesh();
    // write cell types header
    rFileStream << "\nCELL_TYPES " << local_mesh.NumberOfConditions() + local_mesh.NumberOfElements() << "\n";

    // write elements types
    for (const auto& elem : local_mesh.Elements())
    {
        int  cell_type = -1;

        if (elem.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3)
             cell_type = 5; // triangle 2d
        if (elem.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
             cell_type = 9; // quad 2d
        if (elem.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4)
             cell_type = 10; // tetra 3d
        if (elem.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line2D2)
             cell_type = 3; // line 2d or 3d
        if (elem.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line3D2)
             cell_type = 3; // line 2d or 3d
        if (elem.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Point2D)
             cell_type = 1;
        if (elem.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Point3D)
             cell_type = 1;
        if( cell_type < 0)
            KRATOS_ERROR<<"Modelpart contains elements with geometries for which no VTK-output is implemented!"<<std::endl;

        WriteScalarDataToFile( (int)cell_type, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
    }

    // write conditions types
    for (const auto& condition : local_mesh.Conditions())
    {
        int  cell_type = -1;

        if (condition.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3)
             cell_type = 5; // triangle 2d
        if (condition.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
             cell_type = 9; // quad 2d
        if (condition.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4)
             cell_type = 10; // tetra 3d
        if (condition.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line2D2)
             cell_type = 3; // line 2d or 3d
        if (condition.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line3D2)
             cell_type = 3; // line 2d or 3d
        if (condition.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Point2D)
             cell_type = 1;
        if (condition.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Point3D)
             cell_type = 1;
        if( cell_type < 0)
            KRATOS_ERROR<<"Modelpart contains elements with geometries for which no VTK-output is implemented!"<<std::endl;

        WriteScalarDataToFile( (int)cell_type, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
    }


}

void VtkOutput::WriteNodalResultsAsPointData(const ModelPart &rModelPart, std::ofstream& rFileStream)
{
    const auto& local_mesh = rModelPart.GetCommunicator().LocalMesh();
    // write nodal results header
    Parameters nodal_results = this->mrOutputSettings["nodal_solution_step_data_variables"];
    rFileStream << "POINT_DATA " << local_mesh.NumberOfNodes() << "\n";

    for (unsigned int entry = 0; entry < nodal_results.size(); entry++)
    {
        // write nodal results variable header
        std::string nodal_result_name = nodal_results[entry].GetString();
         VtkOutput::WriteDataType data_characteristic; // 0: unknown, 1: Scalar value, 2: 3 component vector
        if (KratosComponents<Variable<double>>::Has(nodal_result_name))
        {
            data_characteristic =  VtkOutput::WriteDataType::VTK_SCALAR;
            rFileStream << "SCALARS " << nodal_result_name << " double"
                       << " 1"
                       << "\n";
            rFileStream << "LOOKUP_TABLE default"
                       << "\n";
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(nodal_result_name))
        {
            data_characteristic =  VtkOutput::WriteDataType::VTK_VECTOR;
            rFileStream << "VECTORS " << nodal_result_name << " double"
                       << "\n";
        }

        // write nodal results
        if (data_characteristic == VtkOutput::WriteDataType::VTK_SCALAR)
        {
            for(const auto& node : local_mesh.Nodes())
            {
                Variable<double> nodal_result_variable = KratosComponents<Variable<double>>::Get(nodal_result_name);
                const auto &nodal_result = node.FastGetSolutionStepValue(nodal_result_variable);
                WriteScalarDataToFile(nodal_result, rFileStream);
                rFileStream <<"\n";
            }
        }
        else if (data_characteristic == VtkOutput::WriteDataType::VTK_VECTOR)
        {
            for(const auto& node : local_mesh.Nodes())
            {
                Variable<array_1d<double, 3>> nodal_result_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(nodal_result_name);
                const auto &nodal_result = node.FastGetSolutionStepValue(nodal_result_variable);
                WriteVectorDataToFile(nodal_result, rFileStream);
                rFileStream << "\n";
            }
        }
    }
}

void VtkOutput::WriteElementData(const ModelPart &rModelPart, std::ofstream& rFileStream)
{
    const auto& local_mesh = rModelPart.GetCommunicator().LocalMesh();
    Parameters element_results = this->mrOutputSettings["element_data_value_variables"];// list of element results

    // write cells header
    if (local_mesh.NumberOfElements() > 0 || local_mesh.NumberOfConditions() > 0)
    {
        rFileStream << "CELL_DATA " << local_mesh.NumberOfElements() + local_mesh.NumberOfConditions()<< "\n";
        for (unsigned int entry = 0; entry < element_results.size(); entry++)
        {

            std::string element_result_name = element_results[entry].GetString();
            VtkOutput::WriteDataType data_characteristic; // 0: unknown, 1: Scalar value, 2: 3 component vector

            if (KratosComponents<Variable<double>>::Has(element_result_name))
            {
                data_characteristic =  VtkOutput::WriteDataType::VTK_SCALAR;
                rFileStream << "SCALARS " << element_result_name << " double"
                           << " 1"
                           << "\n";
                rFileStream << "LOOKUP_TABLE default"
                           << "\n";
            }
            else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(element_result_name))
            {
                data_characteristic =  VtkOutput::WriteDataType::VTK_VECTOR;
                rFileStream << "VECTORS " << element_result_name << " double"
                           << "\n";
            }

            // write nodal results
            if (data_characteristic ==  VtkOutput::WriteDataType::VTK_SCALAR)
            {
                for (const auto& elem : local_mesh.Elements())
                {
                    Variable<double> element_result_variable = KratosComponents<Variable<double>>::Get(element_result_name);
                    const auto &element_result = elem.GetValue(element_result_variable);
                    WriteScalarDataToFile(element_result, rFileStream);
                    rFileStream <<"\n";
                }
                for (const auto& condition : local_mesh.Conditions())
                {
                    Variable<double> element_result_variable = KratosComponents<Variable<double>>::Get(element_result_name);
                    const auto &condition_result = condition.GetValue(element_result_variable);
                    WriteScalarDataToFile(condition_result, rFileStream);
                    rFileStream <<"\n";
                }
            }
            else if (data_characteristic == VtkOutput::WriteDataType::VTK_VECTOR)
            {
                for (const auto& elem : local_mesh.Elements())
                {
                    Variable<array_1d<double, 3>> element_result_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(element_result_name);
                    const auto &element_result = elem.GetValue(element_result_variable);
                    WriteVectorDataToFile(element_result, rFileStream);
                    rFileStream << "\n";
                }

                for (const auto& condition : local_mesh.Conditions())
                {
                    Variable<array_1d<double, 3>> element_result_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(element_result_name);
                    const auto &condition_result = condition.GetValue(element_result_variable);
                    WriteVectorDataToFile(condition_result, rFileStream);
                    rFileStream << "\n";
                }
            }
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
    WriteNodalResultsAsPointData(rModelPart, output_file);
    WriteElementData(rModelPart, output_file);

    output_file.close();
}

void VtkOutput::PrintOutput()
{
    //For whole model part
    WriteModelPart(mrModelPart);
    //For sub model parts
    bool print_sub_model_parts = this->mrOutputSettings["output_sub_model_parts"].GetBool();

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
    int step = rModelPart.GetProcessInfo()[STEP];
    std::string output_file_name = mrOutputSettings["folder_name"].GetString() +"/"+
                                    rModelPart.Name() + "_" + std::to_string(rank) + "_" + std::to_string(step) + ".vtk";
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
        ForceBigEndian((unsigned char *)&rData);
        rFileStream.write((char *)(&rData), sizeof(TData));
    }
}

template <class TData>
void VtkOutput::WriteVectorDataToFile(const TData& rData, std::ofstream& rFileStream)
{

    if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII)
    {
        rFileStream << rData[0] << " ";
        rFileStream << rData[1] << " ";
        rFileStream << rData[2] << " ";
    }
    else if (mFileFormat == VtkOutput::FileFormat::VTK_BINARY)
    {
        ForceBigEndian((unsigned char *)&rData[0]);
        rFileStream.write((char *)(&rData[0]), sizeof(rData[0]));

        ForceBigEndian((unsigned char *)&rData[1]);
        rFileStream.write((char *)(&rData[1]), sizeof(rData[1]));

        ForceBigEndian((unsigned char *)&rData[2]);
        rFileStream.write((char *)(&rData[2]), sizeof(rData[2]));
    }

}

} // namespace Kratos
