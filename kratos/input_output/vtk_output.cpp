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

unsigned int VtkOutput::DetermineVtkCellListSize(const ModelPart &rModelPart)
{
    unsigned int vtk_cell_list_size_elem = 0;
    const auto& local_mesh = rModelPart.GetCommunicator().LocalMesh();

    const auto elements_begin = local_mesh.ElementsBegin();
    const int num_elements = local_mesh.NumberOfElements();
    #pragma omp parallel for reduction(+:vtk_cell_list_size_elem)
    for (int i=0; i<num_elements; ++i)
    {
        auto element_i = elements_begin + i;
        vtk_cell_list_size_elem += element_i->GetGeometry().PointsNumber() +1;
    }

    unsigned int vtk_cell_list_size_cond = 0;
    const auto conditions_begin = local_mesh.ConditionsBegin();
    const int num_conditions = local_mesh.NumberOfConditions();
    #pragma omp parallel for reduction(+:vtk_cell_list_size_cond)
    for (int i=0; i<num_conditions; ++i)
    {
        auto conditions_i = conditions_begin + i;
        vtk_cell_list_size_cond += conditions_i->GetGeometry().PointsNumber() +1;
    }

    return vtk_cell_list_size_cond+vtk_cell_list_size_elem;
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

    WriteConnectivity<Element> (local_mesh.Elements().GetContainer(), rFileStream);
    WriteConnectivity<Condition> (local_mesh.Conditions().GetContainer(), rFileStream);
}

template <typename TEntity>
void VtkOutput::WriteConnectivity(const typename PointerVectorSet<TEntity, IndexedObject>::ContainerType& rContainer, std::ofstream& rFileStream)
{
    // write Conditions
    for (const auto& entity : rContainer)
    {
        const auto& entity_geometry = entity->GetGeometry();
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

void VtkOutput::WriteConditionAndElementTypes(const ModelPart &rModelPart, std::ofstream& rFileStream)
{

    const auto& local_mesh = rModelPart.GetCommunicator().LocalMesh();
    // write cell types header
    rFileStream << "\nCELL_TYPES " << local_mesh.NumberOfConditions() + local_mesh.NumberOfElements() << "\n";
    WriteCellType<Element> (local_mesh.Elements().GetContainer(), rFileStream);
    WriteCellType<Condition> (local_mesh.Conditions().GetContainer(), rFileStream);
}

template <typename TEntity>
void VtkOutput::WriteCellType(const typename PointerVectorSet<TEntity, IndexedObject>::ContainerType& rContainer, std::ofstream& rFileStream)
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
        if( geo_type_vtk_cell_type_map.count(entity->GetGeometry().GetGeometryType()) > 0  )
            cell_type = geo_type_vtk_cell_type_map[ entity->GetGeometry().GetGeometryType() ];
        else
            KRATOS_ERROR<<"Modelpart contains elements or conditions with geometries for which no VTK-output is implemented!"<<std::endl;

        WriteScalarDataToFile( (int)cell_type, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream << "\n";
    }
}


void VtkOutput::WriteNodalResultsAsPointData(const ModelPart &rModelPart, std::ofstream& rFileStream)
{
    const auto& local_mesh = rModelPart.GetCommunicator().LocalMesh();
    // write nodal results header
    Parameters nodal_results = this->mOutputSettings["nodal_solution_step_data_variables"];
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
    Parameters element_results = this->mOutputSettings["element_data_value_variables"];// list of element results

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

    double label = 0;
    std::string output_control = mOutputSettings["output_control_type"].GetString();
    if(output_control == "step")
        label = rModelPart.GetProcessInfo()[STEP];
    else if(output_control == "time")
        label = rModelPart.GetProcessInfo()[TIME];
    else
        KRATOS_ERROR<<"Option for output_control_type : "<<output_control<<" not recognised.!"<<std::endl
            <<"Possible output_control_type options for VTKOutput are :: step and time "<<std::endl;

    // Putting every thing together
    std::string output_file_name = mOutputSettings["folder_name"].GetString() +"/"+
                                    model_part_name + "_" + std::to_string(rank) + "_" + std::to_string(label) + ".vtk";
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
