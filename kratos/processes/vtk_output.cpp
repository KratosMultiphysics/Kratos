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

namespace Kratos
{

VtkOutput::VtkOutput(ModelPart &rModelPart, Parameters rParameters) : mrModelPart(rModelPart), mrOutputSettings(rParameters)
{
    mDefaultPrecision = 7;
    mStep = 0;
    this->mDoneTest = false;
    this->mShouldSwap = false;
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

    for (const auto& elem : rModelPart.Elements())
    {
        vtk_cell_list_size++;
        vtk_cell_list_size += elem.GetGeometry().size();
    }

    for (const auto& condition : rModelPart.Conditions())
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

void VtkOutput::WriteHeader(const ModelPart &rModelPart)
{
    std::string output_file_name = GetOutputFileName(rModelPart);
    std::ofstream output_file;

    output_file.open(output_file_name, std::ios::out | std::ios::binary | std::ios::trunc);
    output_file << "# vtk DataFile Version 4.0"
               << "\n";
    output_file << "vtk output"
               << "\n";
    output_file << "ASCII"
               << "\n";
    output_file << "DATASET UNSTRUCTURED_GRID"
               << "\n";
    output_file.close();
}

void VtkOutput::WriteMesh(const ModelPart &rModelPart)
{
    WriteNodes(rModelPart);
    WriteConditionsAndElements(rModelPart);
    WriteConditionAndElementTypes(rModelPart);
}

void VtkOutput::WriteNodes(const ModelPart &rModelPart)
{
    std::string output_file_name = GetOutputFileName(rModelPart);
    std::ofstream output_file;
    output_file.open(output_file_name, std::ios::out | std::ios::app | std::ios::binary);
    output_file << std::scientific;
    output_file << std::setprecision(mDefaultPrecision);

    // write nodes header
    output_file << "POINTS " << rModelPart.NumberOfNodes() << " double"
               << "\n";

    // write nodes
    for(const auto& node : rModelPart.Nodes())
    {
        const auto& coordinates = node.Coordinates();
        output_file << " " << coordinates(0);
        output_file << " " << coordinates(1);
        output_file << " " << coordinates(2) << "\n";
    }

    output_file.close();
}

void VtkOutput::WriteConditionsAndElements(const ModelPart &rModelPart)
{
    std::string output_file_name = GetOutputFileName(rModelPart);
    std::ofstream output_file;
    output_file.open(output_file_name, std::ios::out | std::ios::app | std::ios::binary);

    // write cells header
    output_file << "CELLS " << rModelPart.NumberOfConditions() + rModelPart.NumberOfElements() << " " << mVtkCellListSize << "\n";

    // write elements
    for (const auto& elem : rModelPart.Elements())
    {
        auto& elem_geometry = elem.GetGeometry();
        const unsigned int numberOfNodes = elem_geometry.size();

        output_file << numberOfNodes;
        for (unsigned int i = 0; i < numberOfNodes; i++)
            output_file << " " << mKratosIdToVtkId[elem_geometry[i].Id()];

        output_file << "\n";
    }

    // write Conditions
    for (const auto& condition : rModelPart.Conditions())
    {
        auto& condition_geometry = condition.GetGeometry();
        const unsigned int numberOfNodes = condition_geometry.size();

        output_file << numberOfNodes;
        for (unsigned int i = 0; i < numberOfNodes; i++)
            output_file << " " << mKratosIdToVtkId[condition_geometry[i].Id()];
        output_file << "\n";
    }

    output_file.close();
}

void VtkOutput::WriteConditionAndElementTypes(const ModelPart &rModelPart)
{
    std::string output_file_name = GetOutputFileName(rModelPart);
    std::ofstream output_file;
    output_file.open(output_file_name, std::ios::out | std::ios::app | std::ios::binary);

    // write cell types header
    output_file << "CELL_TYPES " << rModelPart.NumberOfConditions() + rModelPart.NumberOfElements() << "\n";

    // write elements types
    for (const auto& elem : rModelPart.Elements())
    {
        const unsigned int numberOfNodes = elem.GetGeometry().size();
        unsigned int element_type;

        if (numberOfNodes == 3)
            element_type = 5;
        else if (numberOfNodes == 4)
        {
            if (mrModelPart.GetProcessInfo()[DOMAIN_SIZE] == 2)
                element_type = 9;
            else
                element_type = 10;
        }
        else if (numberOfNodes == 2)
            element_type = 3;
        else if (numberOfNodes == 1)
            element_type = 1;
        else
            KRATOS_THROW_ERROR(std::runtime_error, "Modelpart contains elements with geometries for which no VTK-output is implemented!", "")

        output_file << element_type << "\n";
    }

    // write conditions types
    for (const auto& condition : rModelPart.Conditions())
    {
        const unsigned int numberOfNodes = condition.GetGeometry().size();
        unsigned int element_type;

        if (numberOfNodes == 3)
            element_type = 5;
        else if (numberOfNodes == 4)
        {
            if (mrModelPart.GetProcessInfo()[DOMAIN_SIZE] == 2)
                element_type = 9;
            else
                element_type = 10;
        }
        else if (numberOfNodes == 2)
            element_type = 3;
        else if (numberOfNodes == 1)
            element_type = 1;
        else
            KRATOS_ERROR << "Modelpart contains conditions with geometries for which no VTK-output is implemented!" << std::endl;

        output_file << element_type << "\n";
    }

    output_file.close();
}

void VtkOutput::WriteNodalResultsAsPointData(const ModelPart &rModelPart)
{
    std::string output_file_name = GetOutputFileName(rModelPart);
    std::ofstream output_file;
    output_file.open(output_file_name, std::ios::out | std::ios::app | std::ios::binary);
    // write nodal results header
    Parameters nodalResults = this->mrOutputSettings["nodal_solution_step_data_variables"];
    output_file << "POINT_DATA " << rModelPart.NumberOfNodes() << "\n";

    for (unsigned int entry = 0; entry < nodalResults.size(); entry++)
    {
        // write nodal results variable header
        std::string nodal_result_name = nodalResults[entry].GetString();
        unsigned int dataCharacteristic = 0; // 0: unknown, 1: Scalar value, 2: 3 component vector
        if (KratosComponents<Variable<double>>::Has(nodal_result_name))
        {
            dataCharacteristic = 1;
            output_file << "SCALARS " << nodal_result_name << " double"
                       << " 1"
                       << "\n";
            output_file << "LOOKUP_TABLE default"
                       << "\n";
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(nodal_result_name))
        {
            dataCharacteristic = 2;
            output_file << "VECTORS " << nodal_result_name << " double"
                       << "\n";
        }

        // write nodal results
        output_file << std::scientific;
        output_file << std::setprecision(mDefaultPrecision);
        for(const auto& node : rModelPart.Nodes())
        {
            if (dataCharacteristic == 1)
            {
                Variable<double> nodal_result_variable = KratosComponents<Variable<double>>::Get(nodal_result_name);
                const double &nodal_result = node.FastGetSolutionStepValue(nodal_result_variable);
                output_file << nodal_result << "\n";
            }
            else if (dataCharacteristic == 2)
            {
                Variable<array_1d<double, 3>> nodal_result_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(nodal_result_name);
                const array_1d<double, 3> &nodal_result = node.FastGetSolutionStepValue(nodal_result_variable);
                output_file << nodal_result[0] << " ";
                output_file << nodal_result[1] << " ";
                output_file << nodal_result[2] << "\n";
            }
        }
    }

    output_file.close();
}

void VtkOutput::WriteElementData(const ModelPart &rModelPart)
{
    std::string output_file_name = GetOutputFileName(rModelPart);
    std::ofstream output_file;
    output_file.open(output_file_name, std::ios::out | std::ios::app);
    Parameters elementResults = this->mrOutputSettings["element_data_value_variables"];// list of element results

    // write cells header
    if (rModelPart.NumberOfElements() > 0)
    {
        output_file << "CELL_DATA " << rModelPart.NumberOfElements() << "\n";
        for (unsigned int entry = 0; entry < elementResults.size(); entry++)
        {

            std::string elementResultName = elementResults[entry].GetString();
            unsigned int dataCharacteristic = 0; // 0: unknown, 1: Scalar value, 2: 3 component vector

            if (KratosComponents<Variable<double>>::Has(elementResultName))
            {
                dataCharacteristic = 1;
                output_file << "SCALARS " << elementResultName << " double"
                           << " 1"
                           << "\n";
                output_file << "LOOKUP_TABLE default"
                           << "\n";
            }
            else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(elementResultName))
            {
                dataCharacteristic = 2;
                output_file << "VECTORS " << elementResultName << " double"
                           << "\n";
            }

            // write nodal results
            output_file << std::scientific;
            output_file << std::setprecision(mDefaultPrecision);
            for (const auto& elem : rModelPart.Elements())
            {
                if (dataCharacteristic == 1)
                {
                    Variable<double> elementResultVariable = KratosComponents<Variable<double>>::Get(elementResultName);
                    const double &elementResult = elem.GetValue(elementResultVariable);
                    output_file << elementResult << "\n";
                }
                else if (dataCharacteristic == 2)
                {
                    Variable<array_1d<double, 3>> elementResultVariable = KratosComponents<Variable<array_1d<double, 3>>>::Get(elementResultName);
                    const array_1d<double, 3> &elementResult = elem.GetValue(elementResultVariable);
                    output_file << elementResult[0] << " ";
                    output_file << elementResult[1] << " ";
                    output_file << elementResult[2] << "\n";
                }
            }
        }
        output_file.close();
    }
}

//#############################################For creating vtk files in binary format##########################################################

void VtkOutput::WriteHeaderBinary(const ModelPart &rModelPart)
{
    std::string output_file_name = GetOutputFileName(rModelPart);
    std::ofstream output_file;

    output_file.open(output_file_name, std::ios::out | std::ios::binary);
    output_file << "# vtk DataFile Version 4.0"
               << "\n";
    output_file << "vtk output"
               << "\n";
    output_file << "BINARY"
               << "\n";
    output_file << "DATASET UNSTRUCTURED_GRID"
               << "\n";
    output_file.close();
}

void VtkOutput::WriteMeshBinary(const ModelPart &rModelPart)
{

    WriteNodesBinary(rModelPart);

    WriteConditionsAndElementsBinary(rModelPart);

    WriteConditionAndElementTypesBinary(rModelPart);
}

void VtkOutput::WriteNodesBinary(const ModelPart &rModelPart)
{
    std::string output_file_name = GetOutputFileName(rModelPart);
    std::ofstream output_file;
    output_file.open(output_file_name, std::ios::out | std::ios::app | std::ios::binary);

    // write nodes header
    output_file << "\nPOINTS " << rModelPart.NumberOfNodes() << " double"
               << "\n";

    // write nodes
    for(const auto& node : rModelPart.Nodes())
    {
        const auto& coordinates = node.Coordinates();
        double x_coordinate = coordinates(0);
        double y_coordinate = coordinates(1);
        double z_coordinate = coordinates(2);
        ForceBigEndian((unsigned char *)&coordinates(0));
        output_file.write((char *)(&x_coordinate), sizeof(double));
        ForceBigEndian((unsigned char *)&coordinates(1));
        output_file.write((char *)(&y_coordinate), sizeof(double));
        ForceBigEndian((unsigned char *)&coordinates(2));
        output_file.write((char *)(&z_coordinate), sizeof(double));
    }

    output_file.close();
}

void VtkOutput::WriteConditionsAndElementsBinary(const ModelPart &rModelPart)
{
    std::string output_file_name = GetOutputFileName(rModelPart);
    std::ofstream output_file;
    output_file.open(output_file_name, std::ios::out | std::ios::app | std::ios::binary);

    // write cells header
    output_file << "\nCELLS " << rModelPart.NumberOfConditions() + rModelPart.NumberOfElements() << " " << mVtkCellListSize << "\n";

    // write elements
    for (const auto& elem : rModelPart.Elements())
    {

        const ModelPart::ConditionType::GeometryType &elem_geometry = elem.GetGeometry();

        unsigned int numberOfNodes = elem_geometry.size();

        ForceBigEndian((unsigned char *)&numberOfNodes);

        output_file.write((char *)(&numberOfNodes), sizeof(unsigned int));

        for (unsigned int i = 0; i < elem_geometry.size(); i++)
        {
            int nodenum = mKratosIdToVtkId[elem_geometry[i].Id()];
            ForceBigEndian((unsigned char *)&nodenum);
            output_file.write((char *)(&nodenum), sizeof(int));
        }
    }

    // write Conditions
    for (const auto& condition : rModelPart.Conditions())
    {
        const ModelPart::ConditionType::GeometryType &condition_geometry = condition.GetGeometry();
        unsigned int numberOfNodes = condition_geometry.size();

        ForceBigEndian((unsigned char *)&numberOfNodes);
        output_file.write((char *)(&numberOfNodes), sizeof(unsigned int));

        for (unsigned int i = 0; i < condition_geometry.size(); i++)
        {

            int nodenum = mKratosIdToVtkId[condition_geometry[i].Id()];
            ForceBigEndian((unsigned char *)&nodenum);
            output_file.write((char *)(&nodenum), sizeof(int));
        }
    }

    output_file.close();
}

void VtkOutput::WriteConditionAndElementTypesBinary(const ModelPart &rModelPart)
{
    std::string output_file_name = GetOutputFileName(rModelPart);
    std::ofstream output_file;
    output_file.open(output_file_name, std::ios::out | std::ios::app | std::ios::binary);

    // write cell types header
    output_file << "\nCELL_TYPES " << rModelPart.NumberOfConditions() + rModelPart.NumberOfElements() << "\n";

    // write elements types
    for (const auto& elem : rModelPart.Elements())
    {
        const unsigned int numberOfNodes = elem.GetGeometry().size();
        unsigned int element_type;

        if (numberOfNodes == 3)
            element_type = 5;
        else if (numberOfNodes == 4)
        {
            if (mrModelPart.GetProcessInfo()[DOMAIN_SIZE] == 2)
                element_type = 9;
            else
                element_type = 10;
        }
        else if (numberOfNodes == 2)
            element_type = 3;
        else if (numberOfNodes == 1)
            element_type = 1;
        else
            KRATOS_ERROR << "Modelpart contains elements with geometries for which no VTK-output is implemented!" << std::endl;

        ForceBigEndian((unsigned char *)&element_type);
        output_file.write((char *)(&element_type), sizeof(int));
    }

    // write conditions types
    for (const auto& condition : rModelPart.Conditions())
    {
        const unsigned int numberOfNodes = condition.GetGeometry().size();
        unsigned int element_type;

        if (numberOfNodes == 3)
            element_type = 5;
        else if (numberOfNodes == 4)
        {
            if (mrModelPart.GetProcessInfo()[DOMAIN_SIZE] == 2)
                element_type = 9;
            else
                element_type = 10;
        }
        else if (numberOfNodes == 2)
            element_type = 3;
        else if (numberOfNodes == 1)
            element_type = 1;
        else
            KRATOS_ERROR << "Modelpart contains conditions with geometries for which no VTK-output is implemented!" << std::endl;

        ForceBigEndian((unsigned char *)&element_type);
        output_file.write((char *)(&element_type), sizeof(int));
    }

    output_file.close();
}

void VtkOutput::WriteNodalResultsAsPointDataBinary(const ModelPart &rModelPart)
{
    std::string output_file_name = GetOutputFileName(rModelPart);
    std::ofstream output_file;
    output_file.open(output_file_name, std::ios::out | std::ios::app | std::ios::binary);
    // write nodal results header
    Parameters nodalResults = this->mrOutputSettings["nodal_solution_step_data_variables"];
    output_file << "\nPOINT_DATA " << rModelPart.NumberOfNodes() << "\n";

    for (unsigned int entry = 0; entry < nodalResults.size(); entry++)
    {
        // write nodal results variable header
        std::string nodal_result_name = nodalResults[entry].GetString();
        unsigned int dataCharacteristic = 0; // 0: unknown, 1: Scalar value, 2: 3 DOF global translation vector
        if (KratosComponents<Variable<double>>::Has(nodal_result_name))
        {
            dataCharacteristic = 1;
            output_file << "SCALARS " << nodal_result_name << " double"
                       << "\n";
            output_file << "LOOKUP_TABLE default"
                       << "\n";
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(nodal_result_name))
        {
            dataCharacteristic = 2;
            output_file << "VECTORS " << nodal_result_name << " double"
                       << "\n";
        }

        // write nodal results

        for(const auto& node : rModelPart.Nodes())
        {
            if (dataCharacteristic == 1)
            {
                Variable<double> nodal_result_variable = KratosComponents<Variable<double>>::Get(nodal_result_name);
                const double nodal_result = node.FastGetSolutionStepValue(nodal_result_variable);
                ForceBigEndian((unsigned char *)&nodal_result);
                output_file.write((char *)(&nodal_result), sizeof(double));
            }
            else if (dataCharacteristic == 2)
            {
                Variable<array_1d<double, 3>> nodal_result_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(nodal_result_name);
                const array_1d<double, 3> nodal_result = node.FastGetSolutionStepValue(nodal_result_variable);
                double num1 = nodal_result[0];
                ForceBigEndian((unsigned char *)&num1);
                output_file.write((char *)(&num1), sizeof(double));
                double num2 = nodal_result[1];
                ForceBigEndian((unsigned char *)&num2);
                output_file.write((char *)(&num2), sizeof(double));
                double num3 = nodal_result[2];
                ForceBigEndian((unsigned char *)&num3);
                output_file.write((char *)(&num3), sizeof(double));
            }
        }
    }

    output_file.close();
}

void VtkOutput::WriteElementDataBinary(const ModelPart &rModelPart)
{
    std::string output_file_name = GetOutputFileName(rModelPart);
    std::ofstream output_file;
    output_file.open(output_file_name, std::ios::out | std::ios::app);
    Parameters elementResults = this->mrOutputSettings["element_data_value_variables"];// list of element resultss
    if (rModelPart.NumberOfElements() > 0)
    {
        // write cells header
        output_file << "\nCELL_DATA " << rModelPart.NumberOfElements() << "\n";

        for (unsigned int entry = 0; entry < elementResults.size(); entry++)
        {

            std::string elementResultName = elementResults[entry].GetString();
            unsigned int dataCharacteristic = 0; // 0: unknown, 1: Scalar value, 2: 3 DOF global translation vector

            if (KratosComponents<Variable<double>>::Has(elementResultName))
            {
                dataCharacteristic = 1;
                output_file << "SCALARS " << elementResultName << " double"
                           << "\n";
                output_file << "LOOKUP_TABLE default"
                           << "\n";
            }
            else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(elementResultName))
            {
                dataCharacteristic = 2;
                output_file << "VECTORS " << elementResultName << " double"
                           << "\n";
            }

            // write nodal results

            for (const auto& elem : rModelPart.Elements())
            {
                if (dataCharacteristic == 1)
                {
                    Variable<double> elementResultVariable = KratosComponents<Variable<double>>::Get(elementResultName);
                    double elementResult = elem.GetValue(elementResultVariable);
                    ForceBigEndian((unsigned char *)&elementResult);
                    output_file.write((char *)(&elementResult), sizeof(double));
                }
                else if (dataCharacteristic == 2)
                {
                    Variable<array_1d<double, 3>> elementResultVariable = KratosComponents<Variable<array_1d<double, 3>>>::Get(elementResultName);
                    array_1d<double, 3> elementResult = elem.GetValue(elementResultVariable);
                    double num1 = elementResult[0];
                    ForceBigEndian((unsigned char *)&num1);
                    output_file.write((char *)(&num1), sizeof(double));
                    double num2 = elementResult[1];
                    ForceBigEndian((unsigned char *)&num2);
                    output_file.write((char *)(&num2), sizeof(double));
                    double num3 = elementResult[2];
                    ForceBigEndian((unsigned char *)&num3);
                    output_file.write((char *)(&num3), sizeof(double));
                }
            }
        }
        output_file.close();
    }
}

//#################################################################End of Binary vtk ################################################################

void VtkOutput::WriteModelPart(const ModelPart &rModelPart)
{
    Initialize(rModelPart);
    std::string type = this->mrOutputSettings["file_format"].GetString();
    if (type == "Ascii")
    {
        WriteHeader(rModelPart);
        WriteMesh(rModelPart);
        WriteNodalResultsAsPointData(rModelPart);
        WriteElementData(rModelPart);
    }
    else if (type == "Binary")
    {
        WriteHeaderBinary(rModelPart);
        WriteMeshBinary(rModelPart);
        WriteNodalResultsAsPointDataBinary(rModelPart);
        WriteElementDataBinary(rModelPart);
    }
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
            WriteModelPart(sub_model_part);
        }
    }
    ++mStep;
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

///@}

std::string VtkOutput::GetOutputFileName(const ModelPart &rModelPart)
{
    int rank = 0;
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    std::string outputFilename = mrOutputSettings["folder_name"].GetString() +"/"+
                                    rModelPart.Name() + "_" + std::to_string(rank) + "_" + std::to_string(mStep) + ".vtk";
    return outputFilename;
}

} // namespace Kratos
