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
    std::string outputFileName = GetOutputFileName(rModelPart);
    std::ofstream outputFile;

    outputFile.open(outputFileName, std::ios::out | std::ios::binary | std::ios::trunc);
    outputFile << "# vtk DataFile Version 4.0"
               << "\n";
    outputFile << "vtk output"
               << "\n";
    outputFile << "ASCII"
               << "\n";
    outputFile << "DATASET UNSTRUCTURED_GRID"
               << "\n";
    outputFile.close();
}

void VtkOutput::WriteMesh(const ModelPart &rModelPart)
{
    WriteNodes(rModelPart);
    WriteConditionsAndElements(rModelPart);
    WriteConditionAndElementTypes(rModelPart);
}

void VtkOutput::WriteNodes(const ModelPart &rModelPart)
{
    std::string outputFileName = GetOutputFileName(rModelPart);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);
    outputFile << std::scientific;
    outputFile << std::setprecision(mDefaultPrecision);

    // write nodes header
    outputFile << "POINTS " << rModelPart.NumberOfNodes() << " float"
               << "\n";

    // write nodes
    for(const auto& node : rModelPart.Nodes())
    {
        auto& coordinates = node.Coordinates();
        outputFile << " " << coordinates(0);
        outputFile << " " << coordinates(1);
        outputFile << " " << coordinates(2) << "\n";
    }

    outputFile.close();
}

void VtkOutput::WriteConditionsAndElements(const ModelPart &rModelPart)
{
    std::string outputFileName = GetOutputFileName(rModelPart);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);

    // write cells header
    outputFile << "CELLS " << rModelPart.NumberOfConditions() + rModelPart.NumberOfElements() << " " << mVtkCellListSize << "\n";

    // write elements
    for (const auto& elem : rModelPart.Elements())
    {
        auto& elem_geometry = elem.GetGeometry();
        const unsigned int numberOfNodes = elem_geometry.size();

        outputFile << numberOfNodes;
        for (unsigned int i = 0; i < numberOfNodes; i++)
            outputFile << " " << mKratosIdToVtkId[elem_geometry[i].Id()];

        outputFile << "\n";
    }

    // write Conditions
    for (const auto& condition : rModelPart.Conditions())
    {
        auto& condition_geometry = condition.GetGeometry();
        const unsigned int numberOfNodes = condition_geometry.size();

        outputFile << numberOfNodes;
        for (unsigned int i = 0; i < numberOfNodes; i++)
            outputFile << " " << mKratosIdToVtkId[condition_geometry[i].Id()];
        outputFile << "\n";
    }

    outputFile.close();
}

void VtkOutput::WriteConditionAndElementTypes(const ModelPart &rModelPart)
{
    std::string outputFileName = GetOutputFileName(rModelPart);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);

    // write cell types header
    outputFile << "CELL_TYPES " << rModelPart.NumberOfConditions() + rModelPart.NumberOfElements() << "\n";

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

        outputFile << element_type << "\n";
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

        outputFile << element_type << "\n";
    }

    outputFile.close();
}

void VtkOutput::WriteNodalResultsAsPointData(const ModelPart &rModelPart)
{
    std::string outputFileName = GetOutputFileName(rModelPart);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);
    // write nodal results header
    Parameters nodalResults = this->mrOutputSettings["nodal_solution_step_data_variables"];
    outputFile << "POINT_DATA " << rModelPart.NumberOfNodes() << "\n";

    for (unsigned int entry = 0; entry < nodalResults.size(); entry++)
    {
        // write nodal results variable header
        std::string nodalResultName = nodalResults[entry].GetString();
        unsigned int dataCharacteristic = 0; // 0: unknown, 1: Scalar value, 2: 3 component vector
        if (KratosComponents<Variable<double>>::Has(nodalResultName))
        {
            dataCharacteristic = 1;
            outputFile << "SCALARS " << nodalResultName << " float"
                       << " 1"
                       << "\n";
            outputFile << "LOOKUP_TABLE default"
                       << "\n";
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(nodalResultName))
        {
            dataCharacteristic = 2;
            outputFile << "VECTORS " << nodalResultName << " float"
                       << "\n";
        }

        // write nodal results
        outputFile << std::scientific;
        outputFile << std::setprecision(mDefaultPrecision);
        for(const auto& node : rModelPart.Nodes())
        {
            if (dataCharacteristic == 1)
            {
                Variable<double> nodalResultVariable = KratosComponents<Variable<double>>::Get(nodalResultName);
                const double &nodalResult = node.FastGetSolutionStepValue(nodalResultVariable);
                outputFile << nodalResult << "\n";
            }
            else if (dataCharacteristic == 2)
            {
                Variable<array_1d<double, 3>> nodalResultVariable = KratosComponents<Variable<array_1d<double, 3>>>::Get(nodalResultName);
                const array_1d<double, 3> &nodalResult = node.FastGetSolutionStepValue(nodalResultVariable);
                outputFile << nodalResult[0] << " ";
                outputFile << nodalResult[1] << " ";
                outputFile << nodalResult[2] << "\n";
            }
        }
    }

    outputFile.close();
}

void VtkOutput::WriteElementData(const ModelPart &rModelPart)
{
    std::string outputFileName = GetOutputFileName(rModelPart);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app);
    Parameters elementResults = this->mrOutputSettings["element_data_value_variables"];// list of element results

    // write cells header
    if (rModelPart.NumberOfElements() > 0)
    {
        outputFile << "CELL_DATA " << rModelPart.NumberOfElements() << "\n";
        for (unsigned int entry = 0; entry < elementResults.size(); entry++)
        {

            std::string elementResultName = elementResults[entry].GetString();
            unsigned int dataCharacteristic = 0; // 0: unknown, 1: Scalar value, 2: 3 component vector

            if (KratosComponents<Variable<double>>::Has(elementResultName))
            {
                dataCharacteristic = 1;
                outputFile << "SCALARS " << elementResultName << " float"
                           << " 1"
                           << "\n";
                outputFile << "LOOKUP_TABLE default"
                           << "\n";
            }
            else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(elementResultName))
            {
                dataCharacteristic = 2;
                outputFile << "VECTORS " << elementResultName << " float"
                           << "\n";
            }

            // write nodal results
            outputFile << std::scientific;
            outputFile << std::setprecision(mDefaultPrecision);
            for (const auto& elem : rModelPart.Elements())
            {
                if (dataCharacteristic == 1)
                {
                    Variable<double> elementResultVariable = KratosComponents<Variable<double>>::Get(elementResultName);
                    const double &elementResult = elem.GetValue(elementResultVariable);
                    outputFile << elementResult << "\n";
                }
                else if (dataCharacteristic == 2)
                {
                    Variable<array_1d<double, 3>> elementResultVariable = KratosComponents<Variable<array_1d<double, 3>>>::Get(elementResultName);
                    const array_1d<double, 3> &elementResult = elem.GetValue(elementResultVariable);
                    outputFile << elementResult[0] << " ";
                    outputFile << elementResult[1] << " ";
                    outputFile << elementResult[2] << "\n";
                }
            }
        }
        /*
        outputFile << "SCALARS SPLIT_ELEMENT float 1\nLOOKUP_TABLE default\n";

        // write element results for active
        for (const auto& elem : rModelPart.Elements())
        {
            //outputFile << numberOfNodes;
            bool is_split = elem.GetValue(SPLIT_ELEMENT);
            outputFile << is_split << "\n";
        }
    */
        outputFile.close();
    }
}

//#############################################For creating vtk files in binary format##########################################################

void VtkOutput::WriteHeaderBinary(const ModelPart &rModelPart)
{
    std::string outputFileName = GetOutputFileName(rModelPart);
    std::ofstream outputFile;

    outputFile.open(outputFileName, std::ios::out | std::ios::binary);
    outputFile << "# vtk DataFile Version 4.0"
               << "\n";
    outputFile << "vtk output"
               << "\n";
    outputFile << "BINARY"
               << "\n";
    outputFile << "DATASET UNSTRUCTURED_GRID"
               << "\n";
    outputFile.close();
}

void VtkOutput::WriteMeshBinary(const ModelPart &rModelPart)
{

    WriteNodesBinary(rModelPart);

    WriteConditionsAndElementsBinary(rModelPart);

    WriteConditionAndElementTypesBinary(rModelPart);
}

void VtkOutput::WriteNodesBinary(const ModelPart &rModelPart)
{
    std::string outputFileName = GetOutputFileName(rModelPart);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);

    // write nodes header
    outputFile << "\nPOINTS " << rModelPart.NumberOfNodes() << " float"
               << "\n";

    // write nodes
    for(const auto& node : rModelPart.Nodes())
    {
        float x_coordinate = node.X();
        float y_coordinate = node.Y();
        float z_coordinate = node.Z();
        ForceBigEndian((unsigned char *)&x_coordinate);
        outputFile.write((char *)(&x_coordinate), sizeof(float));
        ForceBigEndian((unsigned char *)&y_coordinate);
        outputFile.write((char *)(&y_coordinate), sizeof(float));
        ForceBigEndian((unsigned char *)&z_coordinate);
        outputFile.write((char *)(&z_coordinate), sizeof(float));
    }

    outputFile.close();
}

void VtkOutput::WriteConditionsAndElementsBinary(const ModelPart &rModelPart)
{
    std::string outputFileName = GetOutputFileName(rModelPart);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);

    // write cells header
    outputFile << "\nCELLS " << rModelPart.NumberOfConditions() + rModelPart.NumberOfElements() << " " << mVtkCellListSize << "\n";

    // write elements
    for (const auto& elem : rModelPart.Elements())
    {

        const ModelPart::ConditionType::GeometryType &elem_geometry = elem.GetGeometry();

        unsigned int numberOfNodes = elem_geometry.size();

        ForceBigEndian((unsigned char *)&numberOfNodes);

        outputFile.write((char *)(&numberOfNodes), sizeof(unsigned int));

        for (unsigned int i = 0; i < elem_geometry.size(); i++)
        {
            int nodenum = mKratosIdToVtkId[elem_geometry[i].Id()];
            ForceBigEndian((unsigned char *)&nodenum);
            outputFile.write((char *)(&nodenum), sizeof(int));
        }
    }

    // write Conditions
    for (const auto& condition : rModelPart.Conditions())
    {
        const ModelPart::ConditionType::GeometryType &condition_geometry = condition.GetGeometry();
        unsigned int numberOfNodes = condition_geometry.size();

        ForceBigEndian((unsigned char *)&numberOfNodes);
        outputFile.write((char *)(&numberOfNodes), sizeof(unsigned int));

        for (unsigned int i = 0; i < condition_geometry.size(); i++)
        {

            int nodenum = mKratosIdToVtkId[condition_geometry[i].Id()];
            ForceBigEndian((unsigned char *)&nodenum);
            outputFile.write((char *)(&nodenum), sizeof(int));
        }
    }

    outputFile.close();
}

void VtkOutput::WriteConditionAndElementTypesBinary(const ModelPart &rModelPart)
{
    std::string outputFileName = GetOutputFileName(rModelPart);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);

    // write cell types header
    outputFile << "\nCELL_TYPES " << rModelPart.NumberOfConditions() + rModelPart.NumberOfElements() << "\n";

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
        outputFile.write((char *)(&element_type), sizeof(int));
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
        outputFile.write((char *)(&element_type), sizeof(int));
    }

    outputFile.close();
}

void VtkOutput::WriteNodalResultsAsPointDataBinary(const ModelPart &rModelPart)
{
    std::string outputFileName = GetOutputFileName(rModelPart);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);
    // write nodal results header
    Parameters nodalResults = this->mrOutputSettings["nodal_solution_step_data_variables"];
    outputFile << "\nPOINT_DATA " << rModelPart.NumberOfNodes() << "\n";

    for (unsigned int entry = 0; entry < nodalResults.size(); entry++)
    {
        // write nodal results variable header
        std::string nodalResultName = nodalResults[entry].GetString();
        unsigned int dataCharacteristic = 0; // 0: unknown, 1: Scalar value, 2: 3 DOF global translation vector
        if (KratosComponents<Variable<double>>::Has(nodalResultName))
        {
            dataCharacteristic = 1;
            outputFile << "SCALARS " << nodalResultName << " float"
                       << "\n";
            outputFile << "LOOKUP_TABLE default"
                       << "\n";
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(nodalResultName))
        {
            dataCharacteristic = 2;
            outputFile << "VECTORS " << nodalResultName << " float"
                       << "\n";
        }

        // write nodal results

        for(const auto& node : rModelPart.Nodes())
        {
            if (dataCharacteristic == 1)
            {
                Variable<double> nodalResultVariable = KratosComponents<Variable<double>>::Get(nodalResultName);
                float nodalResult = node.FastGetSolutionStepValue(nodalResultVariable);
                ForceBigEndian((unsigned char *)&nodalResult);
                outputFile.write((char *)(&nodalResult), sizeof(float));
            }
            else if (dataCharacteristic == 2)
            {
                Variable<array_1d<double, 3>> nodalResultVariable = KratosComponents<Variable<array_1d<double, 3>>>::Get(nodalResultName);
                array_1d<double, 3> nodalResult = node.FastGetSolutionStepValue(nodalResultVariable);
                float num1 = nodalResult[0];
                ForceBigEndian((unsigned char *)&num1);
                outputFile.write((char *)(&num1), sizeof(float));
                float num2 = nodalResult[1]; 
                ForceBigEndian((unsigned char *)&num2);
                outputFile.write((char *)(&num2), sizeof(float));
                float num3 = nodalResult[2];
                ForceBigEndian((unsigned char *)&num3);
                outputFile.write((char *)(&num3), sizeof(float));
            }
        }
    }

    outputFile.close();
}

void VtkOutput::WriteElementDataBinary(const ModelPart &rModelPart)
{
    std::string outputFileName = GetOutputFileName(rModelPart);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app);
    Parameters elementResults = this->mrOutputSettings["element_data_value_variables"];// list of element resultss
    if (rModelPart.NumberOfElements() > 0)
    {
        // write cells header
        outputFile << "\nCELL_DATA " << rModelPart.NumberOfElements() << "\n";

        for (unsigned int entry = 0; entry < elementResults.size(); entry++)
        {

            std::string elementResultName = elementResults[entry].GetString();
            unsigned int dataCharacteristic = 0; // 0: unknown, 1: Scalar value, 2: 3 DOF global translation vector

            if (KratosComponents<Variable<double>>::Has(elementResultName))
            {
                dataCharacteristic = 1;
                outputFile << "SCALARS " << elementResultName << " float"
                           << "\n";
                outputFile << "LOOKUP_TABLE default"
                           << "\n";
            }
            else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(elementResultName))
            {
                dataCharacteristic = 2;
                outputFile << "VECTORS " << elementResultName << " float"
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
                    outputFile.write((char *)(&elementResult), sizeof(float));
                }
                else if (dataCharacteristic == 2)
                {
                    Variable<array_1d<double, 3>> elementResultVariable = KratosComponents<Variable<array_1d<double, 3>>>::Get(elementResultName);
                    array_1d<double, 3> elementResult = elem.GetValue(elementResultVariable);
                    float num1 = elementResult[0];
                    ForceBigEndian((unsigned char *)&num1);
                    outputFile.write((char *)(&num1), sizeof(float));
                    float num2 = elementResult[1];
                    ForceBigEndian((unsigned char *)&num2);
                    outputFile.write((char *)(&num2), sizeof(float));
                    float num3 = elementResult[2];
                    ForceBigEndian((unsigned char *)&num3);
                    outputFile.write((char *)(&num3), sizeof(float));
                }
            }
        }
        outputFile.close();
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
