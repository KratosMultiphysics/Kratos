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
#include "vtk_output_process.h"

namespace Kratos
{

VtkOutputProcess::VtkOutputProcess(ModelPart &model_part, Parameters rParameters) : mrModelPart(model_part), mrOutputSettings(rParameters)
{
    mDefaultPrecision = 7;
    mStep = 0;
    this->mDoneTest = false;
    this->mShouldSwap = false;
}

VtkOutputProcess::~VtkOutputProcess(){};

void VtkOutputProcess::ExecuteInitialize()
{}

void VtkOutputProcess::ExecuteBeforeSolutionLoop()
{}

void VtkOutputProcess::ExecuteInitializeSolutionStep()
{}

void VtkOutputProcess::ExecuteFinalizeSolutionStep()
{
    PrintOutput();
}

void VtkOutputProcess::ExecuteBeforeOutputStep()
{}

void VtkOutputProcess::ExecuteAfterOutputStep()
{}

void VtkOutputProcess::ExecuteFinalize()
{}

int VtkOutputProcess::Check()
{return  0;}



void VtkOutputProcess::CreateMapFromKratosIdToVTKId(ModelPart &model_part)
{
    int vtk_id = 0;

    for (ModelPart::NodeIterator node_i = model_part.NodesBegin(); node_i != model_part.NodesEnd(); ++node_i)
    {
        int KratosId = node_i->Id();
        mKratosIdToVtkId[KratosId] = vtk_id;
        vtk_id++;
    }
}

unsigned int VtkOutputProcess::DetermineVtkCellListSize(ModelPart &model_part)
{
    unsigned int vtk_cell_list_size = 0;

    for (ModelPart::ElementIterator elem_i = model_part.ElementsBegin(); elem_i != model_part.ElementsEnd(); ++elem_i)
    {
        vtk_cell_list_size++;
        vtk_cell_list_size += elem_i->GetGeometry().size();
    }

    for (ModelPart::ConditionIterator condition_i = model_part.ConditionsBegin(); condition_i != model_part.ConditionsEnd(); ++condition_i)
    {
        vtk_cell_list_size++;
        vtk_cell_list_size += condition_i->GetGeometry().size();
    }

    return vtk_cell_list_size;
}

void VtkOutputProcess::Initialize(ModelPart &model_part)
{
    CreateMapFromKratosIdToVTKId(model_part);
    mVtkCellListSize = DetermineVtkCellListSize(model_part);
}

void VtkOutputProcess::WriteHeader(ModelPart &model_part)
{
    std::string outputFileName = GetOutputFileName(model_part);
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

void VtkOutputProcess::WriteMesh(ModelPart &model_part)
{
    WriteNodes(model_part);
    WriteConditionsAndElements(model_part);
    WriteConditionAndElementTypes(model_part);
}

void VtkOutputProcess::WriteNodes(ModelPart &model_part)
{
    std::string outputFileName = GetOutputFileName(model_part);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);
    outputFile << std::scientific;
    outputFile << std::setprecision(mDefaultPrecision);

    // write nodes header
    outputFile << "POINTS " << model_part.NumberOfNodes() << " float"
               << "\n";

    // write nodes
    for (ModelPart::NodeIterator node_i = model_part.NodesBegin(); node_i != model_part.NodesEnd(); ++node_i)
    {
        double x_coordinate = node_i->X();
        double y_coordinate = node_i->Y();
        double z_coordinate = node_i->Z();
        outputFile << " " << x_coordinate;
        outputFile << " " << y_coordinate;
        outputFile << " " << z_coordinate << "\n";
    }

    outputFile.close();
}

void VtkOutputProcess::WriteConditionsAndElements(ModelPart &model_part)
{
    std::string outputFileName = GetOutputFileName(model_part);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);

    // write cells header
    outputFile << "CELLS " << model_part.NumberOfConditions() + model_part.NumberOfElements() << " " << mVtkCellListSize << "\n";

    // write elements
    for (ModelPart::ElementIterator elem_i = model_part.ElementsBegin(); elem_i != model_part.ElementsEnd(); ++elem_i)
    {
        ModelPart::ConditionType::GeometryType &elem_geometry = elem_i->GetGeometry();
        const unsigned int numberOfNodes = elem_geometry.size();

        outputFile << numberOfNodes;
        for (unsigned int i = 0; i < numberOfNodes; i++)
            outputFile << " " << mKratosIdToVtkId[elem_geometry[i].Id()];

        outputFile << "\n";
    }

    // write Conditions
    for (ModelPart::ConditionIterator condition_i = model_part.ConditionsBegin(); condition_i != model_part.ConditionsEnd(); ++condition_i)
    {
        ModelPart::ConditionType::GeometryType &condition_geometry = condition_i->GetGeometry();
        const unsigned int numberOfNodes = condition_geometry.size();

        outputFile << numberOfNodes;
        for (unsigned int i = 0; i < numberOfNodes; i++)
            outputFile << " " << mKratosIdToVtkId[condition_geometry[i].Id()];
        outputFile << "\n";
    }

    outputFile.close();
}

void VtkOutputProcess::WriteConditionAndElementTypes(ModelPart &model_part)
{
    std::string outputFileName = GetOutputFileName(model_part);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);

    // write cell types header
    outputFile << "CELL_TYPES " << model_part.NumberOfConditions() + model_part.NumberOfElements() << "\n";

    // write elements types
    for (ModelPart::ElementIterator elem_i = model_part.ElementsBegin(); elem_i != model_part.ElementsEnd(); ++elem_i)
    {
        const unsigned int numberOfNodes = elem_i->GetGeometry().size();
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
    for (ModelPart::ConditionIterator condition_i = model_part.ConditionsBegin(); condition_i != model_part.ConditionsEnd(); ++condition_i)
    {
        const unsigned int numberOfNodes = condition_i->GetGeometry().size();
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
            KRATOS_THROW_ERROR(std::runtime_error, "Modelpart contains conditions with geometries for which no VTK-output is implemented!", "")

        outputFile << element_type << "\n";
    }

    outputFile.close();
}

void VtkOutputProcess::WriteNodalResultsAsPointData(ModelPart &model_part)
{
    std::string outputFileName = GetOutputFileName(model_part);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);
    // write nodal results header
    Parameters nodalResults = this->mrOutputSettings["nodal_solution_step_data_variables"];
    outputFile << "POINT_DATA " << model_part.NumberOfNodes() << "\n";

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
        for (ModelPart::NodeIterator node_i = model_part.NodesBegin(); node_i != model_part.NodesEnd(); ++node_i)
        {
            if (dataCharacteristic == 1)
            {
                Variable<double> nodalResultVariable = KratosComponents<Variable<double>>::Get(nodalResultName);
                double &nodalResult = node_i->FastGetSolutionStepValue(nodalResultVariable);
                outputFile << nodalResult << "\n";
            }
            else if (dataCharacteristic == 2)
            {
                Variable<array_1d<double, 3>> nodalResultVariable = KratosComponents<Variable<array_1d<double, 3>>>::Get(nodalResultName);
                array_1d<double, 3> &nodalResult = node_i->FastGetSolutionStepValue(nodalResultVariable);
                outputFile << nodalResult[0] << " ";
                outputFile << nodalResult[1] << " ";
                outputFile << nodalResult[2] << "\n";
            }
        }
    }

    outputFile.close();
}

void VtkOutputProcess::WriteElementData(ModelPart &model_part)
{
    std::string outputFileName = GetOutputFileName(model_part);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app);
    Parameters elementResults = this->mrOutputSettings["element_data_value_variables"];// list of element results

    // write cells header
    if (model_part.NumberOfElements() > 0)
    {
        outputFile << "CELL_DATA " << model_part.NumberOfElements() << "\n";
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
            for (ModelPart::ElementIterator elem_i = model_part.ElementsBegin(); elem_i != model_part.ElementsEnd(); ++elem_i)
            {
                if (dataCharacteristic == 1)
                {
                    Variable<double> elementResultVariable = KratosComponents<Variable<double>>::Get(elementResultName);
                    double &elementResult = elem_i->GetValue(elementResultVariable);
                    outputFile << elementResult << "\n";
                }
                else if (dataCharacteristic == 2)
                {
                    Variable<array_1d<double, 3>> elementResultVariable = KratosComponents<Variable<array_1d<double, 3>>>::Get(elementResultName);
                    array_1d<double, 3> &elementResult = elem_i->GetValue(elementResultVariable);
                    outputFile << elementResult[0] << " ";
                    outputFile << elementResult[1] << " ";
                    outputFile << elementResult[2] << "\n";
                }
            }
        }
        /*
        outputFile << "SCALARS SPLIT_ELEMENT float 1\nLOOKUP_TABLE default\n";

        // write element results for active
        for (ModelPart::ElementIterator elem_i = model_part.ElementsBegin(); elem_i != model_part.ElementsEnd(); ++elem_i)
        {
            //outputFile << numberOfNodes;
            bool is_split = elem_i->GetValue(SPLIT_ELEMENT);
            outputFile << is_split << "\n";
        }
    */
        outputFile.close();
    }
}

//#############################################For creating vtk files in binary format##########################################################

void VtkOutputProcess::WriteHeaderBinary(ModelPart &model_part)
{
    std::string outputFileName = GetOutputFileName(model_part);
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

void VtkOutputProcess::WriteMeshBinary(ModelPart &model_part)
{

    WriteNodesBinary(model_part);

    WriteConditionsAndElementsBinary(model_part);

    WriteConditionAndElementTypesBinary(model_part);
}

void VtkOutputProcess::WriteNodesBinary(ModelPart &model_part)
{
    std::string outputFileName = GetOutputFileName(model_part);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);

    // write nodes header
    outputFile << "\nPOINTS " << model_part.NumberOfNodes() << " float"
               << "\n";

    // write nodes
    for (ModelPart::NodeIterator node_i = model_part.NodesBegin(); node_i != model_part.NodesEnd(); ++node_i)
    {
        float x_coordinate = node_i->X();
        float y_coordinate = node_i->Y();
        float z_coordinate = node_i->Z();
        ForceBigEndian((unsigned char *)&x_coordinate);
        outputFile.write((char *)(&x_coordinate), sizeof(float));
        ForceBigEndian((unsigned char *)&y_coordinate);
        outputFile.write((char *)(&y_coordinate), sizeof(float));
        ForceBigEndian((unsigned char *)&z_coordinate);
        outputFile.write((char *)(&z_coordinate), sizeof(float));
    }

    outputFile.close();
}

void VtkOutputProcess::WriteConditionsAndElementsBinary(ModelPart &model_part)
{
    std::string outputFileName = GetOutputFileName(model_part);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);

    // write cells header
    outputFile << "\nCELLS " << model_part.NumberOfConditions() + model_part.NumberOfElements() << " " << mVtkCellListSize << "\n";

    // write elements
    for (ModelPart::ElementIterator elem_i = model_part.ElementsBegin(); elem_i != model_part.ElementsEnd(); ++elem_i)
    {

        ModelPart::ConditionType::GeometryType &elem_geometry = elem_i->GetGeometry();

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
    for (ModelPart::ConditionIterator condition_i = model_part.ConditionsBegin(); condition_i != model_part.ConditionsEnd(); ++condition_i)
    {
        ModelPart::ConditionType::GeometryType &condition_geometry = condition_i->GetGeometry();
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

void VtkOutputProcess::WriteConditionAndElementTypesBinary(ModelPart &model_part)
{
    std::string outputFileName = GetOutputFileName(model_part);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);

    // write cell types header
    outputFile << "\nCELL_TYPES " << model_part.NumberOfConditions() + model_part.NumberOfElements() << "\n";

    // write elements types
    for (ModelPart::ElementIterator elem_i = model_part.ElementsBegin(); elem_i != model_part.ElementsEnd(); ++elem_i)
    {
        const unsigned int numberOfNodes = elem_i->GetGeometry().size();
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

        ForceBigEndian((unsigned char *)&element_type);
        outputFile.write((char *)(&element_type), sizeof(int));
    }

    // write conditions types
    for (ModelPart::ConditionIterator condition_i = model_part.ConditionsBegin(); condition_i != model_part.ConditionsEnd(); ++condition_i)
    {
        const unsigned int numberOfNodes = condition_i->GetGeometry().size();
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
            KRATOS_THROW_ERROR(std::runtime_error, "Modelpart contains conditions with geometries for which no VTK-output is implemented!", "")

        ForceBigEndian((unsigned char *)&element_type);
        outputFile.write((char *)(&element_type), sizeof(int));
    }

    outputFile.close();
}

void VtkOutputProcess::WriteNodalResultsAsPointDataBinary(ModelPart &model_part)
{
    std::string outputFileName = GetOutputFileName(model_part);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);
    // write nodal results header
    Parameters nodalResults = this->mrOutputSettings["nodal_solution_step_data_variables"];
    outputFile << "\nPOINT_DATA " << model_part.NumberOfNodes() << "\n";

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

        for (ModelPart::NodeIterator node_i = model_part.NodesBegin(); node_i != model_part.NodesEnd(); ++node_i)
        {
            if (dataCharacteristic == 1)
            {
                Variable<double> nodalResultVariable = KratosComponents<Variable<double>>::Get(nodalResultName);
                float nodalResult = node_i->FastGetSolutionStepValue(nodalResultVariable);
                ForceBigEndian((unsigned char *)&nodalResult);
                outputFile.write((char *)(&nodalResult), sizeof(float));
            }
            else if (dataCharacteristic == 2)
            {
                Variable<array_1d<double, 3>> nodalResultVariable = KratosComponents<Variable<array_1d<double, 3>>>::Get(nodalResultName);
                array_1d<double, 3> nodalResult = node_i->FastGetSolutionStepValue(nodalResultVariable);
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

void VtkOutputProcess::WriteElementDataBinary(ModelPart &model_part)
{
    std::string outputFileName = GetOutputFileName(model_part);
    std::ofstream outputFile;
    outputFile.open(outputFileName, std::ios::out | std::ios::app);
    Parameters elementResults = this->mrOutputSettings["element_data_value_variables"];// list of element resultss
    if (model_part.NumberOfElements() > 0)
    {
        // write cells header
        outputFile << "\nCELL_DATA " << model_part.NumberOfElements() << "\n";

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

            for (ModelPart::ElementIterator elem_i = model_part.ElementsBegin(); elem_i != model_part.ElementsEnd(); ++elem_i)
            {
                if (dataCharacteristic == 1)
                {
                    Variable<double> elementResultVariable = KratosComponents<Variable<double>>::Get(elementResultName);
                    double elementResult = elem_i->GetValue(elementResultVariable);
                    ForceBigEndian((unsigned char *)&elementResult);
                    outputFile.write((char *)(&elementResult), sizeof(float));
                }
                else if (dataCharacteristic == 2)
                {
                    Variable<array_1d<double, 3>> elementResultVariable = KratosComponents<Variable<array_1d<double, 3>>>::Get(elementResultName);
                    array_1d<double, 3> elementResult = elem_i->GetValue(elementResultVariable);
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

void VtkOutputProcess::PrintOutputModelPart(ModelPart &modelPart)
{
    Initialize(modelPart);
    std::string type = this->mrOutputSettings["file_format"].GetString();
    if (type == "ASCII")
    {
        WriteHeader(modelPart);
        WriteMesh(modelPart);
        WriteNodalResultsAsPointData(modelPart);
        WriteElementData(modelPart);
    }
    else
    {
        Initialize(modelPart);
        WriteHeaderBinary(modelPart);
        WriteMeshBinary(modelPart);
        WriteNodalResultsAsPointDataBinary(modelPart);
        WriteElementDataBinary(modelPart);
    }
}

void VtkOutputProcess::PrintOutput()
{
    //For whole model part
    PrintOutputModelPart(mrModelPart);
    //For sub model parts
    bool print_sub_model_parts = this->mrOutputSettings["output_sub_model_parts"].GetBool();

    if(print_sub_model_parts)
    {
        std::vector<std::string> subModelPartNames = mrModelPart.GetSubModelPartNames();
        for (auto subModelPartName : subModelPartNames)
        {
            ModelPart &subModelPart = mrModelPart.GetSubModelPart(subModelPartName);
            PrintOutputModelPart(subModelPart);
        }
        ++mStep;
    }
}

void VtkOutputProcess::ForceBigEndian(unsigned char *bytes)
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

std::string VtkOutputProcess::GetOutputFileName(ModelPart &model_part)
{
    int rank = 0;
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    std::string outputFilename = mrOutputSettings["folder_name"].GetString() +"/"+
                                    model_part.Name() + "_" + std::to_string(rank) + "_" + std::to_string(mStep) + ".vtk";
    return outputFilename;
}

} // namespace Kratos
