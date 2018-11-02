// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//                   Navaneeth K Narayanan
//
//  Reference :      This class is adapted from applications/ShapeOptimizationapplication/custom_utilities/input_output/vtk_file_io.h
// ==============================================================================

#if !defined(VTK_OUTPUT_H)
#define VTK_OUTPUT_H
// System includes
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip> // for std::setprecision
#include <map>
// For MPI-parallel output
#ifdef KRATOS_USING_MPI
#include "mpi.h"
#endif
#include "includes/deprecated_variables.h"
#include "includes/kratos_parameters.h"
// project includes

namespace Kratos
{
/** \brief Quaternion
	* A simple class that has functionality to write vtk output
	*/
class VtkOutput
{

    typedef ProcessInfo ProcessInfoType;

  public:
    /// Pointer definition of VtkOutput
    KRATOS_CLASS_POINTER_DEFINITION(VtkOutput);

    ///@name Life Cycle
    ///@{

    /**
		Creates a VtkOutput data object
		*/
    VtkOutput(ModelPart &model_part, std::string outPutFileName, Parameters rParameters) : mr_model_part(model_part), mrOutputSettings(rParameters)
    {
        mcaseName = outPutFileName;
        mDefaultPrecision = 7;
        step = 0;
        this->mDoneTest = false;
        this->mShouldSwap = false;
    }
    /// Destructor.
    virtual ~VtkOutput(){};

    ///@}
    ///@name Operators
    ///@{
    std::map<int, int> createMapFromKratosIdToVTKId(ModelPart &model_part)
    {
        std::map<int, int> kratos_id_to_vtk;
        int vtk_id = 0;

        for (ModelPart::NodeIterator node_i = model_part.NodesBegin(); node_i != model_part.NodesEnd(); ++node_i)
        {
            int KratosId = node_i->Id();
            kratos_id_to_vtk[KratosId] = vtk_id;
            vtk_id++;
        }

        return kratos_id_to_vtk;
    }
    std::size_t determineVtkCellListSize(ModelPart &model_part)
    {
        std::size_t vtk_cell_list_size = 0;

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

    void initialize(ModelPart &model_part)
    {
        mKratosIdToVtkId = createMapFromKratosIdToVTKId(model_part);
        mVtkCellListSize = determineVtkCellListSize(model_part);
    }

    void writeHeader(ModelPart &model_part)
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

    void writeMesh(ModelPart &model_part)
    {
        writeNodes(model_part);
        writeConditionsAndElements(model_part);
        writeConditionAndElementTypes(model_part);
    }

    void writeNodes(ModelPart &model_part)
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

    void writeConditionsAndElements(ModelPart &model_part)
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
            const std::size_t numberOfNodes = elem_geometry.size();

            outputFile << numberOfNodes;
            for (std::size_t i = 0; i < numberOfNodes; i++)
                outputFile << " " << mKratosIdToVtkId[elem_geometry[i].Id()];

            outputFile << "\n";
        }

        // write Conditions
        for (ModelPart::ConditionIterator condition_i = model_part.ConditionsBegin(); condition_i != model_part.ConditionsEnd(); ++condition_i)
        {
            ModelPart::ConditionType::GeometryType &condition_geometry = condition_i->GetGeometry();
            const std::size_t numberOfNodes = condition_geometry.size();

            outputFile << numberOfNodes;
            for (std::size_t i = 0; i < numberOfNodes; i++)
                outputFile << " " << mKratosIdToVtkId[condition_geometry[i].Id()];
            outputFile << "\n";
        }

        outputFile.close();
    }

    void writeConditionAndElementTypes(ModelPart &model_part)
    {
        std::string outputFileName = GetOutputFileName(model_part);
        std::ofstream outputFile;
        outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);

        // write cell types header
        outputFile << "CELL_TYPES " << model_part.NumberOfConditions() + model_part.NumberOfElements() << "\n";

        // write elements types
        for (ModelPart::ElementIterator elem_i = model_part.ElementsBegin(); elem_i != model_part.ElementsEnd(); ++elem_i)
        {
            const std::size_t numberOfNodes = elem_i->GetGeometry().size();
            std::size_t element_type;

            if (numberOfNodes == 3)
                element_type = 5;
            else if (numberOfNodes == 4)
            {
                if (mr_model_part.GetProcessInfo()[DOMAIN_SIZE] == 2)
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
            const std::size_t numberOfNodes = condition_i->GetGeometry().size();
            std::size_t element_type;

            if (numberOfNodes == 3)
                element_type = 5;
            else if (numberOfNodes == 4)
            {
                if (mr_model_part.GetProcessInfo()[DOMAIN_SIZE] == 2)
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

    void writeNodalResultsAsPointData(ModelPart &model_part)
    {
        std::string outputFileName = GetOutputFileName(model_part);
        std::ofstream outputFile;
        outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);
        // write nodal results header
        Parameters nodalResults = this->mrOutputSettings["result_file_configuration"]["nodal_results"];
        outputFile << "POINT_DATA " << model_part.NumberOfNodes() << "\n";

        for (std::size_t entry = 0; entry < nodalResults.size(); entry++)
        {
            // write nodal results variable header
            std::string nodalResultName = nodalResults[entry].GetString();
            std::size_t dataCharacteristic = 0; // 0: unknown, 1: Scalar value, 2: 3 DOF global translation vector
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

    void writeElementData(ModelPart &model_part)
    {
        std::string outputFileName = GetOutputFileName(model_part);
        std::ofstream outputFile;
        outputFile.open(outputFileName, std::ios::out | std::ios::app);
        std::vector<std::string> elementResults = {"ELEMENTAL_DISTANCES"}; // list of element results
        // write cells header
        if (model_part.NumberOfElements() > 0)
        {
            outputFile << "CELL_DATA " << model_part.NumberOfElements() << "\n";
            outputFile << "SCALARS ACTIVE float 1\nLOOKUP_TABLE default\n";

            // write element results for active
            for (ModelPart::ElementIterator elem_i = model_part.ElementsBegin(); elem_i != model_part.ElementsEnd(); ++elem_i)
            {
                //outputFile << numberOfNodes;
                if ((elem_i)->IsDefined(ACTIVE))
                {

                    outputFile << elem_i->Is(ACTIVE) << "\n";
                }

                else
                    outputFile << "1\n";
            }

            for (std::size_t entry = 0; entry < elementResults.size(); entry++)
            {

                std::string elementResultName = elementResults[entry];
                std::size_t dataCharacteristic = 0; // 0: unknown, 1: Scalar value, 2: 3 DOF global translation vector

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

        outputFile << "SCALARS SPLIT_ELEMENT float 1\nLOOKUP_TABLE default\n";

        // write element results for active
        for (ModelPart::ElementIterator elem_i = model_part.ElementsBegin(); elem_i != model_part.ElementsEnd(); ++elem_i)
        {
            //outputFile << numberOfNodes;
            bool is_split = elem_i->GetValue(SPLIT_ELEMENT);
            outputFile << is_split << "\n";
        }



            outputFile.close();
        }
    }

    //#############################################For creating vtk files in binary format##########################################################

    void writeHeaderBinary(ModelPart &model_part)
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

    void writeMeshBinary(ModelPart &model_part)
    {

        writeNodesBinary(model_part);

        writeConditionsAndElementsBinary(model_part);

        writeConditionAndElementTypesBinary(model_part);
    }

    void writeNodesBinary(ModelPart &model_part)
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
            force_big_endian((unsigned char *)&x_coordinate);
            outputFile.write((char *)(&x_coordinate), sizeof(float));
            force_big_endian((unsigned char *)&y_coordinate);
            outputFile.write((char *)(&y_coordinate), sizeof(float));
            force_big_endian((unsigned char *)&z_coordinate);
            outputFile.write((char *)(&z_coordinate), sizeof(float));
        }

        outputFile.close();
    }

    void writeConditionsAndElementsBinary(ModelPart &model_part)
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

            std::size_t numberOfNodes = elem_geometry.size();

            force_big_endian((unsigned char *)&numberOfNodes);

            outputFile.write((char *)(&numberOfNodes), sizeof(std::size_t));

            for (std::size_t i = 0; i < elem_geometry.size(); i++)
            {
                int nodenum = mKratosIdToVtkId[elem_geometry[i].Id()];
                force_big_endian((unsigned char *)&nodenum);
                outputFile.write((char *)(&nodenum), sizeof(int));
            }
        }

        // write Conditions
        for (ModelPart::ConditionIterator condition_i = model_part.ConditionsBegin(); condition_i != model_part.ConditionsEnd(); ++condition_i)
        {
            ModelPart::ConditionType::GeometryType &condition_geometry = condition_i->GetGeometry();
            std::size_t numberOfNodes = condition_geometry.size();

            force_big_endian((unsigned char *)&numberOfNodes);
            outputFile.write((char *)(&numberOfNodes), sizeof(std::size_t));

            for (std::size_t i = 0; i < condition_geometry.size(); i++)
            {

                int nodenum = mKratosIdToVtkId[condition_geometry[i].Id()];
                force_big_endian((unsigned char *)&nodenum);
                outputFile.write((char *)(&nodenum), sizeof(int));
            }
        }

        outputFile.close();
    }

    void writeConditionAndElementTypesBinary(ModelPart &model_part)
    {
        std::string outputFileName = GetOutputFileName(model_part);
        std::ofstream outputFile;
        outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);

        // write cell types header
        outputFile << "\nCELL_TYPES " << model_part.NumberOfConditions() + model_part.NumberOfElements() << "\n";

        // write elements types
        for (ModelPart::ElementIterator elem_i = model_part.ElementsBegin(); elem_i != model_part.ElementsEnd(); ++elem_i)
        {
            const std::size_t numberOfNodes = elem_i->GetGeometry().size();
            std::size_t element_type;

            if (numberOfNodes == 3)
                element_type = 5;
            else if (numberOfNodes == 4)
            {
                if (mr_model_part.GetProcessInfo()[DOMAIN_SIZE] == 2)
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

            force_big_endian((unsigned char *)&element_type);
            outputFile.write((char *)(&element_type), sizeof(int));
        }

        // write conditions types
        for (ModelPart::ConditionIterator condition_i = model_part.ConditionsBegin(); condition_i != model_part.ConditionsEnd(); ++condition_i)
        {
            const std::size_t numberOfNodes = condition_i->GetGeometry().size();
            std::size_t element_type;

            if (numberOfNodes == 3)
                element_type = 5;
            else if (numberOfNodes == 4)
            {
                if (mr_model_part.GetProcessInfo()[DOMAIN_SIZE] == 2)
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

            force_big_endian((unsigned char *)&element_type);
            outputFile.write((char *)(&element_type), sizeof(int));
        }

        outputFile.close();
    }

    void writeNodalResultsAsPointDataBinary(ModelPart &model_part)
    {
        std::string outputFileName = GetOutputFileName(model_part);
        std::ofstream outputFile;
        outputFile.open(outputFileName, std::ios::out | std::ios::app | std::ios::binary);
        // write nodal results header
        Parameters nodalResults = this->mrOutputSettings["result_file_configuration"]["nodal_results"];
        outputFile << "\nPOINT_DATA " << model_part.NumberOfNodes() << "\n";

        for (std::size_t entry = 0; entry < nodalResults.size(); entry++)
        {
            // write nodal results variable header
            std::string nodalResultName = nodalResults[entry].GetString();
            std::size_t dataCharacteristic = 0; // 0: unknown, 1: Scalar value, 2: 3 DOF global translation vector
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
                    force_big_endian((unsigned char *)&nodalResult);
                    outputFile.write((char *)(&nodalResult), sizeof(float));
                }
                else if (dataCharacteristic == 2)
                {
                    Variable<array_1d<double, 3>> nodalResultVariable = KratosComponents<Variable<array_1d<double, 3>>>::Get(nodalResultName);
                    array_1d<double, 3> nodalResult = node_i->FastGetSolutionStepValue(nodalResultVariable);
                    float num1 = nodalResult[0];
                    force_big_endian((unsigned char *)&num1);
                    outputFile.write((char *)(&num1), sizeof(float));
                    float num2 = nodalResult[1];
                    force_big_endian((unsigned char *)&num2);
                    outputFile.write((char *)(&num2), sizeof(float));
                    float num3 = nodalResult[2];
                    force_big_endian((unsigned char *)&num3);
                    outputFile.write((char *)(&num3), sizeof(float));
                }
            }
        }

        outputFile.close();
    }

    void writeElementDataBinary(ModelPart &model_part)
    {
        std::string outputFileName = GetOutputFileName(model_part);
        std::ofstream outputFile;
        outputFile.open(outputFileName, std::ios::out | std::ios::app);
        std::vector<std::string> elementResults = {}; //list of element results
        if (model_part.NumberOfElements() > 0)
        {
            // write cells header
            outputFile << "\nCELL_DATA " << model_part.NumberOfElements() << "\n";
            outputFile << "SCALARS ACTIVE float \nLOOKUP_TABLE default\n";

            // write element results for active
            for (ModelPart::ElementIterator elem_i = model_part.ElementsBegin(); elem_i != model_part.ElementsEnd(); ++elem_i)
            {
                //outputFile << numberOfNodes;

                if ((elem_i)->IsDefined(ACTIVE))
                {
                    float is_active = elem_i->Is(ACTIVE);
                    force_big_endian((unsigned char *)&is_active);
                    outputFile.write((char *)(&is_active), sizeof(float));
                }

                else
                {

                    float is_active = 1;
                    force_big_endian((unsigned char *)&is_active);
                    outputFile.write((char *)(&is_active), sizeof(float));
                }
            }

            for (std::size_t entry = 0; entry < elementResults.size(); entry++)
            {

                std::string elementResultName = elementResults[entry];
                std::size_t dataCharacteristic = 0; // 0: unknown, 1: Scalar value, 2: 3 DOF global translation vector

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
                        force_big_endian((unsigned char *)&elementResult);
                        outputFile.write((char *)(&elementResult), sizeof(float));
                    }
                    else if (dataCharacteristic == 2)
                    {
                        Variable<array_1d<double, 3>> elementResultVariable = KratosComponents<Variable<array_1d<double, 3>>>::Get(elementResultName);
                        array_1d<double, 3> elementResult = elem_i->GetValue(elementResultVariable);
                        float num1 = elementResult[0];
                        force_big_endian((unsigned char *)&num1);
                        outputFile.write((char *)(&num1), sizeof(float));
                        float num2 = elementResult[1];
                        force_big_endian((unsigned char *)&num2);
                        outputFile.write((char *)(&num2), sizeof(float));
                        float num3 = elementResult[2];
                        force_big_endian((unsigned char *)&num3);
                        outputFile.write((char *)(&num3), sizeof(float));
                    }
                }
            }

        outputFile << "SCALARS SPLIT_ELEMENT float \nLOOKUP_TABLE default\n";

        // write element results for active
        for (ModelPart::ElementIterator elem_i = model_part.ElementsBegin(); elem_i != model_part.ElementsEnd(); ++elem_i)
        {
            //outputFile << numberOfNodes;
            float is_split = elem_i->GetValue(SPLIT_ELEMENT);
            force_big_endian((unsigned char *)&is_split);
            outputFile.write((char *)(&is_split), sizeof(float));
        }

            outputFile.close();
        }
    }

    //#################################################################End of Binary vtk ################################################################

    void PrintOutputSubModelPart(ModelPart &modelPart)
    {
        initialize(modelPart);
        std::string type = this->mrOutputSettings["result_file_configuration"]["gidpost_flags"]["GiDPostMode"].GetString();
        if (type == "GiD_PostBinary")
        {
            initialize(modelPart);
            writeHeaderBinary(modelPart);
            writeMeshBinary(modelPart);
            writeNodalResultsAsPointDataBinary(modelPart);
            writeElementDataBinary(modelPart);
        }

        else
        {

            writeHeader(modelPart);
            writeMesh(modelPart);
            writeNodalResultsAsPointData(modelPart);
            writeElementData(modelPart);
        }
    }

    void PrintOutput()
    {
        //For whole model part
        PrintOutputSubModelPart(mr_model_part);
        //For sub model parts
        std::vector<std::string> subModelPartNames = mr_model_part.GetSubModelPartNames();
        for (auto subModelPartName : subModelPartNames)
        {
            ModelPart &subModelPart = mr_model_part.GetSubModelPart(subModelPartName);
            PrintOutputSubModelPart(subModelPart);
        }
        ++step;
    }

    void force_big_endian(unsigned char *bytes)
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

    ///@name Access
    ///@{

    ///@

    ///@name Static Operations
    ///@{
    /**
		 * Returns the string containing a detailed description of this object.
		 * @return the string with informations
		 */
    virtual void GetInfo() const
    {
    }

    ///@}
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << " VtkOutput object " << std::endl;
    }

  private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{
    ModelPart &mr_model_part;
    std::string mcaseName;
    std::string mOutputFilename;

    Parameters mrOutputSettings;
    std::size_t mDefaultPrecision;
    std::map<int, int> mKratosIdToVtkId;
    std::size_t mVtkCellListSize;
    std::size_t step;
    bool mDoneTest;
    bool mShouldSwap;

    ///@}

    std::string GetOutputFileName(ModelPart &model_part)
    {
        int rank = 0;
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
        std::string outputFilename = model_part.Name() + "_" + std::to_string(rank) + "_" + std::to_string(step) + ".vtk";
        return outputFilename;
    }

    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer &rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, VtkOutput);
    }

    virtual void load(Serializer &rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, VtkOutput);
    }

    ///@}
};

///@name Input/Output funcitons
///@{

inline std::istream &operator>>(std::istream &rIStream, VtkOutput &rThis);

inline std::ostream &operator<<(std::ostream &rOStream, const VtkOutput &rThis)
{
    return rOStream;
}

///@}

} // namespace Kratos

#endif // CONSTRAINT_SLAVE_H_INCLUDED
