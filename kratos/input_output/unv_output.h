//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//
//

#ifndef KRATOS_UNV_OUTPUT_H_INCLUDED
#define KRATOS_UNV_OUTPUT_H_INCLUDED

/* System includes */
#include <iostream>
#include <fstream>
#include <string>

/* External includes */

/* Project includes */
#include "includes/model_part.h"
#include "includes/exception.h"

namespace Kratos {


/**
 * @brief Provides a tool to write UNV files.
 * 
 * Currently 3 datasets are supported: 
 * 2411 - Node    Dataset
 * 2412 - Element Dataset
 * 2414 - Result  Dataset
 * 
 */
class KRATOS_API(KRATOS_CORE) UnvOutput {
public:

    enum class DatasetID {
        NODES_DATASET               = 2411,
        ELEMENTS_DATASET            = 2412,
        RESULTS_DATASET             = 2414
    };
    
    enum class DatasetLocation {
        DATA_AT_NODES               = 1,
        DATA_AT_ELEMENTS            = 2,
        DATA_AT_NODES_ON_ELEMENTS   = 3,
        DATA_AT_POINTS              = 5,
        DATA_ON_ELEMENTS_AT_NODES   = 6
    };

    enum class ModelType {
        UNKNOWN         = 0,
        STRUCTURAL      = 1,
        HEAT_TRANSFER   = 2,
        FLUID_FLOW      = 3
    };

    enum class AnalysisType {
        UNKNOWN                             = 0,
        STATIC                              = 1,
        NORMAL_MODE                         = 2,
        COMPLEX_EIGENVALUE_FIRST_ORDER      = 3,
        TRANSIENT                           = 4,
        FREQUENCY_RESPONSE                  = 5,
        BUCKLING                            = 6,
        COMPLEX_EIGENVALUE_SECOND_ORDER     = 7,
        STATIC_NON_LINEAR                   = 9,
        CRAIG_BAMPTON_CONSTRAINT_MODES      = 10,
        EQUIVALENT_ATTACHMENT_MODES         = 11,
        EFFECTIVE_MASS_MODES                = 12,
        EFFECTIVE_MASS_MATRIX               = 13,
        EFFECTIVE_MASS_MATRIX_COPY          = 14,   // This record is duplicared intentionally
        DISTRIBUTED_LOAD_LOAD_DISTRIBUTION  = 15,
        DISTRIBUTED_LOAD_ATTACHMENT_MODES   = 16
    };

    // 3_DOF_* Records have their name changed to D3_DOF to avoid conflicts.
    enum class DataCharacteristics {
        UNKNOWN = 0,
        SCALAR = 1,
        D3_DOF_GLOBAL_TRANSLATION_VECTOR = 2,
        D3_DOF_GLOBAL_TRANSLATION_ROTATION_VECTOR = 3,
        SYMMETRIC_GLOBAL_TENSOR = 4,
        STRESS_RESULTANTS = 6
    };

    // This record is to big, I will consider moving these enums to another file.
    // enum class ResultType {

    // };

    enum class DataType {
        INTEGER = 1,
        SINGLE_PRECISION_FLOATING_POINT = 2,
        DOUBLE_PRECISION_FLOATING_POINT = 4,
        SINGLE_PRECISION_COMPLEX = 5,
        DOUBLE_PRECISION_COMPLEX = 6
    };

    KRATOS_CLASS_POINTER_DEFINITION(UnvOutput);

    template <typename Enumeration>
    auto as_integer(Enumeration const value)
        -> typename std::underlying_type<Enumeration>::type
    {
        return static_cast<typename std::underlying_type<Enumeration>::type>(value);
    }

    UnvOutput(Kratos::ModelPart &modelPart, const std::string &outFileWithoutExtension);


    void InitializeOutputFile();

    /**
     * @brief Writes 'mrOutputModelPart' associated mesh.
     * 
     */
    void WriteMesh();

    /**
     * @brief Writes 'mrOutputModelPart' associated nodes.
     * 
     */
    void WriteNodes();

    /**
     * @brief Writes 'mrOutputModelPart' associated conditions.
     * 
     */
    void WriteElements();

    /**
     * @brief Writes a result dataset containing the rVariable value for a given timestep
     * 
     * @param rVariable Kratos Variable to be printed
     * @param timeStep  TimeStep.
     */
    void WriteNodalResults(const Variable<bool>& rVariable, const double timeStep);
    void WriteNodalResults(const Variable<int>& rVariable, const double timeStep);
    void WriteNodalResults(const Variable<double>& rVariable, const double timeStep);
    void WriteNodalResults(const Variable<array_1d<double,3>>& rVariable, const double timeStep);
    void WriteNodalResults(const Variable<Vector>& rVariable, const double timeStep);
    void WriteNodalResults(const Variable<Matrix>& rVariable, const double timeStep);

    /**
     * @brief Returns the type of unv data associated to a Kratos Variable
     * 
     * Vectors and Matrices are not supported at this time.
     * 
     * @return UnvOutput::DataCharacteristics If of the unv data type
     */
    UnvOutput::DataCharacteristics GetDataType(const Variable<bool>&);
    UnvOutput::DataCharacteristics GetDataType(const Variable<int>&);
    UnvOutput::DataCharacteristics GetDataType(const Variable<double>&);
    UnvOutput::DataCharacteristics GetDataType(const Variable<array_1d<double,3>>&);
    UnvOutput::DataCharacteristics GetDataType(const Variable<Vector>&);
    UnvOutput::DataCharacteristics GetDataType(const Variable<Matrix>&);

    /**
     * @brief Writes the variable value for a node.
     * 
     * Vectors and Matrices are not supported at this time.
     * 
     * @param outputFile Output file 
     * @param node       Input node
     * @param rVariable  Variable to print
     */
    void WriteNodalResultValues(std::ofstream &outputFile, const Node& node, const Variable<bool>& rVariable);
    void WriteNodalResultValues(std::ofstream &outputFile, const Node& node, const Variable<int>& rVariable);
    void WriteNodalResultValues(std::ofstream &outputFile, const Node& node, const Variable<double>& rVariable);
    void WriteNodalResultValues(std::ofstream &outputFile, const Node& node, const Variable<array_1d<double,3>>& rVariable);
    void WriteNodalResultValues(std::ofstream &outputFile, const Node& node, const Variable<Vector>& rVariable);
    void WriteNodalResultValues(std::ofstream &outputFile, const Node& node, const Variable<Matrix>& rVariable);

    /**
     * @brief Get the id of the UNV variable name corresponding to rVariable. 1000+ if none found.
     * 
     * @param rVariable Kratos Variable
     * @return int      Id of the unv variable corresponding to rVariable
     */
    template<class TVariablebleType>
    int GetUnvVariableName(const TVariablebleType& rVariable) {
        if(rVariable == VELOCITY)       return 11;
        if(rVariable == TEMPERATURE)    return 5;
        if(rVariable == PRESSURE)       return 117;

        return 1000;
    }

    /**
     * @brief Writes a result dataset using the results in node mode
     * 
     * Fromat: Partially extracted from: http://users.ices.utexas.edu
     * R.  1: unique number of dataset (dataset_label)
     * R.  2: text describing content (dataset_name)
     * R.  3: data belongs to: nodes, elements,...
     *        (dataset_location)
     * R.  4: user-specified text (id_lines_1_to_5[0])
     * R.  5: user-specified text (id_lines_1_to_5[1])
     * R.  6: user-specified text (id_lines_1_to_5[2])
     * R.  7: user-specified text (id_lines_1_to_5[3])
     * R.  8: user-specified text (id_lines_1_to_5[4])
     * R.  9: (model_type) (analysis_type) (data_characteristic) (result_type) (data_type) (nvaldc)
     * R. 10: (design_set_id) (iteration_number) (solution_set_id) (boundary_condition) (load_set) (mode_number) (time_stamp_number) (frequency_number)
     * R. 11: (creation_option) (Unknown)*7
     * R. 12: (time) (frequency) (eigenvalue) (nodal_mass) (viscous_damping_ratio) (hysteretic_damping_ratio)
     * R. 13: (eigenvalue_re) (eigenvalue_im) (modalA_re) (modalA_im) (modalB_re) (modalB_im)
     * 
     * For nodes (Repeat for every node):
     * 
     * R. 14: (node_id)
     * R. 15: (result)*nvaldc
     * 
     * @param rVariable     Variable to be printed 
     * @param numComponents Number of components of the variable 
     * @param timeStep      Current TimeStep
     */
    template<class TVariablebleType>
    void WriteNodalResultRecords(const TVariablebleType& rVariable, const int numComponents, const double timeStep) {
        std::ofstream outputFile;
        outputFile.open(mOutputFileName, std::ios::out | std::ios::app);

        std::string dataSetName = "NodalResults";
        const std::string& dataSetLabel = rVariable.Name();

        outputFile << std::setw(6)  << "-1" << "\n";                                                // Begin block
        outputFile << std::setw(6)  << as_integer(DatasetID::RESULTS_DATASET) << "\n";              // DatasetID

        outputFile << std::setw(10) << dataSetLabel << "\n";                                        // Record 1 - Label
        outputFile << std::setw(6)  << dataSetName << "\n";                                         // Record 2 - Name
        outputFile << std::setw(10) << as_integer(DatasetLocation::DATA_AT_NODES) << "\n";          // Record 3

        // String records, seems like you can put anything you want.
        outputFile << "" << "\n";                                                                   // Record 4
        outputFile << "" << "\n";                                                                   // Record 5
        outputFile << "" << "\n";                                                                   // Record 6
        outputFile << "" << "\n";                                                                   // Record 7
        outputFile << "" << "\n";                                                                   // Record 8
        
        // ModelType, AnalysisType, DataCharacteristic, ResultType, DataType, NumberOfDataValues    // Record 9
        outputFile << std::setw(10) << as_integer(ModelType::STRUCTURAL); 
        outputFile << std::setw(10) << as_integer(AnalysisType::TRANSIENT);
        outputFile << std::setw(10) << as_integer(GetDataType(rVariable));
        outputFile << std::setw(10) << GetUnvVariableName(rVariable);
        outputFile << std::setw(10) << as_integer(DataType::SINGLE_PRECISION_FLOATING_POINT);
        outputFile << std::setw(10) << numComponents; 
        outputFile << "\n";

        // DesignSetId, IterationNumber, SolutionSetId, BoundaryCondition, LoadSet, ModeNumber, TimeStampNumber, FrequencyNumber
        outputFile << std::setw(10) << 0;                                                           // Record 10
        outputFile << std::setw(10) << timeStep;
        outputFile << std::setw(10) << 0;
        outputFile << std::setw(10) << 0;
        outputFile << std::setw(10) << 0;
        outputFile << std::setw(10) << 1;
        outputFile << std::setw(10) << timeStep;
        outputFile << std::setw(10) << 0;
        outputFile << "\n";

        // CreationOption, (Unknown)*7
        outputFile << std::setw(10) << 0;                                                           // Record 11
        outputFile << std::setw(10) << 0;
        outputFile << "\n";

        // Time, Frequency Eigenvalue NodalMass ViscousDampingRatio, HystereticDampingRatio
        outputFile << std::setw(13) << timeStep * 0.1;                                              // Record 12
        outputFile << std::setw(13) << "0.00000E+00";
        outputFile << std::setw(13) << "0.00000E+00";
        outputFile << std::setw(13) << "0.00000E+00";
        outputFile << std::setw(13) << "0.00000E+00";
        outputFile << std::setw(13) << "0.00000E+00";
        outputFile << "\n";

        // Eigenvalue_re, Eigenvalue_im, ModalA_re, ModalA_im, ModalB_re, ModalB_im
        outputFile << std::setw(13) << "0.00000E+00";                                               // Record 13
        outputFile << std::setw(13) << "0.00000E+00";
        outputFile << std::setw(13) << "0.00000E+00";
        outputFile << std::setw(13) << "0.00000E+00";
        outputFile << std::setw(13) << "0.00000E+00";
        outputFile << std::setw(13) << "0.00000E+00";
        outputFile << "\n";

        // Data at nodes:
        for (auto &node_i : mrOutputModelPart.Nodes()) {
            int node_label = node_i.Id();
            outputFile << std::setw(6) << node_label << "\n";                                       // Record 14 - Node Number
            WriteNodalResultValues(outputFile, node_i, rVariable);                                  // Record 15 - NumberOfDataValues' data of the node
        }
        
        outputFile << std::setw(6) << "-1" << "\n";
        outputFile.close();
    }

private:

    Kratos::ModelPart &mrOutputModelPart;
    std::string mOutputFileName;
    
};
}

#endif //KRATOS_UNV_OUTPUT_H_INCLUDED
