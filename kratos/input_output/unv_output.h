//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"

namespace Kratos 
{

///@name Kratos Classes
///@{

/**
 * @ingroup KratosCore
 * @class UnvOutput
 * @brief Provides a tool to write UNV files.
 * @details This class provides a tool to write UNV files, which are used for data exchange between different finite element analysis software. The UNV format is a widely used standard for representing finite element models and results. The UnvOutput class allows users to write mesh data (nodes, elements, and conditions) as well as nodal results (variables) to a UNV file. It supports various types of variables, including scalar, vector, and matrix types, and provides methods to specify the dataset location and characteristics of the data being written. The class is designed to be flexible and extensible, allowing for easy integration into different simulation workflows.
 * Currently 3 datasets are supported: 
 * 2411 - Node    Dataset
 * 2412 - Element Dataset
 * 2414 - Result  Dataset
 * @author Carlos Roig
 */
class KRATOS_API(KRATOS_CORE) UnvOutput 
{
public:
    ///@name Type Definitions
    ///@{

    /**
     * @brief Pointer definition of UnvOutput
     * @details This defines a pointer to UnvOutput
     * - NODES_DATASET: If a dataset is located at nodes, the data will be printed for each node in the dataset.
     * - ELEMENTS_DATASET: If a dataset is located at elements, the data will be printed for each element in the dataset.
     * - RESULTS_DATASET: If a dataset is located at results, the data will be printed for each result in the dataset.
     */
    enum class DatasetID {
        NODES_DATASET               = 2411,
        ELEMENTS_DATASET            = 2412,
        RESULTS_DATASET             = 2414
    };
    
    /**
     * @brief Constructor of the class
     * @details This constructor initializes the UnvOutput class with a reference to a ModelPart and an output file name.
     * - DATA_AT_NODES: If a dataset is located at nodes, the data will be printed for each node in the dataset.
     * - DATA_AT_ELEMENTS: If a dataset is located at elements, the data will be printed for each element in the dataset.
     * - DATA_AT_NODES_ON_ELEMENTS: If a dataset is located at nodes on elements, the data will be printed for each node on each element in the dataset.
     * - DATA_AT_POINTS: If a dataset is located at points, the data will be printed for each point in the dataset.
     * - DATA_ON_ELEMENTS_AT_NODES: If a dataset is located on elements at nodes, the data will be printed for each element at each node in the dataset.
     */
    enum class DatasetLocation {
        DATA_AT_NODES               = 1,
        DATA_AT_ELEMENTS            = 2,
        DATA_AT_NODES_ON_ELEMENTS   = 3,
        DATA_AT_POINTS              = 5,
        DATA_ON_ELEMENTS_AT_NODES   = 6
    };

    /**
     * @brief Destructor of the class
     * @details This destructor cleans up any resources used by the UnvOutput class.
     * UNKNOWN: If the dataset location is unknown, the data will not be printed.
     * STRUCTURAL: If the dataset location is structural, the data will be printed for each structural element in the dataset.
     * HEAT_TRANSFER: If the dataset location is heat transfer, the data will be printed for each heat transfer element in the dataset.
     * FLUID_FLOW: If the dataset location is fluid flow, the data will be printed for each fluid flow element in the dataset.
     */
    enum class ModelType {
        UNKNOWN         = 0,
        STRUCTURAL      = 1,
        HEAT_TRANSFER   = 2,
        FLUID_FLOW      = 3
    };

    /**
     * @brief Enumeration for different types of analysis in UNV output.
     * @details This enumeration defines various types of analysis that can be specified in the UNV output, such as static, transient, and modal analyses.
     * - UNKNOWN: If the analysis type is unknown, the data will not be printed.
     * - STATIC: If the analysis type is static, the data will be printed for each static analysis in the dataset.
     * - NORMAL_MODE: If the analysis type is normal mode, the data will be printed for each normal mode analysis in the dataset.
     * - COMPLEX_EIGENVALUE_FIRST_ORDER: If the analysis type is complex eigenvalue first order, the data will be printed for each complex eigenvalue first order analysis in the dataset.
     * - TRANSIENT: If the analysis type is transient, the data will be printed for each transient analysis in the dataset.
     * - FREQUENCY_RESPONSE: If the analysis type is frequency response, the data will be printed for each frequency response analysis in the dataset.
     * - BUCKLING: If the analysis type is buckling, the data will be printed for each buckling analysis in the dataset.
     * - COMPLEX_EIGENVALUE_SECOND_ORDER: If the analysis type is complex eigenvalue second order, the data will be printed for each complex eigenvalue second order analysis in the dataset.
     * - STATIC_NON_LINEAR: If the analysis type is static non-linear, the data will be printed for each static non-linear analysis in the dataset.
     * - CRAIG_BAMPTON_CONSTRAINT_MODES: If the analysis type is Craig-Bampton constraint modes, the data will be printed for each Craig-Bampton constraint mode analysis in the dataset.
     * - EQUIVALENT_ATTACHMENT_MODES: If the analysis type is equivalent attachment modes, the data will be printed for each equivalent attachment mode analysis in the dataset.
     * - EFFECTIVE_MASS_MODES: If the analysis type is effective mass modes, the data will be printed for each effective mass mode analysis in the dataset.
     * - EFFECTIVE_MASS_MATRIX: If the analysis type is effective mass matrix, the data will be printed for each effective mass matrix analysis in the dataset.
     * - EFFECTIVE_MASS_MATRIX_COPY: If the analysis type is effective mass matrix copy, the data will be printed for each effective mass matrix copy analysis in the dataset.
     * - DISTRIBUTED_LOAD_LOAD_DISTRIBUTION: If the analysis type is distributed load load distribution, the data will be printed for each distributed load load distribution analysis in the dataset.
     * - DISTRIBUTED_LOAD_ATTACHMENT_MODES: If the analysis type is distributed load attachment modes, the data will be printed for each distributed load attachment mode analysis in the dataset.
     */
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

    /**
     * @brief Enumeration for different characteristics of data in UNV output.
     * @details This enumeration defines various characteristics of data that can be specified in the UNV output, such as scalar, vector, and tensor data.
     * - UNKNOWN: If the data characteristic is unknown, the data will not be printed.
     * - SCALAR: If the data characteristic is scalar, the data will be printed as a single value for each node or element in the dataset.
     * - D3_DOF_GLOBAL_TRANSLATION_VECTOR: If the data characteristic is a 3-DOF global translation vector, the data will be printed as a 3-component vector for each node or element in the dataset.
     * - D3_DOF_GLOBAL_TRANSLATION_ROTATION_VECTOR: If the data characteristic is a 3-DOF global translation and rotation vector, the data will be printed as a 6-component vector for each node or element in the dataset.
     * - SYMMETRIC_GLOBAL_TENSOR: If the data characteristic is a symmetric global tensor, the data will be printed as a 6-component vector for each node or element in the dataset.
     * - STRESS_RESULTANTS: If the data characteristic is stress resultants, the data will be printed as a 6-component vector for each node or element in the dataset.
     * @note 3_DOF_* Records have their name changed to D3_DOF to avoid conflicts.
     */
    enum class DataCharacteristics {
        UNKNOWN = 0,
        SCALAR = 1,
        D3_DOF_GLOBAL_TRANSLATION_VECTOR = 2,
        D3_DOF_GLOBAL_TRANSLATION_ROTATION_VECTOR = 3,
        SYMMETRIC_GLOBAL_TENSOR = 4,
        STRESS_RESULTANTS = 6
    };

    /**
     * @brief Enum class representing the data type of the UNV output.
     * @details This enum class defines the data type of the UNV output, which can be one of the following:
     * - INTEGER: The data is of integer type.
     * - SINGLE_PRECISION_FLOATING_POINT: The data is of single precision floating point type.
     * - DOUBLE_PRECISION_FLOATING_POINT: The data is of double precision floating point type.
     * - SINGLE_PRECISION_COMPLEX: The data is of single precision complex type.
     * - DOUBLE_PRECISION_COMPLEX: The data is of double precision complex type.
     */
    enum class DataType {
        INTEGER = 1,
        SINGLE_PRECISION_FLOATING_POINT = 2,
        DOUBLE_PRECISION_FLOATING_POINT = 4,
        SINGLE_PRECISION_COMPLEX = 5,
        DOUBLE_PRECISION_COMPLEX = 6
    };

    /**
     * @brief Enum class representing the write type of the UNV output.
     * @details This enum class defines the write type of the UNV output, which can be one of the following:
     * - HISTORICAL: The data is written as historical data.
     * - NON_HISTORICAL: The data is written as non-historical data.
     * - FLAG: The data is written as a flag.
     */
    enum class WriteType
    {
        HISTORICAL = 0,
        NON_HISTORICAL = 1
    };

    /// Pointer definition of UnvOutput
    KRATOS_CLASS_POINTER_DEFINITION(UnvOutput);

    /**
     * @brief Converts an enumeration value to its underlying integer type.
     * @tparam Enumeration The enumeration type to convert.
     * @param Value The enumeration value to convert.
     * @return The underlying integer value of the enumeration.
     */
    template <typename Enumeration>
    auto as_integer(const Enumeration Value)
        -> typename std::underlying_type<Enumeration>::type
    {
        return static_cast<typename std::underlying_type<Enumeration>::type>(Value);
    }

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief This constructor initializes the UnvOutput class with a reference to a ModelPart and an output file name.
     * @param rModelPart Reference to the ModelPart to be outputted.
     * @param rrOutputFileWithoutExtension Name of the output file without extension.
     */
    UnvOutput(
        Kratos::ModelPart&
        rModelPart, 
        const std::string& rrOutputFileWithoutExtension
        );

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Initializes the output file with the header and the mesh.
     */
    void InitializeOutputFile();

    /**
     * @brief Writes 'mrOutputModelPart' associated mesh.
     * 
     */
    void WriteMesh();

    /**
     * @brief Writes 'mrOutputModelPart' associated nodes.
     * @details This function writes the nodes of the 'mrOutputModelPart' to the output file. It iterates through each node in the model part, retrieves its coordinates, and writes them to the output file in the appropriate format.
     */
    void WriteNodes();

    /**
     * @brief Writes 'mrOutputModelPart' associated conditions.
     * @details This function writes the conditions of the 'mrOutputModelPart' to the output file. It iterates through each condition in the model part, retrieves its geometry, and writes the appropriate data based on the geometry type (e.g., triangles, quadrilaterals, tetrahedra, hexahedra). The function also handles unsupported geometries by throwing an error. 
     */
    void WriteElements();

    /**
     * @brief Writes 'mrOutputModelPart' associated conditions.
     * @details This function writes the conditions of the 'mrOutputModelPart' to the output file. It iterates through each condition in the model part, retrieves its geometry, and writes the appropriate data based on the geometry type (e.g., triangles, quadrilaterals). The function also handles unsupported geometries by throwing an error. 
     */
    void WriteConditions();

    /**
     * @brief Writes a result dataset containing the rVariable value for a given timestep
     * @details This function writes a result dataset to the output file, containing the values of the specified variable (rVariable) for each node in the 'mrOutputModelPart' at the given timestep. The function supports various variable types, including bool, int, double, array_1d<double,3>, Vector, and Matrix. It retrieves the values of the variable for each node and writes them to the output file in the appropriate format.
     * @param rVariable Kratos Variable to be printed
     * @param TimeStep  TimeStep.
     */
    void WriteNodalResults(const Variable<bool>& rVariable, const double TimeStep);
    void WriteNodalResults(const Variable<int>& rVariable, const double TimeStep);
    void WriteNodalResults(const Variable<double>& rVariable, const double TimeStep);
    void WriteNodalResults(const Variable<array_1d<double,3>>& rVariable, const double TimeStep);
    void WriteNodalResults(const Variable<Vector>& rVariable, const double TimeStep);
    void WriteNodalResults(const Variable<Matrix>& rVariable, const double TimeStep);

    /**
     * @brief Writes a result dataset containing the rVariable value for a given timestep (non-historical)
     * @details This function writes a result dataset (non-historical) to the output file, containing the values of the specified variable (rVariable) for each node in the 'mrOutputModelPart' at the given timestep. The function supports various variable types, including bool, int, double, array_1d<double,3>, Vector, and Matrix. It retrieves the values of the variable for each node and writes them to the output file in the appropriate format.
     * @param rVariable Kratos Variable to be printed
     * @param TimeStep  TimeStep.
     */
    void WriteNodalNonHistoricalResults(const Variable<bool>& rVariable, const double TimeStep);
    void WriteNodalNonHistoricalResults(const Variable<int>& rVariable, const double TimeStep);
    void WriteNodalNonHistoricalResults(const Variable<double>& rVariable, const double TimeStep);
    void WriteNodalNonHistoricalResults(const Variable<array_1d<double,3>>& rVariable, const double TimeStep);
    void WriteNodalNonHistoricalResults(const Variable<Vector>& rVariable, const double TimeStep);
    void WriteNodalNonHistoricalResults(const Variable<Matrix>& rVariable, const double TimeStep);

    /**
     * @brief Returns the type of unv data associated to a Kratos Variable
     * @details Vectors and Matrices are not supported at this time.
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
     * @details Vectors and Matrices are not supported at this time.
     * @param rOutputFile Output file 
     * @param rNode       Input node
     * @param rVariable   Variable to print
     */
    void WriteNodalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<bool>& rVariable);
    void WriteNodalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<int>& rVariable);
    void WriteNodalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<double>& rVariable);
    void WriteNodalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<array_1d<double,3>>& rVariable);
    void WriteNodalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<Vector>& rVariable);
    void WriteNodalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<Matrix>& rVariable);

    /**
     * @brief Writes the variable value for a node (non-historical).
     * @details Vectors and Matrices are not supported at this time.
     * @param rOutputFile Output file 
     * @param rNode       Input node
     * @param rVariable   Variable to print
     */
    void WriteNodalNonHistoricalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<bool>& rVariable);
    void WriteNodalNonHistoricalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<int>& rVariable);
    void WriteNodalNonHistoricalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<double>& rVariable);
    void WriteNodalNonHistoricalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<array_1d<double,3>>& rVariable);
    void WriteNodalNonHistoricalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<Vector>& rVariable);
    void WriteNodalNonHistoricalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<Matrix>& rVariable);

    /**
     * @brief Get the id of the UNV variable name corresponding to rVariable. 1000+ if none found.
     * @details This function returns the id of the UNV variable name corresponding to the given Kratos variable (rVariable). If no corresponding UNV variable is found, it returns a value greater than or equal to 1000.
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
     * @details This function writes a result dataset using the results in node mode. The format is partially extracted from: http://users.ices.utexas.edu
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
     * @param rVariable          Variable to be printed 
     * @param NumberOfComponents Number of components of the variable 
     * @param TimeStep           Current TimeStep
     * @tparam TVariablebleType  Type of the variable to be printed
     * @tparam TWriteType        Type of the write operation (historical or non-historical)
     */
    template<class TVariablebleType, WriteType TWriteType = WriteType::HISTORICAL>
    void WriteNodalResultRecords(const TVariablebleType& rVariable, const int NumberOfComponents, const double TimeStep) {
        std::ofstream rOutputFile;
        rOutputFile.open(mOutputFileName, std::ios::out | std::ios::app);

        std::string data_set_name = "NodalResults";
        const std::string& r_data_set_label = rVariable.Name();

        rOutputFile << std::setw(6)  << "-1" << "\n";                                                // Begin block
        rOutputFile << std::setw(6)  << as_integer(DatasetID::RESULTS_DATASET) << "\n";              // DatasetID

        rOutputFile << std::setw(10) << r_data_set_label << "\n";                                        // Record 1 - Label
        rOutputFile << std::setw(6)  << data_set_name << "\n";                                         // Record 2 - Name
        rOutputFile << std::setw(10) << as_integer(DatasetLocation::DATA_AT_NODES) << "\n";          // Record 3

        // String records, seems like you can put anything you want.
        rOutputFile << "\n" << "\n" << "\n" << "\n" << "\n";                                         // Records 4-8
        
        // ModelType, AnalysisType, DataCharacteristic, ResultType, DataType, NumberOfDataValues    // Record 9
        rOutputFile << std::setw(10) << as_integer(ModelType::STRUCTURAL); 
        rOutputFile << std::setw(10) << as_integer(AnalysisType::TRANSIENT);
        rOutputFile << std::setw(10) << as_integer(GetDataType(rVariable));
        rOutputFile << std::setw(10) << GetUnvVariableName(rVariable);
        rOutputFile << std::setw(10) << as_integer(DataType::SINGLE_PRECISION_FLOATING_POINT);
        rOutputFile << std::setw(10) << NumberOfComponents; 
        rOutputFile << "\n";

        // DesignSetId, IterationNumber, SolutionSetId, BoundaryCondition, LoadSet, ModeNumber, TimeStampNumber, FrequencyNumber
        rOutputFile << std::setw(10) << 0;                                                           // Record 10
        rOutputFile << std::setw(10) << TimeStep;
        rOutputFile << std::setw(10) << 0;
        rOutputFile << std::setw(10) << 0;
        rOutputFile << std::setw(10) << 0;
        rOutputFile << std::setw(10) << 1;
        rOutputFile << std::setw(10) << TimeStep;
        rOutputFile << std::setw(10) << 0;
        rOutputFile << "\n";

        // CreationOption, (Unknown)*7
        rOutputFile << std::setw(10) << 0;                                                           // Record 11
        rOutputFile << std::setw(10) << 0;
        rOutputFile << "\n";

        // Time, Frequency Eigenvalue NodalMass ViscousDampingRatio, HystereticDampingRatio
        rOutputFile << std::setw(13) << TimeStep * 0.1;                                              // Record 12
        rOutputFile << std::setw(13) << "0.00000E+00";
        rOutputFile << std::setw(13) << "0.00000E+00";
        rOutputFile << std::setw(13) << "0.00000E+00";
        rOutputFile << std::setw(13) << "0.00000E+00";
        rOutputFile << std::setw(13) << "0.00000E+00";
        rOutputFile << "\n";

        // Eigenvalue_re, Eigenvalue_im, ModalA_re, ModalA_im, ModalB_re, ModalB_im
        rOutputFile << std::setw(13) << "0.00000E+00";                                               // Record 13
        rOutputFile << std::setw(13) << "0.00000E+00";
        rOutputFile << std::setw(13) << "0.00000E+00";
        rOutputFile << std::setw(13) << "0.00000E+00";
        rOutputFile << std::setw(13) << "0.00000E+00";
        rOutputFile << std::setw(13) << "0.00000E+00";
        rOutputFile << "\n";

        // Data at nodes
        int node_label;
        for (auto& r_node : mrOutputModelPart.Nodes()) {
            node_label = r_node.Id();
            rOutputFile << std::setw(6) << node_label << "\n";
            if constexpr (TWriteType == WriteType::HISTORICAL) {                                    // Record 14 - Node Number
                WriteNodalResultValues(rOutputFile, r_node, rVariable);                                  // Record 15 - NumberOfDataValues' data of the node
            } else if constexpr (TWriteType == WriteType::NON_HISTORICAL) {
                WriteNodalNonHistoricalResultValues(rOutputFile, r_node, rVariable);
            } else {
                KRATOS_ERROR << "Unknown WriteType" << std::endl;
            }
        }
        
        rOutputFile << std::setw(6) << "-1" << "\n";
        rOutputFile.close();
    }
    ///@}
private:
    ///@name Member Variables
    ///@{

    Kratos::ModelPart& mrOutputModelPart; /// Reference to the model part to be printed
    std::string mOutputFileName;          /// Name of the output file
    
    ///@}
}; /// class UnvOutput

///@}

} /// namespace Kratos
