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
//  Contributors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/geometry_data.h"
#include "utilities/string_utilities.h"
#include "processes/integration_values_extrapolation_to_nodes_process.h"

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
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) UnvOutput 
{
public:
    ///@name Type Definitions
    ///@{

    // Scalar UNV variables: quantities that are single-valued per node/element
    static inline const std::vector<std::pair<int, std::string>> unv_scalar_variables = {
        {7,   "Strain Energy"},
        {10,  "Kinetic Energy"},
        {13,  "Strain Energy Density"},
        {14,  "Kinetic Energy Density"},
        {15,  "Hydrostatic Pressure"},
        {17,  "Code Check Value"},
        {18,  "Coefficient of Pressure"},
        {21,  "Failure Index for Ply"},
        {22,  "Failure Index for Bonding"},
        {23,  "Reaction Heat Flow"},
        {24,  "Stress Error Density"},
        {28,  "Length"},
        {29,  "Area"},
        {30,  "Volume"},
        {31,  "Mass"},
        {36,  "Strain Energy Error Norm"},
        {38,  "Heat Transfer Coefficient"},
        {40,  "Kinetic Energy Dissipation Rate"},
        {41,  "Strain Energy Error"},
        {42,  "Mass Flow"},
        {44,  "Heat Flow"},
        {45,  "View Factor"},
        {46,  "Heat Load"},
        {50,  "Contact Pressure"},
        {55,  "Infrared Heat Flux"},
        {56,  "Diffuse Solar Heat Flux"},
        {57,  "Collimated Solar Heat Flux"},
        {58,  "Safety Factor"},
        {59,  "Fatigue Damage"},
        {61,  "Fatigue Life"},
        {62,  "Quality Index"},
        {101, "Gap Thickness"},
        {108, "Layered Shear Strain Rate"},
        {113, "Bulk Temperature"},
        {114, "Peak Temperature"},
        {115, "Temperature at Fill"},
        {116, "Mass Density"},
        {118, "Volumetric Shrinkage"},
        {119, "Filling Time"},
        {120, "Ejection Time"},
        {121, "No-flow Time"},
        {122, "Weld Line Meeting Angle"},
        {123, "Weld Line Underflow"},
        {124, "Original Runner Diameter"},
        {125, "Optimized Runner Diameter"},
        {126, "Change in Runner Diameter"},
        {127, "Averaged Layered Cure"},
        {129, "Cure Rate"},
        {130, "Cure Time"},
        {131, "Induction Time"},
        {132, "Temperature at Cure"},
        {133, "Percent Gelation"},
        {138, "Part Ejection Time"},
        {139, "Part Peak Temperature"},
        {140, "Part Average Temperature"},
        {145, "Wall Temperature Convergence"},
        {149, "Line Pressure"},
        {150, "Reynold's Number"},
        {151, "Line Film Coefficient"},
        {152, "Line Temperature"},
        {153, "Line Bulk Temperature"},
        {154, "Mold Temperature"},
        {156, "Rod Heater Temperature"},
        {158, "Original Line Diameter"},
        {159, "Optimized Line Diameter"},
        {160, "Change in Line Diameter"},
        {161, "Air Traps"},
        {162, "Weld Lines"},
        {163, "Injection Growth"},
        {164, "Temp Diff (Celcius)"},
        {165, "Shear Rate"},
        {166, "Viscosity"},
        {167, "Percentage"},
        {168, "Time"},
        {170, "Speed"},
        {171, "Flow Rate"},
        {172, "Thickness Ratio"},
        {301, "Sound Pressure"},
        {302, "Sound Power"},
        {304, "Sound Energy"},
        {305, "Sound Energy Density"},
    };

    // 3-component vector UNV variables: quantities with X, Y, Z components per node/element
    static inline const std::vector<std::pair<int, std::string>> unv_vector3_variables = {
        {16,  "Heat Gradient"},
        {32,  "Constraint Force"},
        {39,  "Temperature Gradient"},
        {43,  "Mass Flux"},
        {49,  "Contact Forces"},
        {105, "Flow Vector at Fill"},
        {106, "Bulk Flow Vector"},
        {107, "Core Displacement"},
        {169, "Flow Direction"},
        {303, "Sound Intensity"},
    };

    /**
     * @brief Pointer definition of UnvOutput
     * @details This defines a pointer to UnvOutput
     * - AUTOMATIC: The entity type is automatically determined based on the model part.
     * - ELEMENTS: The entity type is explicitly set to elements, and the data will be printed for each element in the dataset.
     * - CONDITIONS: The entity type is explicitly set to conditions, and the data will be printed for each condition in the dataset.
     */
    enum class EntityType {
        AUTOMATIC = 0,
        ELEMENTS  = 1,
        CONDITIONS= 2
    }; 

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

    /**
     * @brief Enum class representing the container a result dataset refers to.
     * @details This enum class selects whether a result dataset (2414) stores data at nodes,
     * elements or conditions. It drives the dataset location record and the entity iteration.
     * - NODES: The result is written for each node (DATA_AT_NODES).
     * - ELEMENTS: The result is written for each element (DATA_AT_ELEMENTS).
     * - CONDITIONS: The result is written for each condition (DATA_AT_ELEMENTS).
     */
    enum class ResultLocation
    {
        NODES = 0,
        ELEMENTS = 1,
        CONDITIONS = 2
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

    /**
     * @brief Struct to hold the lists of variables to be outputted in UNV format.
     * @details This struct contains vectors of pointers to different types of variables (bool, int, double, array_1d<double, 3>) that are to be outputted in UNV format. It is used to organize and manage the variables for output operations.
     */
    struct VariablesLists
    {
        std::vector<const Variable<bool>*> mBoolVariables;
        std::vector<const Variable<int>*> mIntVariables;
        std::vector<const Variable<double>*> mDoubleVariables;
        std::vector<const Variable<array_1d<double, 3>>*> mArray1DVariables;
        // std::vector<const Variable<Vector>*> mVectorVariables; // NOTE: Current unsupported
        // std::vector<const Variable<Matrix>*> mMatrixVariables; // NOTE: Current unsupported

        /**
         * @brief Initializes the variable lists based on the provided variable names.
         * @param rVariableNames A vector of variable names to be initialized.
         * @param rUnvVariableKeys A map to store the corresponding UNV variable keys for the initialized variables.
         */
        void Initialize(
            const std::vector<std::string>& rVariableNames,
            std::unordered_map<std::size_t, int>& rUnvVariableKeys
            );
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor of the class (Model + Parameters)
     * @details This constructor initializes the UnvOutput class with a reference to a Model and a set of parameters.
     * @param rModel Reference to the Model to be outputted.
     * @param Settings Parameters for the output process, including output file name and other settings.
     */
    UnvOutput(
        Model& rModel,
        Parameters Settings
        );

    /**
     * @brief This constructor initializes the UnvOutput class with a reference to a ModelPart and an output file name.
     * @param rModelPart Reference to the ModelPart to be outputted.
     * @param rOutputFileWithoutExtension Name of the output file without extension.
     */
    UnvOutput(
        Kratos::ModelPart& rModelPart, 
        const std::string& rOutputFileWithoutExtension
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
     * @details This function writes the mesh of the 'mrOutputModelPart' to the output file. It includes the nodes, elements, and conditions of the model part, along with their associated data. The function ensures that the mesh is properly formatted and ready for output in UNV format.
     */
    void WriteMesh();

    /**
     * @brief Writes the results for the specified variable at the given timestep.
     * @details Uses the mHistoricalVariables and mNonHistoricalVariables to write the results for the specified variable at the given timestep. The function retrieves the values of the variable for each node in the 'mrOutputModelPart' and writes them to the output file in the appropriate format. It supports various variable types, including bool, int, double, array_1d<double,3>, Vector, and Matrix.
     */
    void PrintOutput();

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

    ///@}
private:
    ///@name Member Variables
    ///@{

    /// Multiplier used to derive unique UNV labels for the sub-elements of a decomposed entity.
    static constexpr long long SubElementLabelFactor = 16;

    Kratos::ModelPart& mrOutputModelPart;                  /// Reference to the model part to be printed
    std::string mOutputFileName;                           /// Name of the output file

    VariablesLists mHistoricalVariables;                   /// List of nodal historical variables to be printed
    VariablesLists mNonHistoricalVariables;                /// List of nodal non-historical variables to be printed
    VariablesLists mElementVariables;                      /// List of element (non-historical) variables to be printed
    VariablesLists mConditionVariables;                    /// List of condition (non-historical) variables to be printed
    VariablesLists mGaussPointVariables;                   /// List of variables to be averaged over the integration points and printed on elements

    std::vector<std::pair<const Flags*, std::string>> mNodalFlags;     /// Nodal flags to be printed
    std::vector<std::pair<const Flags*, std::string>> mElementFlags;   /// Element flags to be printed
    std::vector<std::pair<const Flags*, std::string>> mConditionFlags; /// Condition flags to be printed

    std::unordered_map<std::size_t, int> mUnvVariableKeys; /// Map to store the keys of the UNV variables for quick access

    EntityType mEntityType = EntityType::AUTOMATIC;        /// Type of entity to be printed (automatic, elements, or conditions)
    unsigned int mDefaultPrecision = 7;                    /// Precision used when writing floating point result values
    bool mDecomposeQuadraticIntoLinear = true;             /// Flag to decompose quadratic geometries into linear sub-elements (for stricter readers, e.g. Simcenter 3D)
    bool mWriteDeformedConfiguration = false;              /// Flag to indicate if the deformed configuration should be written
    bool mWriteIds = false;                                /// Flag to indicate if the entity ids should be written as result datasets
    bool mOutputSubModelParts = false;                     /// Flag to indicate if the sub model parts should be written as UNV groups
    bool mInitializedOutputFile = false;                   /// Flag to indicate if the output file has been initialized
    bool mMeshWritten = false;                             /// Flag to indicate if the mesh has been written to the output file

    /// Process used to extrapolate the integration point values to the nodes (optional)
    IntegrationValuesExtrapolationToNodesProcess::UniquePointer mpGaussToNodesProcess = nullptr;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Describes how a Kratos geometry maps onto a UNV finite-element descriptor.
     * @details Holds the UNV FE descriptor id (dataset 2412), the number of nodes actually written,
     * whether the entity is a beam (requiring the extra orientation record), and the reorder array
     * that converts the Kratos local node ordering into the UNV convention. Surface geometries carry
     * two descriptor ids: one used when the geometry lives in a 2D working space (plane element) and
     * one for a 3D working space (shell). A value of -1 means the geometry has no descriptor for that
     * working space. When @a Degrade is true, the geometry has no exact UNV representation and the
     * reorder array already drops the surplus (interior) Kratos nodes.
     */
    struct UnvElementDescriptor {
        int fe_descriptor_2d = -1;   /// Descriptor when written in a 2D working space (-1 = N/A)
        int fe_descriptor_3d = -1;   /// Descriptor when written in a 3D working space (-1 = N/A)
        int num_nodes = 0;           /// Number of nodes written to the UNV connectivity record
        bool is_beam = false;        /// Whether the extra beam orientation record must be written
        std::vector<int> reorder;    /// Kratos local indices to emit, in order (empty => identity)
        bool degrade = false;        /// Whether the geometry is represented with a loss of interior nodes
    };

    /**
     * @brief Returns the UNV descriptor associated to a Kratos geometry type.
     * @param Type The Kratos geometry type to look up.
     * @param rOut The descriptor filled on success.
     * @return true if the geometry can be represented in UNV, false otherwise (e.g. pyramids).
     */
    static bool TryGetUnvDescriptor(
        const GeometryData::KratosGeometryType Type,
        UnvElementDescriptor& rOut
        );

    /**
     * @brief Describes how a Kratos geometry is decomposed into linear UNV sub-elements.
     * @details Some readers (e.g. Simcenter 3D) do not robustly support quadratic UNV element
     * descriptors. When decomposition is enabled each geometry is written as one or more linear
     * sub-elements. The @a sub_elements list holds, for each sub-element, the Kratos local node
     * indices to emit. Linear geometries decompose into a single identity sub-element. Serendipity
     * solids without interior nodes (Hexahedra20/27, Prism15) reduce to their linear corner element,
     * dropping the mid-side nodes.
     */
    struct UnvLinearDecomposition {
        int fe_descriptor_2d = -1;                 /// Linear descriptor in a 2D working space (-1 = N/A)
        int fe_descriptor_3d = -1;                 /// Linear descriptor in a 3D working space (-1 = N/A)
        int num_nodes = 0;                         /// Number of nodes per sub-element
        bool is_beam = false;                      /// Whether the sub-elements need the beam orientation record
        bool drops_nodes = false;                  /// Whether mid-side nodes are dropped (serendipity solids)
        std::vector<std::vector<int>> sub_elements; /// Kratos local indices of each linear sub-element
    };

    /**
     * @brief Returns the linear decomposition of a Kratos geometry type.
     * @param Type The Kratos geometry type to decompose.
     * @param rOut The decomposition filled on success.
     * @return true if the geometry can be represented in UNV, false otherwise (e.g. pyramids).
     */
    static bool TryGetLinearDecomposition(
        const GeometryData::KratosGeometryType Type,
        UnvLinearDecomposition& rOut
        );

    /**
     * @brief Returns the number of UNV sub-elements a geometry maps to.
     * @details Returns 1 when decomposition is disabled or the geometry is linear; otherwise the
     * number of linear sub-elements from the decomposition table.
     * @param Type The Kratos geometry type.
     */
    int GetSubElementCount(const GeometryData::KratosGeometryType Type) const;

    /**
     * @brief Returns the UNV element label(s) an entity maps to.
     * @details With decomposition disabled the original id is used. With decomposition enabled the
     * labels are @a id * SubElementLabelFactor + k so that each entity owns a unique, contiguous block
     * of labels shared consistently between the mesh, results and group datasets.
     * @param rEntity The element or condition.
     */
    template<class TEntity>
    std::vector<long long> GetEntityLabels(const TEntity& rEntity) const
    {
        std::vector<long long> labels;
        const long long id = static_cast<long long>(rEntity.Id());
        if (!mDecomposeQuadraticIntoLinear) {
            labels.push_back(id);
        } else {
            const int sub_count = GetSubElementCount(rEntity.GetGeometry().GetGeometryType());
            for (int k = 0; k < sub_count; ++k) {
                labels.push_back(id * SubElementLabelFactor + k);
            }
        }
        return labels;
    }

    /**
     * @brief Writes 'mrOutputModelPart' associated nodes.
     * @details This function writes the nodes of the 'mrOutputModelPart' to the output file. It iterates through each node in the model part, retrieves its coordinates, and writes them to the output file in the appropriate format.
     */
    void WriteNodes();

    /**
     * @brief Writes the elements or conditions of 'mrOutputModelPart' as a UNV 2412 dataset.
     * @details This unified writer iterates over the provided container (elements or conditions),
     * resolves each geometry to its UNV descriptor, and writes the 2412 records: the element header,
     * the optional beam orientation record and the (reordered) connectivity wrapped at 8 ids per line.
     * Geometries without an exact UNV representation are degraded (with a one-time warning) or, if they
     * cannot be represented at all (e.g. pyramids), an error is thrown.
     * @param rContainer The elements or conditions container to be written.
     * @tparam TContainerType The type of the entity container.
     */
    template<typename TContainerType>
    void WriteEntities(const TContainerType& rContainer);

    /**
     * @brief Writes the sub model parts of 'mrOutputModelPart' as UNV groups (dataset 2467).
     * @details Each sub model part is written as a permanent group referencing its elements and
     * conditions (entity type code 8) and, optionally, its nodes (entity type code 7).
     */
    void WriteGroups();

    /**
     * @brief Runs the integration point extrapolation to nodes process, if configured.
     * @details Mirrors VtkOutput: when 'gauss_point_variables_extrapolated_to_nodes' are requested,
     * the associated process is executed so that the extrapolated values are available as
     * non-historical nodal data at print time.
     */
    void PrepareGaussPointResults();

    /**
     * @brief Resolves a list of flag names into the provided flag list.
     * @param rFlagNames The names of the flags to resolve.
     * @param rFlagList The list to populate with the resolved (flag pointer, name) pairs.
     */
    void InitializeFlags(
        const std::vector<std::string>& rFlagNames,
        std::vector<std::pair<const Flags*, std::string>>& rFlagList
        );

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
     * @brief Retrieves the value of a variable stored in an entity (node, element or condition).
     * @details For nodes the value can be read from the historical database (FastGetSolutionStepValue)
     * or the non-historical one (GetValue); elements and conditions only expose non-historical data.
     * @param rEntity The entity to read from.
     * @param rVariable The variable to read.
     * @param IsHistorical Whether to read the historical value (nodes only).
     * @return A const reference to the stored value.
     */
    template<class TEntity, class TVarType>
    static const TVarType& GetEntityValue(
        const TEntity& rEntity,
        const Variable<TVarType>& rVariable,
        const bool IsHistorical
        )
    {
        if constexpr (std::is_same_v<TEntity, Node>) {
            return IsHistorical ? rEntity.FastGetSolutionStepValue(rVariable) : rEntity.GetValue(rVariable);
        } else {
            return rEntity.GetValue(rVariable);
        }
    }

    /**
     * @brief Writes the result value record (Record 15) of a single entity.
     * @details Handles scalar (bool/int/double) and 3-component (array_1d<double,3>) variables.
     * @param rOutputFile Output file.
     * @param rEntity Entity to read from (node, element or condition).
     * @param rVariable Variable to print.
     * @param IsHistorical Whether to read the historical value (nodes only).
     */
    template<class TEntity, class TVarType>
    void WriteEntityResultValues(
        std::ofstream& rOutputFile,
        const TEntity& rEntity,
        const Variable<TVarType>& rVariable,
        const bool IsHistorical
        )
    {
        const int width = static_cast<int>(mDefaultPrecision) + 9;
        if constexpr (std::is_same_v<TVarType, array_1d<double, 3>>) {
            const auto& r_value = GetEntityValue(rEntity, rVariable, IsHistorical);
            rOutputFile << std::setw(width) << r_value[0];
            rOutputFile << std::setw(width) << r_value[1];
            rOutputFile << std::setw(width) << r_value[2] << "\n";
        } else {
            rOutputFile << std::setw(width) << GetEntityValue(rEntity, rVariable, IsHistorical) << "\n";
        }
    }

    /**
     * @brief Writes the common header records (1-13) of a results dataset (2414).
     * @param rOutputFile Output file.
     * @param rLabel Dataset label (Record 1), typically the variable/flag name.
     * @param rName Dataset name (Record 2).
     * @param Location Dataset location code (Record 3).
     * @param DataCharacteristic Data characteristic code (Record 9).
     * @param UnvVariableId UNV result type id (Record 9).
     * @param NumberOfComponents Number of data values per entity (Record 9).
     * @param TimeStep Current time step.
     */
    void WriteResultDatasetHeader(
        std::ofstream& rOutputFile,
        const std::string& rLabel,
        const std::string& rName,
        const int Location,
        const int DataCharacteristic,
        const int UnvVariableId,
        const int NumberOfComponents,
        const double TimeStep
        );

    /**
     * @brief Get the id of the UNV variable name corresponding to rVariable. 1000+ if none found.
     * @details This function retrieves the UNV variable name corresponding to the provided Kratos variable (rVariable). It checks if the variable's key exists in the mUnvVariableKeys map. If found, it returns the associated UNV variable name; otherwise, it returns a value greater than or equal to 1000, indicating that no corresponding UNV variable name was found.
     * @param rVariable The Kratos variable for which to retrieve the UNV variable name
     */
    int GetUnvVariableName(const VariableData& rVariable);

    /**
     * @brief Writes a result dataset (2414) for nodes, elements or conditions.
     * @details This function writes a result dataset. The format is partially extracted from: http://users.ices.utexas.edu
     * R.  1: unique number of dataset (dataset_label)
     * R.  2: text describing content (dataset_name)
     * R.  3: data belongs to: nodes, elements,... (dataset_location)
     * R.  4-8: user-specified text
     * R.  9: (model_type) (analysis_type) (data_characteristic) (result_type) (data_type) (nvaldc)
     * R. 10-13: analysis specific records (load/mode/time/eigenvalue)
     *
     * For every entity (node / element / condition):
     * R. 14: (entity_id) [nvaldc for elements/conditions]
     * R. 15: (result)*nvaldc
     *
     * @param rVariable          Variable to be printed
     * @param NumberOfComponents Number of components of the variable
     * @param TimeStep           Current TimeStep
     * @tparam TVariableType     Type of the variable to be printed
     * @tparam TWriteType        Type of the write operation (historical or non-historical)
     * @tparam TLoc              Location the results refer to (nodes, elements, conditions)
     */
    template<class TVariableType, WriteType TWriteType, ResultLocation TLoc>
    void WriteResultRecords(const TVariableType& rVariable, const int NumberOfComponents, const double TimeStep) {
        std::ofstream rOutputFile;
        rOutputFile.open(mOutputFileName, std::ios::out | std::ios::app);
        rOutputFile << std::scientific << std::setprecision(mDefaultPrecision);

        constexpr DatasetLocation location = (TLoc == ResultLocation::NODES) ? DatasetLocation::DATA_AT_NODES : DatasetLocation::DATA_AT_ELEMENTS;
        const std::string data_set_name = (TLoc == ResultLocation::NODES) ? "NodalResults" : "ElementResults";

        WriteResultDatasetHeader(rOutputFile, rVariable.Name(), data_set_name, as_integer(location),
            as_integer(GetDataType(rVariable)), GetUnvVariableName(rVariable), NumberOfComponents, TimeStep);

        if constexpr (TLoc == ResultLocation::NODES) {
            for (auto& r_node : mrOutputModelPart.Nodes()) {
                rOutputFile << std::setw(10) << static_cast<int>(r_node.Id()) << "\n";                // Record 14 - Node Number
                WriteEntityResultValues(rOutputFile, r_node, rVariable, TWriteType == WriteType::HISTORICAL); // Record 15 - Data values
            }
        } else if constexpr (TLoc == ResultLocation::ELEMENTS) {
            for (auto& r_element : mrOutputModelPart.Elements()) {
                for (const long long label : GetEntityLabels(r_element)) {
                    rOutputFile << std::setw(10) << label << std::setw(10) << NumberOfComponents << "\n";
                    WriteEntityResultValues(rOutputFile, r_element, rVariable, false);
                }
            }
        } else {
            for (auto& r_condition : mrOutputModelPart.Conditions()) {
                for (const long long label : GetEntityLabels(r_condition)) {
                    rOutputFile << std::setw(10) << label << std::setw(10) << NumberOfComponents << "\n";
                    WriteEntityResultValues(rOutputFile, r_condition, rVariable, false);
                }
            }
        }

        rOutputFile << std::setw(6) << "-1" << "\n";
        rOutputFile.close();
    }

    /**
     * @brief Writes a scalar result dataset (2414) holding a flag value (0/1, or -1 if undefined).
     * @param rFlag The flag to be printed.
     * @param rName The name of the flag (used as the dataset label).
     * @param rContainer The container of entities (nodes, elements or conditions) to iterate.
     * @param TimeStep Current TimeStep.
     * @tparam TLoc Location the flags refer to (nodes, elements, conditions).
     * @tparam TContainerType The type of the entity container.
     */
    template<ResultLocation TLoc, class TContainerType>
    void WriteFlagRecords(const Flags& rFlag, const std::string& rName, TContainerType& rContainer, const double TimeStep) {
        std::ofstream rOutputFile;
        rOutputFile.open(mOutputFileName, std::ios::out | std::ios::app);
        rOutputFile << std::scientific << std::setprecision(mDefaultPrecision);

        constexpr DatasetLocation location = (TLoc == ResultLocation::NODES) ? DatasetLocation::DATA_AT_NODES : DatasetLocation::DATA_AT_ELEMENTS;
        WriteResultDatasetHeader(rOutputFile, rName, "Flags", as_integer(location),
            as_integer(DataCharacteristics::SCALAR), 2000, 1, TimeStep);

        const int width = static_cast<int>(mDefaultPrecision) + 9;
        for (auto& r_entity : rContainer) {
            const double value = r_entity.IsDefined(rFlag) ? (r_entity.Is(rFlag) ? 1.0 : 0.0) : -1.0;
            if constexpr (TLoc == ResultLocation::NODES) {
                rOutputFile << std::setw(10) << static_cast<int>(r_entity.Id()) << "\n";
                rOutputFile << std::setw(width) << value << "\n";
            } else {
                for (const long long label : GetEntityLabels(r_entity)) {
                    rOutputFile << std::setw(10) << label << std::setw(10) << 1 << "\n";
                    rOutputFile << std::setw(width) << value << "\n";
                }
            }
        }

        rOutputFile << std::setw(6) << "-1" << "\n";
        rOutputFile.close();
    }

    /**
     * @brief Writes a scalar result dataset (2414) holding the entity id.
     * @param rName The name of the dataset label.
     * @param rContainer The container of entities to iterate.
     * @param TimeStep Current TimeStep.
     * @tparam TLoc Location the ids refer to (nodes, elements, conditions).
     * @tparam TContainerType The type of the entity container.
     */
    template<ResultLocation TLoc, class TContainerType>
    void WriteIdRecords(const std::string& rName, TContainerType& rContainer, const double TimeStep) {
        std::ofstream rOutputFile;
        rOutputFile.open(mOutputFileName, std::ios::out | std::ios::app);
        rOutputFile << std::scientific << std::setprecision(mDefaultPrecision);

        constexpr DatasetLocation location = (TLoc == ResultLocation::NODES) ? DatasetLocation::DATA_AT_NODES : DatasetLocation::DATA_AT_ELEMENTS;
        WriteResultDatasetHeader(rOutputFile, rName, "Ids", as_integer(location),
            as_integer(DataCharacteristics::SCALAR), 2001, 1, TimeStep);

        const int width = static_cast<int>(mDefaultPrecision) + 9;
        for (auto& r_entity : rContainer) {
            const int id = static_cast<int>(r_entity.Id());
            if constexpr (TLoc == ResultLocation::NODES) {
                rOutputFile << std::setw(10) << id << "\n";
                rOutputFile << std::setw(width) << static_cast<double>(id) << "\n";
            } else {
                for (const long long label : GetEntityLabels(r_entity)) {
                    rOutputFile << std::setw(10) << label << std::setw(10) << 1 << "\n";
                    rOutputFile << std::setw(width) << static_cast<double>(id) << "\n";
                }
            }
        }

        rOutputFile << std::setw(6) << "-1" << "\n";
        rOutputFile.close();
    }

    /**
     * @brief Writes an element result dataset (2414) whose value is averaged over the integration points.
     * @details For each element, CalculateOnIntegrationPoints is evaluated and the mean over the
     * integration points is written as the element value (DATA_AT_ELEMENTS).
     * @param rVariable Variable to be averaged and printed.
     * @param NumberOfComponents Number of components of the variable.
     * @param TimeStep Current TimeStep.
     * @tparam TVarType The value type of the variable (double or array_1d<double,3>).
     */
    template<class TVarType>
    void WriteGaussPointElementResults(const Variable<TVarType>& rVariable, const int NumberOfComponents, const double TimeStep) {
        std::ofstream rOutputFile;
        rOutputFile.open(mOutputFileName, std::ios::out | std::ios::app);
        rOutputFile << std::scientific << std::setprecision(mDefaultPrecision);

        WriteResultDatasetHeader(rOutputFile, rVariable.Name(), "GaussPointResults", as_integer(DatasetLocation::DATA_AT_ELEMENTS),
            as_integer(GetDataType(rVariable)), GetUnvVariableName(rVariable), NumberOfComponents, TimeStep);

        const int width = static_cast<int>(mDefaultPrecision) + 9;
        const auto& r_process_info = mrOutputModelPart.GetProcessInfo();
        std::vector<TVarType> gauss_values;
        for (auto& r_element : mrOutputModelPart.Elements()) {
            r_element.CalculateOnIntegrationPoints(rVariable, gauss_values, r_process_info);
            TVarType average = rVariable.Zero();
            const std::size_t number_of_gauss_points = gauss_values.size();
            if (number_of_gauss_points > 0) {
                for (const auto& r_value : gauss_values) {
                    average += r_value;
                }
                average /= static_cast<double>(number_of_gauss_points);
            }
            for (const long long label : GetEntityLabels(r_element)) {
                rOutputFile << std::setw(10) << label << std::setw(10) << NumberOfComponents << "\n";
                if constexpr (std::is_same_v<TVarType, array_1d<double, 3>>) {
                    rOutputFile << std::setw(width) << average[0];
                    rOutputFile << std::setw(width) << average[1];
                    rOutputFile << std::setw(width) << average[2] << "\n";
                } else {
                    rOutputFile << std::setw(width) << average << "\n";
                }
            }
        }

        rOutputFile << std::setw(6) << "-1" << "\n";
        rOutputFile.close();
    }

    /**
     * @brief Returns the default parameters for the UnvOutput class.
     * @details This function returns a Parameters object containing the default parameters for the UnvOutput class. The default parameters include the output path, custom name prefix, custom name postfix, and model part name. These parameters can be used to configure the UnvOutput class when creating an instance of it.
     * @return Parameters object containing the default parameters for the UnvOutput class.
     */
    const Parameters GetDefaultParameters() const;

    ///@}
}; /// class UnvOutput

///@}

} /// namespace Kratos
