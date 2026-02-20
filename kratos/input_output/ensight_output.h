//
//   |   /         |
//   ' /    __| _` | __|  _ \   __|
//   . \   |   (   | |    (   |\__ `
//  _|\_\_|  \__,_|\__|\___/ ____/
//         Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <map>
#include <unordered_map>

// External includes

// Project includes
#include "geometries/bounding_box.h"
#include "includes/kratos_parameters.h"
#include "includes/io.h"
#include "processes/integration_values_extrapolation_to_nodes_process.h"

namespace Kratos
{
///@addtogroup KratosCore

///@name Kratos Classes
///@{

/**
 * @struct PartData
 * @brief Structure representing a part in the EnSight 6/Gold output.
 * @details
 * This structure holds all the necessary information for a single part (typically a SubModelPart)
 * when exporting to the EnSight 6/Gold format. It includes:
 * - The list of nodes belonging to the part.
 * - A mapping from Kratos node IDs to local IDs used in the EnSight output.
 * - The geometrical objects (elements or conditions) associated with the part, organized by type.
 * - The part's unique identifier and name.
 * - A flag indicating whether the part contains elements (true) or conditions (false).
 * @var PartNodes Vector of pointers to the nodes belonging to this part.
 * @var KratosIdToLocalId Mapping from Kratos node IDs to local IDs for EnSight output.
 * @var PartGeometricalObjects Map from geometry type name to vector of pointers to geometrical objects (elements/conditions).
 * @var PartId Unique identifier for the part.
 * @var PartName Name of the part (typically the SubModelPart name).
 * @var PartElements True if the part contains elements, false if it contains conditions.
 */
struct PartData
{
    std::vector<const Node*> PartNodes;                                                  ///< The nodes belonging to this part
    std::unordered_map<std::size_t, std::size_t> KratosIdToLocalId;                      ///< Mapping from Kratos node IDs to local IDs
    std::map<std::string, std::vector<const GeometricalObject*>> PartGeometricalObjects; ///< Geometrical objects (elements/conditions) belonging to this part, organized by type (must be ordered to ensure consistency in output)
    unsigned int PartId = 0;                                                             ///< Unique identifier for the part
    std::string PartName;                                                                ///< Name of the part
    bool PartElements = true;                                                            ///< True if the part contains elements, false if it contains conditions
};

/**
 * @class EnSightOutput
 * @brief A class for exporting simulation results to the EnSight 6/Gold file format.
 * @details This class provides functionality for writing mesh geometry, nodal and elemental results, and transient data to files compatible with the EnSight 6/Gold visualization software.
 * The EnSight Case Gold format is a file structure used for scientific visualization data, particularly in computational fluid dynamics (CFD). It's not a single file, but rather a set of files that work together to describe a complete dataset. The main components are:
 * - Case file (.case or .encas): This is the main file that references all other files. It's a text file that contains metadata such as the format type (ensight 6/gold), geometry file names, variable file names, and time step information.
 * - Geometry file (.geo): This file contains the mesh information. It defines the coordinates of the nodes and how those nodes are connected to form elements (e.g., hexahedrons, tetrahedrons). The format is part-based, with each part having its own local coordinate array.
 * - Variable files (.scl, .vec, .ten, etc.): These files store the actual simulation results.
 *   + .scl: Scalar data (one value per node or element).
 *   + .vec: Vector data (three values per node or element).
 *   + .ten: Tensor data (six or nine values per node or element).
 * It supports both ASCII and binary output formats, and can handle results stored at nodes or extrapolated from integration (Gauss) points.
 * Key features:
 * - Writes EnSight 6/Gold case, geometry, and variable files.
 * - Supports output of historical and non-historical nodal variables.
 * - Handles element variables and extrapolation of integration point results to nodes.
 * - Tracks transient simulation data (time values) for time-dependent output.
 * - Allows customization of output settings via Parameters.
 * - Provides helper methods for data collection and file writing.
 * Usage:
 * 1. Construct with a reference to the ModelPart and optional Parameters for output settings.
 * 2. Call PrintOutput() to write the results for the current simulation step.
 * @see EnSight 7 User Manual, Section 11.1 for file format details. Link: https://fr.scribd.com/document/156786965/Ensight-File-Format-Manual
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) EnSightOutput
    : public IO
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the size type
    using SizeType = std::size_t;

    /// Definition of the index type
    using IndexType = std::size_t;

    /// Definition of the geometry type with given NodeType
    using GeometryType = Geometry<Node>;

    /// Definition of the Gauss point process type
    using GPProcessType = IntegrationValuesExtrapolationToNodesProcess;

    /// Definition of the Gauss point process pointer type
    using GPProcessPointerType = GPProcessType::UniquePointer;

    /// Pointer definition of EnSightOutput
    KRATOS_CLASS_POINTER_DEFINITION(EnSightOutput);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor for EnSightOutput
     * @param rModelPart: The model part to output
     * @param ThisParameters: The parameters for the output
     */
    explicit EnSightOutput(
        ModelPart& rModelPart,
        Parameters ThisParameters = Parameters(R"({})" )
        );

    /// Destructor.
    ~EnSightOutput() override = default;

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief Get the default parameters for the EnSightOutput
    * @return Parameters: The default parameters for the EnSightOutput
    * @details The default parameters include the model part name, file format, output precision,
    */
    static Parameters GetDefaultParameters();

    /**
     * @brief Print the output for the current simulation step
     * @details This method writes the geometry, nodal variables, and elemental variables to files.
     * It handles both historical and non-historical data, and supports transient simulations.
     * @param rOutputFilename: The base name for the output files
     */
    void PrintOutput(const std::string& rOutputFilename = "");

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "EnSightOutput";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "EnSightOutput";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "EnSightOutput Data" << "\n";
        rOStream << "Model Part Name: " << mOutputSettings["model_part_name"].GetString() << "\n";
        rOStream << "Output Settings: " << mOutputSettings << std::endl;
    }

    ///@}
protected:
    ///@name Enum's
    ///@{

    /**
     * @enum EnSightFileFormat
     * @brief Specifies the supported EnSight file formats.
     * @details This enumeration defines the available EnSight file formats for output.
     * - EnSight5: The classic EnSight 5 format.
     * - EnSight6: The classic EnSight 6 format.
     * - EnSightGold: The EnSight Gold format, which supports additional features.
     */
    enum class EnSightFileFormat {
        EnSight5, // TODO: To be added. Very similar to EnSight6
        EnSight6,
        EnSightGold
    };

    /**
     * @enum FileFormat
     * @brief Enumeration for file format types, like ASCII or BINARY.
     * @details This enum is used to specify the format in which the output files will be written.
     */
    enum class FileFormat {
        ASCII,
        BINARY
    };

    /**
     * @enum EntityType
     * @brief Enumeration for entity types used in EnSight output.
     * @details Specifies which entities are processed for output.
     * - UNDEFINED: Entity type not specified.
     * - NODE: Nodal data.
     * - ELEMENT: Elemental data.
     * - CONDITION: Condition data (boundary conditions).
     * - AUTOMATIC: Automatically determine entity type based on ModelPart contents.
     * - NONE: No entity output.
     * - GEOMETRY: Geometrical objects (future extension).
     */
    enum class EntityType {
        UNDEFINED,
        NODE,
        ELEMENT,
        CONDITION,
        AUTOMATIC,
        NONE,
        GEOMETRY // TODO: Implement support for geometry entities
    };

    /**
     * @enum VariableType
     * @brief Enumeration for variable types supported in EnSight 6/Gold output.
     * @details Used to classify variables as scalar, vector, or tensor (symmetric/asymmetric), which determines their output format.
     * - SCALAR: Single value per node/element.
     * - VECTOR: Three values per node/element (e.g., displacement, velocity).
     * - VECTOR_UNDEFINED: Unknown size vector type.
     * - TENSOR_SYMMETRIC: Six values per node/element (symmetric tensor, e.g., stress).
     * - TENSOR_ASYMMETRIC: Nine values per node/element (full tensor).
     * - TENSOR_UNDEFINED: Unknown or unsupported tensor type. We need to manually check the size of assume a certain distribution.
     * - COMPLEX_SCALAR: Two values per node/element (real and imaginary parts). // NOTE: Currently unsupported
     * - COMPLEX_VECTOR: Six values per node/element (real and imaginary parts for each component). // NOTE: Currently unsupported
     * - UNKNOWN: Unknown variable type.
     */
    enum class VariableType {
        SCALAR,
        VECTOR,
        VECTOR_UNDEFINED,
        TENSOR_SYMMETRIC,
        TENSOR_ASYMMETRIC,
        TENSOR_UNDEFINED,
        COMPLEX_SCALAR,
        COMPLEX_VECTOR,
        UNKNOWN
    };

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;                                         /// Reference to the model part containing the simulation data
    Parameters mOutputSettings;                                     /// Parameters for output settings, including file names, variable types, and output intervals
    unsigned int mDefaultPrecision;                                 /// The default precision
    unsigned int mStepLabelPrecision;                               /// The precision for the step label in filenames
    EnSightFileFormat mEnSightFileFormat;                           /// The EnSight file format (6 or Gold)
    FileFormat mFileFormat;                                         /// The file format for output (ASCII or BINARY)
    EntityType mEntityType;                                         /// The entities to be written (nodes, elements, conditions, etc.)
    GPProcessPointerType mpGaussToNodesProcess;                     /// Process for extrapolating integration values to nodes
    std::vector<double> mTimeValues;                                /// Vector to store time values for transient simulations
    std::vector<PartData> mPartDatas;                               /// Vector to store data for each part (TODO: Add a method to check that the model parts are not changed during time steps, and avoid manually defined)
    std::unordered_map<std::string, VariableType> mVariableTypeMap; /// Map to store variable names and their types for quick access

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Helper to determine which entities to write
     * @param rModelPart The ModelPart that is currently outputted (can be a SubModelPart)
     */
    EntityType GetEntityType(const ModelPart& rModelPart) const;

    /**
     * @brief Prepares the variable type map for quick access.
     * @details Populates the mVariableTypeMap with variable names and their corresponding types.
     */
    void PrepareVariableTypeMap();

    /**
     * @brief Prepares the extrapolation of integration point results to nodes.
     * @details Initializes or updates the process that transfers results from Gauss points to nodal values, enabling output of element-based variables at the nodes for visualization.
     */
    void PrepareGaussPointResults();

    /**
     * @brief Writes the EnSight 6/Gold case file.
     * @details Generates the main case file (.case or .encas) that references geometry and variable files, and contains metadata such as format type, file names, and time step information.
     */
    void WriteCaseFile();

    /**
     * @brief Writes the geometry file for EnSight 6/Gold output.
     * @param rFileName The name of the geometry file to write.
     * @details Outputs mesh node coordinates and element connectivity in the EnSight 6/Gold geometry format, supporting both ASCII and binary formats.
     */
    void WriteGeometryFile(const std::string& rFileName);

    /**
     * @brief Writes a nodal variable to an EnSight 6/Gold variable file.
     * @param rVariableName The name of the variable to output.
     * @param IsHistorical Whether the variable is stored historically (true) or non-historically (false).
     * @param rFileName The name of the output variable file.
     * @details Supports scalar, vector, and tensor variables, and writes data for each node in the mesh.
     */
    void WriteNodalVariableToFile(
        const std::string& rVariableName,
        const bool IsHistorical,
        const std::string& rFileName
        );

    /**
     * @brief Writes a nodal flag variable to an EnSight 6/Gold variable file.
     * @param rFlagName The name of the flag to output.
     * @param rFileName The name of the output variable file.
     */
    void WriteNodalFlagToFile(
        const std::string& rFlagName,
        const std::string& rFileName
        );

    /**
     * @brief Writes an elemental/conditional variable to an EnSight 6/Gold variable file.
     * @details Supports scalar, vector, and tensor variables, and writes data for each element in the mesh.
     * @param rVariableName The name of the variable to output.
     * @param rFileName The name of the output variable file.
     * @param IsElement Whether the variable is associated with an element (true) or a node (false).
     */
    void WriteGeometricalVariableToFile(
        const std::string& rVariableName,
        const std::string& rFileName,
        const bool IsElement = true
        );

    /**
     * @brief Writes an elemental/conditional flag variable to an EnSight 6/Gold variable file.
     * @param rFlagName The name of the flag to output.
     * @param rFileName The name of the output variable file.
     * @param IsElement Whether the variable is associated with an element (true) or a node (false).
     */
    void WriteGeometricalFlagToFile(
        const std::string& rFlagName,
        const std::string& rFileName,
        const bool IsElement = true
        );

    /**
     * @brief Writes an elemental/conditional Gauss point variable to an EnSight 6/Gold variable file.
     * @details Supports scalar, vector, and tensor variables, and writes data for each element in the mesh.
     * @param rVariableName The name of the variable to output.
     * @param rFileName The name of the output variable file.
     * @param IsElement Whether the variable is associated with an element (true) or a node (false).
     */
    void WriteGeometricalGaussVariableToFile(
        const std::string& rVariableName,
        const std::string& rFileName,
        const bool IsElement = true
        );

    /**
     * @brief Collects data for a specific part (SubModelPart) of the model.
     * @details Gathers nodes and elements for a part, mapping IDs and organizing elements by type for output.
     * @param rSubModelPart The submodel part to collect data from.
     * @param rPartData Output structure to store collected data for the part.
     */
    void CollectPartData(
        const ModelPart& rSubModelPart,
        PartData& rPartData
        ) const;

    /**
     * @brief Generates the output file name for a given label and extension.
     * @details Combines output path, label, and extension to produce a valid file name for EnSight output.
     * @param rFileLabel The label describing the file (e.g., variable name or geometry).
     * @param rFileExtension The file extension (e.g., ".geo", ".scl").
     * @return The constructed output file name as a string.
     */
    std::string GetOutputFileName(
        const std::string& rFileLabel,
        const std::string& rFileExtension
        ) const;

    /**
     * @brief Writes scalar data to an output file stream.
     * @details Outputs scalar values for nodes or elements in the required EnSight 6/Gold format.
     * @tparam TData The type of the data container (e.g., vector of doubles).
     * @param rFileStream The output file stream.
     * @param rData The scalar data to write.
     * @param EndOfLine Whether to add a newline at the end (default is true).
     * @param AddInitialTabulation Whether to add a tabulation before the data (default is false).
     * @param AddEndTabulation Whether to add a tabulation after the data (default is false).
     */
    template <typename TData>
    void WriteScalarData(
        std::ofstream& rFileStream,
        const TData& rData,
        const bool EndOfLine = true,
        const bool AddInitialTabulation = false,
        const bool AddEndTabulation = false
        ) const;

    /**
     * @brief Writes vector data to an output file stream.
     * @details Outputs vector values (e.g., 3D vectors) for nodes or elements in the required EnSight 6/Gold format.
     * @tparam TData The type of the data container (e.g., vector of arrays).
     * @param rFileStream The output file stream.
     * @param rData The vector data to write.
     * @param EndOfLine Whether to add a newline at the end (default is true).
     * @param AddInitialTabulation Whether to add a tabulation before the data (default is false).
     * @param AddEndTabulation Whether to add a tabulation after the data (default is false).
     */
    template <typename TData>
    void WriteVectorData(
        std::ofstream& rFileStream,
        const TData& rData,
        const bool EndOfLine = true,
        const bool AddInitialTabulation = false,
        const bool AddEndTabulation = false
        ) const;

    /**
     * @brief Writes symmetric tensor data to an output file stream.
     * @details Outputs symmetric tensor values (e.g., stress or strain tensors with 6 components) for nodes or elements in the required EnSight 6/Gold format.
     * @tparam TData The type of the data container (e.g., vector of arrays or matrices).
     * @param rFileStream The output file stream.
     * @param rData The symmetric tensor data to write.
     * @param EndOfLine Whether to add a newline at the end (default is true).
     * @param AddInitialTabulation Whether to add a tabulation before the data (default is false).
     * @param AddEndTabulation Whether to add a tabulation after the data (default is false).
     */
    template <typename TData>
    void WriteSymmetricTensorData(
        std::ofstream& rFileStream,
        const TData& rData,
        const bool EndOfLine = true,
        const bool AddInitialTabulation = false,
        const bool AddEndTabulation = false
        ) const;

    /**
     * @brief Writes a string to an output file stream.
     * @details Helper function for writing text lines to EnSight 6/Gold files.
     * @param rFileStream The output file stream.
     * @param rString The string to write.
     * @param EndOfLine Whether to add a newline at the end (default is true).
     * @param AddInitialTabulation Whether to add a tabulation before the data (default is false).
     * @param AddEndTabulation Whether to add a tabulation after the data (default is false).
     */
    void WriteString(
        std::ofstream& rFileStream,
        const std::string& rString,
        const bool EndOfLine = true,
        const bool AddInitialTabulation = false,
        const bool AddEndTabulation = false
        ) const;

    /**
     * @brief Gets the EnSight name for a given geometry type
     * @param rGeometryType The geometry type to get the name for
     * @return The EnSight name as a string
     */
    std::string GetEnSightName(const GeometryData::KratosGeometryType& rGeometryType) const;

    /**
     * @brief Updates the part data for all model parts
     */
    void UpdatePartData();

    /**
     * @brief Determines the EnSight variable type for a given variable name.
     * @details This method inspects the variable (by name) in the provided ModelPart (or the main ModelPart if not specified),
     * and returns the corresponding VariableType enum value for EnSight output. It considers the entity type (node, element, etc.)
     * and whether the variable is historical or non-historical. The result is used to select the correct output format (scalar, vector, tensor).
     * @param rVariableName The name of the variable to check.
     * @param pModelPart Optional pointer to the ModelPart to inspect (defaults to main ModelPart).
     * @param Entity The entity type (node, element, etc.) for which the variable is defined.
     * @param IsHistorical Whether the variable is historical (true) or non-historical (false).
     * @return The corresponding VariableType enum value.
     */
    VariableType GetVariableType(
        const std::string& rVariableName,
        ModelPart* pModelPart = nullptr,
        const EntityType Entity = EntityType::UNDEFINED,
        const bool IsHistorical = true
        ) const;

    /**
     * @brief Returns the file extension corresponding to the given variable type for EnSight output.
     * @details This function maps a VariableType to its associated file extension used in EnSight format.
     * The mapping is as follows:
     * - SCALAR: ".scl"
     * - VECTOR: ".vec"
     * - VECTOR_UNDEFINED: (TODO: to be verified)
     * - TENSOR_SYMMETRIC: ".ten"
     * - TENSOR_ASYMMETRIC: ".ten" (TODO: to be verified)
     * - TENSOR_UNDEFINED: ".ten" (TODO: to be verified)
     * - COMPLEX_SCALAR: ".cplx" (TODO: currently unsupported)
     * - COMPLEX_VECTOR: ".cplx" (TODO: currently unsupported)
     * @param Type The variable type for which the file extension is requested.
     * @return The corresponding file extension as a std::string.
     * @throws If the variable type is unknown, an error is thrown.
     */
    std::string GetExtensionFile(const VariableType Type) const;

    /**
     * @brief Returns the type label corresponding to the given variable type for EnSight output.
     * @details This function maps a VariableType to its associated type label used in EnSight format.
     * The mapping is as follows:
     * - SCALAR: scalar
     * - VECTOR: vector
     * - VECTOR_UNDEFINED: (TODO: to be verified)
     * - TENSOR_SYMMETRIC: symmetric tensor
     * - TENSOR_ASYMMETRIC: asymmetric tensor (TODO: to be verified)
     * - TENSOR_UNDEFINED: undefined tensor (TODO: to be verified)
     * - COMPLEX_SCALAR: complex (TODO: currently unsupported)
     * - COMPLEX_VECTOR: complex vector (TODO: currently unsupported)
     * @param Type The variable type for which the file extension is requested.
     * @return The corresponding type label
     * @throws If the variable type is unknown, an error is thrown.
     */
    std::string GetTypeLabel(const VariableType Type) const;

    /**
     * @brief Computes the bounding box of the given ModelPart, optionally scaled by a coefficient.
     * @details This function calculates the axis-aligned bounding box that encloses all points in the specified ModelPart.
     * The resulting bounding box can be scaled by the provided coefficient.
     * @param rModelPart The ModelPart for which the bounding box is to be computed.
     * @param Coefficient A scaling factor to apply to the computed bounding box.
     * @return BoundingBox<Point> The computed (and possibly scaled) bounding box.
     */
    BoundingBox<Point> ComputeBoundingBox(
        const ModelPart& rModelPart,
        const double Coefficient
        );

    /**
     * @brief Retrieves the connectivity information for a given geometrical object.
     * @details This function fills the provided vector with the connectivity indices that define how the nodes of the given geometrical object are connected.
     * @param rGeometricalObject The geometrical object whose connectivity is to be retrieved.
     * @param rConnectivity Output vector that will be populated with the connectivity indices.
     */
    void GetGeometryConnectivity(
        const GeometricalObject& rGeometricalObject,
        std::vector<std::size_t>& rConnectivity
        );

    /**
     * @brief Computes the average value of a variable over all integration points of a given element.
     * @details This function calculates the average of the specified variable evaluated at each integration point of the provided element.
     * @tparam TData The data type of the variable to be averaged.
     * @param rElement The element over which the average is computed.
     * @param rVariable The variable to be averaged.
     * @return The average value of the variable over all integration points of the element.
     */
    template<typename TData>
    TData GetAverageIntegrationValue(
        const Element& rElement,
        const Variable<TData>& rVariable
        );

    ///@}
};// Class EnSightOutput

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                EnSightOutput& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const EnSightOutput& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos