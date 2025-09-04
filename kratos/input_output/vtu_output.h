//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <string>
#include <variant>
#include <unordered_map>
#include <map>

// External includes

// Project includes
#include "containers/flags.h"
#include "containers/nd_data.h"
#include "containers/variable.h"
#include "expression/container_expression.h"
#include "includes/define.h"
#include "includes/global_variables.h"
#include "includes/io.h"
#include "includes/model_part.h"
#include "tensor_adaptors/tensor_adaptor.h"

/**
 * @class VtuOutput
 * @brief Handles VTU (Visualization Toolkit Unstructured grid) output for Kratos ModelParts.
 * @author Suneth Warnakulasuriya
 * @ingroup InputOutput
 *
 * This class provides functionality to export simulation data from a Kratos ModelPart to VTU files,
 * supporting both ASCII and binary formats. It allows users to register @ref Variable, @ref Flags, @ref ContainerExpression,
 * and @ref TensorAdaptor for output, and supports writing data on nodes, elements, conditions, and integration points.
 * The output can be configured to include submodel parts and supports parallel execution (MPI).
 *
 * @section TypeDefinitions Type Definitions
 * - IndexType: Alias for std::size_t.
 * - FieldPointerType: Variant type for supported field data pointers.
 * - CellContainerPointerType: Variant type for supported cell container pointers.
 * - SupportedVariablePointerType: Variant type for supported @ref Variable pointers.
 * - SupportedContainerExpressionPointerType: Variant type for supported @ref ContainerExpression pointers.
 * - SupportedTensorAdaptorPointerType: Variant type for supported @ref TensorAdaptor pointers.
 * - IndicesMap: Unordered map for index mapping between Kratos node ids and Vtk point indices.
 * - DataMap: Unordered map for @ref Globals::DataLocation and map of data field name and type of the data field.
 *
 * @section VtuOutput_Enums Enums
 * - WriterFormat: Specifies output format (ASCII or BINARY).
 *
 * @section VtuOutput_Structs Structs
 * - UnstructuredGridData: Holds data for an unstructured grid, including points, cells, and associated fields.
 *
 * @section VtuOutput_PrivateMembers Private Member Variables
 * - mrModelPart: Reference to the ModelPart.
 * - mIsInitialConfiguration: Flag for initial configuration output.
 * - mEchoLevel: Echo level for logging.
 * - mOutputFormat: Output format (ASCII/BINARY).
 * - mPrecision: Precision for ASCII output.
 * - mFlags: Registered flags for output.
 * - mVariables: Registered variables for output.
 * - mIntegrationPointVariables: Registered integration point variables.
 * - mUnstructuredGridDataList: List of unstructured grid data.
 * - mTimeStepList: List of time steps.
 *
 * @section VtuOutput_PrivateOperations Private Operations
 * - WriteUnstructuredGridData: Writes unstructured grid data to file.
 * - WriteIntegrationPointData: Writes integration point data to file.
 */
namespace Kratos {
class KRATOS_API(KRATOS_CORE) VtuOutput : public IO
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using FieldPointerType = std::variant<
                                    NDData<unsigned char>::Pointer,
                                    NDData<bool>::Pointer,
                                    NDData<int>::Pointer,
                                    NDData<double>::Pointer
                                >;

    using CellContainerPointerType = std::variant<
                                        ModelPart::ConditionsContainerType::Pointer,
                                        ModelPart::ElementsContainerType::Pointer
                                    >;

    using SupportedVariablePointerType = std::variant<
                                                Variable<int> const *,
                                                Variable<double> const *,
                                                Variable<array_1d<double, 3>> const *,
                                                Variable<array_1d<double, 4>> const *,
                                                Variable<array_1d<double, 6>> const *,
                                                Variable<array_1d<double, 9>> const *,
                                                Variable<Vector> const *,
                                                Variable<Matrix> const *
                                            >;

    using SupportedContainerExpressionPointerType = std::variant<
                                                        ContainerExpression<ModelPart::NodesContainerType>::Pointer,
                                                        ContainerExpression<ModelPart::ConditionsContainerType>::Pointer,
                                                        ContainerExpression<ModelPart::ElementsContainerType>::Pointer
                                                    >;

    using SupportedTensorAdaptorPointerType = std::variant<
                                                    TensorAdaptor<bool>::Pointer,
                                                    TensorAdaptor<int>::Pointer,
                                                    TensorAdaptor<double>::Pointer
                                                >;

    using IndicesMap = std::unordered_map<IndexType, IndexType>;

    template<class T>
    using DataMap = std::unordered_map<Globals::DataLocation, std::map<std::string, T>>;

    KRATOS_CLASS_POINTER_DEFINITION(VtuOutput);

    ///@}
    ///@name Public enums
    ///@{

    /// Enumerations for the output writer format.
    enum WriterFormat
    {
        ASCII,  /// ASCII format.
        BINARY  /// Binary format.
    };

    struct UnstructuredGridData
    {
        const bool                                    UsePointsForDataFieldOutput;  // If true, then the mpPoints represents a container in the model part, otherwise, mpPoints refers to a container which is on the fly generated.
        ModelPart * const                             mpModelPart;                  // Model part associated with the unstructured grid data.
        const ModelPart::NodesContainerType::Pointer  mpPoints;                     // Points to be used in the unstructured grid.
        const std::optional<CellContainerPointerType> mpCells;                      // Cells to be used in the unstructured grid.
        std::map<std::string, FieldPointerType>       mPointFields;                 // Point data fields such as expressions or tensor adaptors.
        std::map<std::string, FieldPointerType>       mCellFields;                  // Cell data fields such as expressions or tensor adaptors.
    };

    ///@}
    ///@name Life cycle
    ///@{

    /**
     * @brief Construct a new Vtu Output
     *
     * @param rModelPart                Model part to be used in the vtu output.
     * @param IsInitialConfiguration    Which nodal positions to be written out.
     * @param OutputFormat              The format of the output.
     * @param Precision                 Precision of the double values to be used when writing the doubles as ASCII.
     * @param OutputSubModelParts       To consider all the submodel parts recursively for output.
     * @param EchoLevel                 Echo level for to print information.
     */
    VtuOutput(
        ModelPart& rModelPart,
        const bool IsInitialConfiguration = true,
        const WriterFormat OutputFormat = WriterFormat::BINARY,
        const IndexType Precision = 9,
        const bool OutputSubModelParts = false,
        const IndexType EchoLevel = 0);

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Adds a flag to the VTU output.
     *
     * This method registers a flag variable to be included in the VTU output file.
     * The flag will be associated with the specified data location (e.g., node, element, condition).
     * This allows the Flag's values to be written to the output file when @ref PringOutput is called.
     *
     * @throws std::runtime_error if @ref AddFlag is called after @ref PrintOutput.
     * @throws std::runtime_error if tries to add the same @p rFlagName more than ones.
     * @throws std::runtime_error if there already exist a field name same as @p rFlagName in one of the compatible data locations with @p DataLocation .
     * @throws std::runtime_error if @p DataLocation is not NodeNonHistorical, Condition, or Element.
     *
     * @param rFlagName The name of the flag to be added.
     * @param rFlagVariable The flag variable to be added.
     * @param DataLocation The location where the flag data is stored (e.g., NodeNonHistorical, Condition, Element).
     *
     * @see @ref FlagTensorAdaptor
     */
    void AddFlag(
        const std::string& rFlagName,
        const Flags& rFlagVariable,
        Globals::DataLocation DataLocation);

    /**
     * @brief Adds a variable to the output for a specified data location.
     *
     * Registers the given variable to be included in the VTU output at the specified
     * data location (e.g., NodeHistorical, NodeNonHistorical, Condition, or Element). This allows the variable's values
     * to be written to the output file when @ref PringOutput is called.
     *
     * @throws std::runtime_error if @ref AddVariable is called after @ref PrintOutput.
     * @throws std::runtime_error if tries to add the same @p pVariable more than ones.
     * @throws std::runtime_error if there already exist a field name same as @p pVariable in one of the compatible data locations with @p DataLocation .
     * @throws std::runtime_error if @p DataLocation is not NodeHistorical, NodeNonHistorical, Condition, or Element.
     *
     * @param pVariable Pointer to the supported variable to be added.
     * @param DataLocation The location in the model where the variable is defined (e.g., NodeHistorical, NodeNonHistorical, Condition, or Element).
     *
     * @see @ref VariableTensorAdaptor
     * @see @ref HistoricalVariableTensorAdaptor
     */
    void AddVariable(
        SupportedVariablePointerType pVariable,
        Globals::DataLocation DataLocation);

    /**
     * @brief Adds a variable to be output at integration points.
     *
     * This method registers a variable for output at the integration points of the mesh elements or conditions.
     * This allows the entities variable's integration point values to be written to the output file when @ref PringOutput is called.
     *
     * @throws std::runtime_error if @ref AddIntegrationPointVariable is called after @ref PrintOutput.
     * @throws std::runtime_error if tries to add the same @p pVariable more than ones.
     * @throws std::runtime_error if there already exist a field name same as @p pVariable in one of the compatible data locations with @p DataLocation .
     * @throws std::runtime_error if @p DataLocation is not Condition or Element.
     *
     * @param pVariable Pointer to the supported variable to be added for output.
     * @param DataLocation Specifies the location (e.g., Condition or Element) where the variable is used to calculate integration point values.
     */
    void AddIntegrationPointVariable(
        SupportedVariablePointerType pVariable,
        Globals::DataLocation DataLocation);

    /**
     * @brief Adds a container expression to the output with a specified name.
     *
     * This method registers a container expression, identified by a unique name,
     * to be included in the VTU output. This copies the internal data of the container
     * expression, therefore at the @ref PrintOutput call, the values at the point of calling
     * @ref AddContainerExpression will be written.
     *
     * @throws std::runtime_error if @ref AddContainerExpression is called after @ref PrintOutput.
     * @throws std::runtime_error if tries to add the same @p rExpressionName more than ones.
     * @throws std::runtime_error if there already exist a field name same as @p rExpressionName in one of the compatible data locations represented by the container of @p pContainerExpression .
     * @throws std::runtime_error if @p pContainerExpression does not corresponds to any of the containers which is written by this @ref VtuOutput .
     *
     * @param rExpressionName The name to associate with the container expression.
     * @param pContainerExpression Pointer to the container expression to be added.
     */
    void AddContainerExpression(
        const std::string& rExpressionName,
        SupportedContainerExpressionPointerType pContainerExpression);

    /**
     * @brief Adds a tensor adaptor to the output.
     *
     * Registers a @ref TensorAdaptor with the specified name, and copies the internal data of the @ref TensorAdaptor
     * to internal storage. When @ref PrintOutput is called, the values of the @ref TensorAdaptor at the point of addition will be written.
     *
     * @throws std::runtime_error if @ref AddTensorAdaptor is called after @ref PrintOutput.
     * @throws std::runtime_error if tries to add the same @p rTensorAdaptorName more than ones.
     * @throws std::runtime_error if there already exist a field name same as @p rTensorAdaptorName in one of the compatible data locations represented by the container of @p pTensorAdaptor .
     * @throws std::runtime_error if @p pTensorAdaptor does not corresponds to any of the containers which is written by this @ref VtuOutput .
     *
     * @param rTensorAdaptorName The name to associate with the tensor adaptor.
     * @param pTensorAdaptor Pointer to the tensor adaptor to be added.
     */
    void AddTensorAdaptor(
        const std::string& rTensorAdaptorName,
        SupportedTensorAdaptorPointerType pTensorAdaptor);

    /**
     * @brief Updates the container expression associated with the given expression name.
     *
     * This method assigns the provided container expression pointer to the specified expression name.
     * It is typically used to update or replace the expression stored in the data fields of @ref VtuOutput.
     *
     * @throws std::runtime_error if there is no field name same as @p rExpressionName in one of the compatible data locations represented by the container of @p pContainerExpression .
     * @throws std::runtime_error if @p pContainerExpression does not corresponds to any of the containers which is written by this @ref VtuOutput .
     *
     * @param rExpressionName The name of the expression to be updated.
     * @param pContainerExpression Pointer to the container expression to associate with the given name.
     */
    void UpdateContainerExpression(
        const std::string& rExpressionName,
        SupportedContainerExpressionPointerType pContainerExpression);

    /**
     * @brief Updates the tensor adaptor associated with the given name.
     *
     * This method assigns a new tensor adaptor to the specified name, allowing
     * for dynamic changes to the tensor adaptor used in the @ref VtuOutput process.
     *
     * @throws std::runtime_error if there is no field name same as @p rTensorAdaptorName in one of the compatible data locations represented by the container of @p pTensorAdaptor .
     * @throws std::runtime_error if @p pTensorAdaptor does not corresponds to any of the containers which is written by this @ref VtuOutput .
     *
     * @param rTensorAdaptorName The name identifying the tensor adaptor to update.
     * @param pTensorAdaptor Pointer to the new tensor adaptor to be associated.
     */
    void UpdateTensorAdaptor(
        const std::string& rTensorAdaptorName,
        SupportedTensorAdaptorPointerType pTensorAdaptor);

    /**
    * @brief Returns the model part.
    * @return The constant reference to the model part.
    */
    const ModelPart& GetModelPart() const;

    /**
     * @brief Prints the vtu output data to a file with the specified filename prefix.
     *
     * Will create the folder with the provided @p rOutputFilenamePrefix . In this folder, following files will be created.
     * - If @p OutputSubModelParts is true, then it will create one vtu file per timestep, per model part's container (in MPI, per rank as well).
     * - If this is called in MPI, then this will create one pvtu file per timestep, per model part's container (will be having information of all the corresponding rank vtu files)
     *
     * Finally it will create one pvd ( @p rOutputFilenamePrefix.pvd ) file for all the model parts and all the timesteps linking
     * - In MPI the pvtu files created.
     * - In shared memory parallelism the vtu files created.
     *
     * @param rOutputFilenamePrefix The prefix to use for the output filename.
     */
    void PrintOutput(const std::string& rOutputFilenamePrefix);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    ModelPart& mrModelPart;

    const bool mIsInitialConfiguration;

    const IndexType mEchoLevel;

    const WriterFormat mOutputFormat;

    const IndexType mPrecision;

    DataMap<Flags const *> mFlags;

    DataMap<SupportedVariablePointerType> mVariables;

    DataMap<SupportedVariablePointerType> mIntegrationPointVariables;

    std::vector<UnstructuredGridData> mUnstructuredGridDataList;

    std::vector<double> mTimeStepList;

    ///@}
    ///@name Private operations
    ///@{

    template<class TXmlDataElementWrapper>
    std::pair<std::string, std::string> WriteUnstructuredGridData(
        const std::string& rOutputPrefix,
        UnstructuredGridData& rUnstructuredGridData,
        TXmlDataElementWrapper& rXmlDataElementWrapper) const;

    template<class TXmlDataElementWrapper>
    std::pair<std::string, std::string> WriteIntegrationPointData(
        const std::string& rOutputPrefix,
        UnstructuredGridData& rUnstructuredGridData,
        TXmlDataElementWrapper& rXmlDataElementWrapper) const;

    ///@}
};
} // namespace Kratos
