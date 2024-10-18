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
#include "includes/define.h"
#include "includes/io.h"
#include "includes/model_part.h"
#include "containers/variable.h"
#include "containers/flags.h"
#include "expression/container_expression.h"

namespace Kratos {
/**
 * @class VtuOutput
 * @brief Class to output Kratos Flags, Variables and ContainerExpressions
 *        to vtu. Supports both shared and distributed memory architectures.
 *
 * @details This class does not create or destroy any folder structures, hence the output
 * file name prefix should have a valid parent directory.
 * @author Suneth Warnakulasuriya
 */
class KRATOS_API(KRATOS_CORE) VtuOutput : public IO
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using SupportedVariables = std::variant<
                                    const Variable<int>*,
                                    const Variable<double>*,
                                    const Variable<array_1d<double, 3>>*,
                                    const Variable<array_1d<double, 4>>*,
                                    const Variable<array_1d<double, 6>>*,
                                    const Variable<array_1d<double, 9>>*>;

    using SupportedCellContainerExpressions = std::variant<
                                                ContainerExpression<ModelPart::ConditionsContainerType>::Pointer,
                                                ContainerExpression<ModelPart::ElementsContainerType>::Pointer>;

    KRATOS_CLASS_POINTER_DEFINITION(VtuOutput);

    KRATOS_DEFINE_LOCAL_FLAG( NODES );
    KRATOS_DEFINE_LOCAL_FLAG( CONDITIONS );
    KRATOS_DEFINE_LOCAL_FLAG( ELEMENTS );

    ///@}
    ///@name Public enums
    ///@{

    /// Enumerations for the output writer format.
    enum WriterFormat
    {
        ASCII,  /// ASCII format.
        BINARY  /// Binary format.
    };

    ///@}
    ///@name Life cycle
    ///@{

    /**
     * @brief Construct a new Vtu Output IO
     * @details Constructs a new VtuOuput IO instance with the given parameters.
     * @param rModelPart                Model part to be used.
     * @param IsInitialConfiguration    If true, the initial configuration is written.
     * @param OutputFormat              Output format. Either ASCII or BINARY supported.
     * @param Precision                 Precision of the double output.
     */
    VtuOutput(
        ModelPart& rModelPart,
        const bool IsInitialConfiguration = true,
        const WriterFormat OutputFormat = WriterFormat::BINARY,
        const IndexType Precision = 9);

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Adds historical variables to the output.
     *
     * @tparam TDataType
     * @param rVariable
     */
    template<class TDataType>
    void AddHistoricalVariable(const Variable<TDataType>& rVariable);

    /**
     * @brief Adds non historical variables to the output.
     *
     * @tparam TDataType
     * @param rVariable         Variable to be added.
     * @param rEntityFlags      Considered container for the variable. Either NODES, CONDITIONS or ELEMENTS
     */
    template<class TDataType>
    void AddNonHistoricalVariable(
        const Variable<TDataType>& rVariable,
        const Flags& rEntityFlags);

    /**
     * @brief Adds flag output.
     *
     * @param rFlagName         Flag name.
     * @param rFlagVariable     Variable to be added.
     * @param rEntityFlags      Considered container for the variable. Either NODES, CONDITIONS or ELEMENTS
     */
    void AddFlagVariable(
        const std::string& rFlagName,
        const Flags& rFlagVariable,
        const Flags& rEntityFlags);

    /**
     * @brief Adds container expressions to the vtu output.
     *
     * This adds container expressions to the output. Proper care should be taken when updating ContainerExpressions because
     * In python, when a container expression is assigned with a new container expression, it does not call the assignment operator.
     * Hence, the new expression takes place. Therefore, when container expressions required to be outputted, then it is best
     * to always clear the existing container expressions and add the new ones. Otherwise, the vtu output may be writing not the
     * latest container expression.
     *
     * @tparam TContainerType
     * @param rExpressionName           Name for the container expression.
     * @param pContainerExpression      Container expression.
     */
    template <class TContainerType>
    void AddContainerExpression(
        const std::string& rExpressionName,
        const typename ContainerExpression<TContainerType>::Pointer pContainerExpression);

    /**
    * @brief Clears the historical variables.
    */
    void ClearHistoricalVariables();

    /**
    * @brief Clears the nodal non-historical variables.
    */
    void ClearNodalNonHistoricalVariables();

    /**
    * @brief Clears the cell non-historical variables.
    */
    void ClearCellNonHistoricalVariables();

    /**
    * @brief Clears the nodal flags.
    */
    void ClearNodalFlags();

    /**
    * @brief Clears the cell flags.
    */
    void ClearCellFlags();

    /**
    * @brief Clears the nodal container expressions.
    */
    void ClearNodalContainerExpressions();

    /**
    * @brief Clears the cell container expressions.
    */
    void ClearCellContainerExpressions();

    /**
    * @brief Returns the model part.
    * @return The constant reference to the model part.
    */
    const ModelPart& GetModelPart() const;

    /**
     * @brief Writes the Vtu file.
     * @details This writes the final vtu file. If this is a transient output, then @ref rOutputFilenamePrefix
     * should have indication of the step, otherwise this will overwrite the same file.
     * In the MPI case, this will create one .vtu file per rank, and a .pvtu file to combine them.
     * @param rOutputFilenamePrefix         Output file name prefix.
     */
    void PrintOutput(const std::string& rOutputFilenamePrefix);

    ///@}
private:
    ///@name Private member variables
    ///@{

    ModelPart& mrModelPart; /// Reference to the model part.

    const bool mIsInitialConfiguration; /// Flag indicating if it is the initial configuration.

    const WriterFormat mOutputFormat; /// The output format for writing the model part.

    const IndexType mPrecision; /// The precision used for writing floating-point values.

    bool mIsConditionsConsidered; /// Flag indicating if conditions are considered.

    bool mIsElementsConsidered; /// Flag indicating if elements are considered.

    // TODO: In the future study to replace the std::map
    // TODO: Study replace string, expensive, with hashes or keys

    std::unordered_map<IndexType, IndexType> mKratosVtuIndicesMap; /// Map to store Kratos VTU indices.

    std::map<std::string, SupportedVariables> mHistoricalVariablesMap; /// Map to store supported historical variables.

    std::map<std::string, SupportedVariables> mNonHistoricalNodalVariablesMap; /// Map to store supported non-historical nodal variables.

    std::map<std::string, SupportedVariables> mNonHistoricalCellVariablesMap; /// Map to store supported non-historical cell variables.

    std::map<std::string, const Flags*> mNodalFlagsMap; /// Map to store nodal flags.

    std::map<std::string, const Flags*> mCellFlagsMap; /// Map to store cell flags.

    std::map<std::string, ContainerExpression<ModelPart::NodesContainerType>::Pointer> mPointContainerExpressionsMap; /// Map to store point container expressions.

    std::map<std::string, SupportedCellContainerExpressions> mCellContainerExpressionsMap; /// Map to store supported cell container expressions.

    ///@}
    ///@name Private operations
    ///@{

    /**
    * @brief Writes the model part to a file.
    * @param rOutputFileNamePrefix The output file name prefix.
    * @param rModelPart            The model part to write.
    */
    void PrintModelPart(
        const std::string& rOutputFileNamePrefix,
        ModelPart& rModelPart) const;

    ///@}
};
} // namespace Kratos
