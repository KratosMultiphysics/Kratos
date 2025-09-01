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
        const bool                                    UsePointsForDataFieldOutput;
        ModelPart * const                             mpModelPart;
        const ModelPart::NodesContainerType::Pointer  mpPoints;
        const std::optional<CellContainerPointerType> mpCells;
        std::map<std::string, FieldPointerType>       mPointFields;
        std::map<std::string, FieldPointerType>       mCellFields;
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
        const IndexType Precision = 9,
        const bool OutputSubModelParts = false,
        const IndexType EchoLevel = 0);

    ///@}
    ///@name Public operations
    ///@{

    void AddFlag(
        const std::string& rFlagName,
        const Flags& rFlagVariable,
        Globals::DataLocation DataLocation);

    void AddVariable(
        SupportedVariablePointerType pVariable,
        Globals::DataLocation DataLocation);

    void AddIntegrationPointVariable(
        SupportedVariablePointerType pVariable,
        Globals::DataLocation DataLocation);

    void AddContainerExpression(
        const std::string& rExpressionName,
        SupportedContainerExpressionPointerType pContainerExpression);

    void AddTensorAdaptor(
        const std::string& rTensorAdaptorName,
        SupportedTensorAdaptorPointerType pTensorAdaptor);

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

    ModelPart& mrModelPart; /// Reference to the model part.

    const bool mIsInitialConfiguration; /// Flag indicating if it is the initial configuration.

    const IndexType mEchoLevel;

    const WriterFormat mOutputFormat; /// The output format for writing the model part.

    const IndexType mPrecision; /// The precision used for writing floating-point values.

    DataMap<Flags const *> mFlags;

    DataMap<SupportedVariablePointerType> mVariables;

    DataMap<SupportedVariablePointerType> mIntegrationPointVariables;

    std::vector<UnstructuredGridData> mListOfUnstructuredGridData;

    std::vector<std::pair<IndexType, double>> mStepInfo;

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
