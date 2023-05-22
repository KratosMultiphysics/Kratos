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

// External includes

// Project includes
#include "includes/io.h"
#include "includes/model_part.h"
#include "containers/variable.h"
#include "containers/flags.h"
#include "containers/container_expression/container_expression.h"


namespace Kratos {
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
    ///@name Public classes
    ///@{

    struct VariableComparator
    {
        bool operator()(
            const SupportedVariables& rFirst,
            const SupportedVariables& rSecond) const
        {
            return std::visit([](auto pFirst, auto pSecond) {
                return pFirst->Key() < pSecond->Key();
            }, rFirst, rSecond);
        }
    };

    enum WriterFormat
    {
        ASCII,
        BINARY
    };

    ///@}
    ///@name Life cycle
    ///@{

    VtuOutput(
        ModelPart& rModelPart,
        const bool IsInitialConfiguration = true,
        const WriterFormat OutputFormat = WriterFormat::BINARY,
        const IndexType Precision = 9);

    ///@}
    ///@name Public operations
    ///@{

    template<class TDataType>
    void AddHistoricalVariable(const Variable<TDataType>& rVariable);

    template<class TDataType>
    void AddNonHistoricalVariable(
        const Variable<TDataType>& rVariable,
        const Flags& rEntityFlags);

    void AddFlagVariable(
        const std::string& rFlagName,
        const Flags& rFlagVariable,
        const Flags& rEntityFlags);

    template <class TContainerType>
    void AddContainerExpression(
        const std::string& rExpressionName,
        const typename ContainerExpression<TContainerType>::Pointer pContainerExpression);

    void ClearHistoricalVariables();

    void ClearNodalNonHistoricalVariables();

    void ClearCellNonHistoricalVariables();

    void ClearNodalFlags();

    void ClearCellFlags();

    void ClearNodalContainerExpressions();

    void ClearCellContainerExpressions();

    const ModelPart& GetModelPart() const;

    void PrintOutput(const std::string& rOutputFilenamePrefix);

    ///@}

private:
    ///@name Private member variables
    ///@{

    ModelPart& mrModelPart;

    const bool mIsInitialConfiguration;

    const WriterFormat mOutputFormat;

    const IndexType mPrecision;

    bool mIsConditionsConsidered;

    bool mIsElementsConsidered;

    std::unordered_map<IndexType, IndexType> mKratosVtuIndicesMap;

    std::unordered_map<std::string, SupportedVariables> mHistoricalVariablesMap;

    std::unordered_map<std::string, SupportedVariables> mNonHistoricalNodalVariablesMap;

    std::unordered_map<std::string, SupportedVariables> mNonHistoricalCellVariablesMap;

    std::unordered_map<std::string, const Flags*> mNodalFlagsMap;

    std::unordered_map<std::string, const Flags*> mCellFlagsMap;

    std::unordered_map<std::string, ContainerExpression<ModelPart::NodesContainerType>::Pointer> mPointContainerExpressionsMap;

    std::unordered_map<std::string, SupportedCellContainerExpressions> mCellContainerExpressionsMap;

    ///@}
    ///@name Private operations
    ///@{

    void WriteModelPart(
        const std::string& rOutputFileNamePrefix,
        ModelPart& rModelPart) const;

    ///@}
};
} // namespace Kratos