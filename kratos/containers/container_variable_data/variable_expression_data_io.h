//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <vector>

// Project includes
#include "containers/container_variable_data/expressions/expression.h"
#include "containers/container_variable_data/expressions/literal/literal_flat_expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

// Dummy generic class to be specialized later

template<class TDataType>
class VariableExpressionDataIO
{
public:
    ///@name Type definitions
    ///@{

    using Pointer = std::shared_ptr<VariableExpressionDataIO<TDataType>>;

    using IndexType = std::size_t;

    ///@}
    ///@name Life cycle
    ///@{

    VariableExpressionDataIO(const TDataType& SampleValue);

    VariableExpressionDataIO(const std::vector<IndexType>& rShape);

    ///@}
    ///@name Public operations
    ///@{

    static Pointer Create(const TDataType& SampleValue);

    static Pointer Create(const std::vector<IndexType>& rShape);

    void Assign(
        TDataType& rOutput,
        const Expression& rExpression,
        const IndexType EntityIndex) const;

    void Read(
        LiteralFlatExpression& rExpression,
        const IndexType EntityIndex,
        const TDataType& Value) const;

    const std::vector<IndexType> GetShape() const { return mShape; }

    ///@}
private:
    ///@name Private members
    ///@{

    std::vector<IndexType> mShape;

    ///@}
};

} // namespace Kratos