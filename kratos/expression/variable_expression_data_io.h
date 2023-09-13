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
#include <vector>
#include <type_traits>

// Project includes
#include "expression/expression.h"
#include "expression/literal_flat_expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

// Dummy generic class to be specialized later

template<class TDataType>
class KRATOS_API(KRATOS_CORE) VariableExpressionDataIO
{
public:
    ///@name Type definitions
    ///@{

    using Pointer = std::shared_ptr<VariableExpressionDataIO<TDataType>>;

    using IndexType = std::size_t;

    using RawType = std::conditional_t<std::is_same_v<TDataType, int>, int, double>;

    using RawLiteralFlatExpression = LiteralFlatExpression<RawType>;

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
        RawLiteralFlatExpression& rExpression,
        const IndexType EntityIndex,
        const TDataType& Value) const;

    const std::vector<IndexType> GetItemShape() const { return mShape; }

    ///@}
private:
    ///@name Private members
    ///@{

    std::vector<IndexType> mShape;

    ///@}
};

} // namespace Kratos