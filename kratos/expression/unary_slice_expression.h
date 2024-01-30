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

// Project includes
#include "expression/expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Unary slice expression used to represent a slicing lazy expression of a given input expression.
 *
 * @details This expression slices the input expression's entity values with an offset and stride (length of components).
 *
 */
class KRATOS_API(KRATOS_CORE) UnarySliceExpression : public Expression {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Life cycle
    ///@{

    UnarySliceExpression(
        Expression::ConstPointer pExpression,
        const IndexType Offset,
        const IndexType Stride);

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(
        Expression::ConstPointer pExpression,
        const IndexType Offset,
        const IndexType Stride);

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType EntityDataBeginIndex,
        const IndexType ComponentIndex) const override;

    const std::vector<IndexType> GetItemShape() const override;

    IndexType GetMaxDepth() const override;

    std::string Info() const override;

    ///@}
protected:
    ///@name Private member variables
    ///@{

    const Expression::ConstPointer mpSourceExpression;

    const IndexType mOffset;

    const IndexType mStride;

    const IndexType mSourceStride;

    ///@}
};

} // namespace Kratos