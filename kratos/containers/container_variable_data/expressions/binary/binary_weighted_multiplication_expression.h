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
#include <string>

// Project includes
#include "containers/container_variable_data/expressions/expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Expression to compute weighted multiplication.
 *
 */
class BinaryWeightedMultiplicationExpression : public Expression {
public:
    ///@name Life cycle
    ///@{

    BinaryWeightedMultiplicationExpression(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    BinaryWeightedMultiplicationExpression(const BinaryWeightedMultiplicationExpression& rOther) = delete;

    ~BinaryWeightedMultiplicationExpression() override = default;

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(
        Expression::Pointer pLeft,
        Expression::Pointer pRight);

    double Evaluate(
        const IndexType EntityDataBeginIndex,
        const IndexType ComponentIndex) const override;

    const std::vector<IndexType> GetShape() const override;

    std::string Info() const override;

    ///@}
protected:
    ///@name Private member variables
    ///@{

    const Expression::Pointer mpLeft;

    const Expression::Pointer mpRight;

    const IndexType mLocalSize;

    ///@}
};

} // namespace Kratos