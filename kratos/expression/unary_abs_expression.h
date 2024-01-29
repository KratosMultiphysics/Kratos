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
 * @brief Unary abs expression used to get the absolute values of a given input expression.
 */
class KRATOS_API(KRATOS_CORE) UnaryAbsExpression : public Expression {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Life cycle
    ///@{

    UnaryAbsExpression(Expression::ConstPointer pExpression);

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(Expression::ConstPointer pExpression);

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

    ///@}
};

} // namespace Kratos