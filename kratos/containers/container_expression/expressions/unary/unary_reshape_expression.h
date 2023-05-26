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
#include "containers/container_expression/expressions/expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) UnaryReshapeExpression : public Expression {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Life cycle
    ///@{

    UnaryReshapeExpression(
        Expression::Pointer pExpression,
        const std::vector<IndexType>& rShape);

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(
        Expression::Pointer pExpression,
        const std::vector<IndexType>& rShape);

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType EntityDataBeginIndex,
        const IndexType ComponentIndex) const override;

    const std::vector<IndexType> GetItemShape() const override;

    std::string Info() const override;

    ///@}
protected:
    ///@name Private member variables
    ///@{

    const Expression::Pointer mpSourceExpression;

    const std::vector<IndexType> mShape;

    ///@}
};

} // namespace Kratos