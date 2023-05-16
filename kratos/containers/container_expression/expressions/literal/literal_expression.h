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
#include <vector>

// Project includes
#include "containers/container_expression/expressions/expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Expression to hold a literal value
 *
 */
template <class TDataType>
class KRATOS_API(KRATOS_CORE) LiteralExpression : public Expression {
public:
    ///@name Life cycle
    ///@{

    LiteralExpression(const TDataType& Value);

    ///@}
    ///@name Public operations
    ///@{

    static Expression::Pointer Create(const TDataType& Value);

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType EntityDataBeginIndex,
        const IndexType ComponentIndex) const override;

    const std::vector<IndexType> GetShape() const override;

    std::string Info() const override;

    ///@}
private:
    ///@name Private member variables
    ///@{

    const TDataType mValue;

    std::vector<IndexType> mShape;

    ///@}
};

} // namespace Kratos