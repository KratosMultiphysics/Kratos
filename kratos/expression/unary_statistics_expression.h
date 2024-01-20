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
#include <limits>

// Project includes
#include "expression/expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

namespace UnaryStatisticOperations
{
    struct Min { static inline constexpr double mInitial = std::numeric_limits<double>::max();  static inline constexpr double Evaluate(const double V1, const double V2) { return std::min(V1, V2); } };
    struct Max { static inline constexpr double mInitial = std::numeric_limits<double>::lowest(); static inline constexpr double Evaluate(const double V1, const double V2) { return std::max(V1, V2); } };
    struct Sum { static inline constexpr double mInitial = 0.0; static inline constexpr double Evaluate(const double V1, const double V2) { return V1 + V2; } };
}

/**
 * @brief Used to create lazy expression to get statistical quantities of entity values.
 *
 */

template<class TOperationType>
class UnaryStatisticsExpression : public Expression {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Life cycle
    ///@{

    UnaryStatisticsExpression(Expression::ConstPointer pExpression);

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

private:
    ///@name Private member variables
    ///@{

    const IndexType mItemComponentCount;

    const Expression::ConstPointer mpSourceExpression;

    ///@}
};

} // namespace Kratos