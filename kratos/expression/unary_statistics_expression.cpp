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

// System includes
#include <type_traits>

// Project includes

// Include base h
#include "unary_statistics_expression.h"

namespace Kratos {

template<class TOperationType>
UnaryStatisticsExpression<TOperationType>::UnaryStatisticsExpression(Expression::ConstPointer pExpression)
    : Expression(pExpression->NumberOfEntities()),
      mItemComponentCount(pExpression->GetItemComponentCount()),
      mpSourceExpression(pExpression)
{
}

template<class TOperationType>
Expression::Pointer UnaryStatisticsExpression<TOperationType>::Create(Expression::ConstPointer pExpression)
{
    return Kratos::make_intrusive<UnaryStatisticsExpression>(std::move(pExpression));
}

template<class TOperationType>
double UnaryStatisticsExpression<TOperationType>::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    const auto source_entity_begin_index = EntityIndex * mItemComponentCount;
    double norm_value = TOperationType::mInitial;
    for (IndexType ComponentIndex = 0; ComponentIndex < mItemComponentCount; ++ComponentIndex) {
        norm_value = TOperationType::Evaluate(norm_value, mpSourceExpression->Evaluate(EntityIndex, source_entity_begin_index, ComponentIndex));
    }
    return norm_value;
}

template<class TOperationType>
const std::vector<std::size_t> UnaryStatisticsExpression<TOperationType>::GetItemShape() const
{
    return {};
}

template<class TOperationType>
std::size_t UnaryStatisticsExpression<TOperationType>::GetMaxDepth() const
{
    return mpSourceExpression->GetMaxDepth() + 1;
}

template<class TOperationType>
std::string UnaryStatisticsExpression<TOperationType>::Info() const
{
    std::stringstream msg;
    if constexpr(std::is_same_v<TOperationType, UnaryStatisticOperations::Min>) {
        msg << "min(" << mpSourceExpression << ")";
    } else if constexpr(std::is_same_v<TOperationType, UnaryStatisticOperations::Max>) {
        msg << "max(" << mpSourceExpression << ")";
    } else if constexpr(std::is_same_v<TOperationType, UnaryStatisticOperations::Sum>) {
        msg << "sum(" << mpSourceExpression << ")";
    } else {
        static_assert(!std::is_same_v<TOperationType, TOperationType>, "Unsupported operation type.");
    }
    return msg.str();
}

// template instantiations
template class UnaryStatisticsExpression<UnaryStatisticOperations::Min>;
template class UnaryStatisticsExpression<UnaryStatisticOperations::Max>;
template class UnaryStatisticsExpression<UnaryStatisticOperations::Sum>;

} // namespace Kratos