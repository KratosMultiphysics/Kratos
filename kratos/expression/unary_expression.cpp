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
#include <sstream>
#include <type_traits>

// Project includes

// Include base h
#include "unary_expression.h"

namespace Kratos {

template <class TOperationType>
UnaryExpression<TOperationType>::UnaryExpression(Expression::ConstPointer pExpression)
    : Expression(pExpression->NumberOfEntities()),
      mpSourceExpression(pExpression)
{
}

template <class TOperationType>
Expression::Pointer UnaryExpression<TOperationType>::Create(Expression::ConstPointer pExpression)
{
    return Kratos::make_intrusive<UnaryExpression<TOperationType>>(std::move(pExpression));
}

template <class TOperationType>
double UnaryExpression<TOperationType>::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return TOperationType::Evaluate(mpSourceExpression->Evaluate(EntityIndex, EntityDataBeginIndex, ComponentIndex));
}

template <class TOperationType>
const std::vector<std::size_t> UnaryExpression<TOperationType>::GetItemShape() const
{
    return mpSourceExpression->GetItemShape();
}

template <class TOperationType>
std::size_t UnaryExpression<TOperationType>::GetMaxDepth() const
{
    return mpSourceExpression->GetMaxDepth() + 1;
}

template <class TOperationType>
std::string UnaryExpression<TOperationType>::Info() const
{
    std::stringstream msg;

    if constexpr(std::is_same_v<TOperationType, UnaryOperations::Absolute>) {
        msg << "Abs";
    } else if constexpr(std::is_same_v<TOperationType, UnaryOperations::Logarithmic>) {
        msg << "Log";
    } else {
        static_assert(!std::is_same_v<TOperationType, TOperationType>, "Unsupported unary operation type.");
    }

    msg << "(" << mpSourceExpression << ")";
    return msg.str();
}

// template instantiations
template class UnaryExpression<UnaryOperations::Absolute>;
template class UnaryExpression<UnaryOperations::Logarithm>;

} // namespace Kratos