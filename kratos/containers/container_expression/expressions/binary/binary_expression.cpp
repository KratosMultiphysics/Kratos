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
#include <vector>

// Project includes

// Include base h
#include "binary_expression.h"

namespace Kratos {

namespace ExpressionHelperUtilities {
using IndexType = std::size_t;

std::string GetShape(const std::vector<IndexType>& rShape)
{
    std::stringstream msg;
    msg << "[";
    for (const auto v : rShape) {
        msg << v << ",";
    }
    if (rShape.size() != 0) {
        msg << "\b";
    }

    msg << "]";
    return msg.str();
}

} // namespace ExpressionHelperUtilities

template <class TOperationType>
BinaryExpression<TOperationType>::BinaryExpression(
    Expression::Pointer pLeft,
    Expression::Pointer pRight)
    : mpLeft(pLeft),
      mpRight(pRight)
{
    const auto& r_left_shape = mpLeft->GetShape();
    const auto& r_right_shape = mpRight->GetShape();

    KRATOS_ERROR_IF_NOT(r_left_shape == r_right_shape || r_right_shape.size() == 0)
        << "Binary operation should have equal shape in left and right side "
           "expressions or right hand side should be a expression with shape {}. ["
        << "lhs shape = " << ExpressionHelperUtilities::GetShape(r_left_shape) << ", "
        << "rhs shape = " << ExpressionHelperUtilities::GetShape(r_right_shape) << " ].\n"
        << "Expression:\n"
        << "lhs = " << *mpLeft << "\n"
        << "rhs = " << *mpRight << "\n";

    KRATOS_ERROR_IF_NOT(mpLeft.get())
        << "Binary operation is provided with uninitialized left hand side "
           "expression.\n";

    KRATOS_ERROR_IF_NOT(mpRight.get())
        << "Binary operation is provided with uninitialized right hand side "
           "expression.\n";
}

template <class TOperationType>
Expression::Pointer BinaryExpression<TOperationType>::Create(
    Expression::Pointer pLeft,
    Expression::Pointer pRight)
{
    return Kratos::make_intrusive<BinaryExpression<TOperationType>>(std::move(pLeft), std::move(pRight));
}

template <class TOperationType>
double BinaryExpression<TOperationType>::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return TOperationType::Evaluate(
        mpLeft->Evaluate(EntityIndex, EntityDataBeginIndex, ComponentIndex),
        mpRight->Evaluate(EntityIndex, EntityDataBeginIndex, ComponentIndex));
}

template <class TOperationType>
const std::vector<std::size_t> BinaryExpression<TOperationType>::GetShape() const
{
    return this->mpLeft->GetShape();
}

template <class TOperationType>
std::string BinaryExpression<TOperationType>::Info() const
{
    std::stringstream msg;
    msg << "(" << this->mpLeft->Info();

    if constexpr(std::is_same_v<TOperationType, BinaryOperations::Addition>) {
        msg << "+";
    } else if constexpr(std::is_same_v<TOperationType, BinaryOperations::Substraction>) {
        msg << "-";
    } else if constexpr(std::is_same_v<TOperationType, BinaryOperations::Multiplication>) {
        msg << "*";
    } else if constexpr(std::is_same_v<TOperationType, BinaryOperations::Division>) {
        msg << "/";
    } else if constexpr(std::is_same_v<TOperationType, BinaryOperations::Power>) {
        msg << "^";
    }

    msg << this->mpRight->Info() << ")";
    return msg.str();
}

// template instantiations
template class BinaryExpression<BinaryOperations::Addition>;
template class BinaryExpression<BinaryOperations::Substraction>;
template class BinaryExpression<BinaryOperations::Multiplication>;
template class BinaryExpression<BinaryOperations::Division>;
template class BinaryExpression<BinaryOperations::Power>;

} // namespace Kratos