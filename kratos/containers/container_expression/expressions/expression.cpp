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
#include <numeric>

// Project includes

// Include base h
#include "expression.h"
#include "containers/container_expression/expressions/literal/literal_expression.h"
#include "containers/container_expression/expressions/binary/binary_expression.h"

namespace Kratos {

Expression::ExpressionIterator::ExpressionIterator()
    : mpExpression(nullptr),
      mEntityIndex(0),
      mEntityDataBeginIndex(0),
      mItemComponentIndex(0),
      mItemComponentCount(0)
{
    KRATOS_ERROR << "The default construction of ExpressionIterator is not allowed.\n";
}

Expression::ExpressionIterator::ExpressionIterator(Expression::Pointer pExpression)
    : mpExpression(pExpression),
      mEntityIndex(0),
      mEntityDataBeginIndex(0),
      mItemComponentIndex(0),
      mItemComponentCount(pExpression->GetItemComponentCount())
{
}

Expression::ExpressionIterator::ExpressionIterator(const ExpressionIterator& rOther)
    : mpExpression(rOther.mpExpression),
      mEntityIndex(rOther.mEntityIndex),
      mEntityDataBeginIndex(rOther.mEntityDataBeginIndex),
      mItemComponentIndex(rOther.mItemComponentIndex),
      mItemComponentCount(rOther.mItemComponentCount)
{
}

Expression::Pointer Expression::ExpressionIterator::GetExpression() const
{
    return mpExpression;
}

double Expression::ExpressionIterator::operator*() const
{
    return mpExpression->Evaluate(mEntityIndex, mEntityDataBeginIndex, mItemComponentIndex);
}

bool Expression::ExpressionIterator::operator==(const ExpressionIterator& rOther) const
{
    return (
        mpExpression.get() == rOther.mpExpression.get() &&
        mpExpression.get() != nullptr &&
        mEntityIndex == rOther.mEntityIndex &&
        mItemComponentIndex == rOther.mItemComponentIndex);
}

bool Expression::ExpressionIterator::operator!=(const ExpressionIterator& rOther) const
{
    return !this->operator==(rOther);
}

Expression::ExpressionIterator& Expression::ExpressionIterator::operator=(const ExpressionIterator& rOther)
{
    mpExpression = rOther.mpExpression;
    mEntityIndex = rOther.mEntityIndex;
    mEntityDataBeginIndex = rOther.mEntityDataBeginIndex;
    mItemComponentIndex = rOther.mItemComponentIndex;
    mItemComponentCount = rOther.mItemComponentCount;
    return *this;
}

Expression::ExpressionIterator& Expression::ExpressionIterator::operator++()
{
    ++mItemComponentIndex;
    if (mItemComponentIndex == mItemComponentCount) {
        mItemComponentIndex = 0;
        ++mEntityIndex;
        mEntityDataBeginIndex = mEntityIndex * mItemComponentCount;
    }
    return *this;
}

Expression::ExpressionIterator Expression::ExpressionIterator::operator++(int)
{
    ExpressionIterator temp = *this;
    ++*this;
    return temp;
}

std::size_t Expression::GetItemComponentCount() const
{
    const auto& r_shape = this->GetItemShape();
    return std::accumulate(
        r_shape.begin(),
        r_shape.end(), 1UL,
        [](const auto V1, const auto V2) { return V1 * V2; });
}

std::size_t Expression::size() const
{
    return this->NumberOfEntities() * this->GetItemComponentCount();
}

Expression::ExpressionIterator Expression::begin() const
{
    return ExpressionIterator(this);
}

Expression::ExpressionIterator Expression::end() const
{
    ExpressionIterator result(this);
    result.mEntityIndex = this->NumberOfEntities();
    result.mEntityDataBeginIndex = result.mEntityIndex * result.mItemComponentCount;
    return result;
}

Expression::ExpressionIterator Expression::cbegin() const
{
    return begin();
}

Expression::ExpressionIterator Expression::cend() const
{
    return end();
}

#define KRATOS_DEFINE_BINARY_EXPRESSION_OPERATOR(OPERATOR_NAME, OPERATOR_CLASS)              \
    Expression::Pointer OPERATOR_NAME(const Expression::Pointer& rpLeft, const double Right) \
    {                                                                                        \
        return BinaryExpression<OPERATOR_CLASS>::Create(                                     \
            rpLeft, LiteralExpression<double>::Create(Right, rpLeft->NumberOfEntities()));   \
    }                                                                                        \
                                                                                             \
    Expression::Pointer OPERATOR_NAME(const double Left, const Expression::Pointer& rpRight) \
    {                                                                                        \
        return BinaryExpression<OPERATOR_CLASS>::Create(                                     \
            LiteralExpression<double>::Create(Left, rpRight->NumberOfEntities()), rpRight);  \
    }                                                                                        \
                                                                                             \
    Expression::Pointer OPERATOR_NAME(const Expression::Pointer& rpLeft,                     \
                                      const Expression::Pointer& rpRight)                    \
    {                                                                                        \
        KRATOS_ERROR_IF_NOT(                                                                 \
            rpLeft->NumberOfEntities() * rpLeft->GetItemComponentCount() ==                  \
            rpRight->NumberOfEntities() * rpRight->GetItemComponentCount())                  \
            << "Operand size mismatch in binary operator: " << #OPERATOR_NAME << "!\n"       \
            << "Left operand: " << *rpLeft << '\n'                                           \
            << "Right operand: " << *rpRight;                                                \
        return BinaryExpression<OPERATOR_CLASS>::Create(rpLeft, rpRight);                    \
    }

KRATOS_DEFINE_BINARY_EXPRESSION_OPERATOR(operator+, BinaryOperations::Addition)

KRATOS_DEFINE_BINARY_EXPRESSION_OPERATOR(operator-, BinaryOperations::Substraction)

KRATOS_DEFINE_BINARY_EXPRESSION_OPERATOR(operator*, BinaryOperations::Multiplication)

KRATOS_DEFINE_BINARY_EXPRESSION_OPERATOR(operator/, BinaryOperations::Division)

KRATOS_DEFINE_BINARY_EXPRESSION_OPERATOR(Pow, BinaryOperations::Power)

#undef KRATOS_DEFINE_BINARY_EXPRESSION_OPERATOR

Expression::Pointer Scale(const Expression::Pointer& rpLeft, const Expression::Pointer& rpRight)
{
    KRATOS_ERROR_IF_NOT(rpLeft->NumberOfEntities() == rpRight->NumberOfEntities())
        << "Operand size mismatch in binary operator: Scale!\n"
        << "Left operand: " << *rpLeft << '\n'
        << "Right operand: " << *rpRight;
    KRATOS_ERROR_IF_NOT(rpLeft->GetItemComponentCount() == rpRight->GetItemComponentCount() || rpRight->GetItemComponentCount() == 1)
        << "Incompatible operand item component sizes in binary operator: Scale!\n"
        << "Left operand: " << *rpLeft << '\n'
        << "Right operand: " << *rpRight;
    return BinaryExpression<BinaryOperations::Multiplication>::Create(
        rpLeft,
        rpRight
    );
}

} // namespace Kratos