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

// System includes
#include <cmath>
#include <vector>

// Project includes
#include "includes/define.h"

// Application includes

// Include base h
#include "expressions.h"

namespace Kratos {

void intrusive_ptr_add_ref(const Expression* x)
{
    x->mReferenceCounter.fetch_add(1, std::memory_order_relaxed);
}

void intrusive_ptr_release(const Expression* x)
{
    if (x->mReferenceCounter.fetch_sub(1, std::memory_order_release) == 1) {
        std::atomic_thread_fence(std::memory_order_acquire);
        delete x;
    }
}

LiteralDoubleExpression::LiteralDoubleExpression(const double Value)
    : mValue(Value)
{
}

Expression::Pointer LiteralDoubleExpression::Create(const double Value)
{
    return Kratos::make_intrusive<LiteralDoubleExpression>(Value);
}

double LiteralDoubleExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return mValue;
}

std::size_t LiteralDoubleExpression::GetDimension() const
{
    return 1;
}

LiteralArray3Expression::LiteralArray3Expression(
    const array_1d<double, 3>& Value,
    const IndexType Dimension)
    : mValue(Value),
      mDimension(Dimension)
{
}

Expression::Pointer LiteralArray3Expression::Create(
    const array_1d<double, 3>& Value,
    const IndexType Dimension)
{
    return Kratos::make_intrusive<LiteralArray3Expression>(Value, Dimension);
}

double LiteralArray3Expression::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return mValue[ComponentIndex];
}

std::size_t LiteralArray3Expression::GetDimension() const
{
    return mDimension;
}

LiteralVectorExpression::LiteralVectorExpression(
    Kratos::shared_ptr<const Vector> pValue,
    const IndexType Dimension)
    : mpValue(pValue),
      mDimension(Dimension)
{
}

Expression::Pointer LiteralVectorExpression::Create(
    Kratos::shared_ptr<Vector> pValue,
    const IndexType Dimension)
{
    return Kratos::make_intrusive<LiteralVectorExpression>(pValue, Dimension);
}

double LiteralVectorExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return (*mpValue)[EntityIndex * mDimension + ComponentIndex];
}

std::size_t LiteralVectorExpression::GetDimension() const
{
    return mDimension;
}

BinaryAddExpression::BinaryAddExpression(
    Expression::Pointer pLeft,
    Expression::Pointer pRight)
    : mpLeft(pLeft),
      mpRight(pRight)
{
    KRATOS_ERROR_IF_NOT(mpLeft->GetDimension() == mpRight->GetDimension() || mpRight->GetDimension() == 1)
        << "Addition should have equal dimensions in left and right side expressions or right hand side should be with dimension = 1. [ Left expresion dimension = "
        << mpLeft->GetDimension() << ", Right expresion dimensions = " << mpRight->GetDimension() << " ].\n";

    KRATOS_ERROR_IF_NOT(mpLeft.get())
        << "Addition is provided with uninitialized left hand side "
           "expression.\n";

    KRATOS_ERROR_IF_NOT(mpRight.get())
        << "Addition is provided with uninitialized right hand side "
           "expression.\n";
}

Expression::Pointer BinaryAddExpression::Create(
    Expression::Pointer pLeft,
    Expression::Pointer pRight)
{
    return Kratos::make_intrusive<BinaryAddExpression>(pLeft, pRight);
}

double BinaryAddExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return mpLeft->Evaluate(EntityIndex, ComponentIndex) + mpRight->Evaluate(EntityIndex, ComponentIndex);
}

std::size_t BinaryAddExpression::GetDimension() const
{
    return this->mpLeft->GetDimension();
}

BinarySubstractExpression::BinarySubstractExpression(
    Expression::Pointer pLeft,
    Expression::Pointer pRight)
    : mpLeft(pLeft),
      mpRight(pRight)
{
    KRATOS_ERROR_IF_NOT(mpLeft->GetDimension() == mpRight->GetDimension() || mpRight->GetDimension() == 1)
        << "Substraction should have equal dimensions in left and right side expressions or right hand side should be with dimension = 1. [ Left expresion dimension = "
        << mpLeft->GetDimension() << ", Right expresion dimensions = " << mpRight->GetDimension() << " ].\n";

    KRATOS_ERROR_IF_NOT(mpLeft.get())
        << "Substraction is provided with uninitialized left hand side "
           "expression.\n";

    KRATOS_ERROR_IF_NOT(mpRight.get())
        << "Substraction is provided with uninitialized right hand side "
           "expression.\n";
}

Expression::Pointer BinarySubstractExpression::Create(
    Expression::Pointer pLeft,
    Expression::Pointer pRight)
{
    return Kratos::make_intrusive<BinarySubstractExpression>(pLeft, pRight);
}

double BinarySubstractExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return mpLeft->Evaluate(EntityIndex, ComponentIndex) - mpRight->Evaluate(EntityIndex, ComponentIndex);
}

std::size_t BinarySubstractExpression::GetDimension() const
{
    return this->mpLeft->GetDimension();
}

BinaryMultiplyExpression::BinaryMultiplyExpression(
    Expression::Pointer pLeft,
    Expression::Pointer pRight)
    : mpLeft(pLeft),
      mpRight(pRight)
{
    KRATOS_ERROR_IF_NOT(mpLeft->GetDimension() == mpRight->GetDimension() || mpRight->GetDimension() == 1)
        << "Multiplication should have equal dimensions in left and right side expressions or right hand side should be with dimension = 1. [ Left expresion dimension = "
        << mpLeft->GetDimension() << ", Right expresion dimensions = " << mpRight->GetDimension() << " ].\n";

    KRATOS_ERROR_IF_NOT(mpLeft.get())
        << "Multiplication is provided with uninitialized left hand side "
           "expression.\n";

    KRATOS_ERROR_IF_NOT(mpRight.get())
        << "Multiplication is provided with uninitialized right hand side "
           "expression.\n";
}

Expression::Pointer BinaryMultiplyExpression::Create(
    Expression::Pointer pLeft,
    Expression::Pointer pRight)
{
    return Kratos::make_intrusive<BinaryMultiplyExpression>(pLeft, pRight);
}

double BinaryMultiplyExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return mpLeft->Evaluate(EntityIndex, ComponentIndex) * mpRight->Evaluate(EntityIndex, ComponentIndex % this->mpRight->GetDimension());
}

std::size_t BinaryMultiplyExpression::GetDimension() const
{
    return this->mpLeft->GetDimension();
}

BinaryDivideExpression::BinaryDivideExpression(
    Expression::Pointer pLeft,
    Expression::Pointer pRight)
    : mpLeft(pLeft),
      mpRight(pRight)
{
    KRATOS_ERROR_IF_NOT(mpLeft->GetDimension() == mpRight->GetDimension() || mpRight->GetDimension() == 1)
        << "Division should have equal dimensions in left and right side expressions or right hand side should be with dimension = 1. [ Left expresion dimension = "
        << mpLeft->GetDimension() << ", Right expresion dimensions = " << mpRight->GetDimension() << " ].\n";

    KRATOS_ERROR_IF_NOT(mpLeft.get())
        << "Division is provided with uninitialized left hand side "
           "expression.\n";

    KRATOS_ERROR_IF_NOT(mpRight.get())
        << "Division is provided with uninitialized right hand side "
           "expression.\n";
}

Expression::Pointer BinaryDivideExpression::Create(
    Expression::Pointer pLeft,
    Expression::Pointer pRight)
{
    return Kratos::make_intrusive<BinaryDivideExpression>(pLeft, pRight);
}

double BinaryDivideExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return mpLeft->Evaluate(EntityIndex, ComponentIndex) / mpRight->Evaluate(EntityIndex, ComponentIndex % this->mpRight->GetDimension());
}

std::size_t BinaryDivideExpression::GetDimension() const
{
    return this->mpLeft->GetDimension();
}

BinaryPowerExpression::BinaryPowerExpression(
    Expression::Pointer pLeft,
    Expression::Pointer pRight)
    : mpLeft(pLeft),
      mpRight(pRight)
{
    KRATOS_ERROR_IF_NOT(mpLeft->GetDimension() == mpRight->GetDimension() || mpRight->GetDimension() == 1)
        << "Power should have equal dimensions in left and right side expressions or right hand side should be with dimension = 1. [ Left expresion dimension = "
        << mpLeft->GetDimension() << ", Right expresion dimensions = " << mpRight->GetDimension() << " ].\n";

    KRATOS_ERROR_IF_NOT(mpLeft.get())
        << "Power is provided with uninitialized left hand side "
           "expression.\n";

    KRATOS_ERROR_IF_NOT(mpRight.get())
        << "Power is provided with uninitialized right hand side "
           "expression.\n";
}

Expression::Pointer BinaryPowerExpression::Create(
    Expression::Pointer pLeft,
    Expression::Pointer pRight)
{
    return Kratos::make_intrusive<BinaryPowerExpression>(pLeft, pRight);
}


double BinaryPowerExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return std::pow(mpLeft->Evaluate(EntityIndex, ComponentIndex), mpRight->Evaluate(EntityIndex, ComponentIndex));
}

std::size_t BinaryPowerExpression::GetDimension() const
{
    return this->mpLeft->GetDimension();
}

}