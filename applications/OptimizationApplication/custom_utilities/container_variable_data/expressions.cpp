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
#include <sstream>
#include <numeric>
#include <functional>

// Project includes
#include "includes/define.h"

// Application includes

// Include base h
#include "expressions.h"

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
    msg << "\b]";
    return msg.str();
}

} // namespace ExpressionHelperUtilities

std::size_t Expression::GetLocalSize() const
{
    const auto& r_shape = this->GetShape();
    return std::transform_reduce(r_shape.begin(), r_shape.end(), 1UL, std::multiplies{},
                                 [](const auto v) { return v; });
}

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

const std::vector<std::size_t> LiteralDoubleExpression::GetShape() const
{
    return {};
}

std::string LiteralDoubleExpression::Info() const
{
    std::stringstream msg;
    msg << mValue;
    return msg.str();
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

const std::vector<std::size_t> LiteralArray3Expression::GetShape() const
{
    return {mDimension};
}

std::string LiteralArray3Expression::Info() const
{
    std::stringstream msg;
    msg << mValue;
    return msg.str();
}

LiteralVectorExpression::LiteralVectorExpression(
    Kratos::shared_ptr<const Vector> pValue,
    const std::vector<IndexType>& rShape)
    : mpValue(pValue),
      mShape(rShape)
{
    mLocalSize = this->GetLocalSize();
}

Expression::Pointer LiteralVectorExpression::Create(
    Kratos::shared_ptr<Vector> pValue,
    const std::vector<IndexType>& rShape)
{
    return Kratos::make_intrusive<LiteralVectorExpression>(pValue, rShape);
}

double LiteralVectorExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return (*mpValue)[EntityIndex * mLocalSize + ComponentIndex];
}

const std::vector<std::size_t> LiteralVectorExpression::GetShape() const
{
    return mShape;
}

std::string LiteralVectorExpression::Info() const
{
    std::stringstream msg;
    msg << "v" << mShape;
    return msg.str();
}

BinaryExpression::BinaryExpression(
    Expression::Pointer pLeft,
    Expression::Pointer pRight)
    : mpLeft(pLeft),
      mpRight(pRight)
{
    const auto& r_left_shape = mpLeft->GetShape();
    const auto& r_right_shape = mpRight->GetShape();

    KRATOS_ERROR_IF_NOT(r_left_shape == r_right_shape || (typeid(*mpRight) == typeid(LiteralDoubleExpression(0))))
        << "Binary operation should have equal shape in left and right side expressions or right hand side should be a scalar. ["
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

const std::vector<std::size_t> BinaryExpression::GetShape() const
{
    return this->mpLeft->GetShape();
}

std::string BinaryExpression::Info() const
{
    std::stringstream msg;
    msg << "(" << this->mpLeft->Info() << this->Operation()
        << this->mpRight->Info() << ")";
    return msg.str();
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
    return mpLeft->Evaluate(EntityIndex, ComponentIndex) * mpRight->Evaluate(EntityIndex, ComponentIndex);
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
    return mpLeft->Evaluate(EntityIndex, ComponentIndex) / mpRight->Evaluate(EntityIndex, ComponentIndex);
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

BinaryWeightedMultiplicationExpression::BinaryWeightedMultiplicationExpression(
    Expression::Pointer pLeft,
    Expression::Pointer pRight)
    : mpLeft(pLeft),
      mpRight(pRight)
{
    const auto& r_left_shape = mpLeft->GetShape();
    const auto& r_right_shape = mpRight->GetShape();

    KRATOS_ERROR_IF_NOT(mpRight->GetLocalSize() == 1 && typeid(*mpRight) == typeid(LiteralVectorExpression(nullptr, {})))
        << "Binary weighted multiplication should have right hand with a vector having shape of [1]. ["
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

Expression::Pointer BinaryWeightedMultiplicationExpression::Create(
    Expression::Pointer pLeft,
    Expression::Pointer pRight)
{
    return Kratos::make_intrusive<BinaryWeightedMultiplicationExpression>(pLeft, pRight);
}

double BinaryWeightedMultiplicationExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return mpLeft->Evaluate(EntityIndex, ComponentIndex) * mpRight->Evaluate(EntityIndex, 0);
}

const std::vector<std::size_t> BinaryWeightedMultiplicationExpression::GetShape() const
{
    return this->mpLeft->GetShape();
}

std::string BinaryWeightedMultiplicationExpression::Info() const
{
    std::stringstream msg;
    msg << "(" << this->mpLeft->Info() << ":"
        << this->mpRight->Info() << ")";
    return msg.str();
}

}