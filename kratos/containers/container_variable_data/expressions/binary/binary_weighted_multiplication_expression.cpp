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
#include <sstream>
#include <vector>

// Project includes
#include "containers/container_variable_data/expressions/literal/literal_flat_expression.h"

// Include base h
#include "binary_weighted_multiplication_expression.h"

namespace Kratos {

namespace BinaryExpressionHelperUtilities {
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

} // namespace BinaryExpressionHelperUtilities

BinaryWeightedMultiplicationExpression::BinaryWeightedMultiplicationExpression(
    Expression::Pointer pLeft,
    Expression::Pointer pRight)
    : mpLeft(pLeft),
      mpRight(pRight),
      mLocalSize(mpLeft->GetLocalSize())
{
    const auto& r_left_shape = mpLeft->GetShape();
    const auto& r_right_shape = mpRight->GetShape();

    KRATOS_ERROR_IF_NOT(mpRight->GetLocalSize() == 1 && r_right_shape.size() == 0)
        << "Binary weighted multiplication should have right hand with a "
           "vector/scalar having shape of []. ["
        << "lhs shape = " << BinaryExpressionHelperUtilities::GetShape(r_left_shape) << ", "
        << "rhs shape = " << BinaryExpressionHelperUtilities::GetShape(r_right_shape) << " ].\n"
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
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return mpLeft->Evaluate(EntityDataBeginIndex, ComponentIndex) * mpRight->Evaluate(EntityDataBeginIndex / mLocalSize, 0);
}

const std::vector<std::size_t> BinaryWeightedMultiplicationExpression::GetShape() const
{
    return this->mpLeft->GetShape();
}

std::string BinaryWeightedMultiplicationExpression::Info() const
{
    std::stringstream msg;
    msg << "(" << this->mpLeft->Info() << ":" << this->mpRight->Info() << ")";
    return msg.str();
}

} // namespace Kratos