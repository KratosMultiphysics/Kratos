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
#include <cmath>

// Project includes

// Include base h
#include "unary_abs_expression.h"

namespace Kratos {

UnaryAbsExpression::UnaryAbsExpression(Expression::ConstPointer pExpression)
    : Expression(pExpression->NumberOfEntities()),
      mpSourceExpression(pExpression)
{
}

Expression::Pointer UnaryAbsExpression::Create(Expression::ConstPointer pExpression)
{
    return Kratos::make_intrusive<UnaryAbsExpression>(std::move(pExpression));
}

double UnaryAbsExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return std::abs(mpSourceExpression->Evaluate(EntityIndex, EntityDataBeginIndex, ComponentIndex));
}

const std::vector<std::size_t> UnaryAbsExpression::GetItemShape() const
{
    return mpSourceExpression->GetItemShape();
}

std::size_t UnaryAbsExpression::GetMaxDepth() const
{
    return mpSourceExpression->GetMaxDepth() + 1;
}

std::string UnaryAbsExpression::Info() const
{
    std::stringstream msg;
    msg << "|" << mpSourceExpression << "|";
    return msg.str();
}

} // namespace Kratos