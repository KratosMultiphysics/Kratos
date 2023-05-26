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

// Project includes

// Include base h
#include "unary_reshape_expression.h"

namespace Kratos {

UnaryReshapeExpression::UnaryReshapeExpression(
    Expression::Pointer pExpression,
    const std::vector<IndexType>& rShape)
    : Expression(pExpression->NumberOfEntities()),
      mpSourceExpression(pExpression),
      mShape(rShape)
{
    KRATOS_ERROR_IF_NOT(this->GetItemComponentCount() == mpSourceExpression->GetItemComponentCount())
        << "The source and destination should have same stride. [ source stride = "
        << mpSourceExpression->GetItemComponentCount() << ", destination stride = "
        << this->GetItemComponentCount() << " ].\n"
        << "Source expression:\n"
        << *mpSourceExpression << "\n";
}

Expression::Pointer UnaryReshapeExpression::Create(
    Expression::Pointer pExpression,
    const std::vector<IndexType>& rShape)
{
    return Kratos::make_intrusive<UnaryReshapeExpression>(std::move(pExpression), std::move(rShape));
}

double UnaryReshapeExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return mpSourceExpression->Evaluate(EntityIndex, EntityDataBeginIndex, ComponentIndex);
}

const std::vector<std::size_t> UnaryReshapeExpression::GetItemShape() const
{
    return mShape;
}

std::string UnaryReshapeExpression::Info() const
{
    std::stringstream msg;
    msg << mpSourceExpression << "-> [";
    for (const auto v : mShape) {
        msg << " " << v;
    }
    msg << " ]";
    return msg.str();
}

} // namespace Kratos