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