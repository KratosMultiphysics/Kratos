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
#include "unary_combine_expression.h"

namespace Kratos {

double UnaryCombineExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    IndexType current_starting = 0;

    for (IndexType i = 0; i < mSourceExpressions.size(); ++i) {
        const IndexType current_stride = mStrides[i];
        if (ComponentIndex < current_stride) {
            return mSourceExpressions[i]->Evaluate(EntityIndex, EntityIndex * (current_stride - current_starting), ComponentIndex - current_starting);
        }
        current_starting = current_stride;
    }

    KRATOS_ERROR << "Component index is greater than the current expressions stride.\n";

    return 0.0;
}

const std::vector<std::size_t> UnaryCombineExpression::GetItemShape() const
{
    return {mStrides.back()};
}

std::string UnaryCombineExpression::Info() const
{
    std::stringstream msg;

    msg << "CombinedExpression: ";

    for (const auto& p_expression : mSourceExpressions) {
        msg << "\n\t" << *p_expression;
    }

    return msg.str();
}

} // namespace Kratos