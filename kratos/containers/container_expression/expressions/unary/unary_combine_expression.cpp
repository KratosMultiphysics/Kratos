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

namespace UnaryCombineExpressionHelperUtilities
{

using const_iterator = std::vector<Expression::Pointer>::const_iterator;

std::size_t GetNumberOfEntities(
    const_iterator Begin,
    const_iterator End)
{
    if (Begin != End) {
        return (*Begin)->NumberOfEntities();
    } else {
        KRATOS_ERROR << "Trying to create a UnaryCombineExpression with empty "
                        "expressions list.\n";
        return 0;
    }
}

} // namespace UnaryCombineExpressionHelperUtilities

UnaryCombineExpression::UnaryCombineExpression(
    const_iterator Begin,
    const_iterator End)
    : Expression(UnaryCombineExpressionHelperUtilities::GetNumberOfEntities(Begin, End)),
      mSourceExpressions(Begin, End)

{
    mStrides.resize(mSourceExpressions.size());

    // every expression should have same number of entities
    IndexType local_stride = 0;
    for (IndexType i = 0; i < mSourceExpressions.size(); ++i) {
        const auto& p_expression = mSourceExpressions[i];

        KRATOS_ERROR_IF_NOT(p_expression->NumberOfEntities() == NumberOfEntities())
            << "Expression number of entities mismatch. [ required number of entities = "
            << NumberOfEntities() << ", found number of entities = "
            << p_expression->NumberOfEntities() << " ].\n"
            << "Expressions:\n"
            << "Reference = " << *mSourceExpressions[0] << "\n"
            << "Current   = " << p_expression << "\n";

        local_stride += p_expression->GetItemComponentCount();
        mStrides[i] = local_stride;
    }
}

Expression::Pointer UnaryCombineExpression::Create(
    const_iterator Begin,
    const_iterator End)
{
    return Kratos::make_intrusive<UnaryCombineExpression>(Begin, End);
}

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