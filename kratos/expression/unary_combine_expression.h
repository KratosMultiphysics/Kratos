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

#pragma once

// System includes
#include <string>
#include <vector>

// Project includes
#include "expression/expression.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Use to create an expression which combines given input expressions list in the order given.
 *
 * @details Instances of this expression combines given list of expression in the given order. This combines the entity values wise,
 *          and not append one expressions all the values after another.
 *
 * @note All the expressions should have the same number of entities.
 *
 */
class UnaryCombineExpression : public Expression {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Life cycle
    ///@{

    template<class TIteratorType>
    UnaryCombineExpression(
        TIteratorType Begin,
        TIteratorType End)
        : Expression(Begin != End ? (*Begin)->NumberOfEntities() : 0),
          mSourceExpressions(Begin, End)

    {
        mStrides.resize(mSourceExpressions.size());

        // every expression should have same number of entities
        IndexType local_stride = 0;
        for (IndexType i = 0; i < mSourceExpressions.size(); ++i) {
            const auto& p_expression = mSourceExpressions[i];

            KRATOS_ERROR_IF_NOT(p_expression->NumberOfEntities() == this->NumberOfEntities())
                << "Expression number of entities mismatch. [ required number of entities = "
                << NumberOfEntities() << ", found number of entities = "
                << p_expression->NumberOfEntities() << " ].\n"
                << "Expressions:\n"
                << "Reference = " << *mSourceExpressions[0] << "\n"
                << "Current   = " << p_expression << "\n";

            local_stride += p_expression->GetItemComponentCount();
            mStrides[i] = local_stride;
        }

        // Corner case: empty expression range provided
        KRATOS_ERROR_IF(this->GetItemComponentCount() == 0)
            << "No expressions were given.\n";
    }

    ///@}
    ///@name Public operations
    ///@{

    template<class TIteratorType>
    static Expression::Pointer Create(
        TIteratorType Begin,
        TIteratorType End)
    {
        return Kratos::make_intrusive<UnaryCombineExpression>(Begin, End);
    }

    double Evaluate(
        const IndexType EntityIndex,
        const IndexType EntityDataBeginIndex,
        const IndexType ComponentIndex) const override
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

    const std::vector<IndexType> GetItemShape() const override
    {
        return mStrides.size() == 1 && mStrides.back() == 1 ? std::vector<IndexType> {} : std::vector<IndexType> {mStrides.back()};
    }

    IndexType GetMaxDepth() const override
    {
        IndexType max_depth = 0;
        for (const auto& p_expression : mSourceExpressions) {
            max_depth = std::max(max_depth, p_expression->GetMaxDepth());
        }
        return max_depth + 1;
    }

    std::string Info() const override
    {
        std::stringstream msg;

        msg << "CombinedExpression: ";

        for (const auto& p_expression : mSourceExpressions) {
            msg << "\n\t" << *p_expression;
        }

        return msg.str();
    }

    ///@}
protected:
    ///@name Private member variables
    ///@{

    const std::vector<Expression::ConstPointer> mSourceExpressions;

    std::vector<IndexType> mStrides;

    ///@}
};

} // namespace Kratos