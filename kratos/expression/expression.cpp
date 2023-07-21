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
#include <algorithm>
#include <vector>

// Project includes
#include "expression/literal_expression.h"
#include "expression/literal_flat_expression.h"
#include "includes/ublas_interface.h"
#include "utilities/parallel_utilities.h"

// Include base h
#include "expression.h"


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

Expression::ExpressionIterator::ExpressionIterator(Expression::ConstPointer pExpression)
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

Expression::ConstPointer Expression::ExpressionIterator::GetExpression() const
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

Expression::Expression(const IndexType NumberOfEntities)
    : mNumberOfEntities(NumberOfEntities)
{
}

Expression::ConstPointer Expression::GetShrinkedExpression() const
{
    const IndexType number_of_entities = this->NumberOfEntities();
    const IndexType number_of_components = this->GetItemComponentCount();
    const auto& r_shape = this->GetItemShape();

    // generate the unique utilized expressions list
    std::set<Expression::ConstPointer> expressions;
    this->FillUtilizedExpressions(expressions);

    if (std::any_of(expressions.begin(), expressions.end(), [](const auto& pExpression) {
        return dynamic_cast<LiteralFlatExpression<char> const*>(pExpression.get()) ||
               dynamic_cast<LiteralFlatExpression<int> const*>(pExpression.get()) ||
               dynamic_cast<LiteralFlatExpression<double> const*>(pExpression.get()); })) {

        auto p_expression = LiteralFlatExpression<double>::Create(number_of_entities, r_shape);

        IndexPartition<IndexType>(number_of_entities).for_each([&](const auto Index) {
            const IndexType entity_data_begin_index = Index * number_of_components;
            for (IndexType i = 0; i < number_of_components; ++i) {
                p_expression->SetData(entity_data_begin_index, i, this->Evaluate(Index, entity_data_begin_index, i));
            }
        });

        return p_expression;
    } else {
        if (r_shape.size() == 0) {
            return LiteralExpression<double>::Create(this->Evaluate(0, 0, 0), number_of_entities);
        } else if (r_shape.size() == 1) {
            std::vector<IndexType> indices;
            indices.resize(r_shape[0]);
            std::iota(indices.begin(), indices.end(), 0);
            Vector value(r_shape[0]);
            std::transform(indices.begin(), indices.end(), value.begin(), [&](const auto Index) { return this->Evaluate(0, 0, Index);});
            return LiteralExpression<Vector>::Create(value, number_of_entities);
        } else if (r_shape.size() == 2) {
            std::vector<IndexType> indices;
            indices.resize(r_shape[0] * r_shape[1]);
            std::iota(indices.begin(), indices.end(), 0);
            Matrix value(r_shape[0], r_shape[1]);
            std::transform(indices.begin(), indices.end(), value.data().begin(), [&](const auto Index) { return this->Evaluate(0, 0, Index);});
            return LiteralExpression<Matrix>::Create(value, number_of_entities);
        } else {
            KRATOS_ERROR << "Unsupported shape found in expression: \"" << *this << "\". The maximum dimension supported is 2. [ shape = " << this->GetItemShape() << " ].\n";
            return LiteralExpression<double>::Create(this->Evaluate(0, 0, 0), number_of_entities);
        }
    }
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


} // namespace Kratos
