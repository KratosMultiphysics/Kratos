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

// Project includes

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
