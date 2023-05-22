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

Expression::ExpressionIterator::ExpressionIterator(Expression::Pointer pExpression)
    : mpExpression(pExpression),
      mEntityIndex(0),
      mEntityDataBeginIndex(0),
      mComponentIndex(0),
      mFlattenedShapeSize(pExpression->GetFlattenedShapeSize())
{
}

Expression::ExpressionIterator::ExpressionIterator(const ExpressionIterator& rOther)
    : mpExpression(rOther.mpExpression),
      mEntityIndex(rOther.mEntityIndex),
      mEntityDataBeginIndex(rOther.mEntityDataBeginIndex),
      mComponentIndex(rOther.mComponentIndex),
      mFlattenedShapeSize(rOther.mFlattenedShapeSize)
{
}

Expression::Pointer Expression::ExpressionIterator::GetExpression() const
{
    return mpExpression;
}

double Expression::ExpressionIterator::operator*() const
{
    return mpExpression->Evaluate(mEntityIndex, mEntityDataBeginIndex, mComponentIndex);
}

bool Expression::ExpressionIterator::operator==(const ExpressionIterator& rOther) const
{
    return (&*mpExpression == &*rOther.mpExpression) && (mEntityIndex == rOther.mEntityIndex) && (mComponentIndex == rOther.mComponentIndex);
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
    mComponentIndex = rOther.mComponentIndex;
    mFlattenedShapeSize = rOther.mFlattenedShapeSize;
    return *this;
}

Expression::ExpressionIterator& Expression::ExpressionIterator::operator++()
{
    ++mComponentIndex;
    if (mComponentIndex == mFlattenedShapeSize) {
        mComponentIndex = 0;
        ++mEntityIndex;
        mEntityDataBeginIndex = mEntityIndex * mFlattenedShapeSize;
    }
    return *this;
}

Expression::ExpressionIterator Expression::ExpressionIterator::operator++(int)
{
    ExpressionIterator temp = *this;
    ++*this;
    return temp;
}

std::size_t Expression::GetFlattenedShapeSize() const
{
    const auto& r_shape = this->GetShape();
    return std::accumulate(
        r_shape.begin(),
        r_shape.end(), 1UL,
        [](const auto V1, const auto V2) { return V1 * V2; });
}

Expression::ExpressionIterator Expression::cbegin() const
{
    return ExpressionIterator(this);
}

Expression::ExpressionIterator Expression::cend() const
{
    ExpressionIterator result(this);
    result.mEntityIndex = this->NumberOfEntities();
    result.mEntityDataBeginIndex = result.mEntityIndex * result.mFlattenedShapeSize;
    return result;
}

} // namespace Kratos