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
#include <vector>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"

// Include base h
#include "literal_expression.h"

namespace Kratos {

template <class TDataType>
std::string LiteralExpression<TDataType>::Info() const
{
    std::stringstream msg;
    msg << mValue;
    return msg.str();
}

template <class TDataType>
const std::vector<std::size_t> LiteralExpression<TDataType>::GetItemShape() const
{
    return mShape;
}

template <class TDataType>
Expression::Pointer LiteralExpression<TDataType>::Create(
    const TDataType& Value,
    const IndexType NumberOfEntities)
{
    return Kratos::make_intrusive<LiteralExpression<TDataType>>(Value, NumberOfEntities);
}

template <class TDataType>
LiteralExpression<TDataType>::LiteralExpression(
    const TDataType& rValue,
    const IndexType NumberOfEntities)
    : Expression(NumberOfEntities),
      mValue(rValue),
      mShape({rValue.size()})
{

}

template <class TDataType>
double LiteralExpression<TDataType>::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return mValue[ComponentIndex];
}

template <>
LiteralExpression<int>::LiteralExpression(
    const int& Value,
    const IndexType NumberOfEntities)
    : Expression(NumberOfEntities),
      mValue(Value),
      mShape({})
{
}

template <>
LiteralExpression<double>::LiteralExpression(
    const double& Value,
    const IndexType NumberOfEntities)
    : Expression(NumberOfEntities),
      mValue(Value),
      mShape({})
{
}

template <>
double LiteralExpression<int>::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return static_cast<double>(mValue);
}

template <>
double LiteralExpression<double>::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return mValue;
}

template <>
LiteralExpression<Matrix>::LiteralExpression(
    const Matrix& Value,
    const IndexType NumberOfEntities)
    : Expression(NumberOfEntities),
      mValue(Value),
      mShape({Value.size1(), Value.size2()})
{
}

template <>
double LiteralExpression<Matrix>::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return mValue.data()[ComponentIndex];
}

// template instantiations
template class LiteralExpression<int>;
template class LiteralExpression<double>;
template class LiteralExpression<array_1d<double, 3>>;
template class LiteralExpression<array_1d<double, 4>>;
template class LiteralExpression<array_1d<double, 6>>;
template class LiteralExpression<array_1d<double, 9>>;
template class LiteralExpression<Vector>;
template class LiteralExpression<Matrix>;

} // namespace Kratos