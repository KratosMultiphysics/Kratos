//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes
#include <sstream>
#include <vector>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"

// Application includes

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
const std::vector<std::size_t> LiteralExpression<TDataType>::GetShape() const
{
    return mShape;
}

template <class TDataType>
Expression::Pointer LiteralExpression<TDataType>::Create(const TDataType& Value)
{
    return Kratos::make_intrusive<LiteralExpression<TDataType>>(Value);
}

template <class TDataType>
LiteralExpression<TDataType>::LiteralExpression(const TDataType& rValue)
    : mValue(rValue),
      mShape({rValue.size()})
{

}

template <class TDataType>
double LiteralExpression<TDataType>::Evaluate(
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return mValue[ComponentIndex];
}

template <>
LiteralExpression<double>::LiteralExpression(const double& Value)
    : mValue(Value),
      mShape({})
{
}

template <>
double LiteralExpression<double>::Evaluate(
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return mValue;
}

// template instantiations
template class LiteralExpression<double>;
template class LiteralExpression<array_1d<double, 3>>;

} // namespace Kratos