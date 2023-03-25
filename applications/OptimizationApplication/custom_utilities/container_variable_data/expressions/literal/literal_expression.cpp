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
LiteralExpression<TDataType>::LiteralExpression(const TDataType& Value)
    : mValue(Value)
{
}

template <class TDataType>
Expression::Pointer LiteralExpression<TDataType>::Create(const TDataType& Value)
{
    return Kratos::make_intrusive<LiteralExpression<TDataType>>(Value);
}

template <class TDataType>
std::string LiteralExpression<TDataType>::Info() const
{
    std::stringstream msg;
    msg << mValue;
    return msg.str();
}

template <>
double LiteralExpression<double>::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return mValue;
}

template <>
const std::vector<std::size_t> LiteralExpression<double>::GetShape() const
{
    return {};
}

template <>
double LiteralExpression<array_1d<double, 3>>::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return mValue[ComponentIndex];
}

template <>
const std::vector<std::size_t> LiteralExpression<array_1d<double, 3>>::GetShape() const
{
    return {3};
}

// template instantiations
template class LiteralExpression<double>;
template class LiteralExpression<array_1d<double, 3>>;

} // namespace Kratos