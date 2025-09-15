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
#include "literal_flat_expression.h"

namespace Kratos {

template<class TRawDataType>
LiteralFlatExpression<TRawDataType>::LiteralFlatExpression(
    const IndexType NumberOfEntities,
    const std::vector<IndexType>& rShape)
    : Expression(NumberOfEntities),
      mShape(rShape),
      mData(NumberOfEntities * this->GetItemComponentCount())
{
}

template<class TRawDataType>
LiteralFlatExpression<TRawDataType>::LiteralFlatExpression(
    TRawDataType* pDataBegin,
    const IndexType NumberOfEntities,
    const std::vector<IndexType>& rShape)
    : Expression(NumberOfEntities),
      mShape(rShape),
      mData(pDataBegin, NumberOfEntities * this->GetItemComponentCount())
{
}

template<class TRawDataType>
typename LiteralFlatExpression<TRawDataType>::Pointer LiteralFlatExpression<TRawDataType>::Create(
    const IndexType NumberOfEntities,
    const std::vector<IndexType>& rShape)
{
    if (rShape.size() == 0) {
        return Kratos::make_intrusive<LiteralScalarFlatExpression<TRawDataType>>(NumberOfEntities, rShape);
    } else {
        return Kratos::make_intrusive<LiteralNonScalarFlatExpression<TRawDataType>>(NumberOfEntities, rShape);
    }
}

template<class TRawDataType>
typename LiteralFlatExpression<TRawDataType>::Pointer LiteralFlatExpression<TRawDataType>::Create(
    TRawDataType* pDataBegin,
    const IndexType NumberOfEntities,
    const std::vector<IndexType>& rShape)
{
    if (rShape.size() == 0) {
        return Kratos::make_intrusive<LiteralScalarFlatExpression<TRawDataType>>(pDataBegin, NumberOfEntities, rShape);
    } else {
        return Kratos::make_intrusive<LiteralNonScalarFlatExpression<TRawDataType>>(pDataBegin, NumberOfEntities, rShape);
    }
}

template<class TRawDataType>
void LiteralFlatExpression<TRawDataType>::SetData(
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex,
    const TRawDataType Value)
{
    *(mData.begin() + EntityDataBeginIndex + ComponentIndex) = Value;
}

template<class TRawDataType>
const std::vector<std::size_t> LiteralFlatExpression<TRawDataType>::GetItemShape() const
{
    return mShape;
}

template<>
std::string LiteralFlatExpression<char>::Info() const
{
    std::stringstream msg;
    msg << "CharVec" << mShape;
    return msg.str();
}

template<>
std::string LiteralFlatExpression<int>::Info() const
{
    std::stringstream msg;
    msg << "IntVec" << mShape;
    return msg.str();
}

template<>
std::string LiteralFlatExpression<double>::Info() const
{
    std::stringstream msg;
    msg << "DoubleVec" << mShape;
    return msg.str();
}

template<class TRawDataType>
double LiteralScalarFlatExpression<TRawDataType>::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return *(this->mData.cbegin() + EntityIndex);
}

template<class TRawDataType>
double LiteralNonScalarFlatExpression<TRawDataType>::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return *(this->mData.cbegin() + EntityDataBeginIndex + ComponentIndex);
}

// template instantiations
template class LiteralFlatExpression<char>;
template class LiteralScalarFlatExpression<char>;
template class LiteralNonScalarFlatExpression<char>;

template class LiteralFlatExpression<int>;
template class LiteralScalarFlatExpression<int>;
template class LiteralNonScalarFlatExpression<int>;

template class LiteralFlatExpression<double>;
template class LiteralScalarFlatExpression<double>;
template class LiteralNonScalarFlatExpression<double>;
} // namespace Kratos