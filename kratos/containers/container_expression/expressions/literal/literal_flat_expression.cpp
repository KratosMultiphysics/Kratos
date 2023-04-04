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

LiteralFlatExpression::LiteralFlatExpression(
    const IndexType NumberOfEntities,
    const std::vector<IndexType>& rShape)
    : mShape(rShape),
      mData(NumberOfEntities * this->GetFlattenedSize())
{
}

LiteralFlatExpression::LiteralFlatExpression(
    double* pDataBegin,
    const IndexType NumberOfEntities,
    const std::vector<IndexType>& rShape)
    : mShape(rShape),
      mData(pDataBegin)
{
}

LiteralFlatExpression::Pointer LiteralFlatExpression::Create(
    const IndexType NumberOfEntities,
    const std::vector<IndexType>& rShape)
{
    if (rShape.size() == 0) {
        return Kratos::make_intrusive<LiteralScalarFlatExpression>(NumberOfEntities, rShape);
    } else {
        return Kratos::make_intrusive<LiteralNonScalarFlatExpression>(NumberOfEntities, rShape);
    }
}

LiteralFlatExpression::Pointer LiteralFlatExpression::Create(
    double* pDataBegin,
    const IndexType NumberOfEntities,
    const std::vector<IndexType>& rShape)
{
    if (rShape.size() == 0) {
        return Kratos::make_intrusive<LiteralScalarFlatExpression>(pDataBegin, NumberOfEntities, rShape);
    } else {
        return Kratos::make_intrusive<LiteralNonScalarFlatExpression>(pDataBegin, NumberOfEntities, rShape);
    }
}

void LiteralFlatExpression::SetData(
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex,
    const double Value)
{
    mData[EntityDataBeginIndex + ComponentIndex] = Value;
}

const std::vector<std::size_t> LiteralFlatExpression::GetShape() const
{
    return mShape;
}

std::string LiteralFlatExpression::Info() const
{
    std::stringstream msg;
    msg << "Vec" << mShape;
    return msg.str();
}

double LiteralScalarFlatExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return mData[EntityIndex];
}

double LiteralNonScalarFlatExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return mData[EntityDataBeginIndex + ComponentIndex];
}

} // namespace Kratos