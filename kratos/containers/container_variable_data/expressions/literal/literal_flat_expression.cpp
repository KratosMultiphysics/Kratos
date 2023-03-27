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
    : mShape(rShape)
{
    mData.resize(NumberOfEntities * this->GetLocalSize());
}

LiteralFlatExpression::Pointer LiteralFlatExpression::Create(
    const IndexType NumberOfEntities,
    const std::vector<IndexType>& rShape)
{
    return Kratos::make_intrusive<LiteralFlatExpression>(NumberOfEntities, rShape);
}

void LiteralFlatExpression::SetData(
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex,
    const double Value)
{
    mData[EntityDataBeginIndex + ComponentIndex] = Value;
}

double LiteralFlatExpression::Evaluate(
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return mData[EntityDataBeginIndex + ComponentIndex];
}

const std::vector<std::size_t> LiteralFlatExpression::GetShape() const
{
    return mShape;
}

bool LiteralFlatExpression::IsScalar() const
{
    return false;
}

std::string LiteralFlatExpression::Info() const
{
    std::stringstream msg;
    msg << "v" << mShape;
    return msg.str();
}

} // namespace Kratos