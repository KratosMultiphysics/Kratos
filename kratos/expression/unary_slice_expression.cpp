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
#include "unary_slice_expression.h"

namespace Kratos {

UnarySliceExpression::UnarySliceExpression(
    Expression::ConstPointer pExpression,
    const IndexType Offset,
    const IndexType Stride)
    : Expression(pExpression->NumberOfEntities()),
      mpSourceExpression(pExpression),
      mOffset(Offset),
      mStride(Stride),
      mSourceStride(mpSourceExpression->GetItemComponentCount())
{
    KRATOS_ERROR_IF(Offset + Stride > mSourceStride)
        << "Offset + stride is greater than the maximum stride of the given expression. [ Offset = "
        << Offset << ", stride = " << Stride << ", expression stride = "
        << mSourceStride << " ].\n" << "Source expression:\n"
        << *mpSourceExpression << "\n";
}

Expression::Pointer UnarySliceExpression::Create(
    Expression::ConstPointer pExpression,
    const IndexType Offset,
    const IndexType Stride)
{
    return Kratos::make_intrusive<UnarySliceExpression>(std::move(pExpression), Offset, Stride);
}

double UnarySliceExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType EntityDataBeginIndex,
    const IndexType ComponentIndex) const
{
    return mpSourceExpression->Evaluate(EntityIndex, EntityIndex * mSourceStride, ComponentIndex + mOffset);
}

const std::vector<std::size_t> UnarySliceExpression::GetItemShape() const
{
    // Following will always convert shape = [1] to shape = [].
    // Hence, if slicing is used, the output for scalars always will be
    // shape = [] irrespective of the input scalar shape. If they still want to
    // have the shape = [1], then UnaryReshapeExpression may be used.
    if (mStride == 1) {
        return {};
    } else {
        return {mStride};
    }
}

std::size_t UnarySliceExpression::GetMaxDepth() const
{
    return mpSourceExpression->GetMaxDepth() + 1;
}

std::string UnarySliceExpression::Info() const
{
    std::stringstream msg;
    msg << mpSourceExpression << "{" << mOffset << ":" << mStride << " }";
    return msg.str();
}

} // namespace Kratos
