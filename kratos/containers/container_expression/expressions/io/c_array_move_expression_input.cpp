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

// Project includes
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"

// Include base h
#include "c_array_move_expression_input.h"

namespace Kratos {

template<class TRawDataType>
CArrayMoveExpressionInput::CArrayMoveExpressionInput(
    TRawDataType* pBegin,
    const int NumberOfEntities,
    int const* pShapeBegin,
    const int ShapeSize)
    : mpCArray(pBegin),
      mNumberOfEntities(NumberOfEntities),
      mpShapeBegin(pShapeBegin),
      mShapeSize(ShapeSize)
{

}

Expression::Pointer CArrayMoveExpressionInput::Execute() const
{
    KRATOS_TRY

    // Convert int indices to IndexType
    std::vector<IndexType> shape(mShapeSize);
    std::copy(mpShapeBegin, mpShapeBegin + mShapeSize, shape.begin());

    return std::visit([&](auto pBegin){
        using data_type = std::remove_const_t<std::remove_pointer_t<decltype(pBegin)>>;

        auto p_expression = LiteralFlatExpression<data_type>::Create(pBegin, mNumberOfEntities, shape);

        return Kratos::intrusive_ptr<Expression>(&*p_expression);

    }, mpCArray);

    KRATOS_CATCH("");
}

// template instantiations
template CArrayMoveExpressionInput::CArrayMoveExpressionInput(int*, const int, int const*, const int);
template CArrayMoveExpressionInput::CArrayMoveExpressionInput(double*, const int, int const*, const int);

} // namespace Kratos
