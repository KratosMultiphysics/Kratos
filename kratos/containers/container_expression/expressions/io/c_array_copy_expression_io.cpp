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
#include <type_traits>

// Project includes
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"
#include "utilities/parallel_utilities.h"

// Include base h
#include "c_array_copy_expression_io.h"


namespace Kratos {

template<class TRawDataType>
CArrayCopyExpressionInput::CArrayCopyExpressionInput(
    TRawDataType const* pBegin,
    const int NumberOfEntities,
    int const* pShapeBegin,
    const int ShapeSize)
    : mpCArray(pBegin),
      mNumberOfEntities(NumberOfEntities),
      mpShapeBegin(pShapeBegin),
      mShapeSize(ShapeSize)
{

}

Expression::Pointer CArrayCopyExpressionInput::Execute() const
{
    KRATOS_TRY

    // Convert int indices to IndexType
    std::vector<IndexType> shape(mShapeSize);
    std::copy(mpShapeBegin, mpShapeBegin + mShapeSize, shape.begin());

    return std::visit([&](auto pBegin){
        using data_type = std::remove_const_t<std::remove_pointer_t<decltype(pBegin)>>;

        auto p_expression = LiteralFlatExpression<data_type>::Create(mNumberOfEntities, shape);
        data_type* data_itr = p_expression->begin();

        const IndexType flattened_size = p_expression->GetItemComponentCount();
        const IndexType total_size = mNumberOfEntities * flattened_size;

        IndexPartition<IndexType>(total_size).for_each([pBegin, data_itr](const IndexType Index) {
            data_itr[Index] = pBegin[Index];
        });

        return Kratos::intrusive_ptr<Expression>(&*p_expression);

    }, mpCArray);

    KRATOS_CATCH("");
}

template<class TRawDataType>
CArrayCopyExpressionOutput::CArrayCopyExpressionOutput(
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

void CArrayCopyExpressionOutput::Execute(const Expression& rExpression)
{
    KRATOS_TRY

    std::vector<IndexType> shape(mShapeSize);
    std::copy(mpShapeBegin, mpShapeBegin + mShapeSize, shape.begin());

    KRATOS_ERROR_IF_NOT(shape == rExpression.GetItemShape())
        << "Shape mismatch. [ Requested shape  = " << shape
        << ", available shape = " << rExpression.GetItemShape() << " ].\n";

    const IndexType flattened_size = rExpression.GetItemComponentCount();

    std::visit([&](auto pBegin){
        using data_type = std::remove_const_t<std::remove_pointer_t<decltype(pBegin)>>;

        if constexpr(std::is_same_v<data_type, int>) {
            IndexPartition<IndexType>(mNumberOfEntities).for_each([pBegin, flattened_size, &rExpression](const IndexType EntityIndex) {
                const IndexType entity_data_begin_index = EntityIndex * flattened_size;
                int* p_input_data_begin = pBegin + entity_data_begin_index;
                for (IndexType i = 0; i < flattened_size; ++i) {
                    *(p_input_data_begin+i) = static_cast<int>(rExpression.Evaluate(EntityIndex, entity_data_begin_index, i));
                }
            });
        } else {
            IndexPartition<IndexType>(mNumberOfEntities).for_each([pBegin, flattened_size, &rExpression](const IndexType EntityIndex) {
                const IndexType entity_data_begin_index = EntityIndex * flattened_size;
                double* p_input_data_begin = pBegin + entity_data_begin_index;
                for (IndexType i = 0; i < flattened_size; ++i) {
                    *(p_input_data_begin+i) = rExpression.Evaluate(EntityIndex, entity_data_begin_index, i);
                }
            });
        }
    }, mpCArray);

    KRATOS_CATCH("");
}

// template instantiations
template CArrayCopyExpressionInput::CArrayCopyExpressionInput(int const*, const int, int const*, const int);
template CArrayCopyExpressionInput::CArrayCopyExpressionInput(double const*, const int, int const*, const int);

template CArrayCopyExpressionOutput::CArrayCopyExpressionOutput(int*, const int, int const*, const int);
template CArrayCopyExpressionOutput::CArrayCopyExpressionOutput(double*, const int, int const*, const int);

} // namespace Kratos
