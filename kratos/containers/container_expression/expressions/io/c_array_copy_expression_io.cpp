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
#include "includes/define.h"
#include "utilities/parallel_utilities.h"

// Include base h
#include "c_array_copy_expression_io.h"


namespace Kratos {

template<class TRawDataType>
CArrayExpressionInput::CArrayExpressionInput(
    TRawDataType const* pBegin,
    const int NumberOfEntities,
    int const* pShapeBegin,
    const int ShapeSize)
    : mpCArray(pBegin),
      mNumberOfEntities(NumberOfEntities),
      mShape(pShapeBegin, pShapeBegin + ShapeSize)
{

}

Expression::Pointer CArrayExpressionInput::Execute() const
{
    KRATOS_TRY

    return std::visit([&](auto pBegin){
        using data_type = std::remove_const_t<std::remove_pointer_t<decltype(pBegin)>>;

        auto p_expression = LiteralFlatExpression<data_type>::Create(mNumberOfEntities, mShape);
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
CArrayExpressionOutput::CArrayExpressionOutput(
    TRawDataType* pBegin,
    const int Size)
    : mpCArray(pBegin),
      mSize(Size)
{

}

void CArrayExpressionOutput::Execute(const Expression& rExpression)
{
    KRATOS_TRY

    const IndexType flattened_size = rExpression.GetItemComponentCount();
    const IndexType number_of_entities = rExpression.NumberOfEntities();

    KRATOS_ERROR_IF_NOT(number_of_entities * flattened_size == static_cast<unsigned int>(mSize))
        << "Shape mismatch. [ Requested size  = " << number_of_entities * flattened_size
        << ", available size = " << mSize << " ].\n";

    std::visit([&](auto pBegin){
        using data_type = std::remove_const_t<std::remove_pointer_t<decltype(pBegin)>>;

        if constexpr(std::is_same_v<data_type, int>) {
            IndexPartition<IndexType>(number_of_entities).for_each([pBegin, flattened_size, &rExpression](const IndexType EntityIndex) {
                const IndexType entity_data_begin_index = EntityIndex * flattened_size;
                int* p_input_data_begin = pBegin + entity_data_begin_index;
                for (IndexType i = 0; i < flattened_size; ++i) {
                    *(p_input_data_begin+i) = static_cast<int>(rExpression.Evaluate(EntityIndex, entity_data_begin_index, i));
                }
            });
        } else {
            IndexPartition<IndexType>(number_of_entities).for_each([pBegin, flattened_size, &rExpression](const IndexType EntityIndex) {
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
template KRATOS_API(KRATOS_CORE) CArrayExpressionInput::CArrayExpressionInput(int const*, const int, int const*, const int);
template KRATOS_API(KRATOS_CORE) CArrayExpressionInput::CArrayExpressionInput(double const*, const int, int const*, const int);

template KRATOS_API(KRATOS_CORE) CArrayExpressionOutput::CArrayExpressionOutput(int*, const int);
template KRATOS_API(KRATOS_CORE) CArrayExpressionOutput::CArrayExpressionOutput(double*, const int);

} // namespace Kratos
