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
#include "expression/literal_flat_expression.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

// Include base h
#include "c_array_expression_io.h"

namespace Kratos {

template<class TRawDataType>
CArrayExpressionIO::Input::Input(
    TRawDataType const* pBegin,
    const int NumberOfEntities,
    int const* pShapeBegin,
    const int ShapeSize)
    : mpCArray(pBegin),
      mNumberOfEntities(NumberOfEntities),
      mShape(pShapeBegin, pShapeBegin + ShapeSize)
{

}

Expression::Pointer CArrayExpressionIO::Input::Execute() const
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
CArrayExpressionIO::MoveInput::MoveInput(
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

Expression::Pointer CArrayExpressionIO::MoveInput::Execute() const
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

template<class TRawDataType>
CArrayExpressionIO::Output::Output(
    TRawDataType* pBegin,
    const int Size)
    : mpCArray(pBegin),
      mSize(Size)
{

}

void CArrayExpressionIO::Output::Execute(const Expression& rExpression)
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

template<class TRawDataType, class TContainerType, MeshType TMeshType>
void CArrayExpressionIO::Read(
    ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
    TRawDataType const* pBegin,
    const int NumberOfEntities,
    int const* pShapeBegin,
    const int ShapeSize)
{
    rContainerExpression.SetExpression(
        Input(pBegin, NumberOfEntities, pShapeBegin, ShapeSize).Execute());
}

template<class TContainerType, MeshType TMeshType>
void CArrayExpressionIO::Read(
    ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
    const Vector& rValues,
    const std::vector<IndexType>& rShape)
{
    const IndexType number_of_entities = rContainerExpression.GetContainer().size();
    auto p_flat_expression = LiteralFlatExpression<double>::Create(number_of_entities, rShape);
    const IndexType flat_size = p_flat_expression->GetItemComponentCount() * number_of_entities;

    KRATOS_ERROR_IF_NOT(flat_size == rValues.size())
        << "Vector does not contain required number of values for the given container and the shape. [ Vector size = "
        << rValues.size() << ", shape = " << rShape << " number of entities = "
        << number_of_entities << ", required number of values = "
        << flat_size << " ].\n";

    IndexPartition<IndexType>(flat_size).for_each([&p_flat_expression, &rValues](const IndexType Index) {
        *(p_flat_expression->begin() + Index) = rValues[Index];
    });

    rContainerExpression.SetExpression(p_flat_expression);
}

template<class TRawDataType, class TContainerType, MeshType TMeshType>
void CArrayExpressionIO::Move(
    ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
    TRawDataType* pBegin,
    const int NumberOfEntities,
    int const* pShapeBegin,
    const int ShapeSize)
{
    rContainerExpression.SetExpression(
        MoveInput(pBegin, NumberOfEntities, pShapeBegin, ShapeSize).Execute());
}

template<class TContainerType, MeshType TMeshType>
void CArrayExpressionIO::Move(
    ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
    Vector& rValues,
    const std::vector<IndexType>& rShape)
{
    const IndexType number_of_entities = rContainerExpression.GetContainer().size();
    auto p_flat_expression = LiteralFlatExpression<double>::Create(rValues.data().begin(), number_of_entities, rShape);
    const IndexType flat_size = p_flat_expression->GetItemComponentCount() * number_of_entities;

    KRATOS_ERROR_IF_NOT(flat_size == rValues.size())
        << "Vector does not contain required number of values for the given container and the shape. [ Vector size = "
        << rValues.size() << ", shape = " << rShape << " number of entities = "
        << number_of_entities << ", required number of values = "
        << flat_size << " ].\n";

    rContainerExpression.SetExpression(p_flat_expression);
}

template<class TRawDataType, class TContainerType, MeshType TMeshType>
void CArrayExpressionIO::Write(
    const ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
    TRawDataType* pBegin,
    const int mSize)
{
    Output(pBegin, mSize).Execute(rContainerExpression.GetExpression());
}

template<class TContainerType, MeshType TMeshType>
void CArrayExpressionIO::Write(
    const ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
    Vector& rValues)
{
    const IndexType flat_size = rContainerExpression.GetItemComponentCount() * rContainerExpression.GetContainer().size();
    if (rValues.size() != flat_size) {
        rValues.resize(flat_size, false);
    }

    Output(rValues.data().begin(), rValues.size()).Execute(rContainerExpression.GetExpression());
}

#define KRATOS_INSTANTIATE_C_ARRAY_EXPRESSION_IO_METHODS(RAW_DATA_TYPE, CONTAINER_TYPE, MESH_TYPE)                                                                      \
    template void KRATOS_API(KRATOS_CORE) CArrayExpressionIO::Read(ContainerExpression<CONTAINER_TYPE, MESH_TYPE>&, RAW_DATA_TYPE const*, const int, int const*, const int); \
    template void KRATOS_API(KRATOS_CORE) CArrayExpressionIO::Move(ContainerExpression<CONTAINER_TYPE, MESH_TYPE>&, RAW_DATA_TYPE*, const int, int const*, const int);       \
    template void KRATOS_API(KRATOS_CORE) CArrayExpressionIO::Write(const ContainerExpression<CONTAINER_TYPE, MESH_TYPE>&, RAW_DATA_TYPE*, const int);

#define KRATOS_INSTANTIATE_C_ARRAY_EXPRESSION_VECTOR_AND_OTHER_IO_METHODS(CONTAINER_TYPE, MESH_TYPE)                                                                          \
    template void KRATOS_API(KRATOS_CORE) CArrayExpressionIO::Read(ContainerExpression<CONTAINER_TYPE, MESH_TYPE>&, const Vector&, const std::vector<IndexType>&);  \
    template void KRATOS_API(KRATOS_CORE) CArrayExpressionIO::Move(ContainerExpression<CONTAINER_TYPE, MESH_TYPE>&, Vector&, const std::vector<IndexType>&);        \
    template void KRATOS_API(KRATOS_CORE) CArrayExpressionIO::Write(const ContainerExpression<CONTAINER_TYPE, MESH_TYPE>&, Vector&);                                \
    KRATOS_INSTANTIATE_C_ARRAY_EXPRESSION_IO_METHODS(int, CONTAINER_TYPE, MESH_TYPE)                                                                                \
    KRATOS_INSTANTIATE_C_ARRAY_EXPRESSION_IO_METHODS(double, CONTAINER_TYPE, MESH_TYPE)

#define KRATOS_INSTANTIATE_C_ARRAY_EXPRESSION_IO_METHODS_1(MESH_TYPE)                                                 \
    KRATOS_INSTANTIATE_C_ARRAY_EXPRESSION_VECTOR_AND_OTHER_IO_METHODS(ModelPart::NodesContainerType, MESH_TYPE)       \
    KRATOS_INSTANTIATE_C_ARRAY_EXPRESSION_VECTOR_AND_OTHER_IO_METHODS(ModelPart::ConditionsContainerType, MESH_TYPE)  \
    KRATOS_INSTANTIATE_C_ARRAY_EXPRESSION_VECTOR_AND_OTHER_IO_METHODS(ModelPart::ElementsContainerType, MESH_TYPE)

// template instantiations
template KRATOS_API(KRATOS_CORE) CArrayExpressionIO::Input::Input(int const*, const int, int const*, const int);
template KRATOS_API(KRATOS_CORE) CArrayExpressionIO::Input::Input(double const*, const int, int const*, const int);

template KRATOS_API(KRATOS_CORE) CArrayExpressionIO::MoveInput::MoveInput(int*, const int, int const*, const int);
template KRATOS_API(KRATOS_CORE) CArrayExpressionIO::MoveInput::MoveInput(double*, const int, int const*, const int);

template KRATOS_API(KRATOS_CORE) CArrayExpressionIO::Output::Output(int*, const int);
template KRATOS_API(KRATOS_CORE) CArrayExpressionIO::Output::Output(double*, const int);

KRATOS_INSTANTIATE_C_ARRAY_EXPRESSION_IO_METHODS_1(MeshType::Local)
KRATOS_INSTANTIATE_C_ARRAY_EXPRESSION_IO_METHODS_1(MeshType::Ghost)
KRATOS_INSTANTIATE_C_ARRAY_EXPRESSION_IO_METHODS_1(MeshType::Interface)

#undef KRATOS_INSTANTIATE_C_ARRAY_EXPRESSION_IO_METHODS
#undef KRATOS_INSTANTIATE_C_ARRAY_EXPRESSION_IO_METHODS_1
#undef KRATOS_INSTANTIATE_C_ARRAY_EXPRESSION_VECTOR_AND_OTHER_IO_METHODS

} // namespace Kratos
