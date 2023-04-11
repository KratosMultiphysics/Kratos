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
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "containers/container_expression/expressions/literal/literal_expression.h"
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"

// Include base h
#include "container_expression.h"

namespace Kratos {

template <class TContainerType>
ContainerExpression<TContainerType>::ContainerExpression(ModelPart& rModelPart)
    : mpExpression(),
      mpModelPart(&rModelPart)
{
}

template <class TContainerType>
ContainerExpression<TContainerType>::ContainerExpression(
    const ContainerExpression& rOther)
    : mpExpression(rOther.mpExpression),
      mpModelPart(rOther.mpModelPart)
{
}
template <class TContainerType>
void ContainerExpression<TContainerType>::CopyFrom(
    const ContainerExpression<TContainerType>& rOther)
{
    KRATOS_ERROR_IF(this->GetContainer().size() != rOther.GetContainer().size())
        << "Mismatching model parts found with different number of entities in CopyFrom operation.\n"
        << this->mpModelPart->FullName()
        << ", origin model part name: " << rOther.GetModelPart().FullName() << " ].\n";

    mpExpression = rOther.mpExpression;
}

template <class TContainerType>
void ContainerExpression<TContainerType>::Read(
    double const* pBegin,
    const int NumberOfEntities,
    int const* pShapeBegin,
    const int ShapeSize)
{
    KRATOS_TRY

    const IndexType number_of_entities = NumberOfEntities;

    KRATOS_ERROR_IF_NOT(number_of_entities == this->GetContainer().size())
        << "Number of entities does not match with the local container size. [ "
           "NumberOfEntities = "
        << NumberOfEntities
        << ", local container size = " << this->GetContainer().size() << " ].\n";

    // Convert int indices to IndexType
    std::vector<IndexType> shape(ShapeSize);
    std::copy(pShapeBegin,
              pShapeBegin + ShapeSize,
              shape.begin());

    auto p_expression = LiteralFlatExpression::Create(number_of_entities, shape);
    this->mpExpression = p_expression;

    const IndexType flattened_size = this->GetExpression().GetFlattenedSize();

    IndexPartition<IndexType>(number_of_entities).for_each([pBegin, flattened_size, &p_expression](const IndexType EntityIndex) {
        const IndexType entity_data_begin_index = EntityIndex * flattened_size;
        double const* p_input_data_begin = pBegin + entity_data_begin_index;
        for (IndexType i = 0; i < flattened_size; ++i) {
            p_expression->SetData(entity_data_begin_index, i, *(p_input_data_begin+i));
        }
    });

    KRATOS_CATCH("");
}

template <class TContainerType>
void ContainerExpression<TContainerType>::MoveFrom(
    double* pBegin,
    const int NumberOfEntities,
    int const* pShapeBegin,
    const int ShapeSize)
{
    KRATOS_TRY

    const IndexType number_of_entities = NumberOfEntities;

    KRATOS_ERROR_IF_NOT(number_of_entities == this->GetContainer().size())
        << "Number of entities does not match with the local container size. [ "
           "NumberOfEntities = "
        << NumberOfEntities
        << ", local container size = " << this->GetContainer().size() << " ].\n";

    std::vector<IndexType> shape(ShapeSize);
    std::copy(pShapeBegin, pShapeBegin + ShapeSize, shape.begin());

    auto p_expression = LiteralFlatExpression::Create(pBegin, number_of_entities, shape);
    this->mpExpression = p_expression;

    KRATOS_CATCH("");
}

template <class TContainerType>
void ContainerExpression<TContainerType>::Evaluate(
    double* pBegin,
    const int NumberOfEntities,
    int const* pShapeBegin,
    const int ShapeSize) const
{
    KRATOS_TRY

    const IndexType number_of_entities = NumberOfEntities;

    KRATOS_ERROR_IF_NOT(number_of_entities == this->GetContainer().size())
        << "Number of entities does not match with the local container size. [ "
           "NumberOfEntities = "
        << NumberOfEntities
        << ", local container size = " << this->GetContainer().size() << " ].\n";

    const auto& r_expression = this->GetExpression();

    std::vector<IndexType> shape(ShapeSize);
    std::copy(pShapeBegin, pShapeBegin + ShapeSize, shape.begin());

    KRATOS_ERROR_IF_NOT(shape == r_expression.GetShape())
        << "Shape mismatch. [ Requested shape  = " << shape
        << ", available shape = " << r_expression.GetShape() << " ].\n";

    const IndexType flattened_size = r_expression.GetFlattenedSize();

    IndexPartition<IndexType>(number_of_entities).for_each([pBegin, flattened_size, &r_expression](const IndexType EntityIndex) {
        const IndexType entity_data_begin_index = EntityIndex * flattened_size;
        double* p_input_data_begin = pBegin + entity_data_begin_index;
        for (IndexType i = 0; i < flattened_size; ++i) {
            *(p_input_data_begin+i) = r_expression.Evaluate(EntityIndex, entity_data_begin_index, i);
        }
    });

    KRATOS_CATCH("");
}

template <class TContainerType>
void ContainerExpression<TContainerType>::SetDataToZero()
{
    mpExpression = LiteralExpression<double>::Create(0.0);
}

template <class TContainerType>
void ContainerExpression<TContainerType>::SetExpression(Expression::Pointer pExpression)
{
    this->mpExpression = pExpression;
}

template <class TContainerType>
const Expression& ContainerExpression<TContainerType>::GetExpression() const
{
    return *(*mpExpression);
}


template <class TContainerType>
const Expression::Pointer ContainerExpression<TContainerType>::pGetExpression() const
{
    return *mpExpression;
}

template <class TContainerType>
const std::vector<std::size_t> ContainerExpression<TContainerType>::GetShape() const
{
    return this->GetExpression().GetShape();
}

template <class TContainerType>
std::size_t ContainerExpression<TContainerType>::GetFlattenedSize() const
{
    return this->GetExpression().GetFlattenedSize();
}

template <class TContainerType>
ModelPart& ContainerExpression<TContainerType>::GetModelPart()
{
    return *mpModelPart;
}

template <class TContainerType>
const ModelPart& ContainerExpression<TContainerType>::GetModelPart() const
{
    return *mpModelPart;
}

template <>
ModelPart::NodesContainerType& ContainerExpression<ModelPart::NodesContainerType>::GetContainer()
{
    return mpModelPart->GetCommunicator().LocalMesh().Nodes();
}

template <>
ModelPart::ConditionsContainerType& ContainerExpression<ModelPart::ConditionsContainerType>::GetContainer()
{
    return mpModelPart->GetCommunicator().LocalMesh().Conditions();
}

template <>
ModelPart::ElementsContainerType& ContainerExpression<ModelPart::ElementsContainerType>::GetContainer()
{
    return mpModelPart->GetCommunicator().LocalMesh().Elements();
}

template <>
const ModelPart::NodesContainerType& ContainerExpression<ModelPart::NodesContainerType>::GetContainer() const
{
    return mpModelPart->GetCommunicator().LocalMesh().Nodes();
}

template <>
const ModelPart::ConditionsContainerType& ContainerExpression<ModelPart::ConditionsContainerType>::GetContainer() const
{
    return mpModelPart->GetCommunicator().LocalMesh().Conditions();
}

template <>
const ModelPart::ElementsContainerType& ContainerExpression<ModelPart::ElementsContainerType>::GetContainer() const
{
    return mpModelPart->GetCommunicator().LocalMesh().Elements();
}

template <class TContainerType>
std::string ContainerExpression<TContainerType>::Info() const
{
    std::stringstream msg;

    msg << "VariableDataInfo: [ ModelPartName: "
        << this->GetModelPart().FullName()
        << ", Number of entities: " << this->GetContainer().size();

    if (mpExpression.has_value()) {
        msg << ", DataDimension: " << this->GetFlattenedSize();
        msg << ", Expression: " << this->GetExpression();
    } else {
        msg << ", Expression: not initialized";
    }

    msg << " ].\n";

    return msg.str();
}

template <class TContainerType>
std::string ContainerExpression<TContainerType>::PrintData() const
{
    std::stringstream msg;
    msg << this->Info();
    return msg.str();
}

// template instantiations
template class ContainerExpression<ModelPart::NodesContainerType>;
template class ContainerExpression<ModelPart::ConditionsContainerType>;
template class ContainerExpression<ModelPart::ElementsContainerType>;

} // namespace Kratos