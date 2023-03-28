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
#include "containers/container_variable_data/expressions/literal/literal_expression.h"
#include "containers/container_variable_data/expressions/literal/literal_flat_expression.h"

// Include base h
#include "container_variable_data.h"

namespace Kratos {

template <class TContainerType>
ContainerVariableData<TContainerType>::ContainerVariableData(ModelPart& rModelPart)
    : mpModelPart(&rModelPart)
{
}

template <class TContainerType>
ContainerVariableData<TContainerType>::ContainerVariableData(
    const ContainerVariableData& rOther)
    : mpExpression(rOther.mpExpression),
      mpModelPart(rOther.mpModelPart)
{
}
template <class TContainerType>
void ContainerVariableData<TContainerType>::CopyFrom(
    const ContainerVariableData<TContainerType>& rOther)
{
    KRATOS_ERROR_IF(this->GetContainer().size() != rOther.GetContainer().size())
        << "Mismatching model parts found with different number of entities in CopyFrom operation.\n"
        << this->mpModelPart->FullName()
        << ", origin model part name: " << rOther.GetModelPart().FullName() << " ].\n";

    mpExpression = rOther.mpExpression;
}

template <class TContainerType>
void ContainerVariableData<TContainerType>::Read(
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

    std::vector<IndexType> shape(ShapeSize);
    for (int i = 0; i < ShapeSize; ++i) {
        shape[i] = *(pShapeBegin++);
    }

    auto p_expression = LiteralFlatExpression::Create(number_of_entities, shape);
    this->mpExpression = p_expression;

    const IndexType local_size = this->GetExpression().GetFlattenedSize();

    IndexPartition<IndexType>(number_of_entities).for_each([pBegin, local_size, &p_expression](const IndexType EntityIndex) {
        const IndexType entity_data_begin_index = EntityIndex * local_size;
        double const* p_input_data_begin = pBegin + entity_data_begin_index;
        for (IndexType i = 0; i < local_size; ++i) {
            p_expression->SetData(entity_data_begin_index, i, *(p_input_data_begin+i));
        }
    });

    KRATOS_CATCH("");
}

template <class TContainerType>
void ContainerVariableData<TContainerType>::MoveFrom(
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
    for (int i = 0; i < ShapeSize; ++i) {
        shape[i] = *(pShapeBegin++);
    }

    auto p_expression = LiteralFlatExpression::Create(pBegin, number_of_entities, shape);
    this->mpExpression = p_expression;

    KRATOS_CATCH("");
}

template <class TContainerType>
void ContainerVariableData<TContainerType>::Evaluate(
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
    for (int i = 0; i < ShapeSize; ++i) {
        shape[i] = *(pShapeBegin++);
    }

    KRATOS_ERROR_IF_NOT(shape == r_expression.GetShape())
        << "Shape mismatch. [ Requested shape  = " << shape
        << ", available shape = " << r_expression.GetShape() << " ].\n";

    const IndexType local_size = r_expression.GetFlattenedSize();

    IndexPartition<IndexType>(number_of_entities).for_each([pBegin, local_size, &r_expression](const IndexType EntityIndex) {
        const IndexType entity_data_begin_index = EntityIndex * local_size;
        double* p_input_data_begin = pBegin + entity_data_begin_index;
        for (IndexType i = 0; i < local_size; ++i) {
            *(p_input_data_begin+i) = r_expression.Evaluate(EntityIndex, entity_data_begin_index, i);
        }
    });

    KRATOS_CATCH("");
}

template <class TContainerType>
void ContainerVariableData<TContainerType>::SetDataToZero()
{
    mpExpression = LiteralExpression<double>::Create(0.0);
}

template <class TContainerType>
void ContainerVariableData<TContainerType>::SetExpression(Expression::Pointer pExpression)
{
    this->mpExpression = pExpression;
}

template <class TContainerType>
const Expression& ContainerVariableData<TContainerType>::GetExpression() const
{
    return *(*mpExpression);
}


template <class TContainerType>
const Expression::Pointer ContainerVariableData<TContainerType>::pGetExpression() const
{
    return *mpExpression;
}

template <class TContainerType>
const std::vector<std::size_t> ContainerVariableData<TContainerType>::GetShape() const
{
    return this->GetExpression().GetShape();
}

template <class TContainerType>
std::size_t ContainerVariableData<TContainerType>::GetFlattenedSize() const
{
    return this->GetExpression().GetFlattenedSize();
}

template <class TContainerType>
ModelPart& ContainerVariableData<TContainerType>::GetModelPart()
{
    return *mpModelPart;
}

template <class TContainerType>
const ModelPart& ContainerVariableData<TContainerType>::GetModelPart() const
{
    return *mpModelPart;
}

template <>
ModelPart::NodesContainerType& ContainerVariableData<ModelPart::NodesContainerType>::GetContainer()
{
    return mpModelPart->GetCommunicator().LocalMesh().Nodes();
}

template <>
ModelPart::ConditionsContainerType& ContainerVariableData<ModelPart::ConditionsContainerType>::GetContainer()
{
    return mpModelPart->GetCommunicator().LocalMesh().Conditions();
}

template <>
ModelPart::ElementsContainerType& ContainerVariableData<ModelPart::ElementsContainerType>::GetContainer()
{
    return mpModelPart->GetCommunicator().LocalMesh().Elements();
}

template <>
const ModelPart::NodesContainerType& ContainerVariableData<ModelPart::NodesContainerType>::GetContainer() const
{
    return mpModelPart->GetCommunicator().LocalMesh().Nodes();
}

template <>
const ModelPart::ConditionsContainerType& ContainerVariableData<ModelPart::ConditionsContainerType>::GetContainer() const
{
    return mpModelPart->GetCommunicator().LocalMesh().Conditions();
}

template <>
const ModelPart::ElementsContainerType& ContainerVariableData<ModelPart::ElementsContainerType>::GetContainer() const
{
    return mpModelPart->GetCommunicator().LocalMesh().Elements();
}

template <class TContainerType>
std::string ContainerVariableData<TContainerType>::Info() const
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
std::string ContainerVariableData<TContainerType>::PrintData() const
{
    std::stringstream msg;
    msg << this->Info();
    return msg.str();
}

// template instantiations
template class ContainerVariableData<ModelPart::NodesContainerType>;
template class ContainerVariableData<ModelPart::ConditionsContainerType>;
template class ContainerVariableData<ModelPart::ElementsContainerType>;

} // namespace Kratos