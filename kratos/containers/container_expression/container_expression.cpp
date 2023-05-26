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

namespace ContainerExpressionHelperUtilities
{
template <class TMeshType>
ModelPart::MeshType& GetMesh(ModelPart& rModelPart);

template<>
ModelPart::MeshType& GetMesh<MeshType::Local>(ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh();
}

template<>
ModelPart::MeshType& GetMesh<MeshType::Ghost>(ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().GhostMesh();
}

template<>
ModelPart::MeshType& GetMesh<MeshType::Interface>(ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().InterfaceMesh();
}

template <class TMeshType>
const ModelPart::MeshType& GetMesh(const ModelPart& rModelPart);

template<>
const ModelPart::MeshType& GetMesh<MeshType::Local>(const ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh();
}

template<>
const ModelPart::MeshType& GetMesh<MeshType::Ghost>(const ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().GhostMesh();
}

template<>
const ModelPart::MeshType& GetMesh<MeshType::Interface>(const ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().InterfaceMesh();
}

template <class TContainerType>
TContainerType& GetContainer(ModelPart::MeshType &rMesh);

template<>
ModelPart::NodesContainerType& GetContainer<ModelPart::NodesContainerType>(ModelPart::MeshType& rMesh)
{
    return rMesh.Nodes();
}

template<>
ModelPart::ConditionsContainerType& GetContainer<ModelPart::ConditionsContainerType>(ModelPart::MeshType& rMesh)
{
    return rMesh.Conditions();
}

template<>
ModelPart::ElementsContainerType& GetContainer<ModelPart::ElementsContainerType>(ModelPart::MeshType& rMesh)
{
    return rMesh.Elements();
}

template <class TContainerType>
const TContainerType& GetContainer(const ModelPart::MeshType &rMesh);

template<>
const ModelPart::NodesContainerType& GetContainer<ModelPart::NodesContainerType>(const ModelPart::MeshType& rMesh)
{
    return rMesh.Nodes();
}

template<>
const ModelPart::ConditionsContainerType& GetContainer<ModelPart::ConditionsContainerType>(const ModelPart::MeshType& rMesh)
{
    return rMesh.Conditions();
}

template<>
const ModelPart::ElementsContainerType& GetContainer<ModelPart::ElementsContainerType>(const ModelPart::MeshType& rMesh)
{
    return rMesh.Elements();
}

template <class TRawDataType, class TContainerExpressionType>
void Read(
    TContainerExpressionType& rContainerExpression,
    TRawDataType const* pBegin,
    const int NumberOfEntities,
    int const* pShapeBegin,
    const int ShapeSize)
{
    KRATOS_TRY

    const IndexType number_of_entities = NumberOfEntities;

    KRATOS_ERROR_IF_NOT(number_of_entities == rContainerExpression.GetContainer().size())
        << "Number of entities does not match with the local container size. [ "
           "NumberOfEntities = "
        << NumberOfEntities
        << ", local container size = " << rContainerExpression.GetContainer().size() << " ].\n";

    // Convert int indices to IndexType
    std::vector<IndexType> shape(ShapeSize);
    std::copy(pShapeBegin,
              pShapeBegin + ShapeSize,
              shape.begin());

    auto p_expression = LiteralFlatExpression<TRawDataType>::Create(number_of_entities, shape);
    rContainerExpression.SetExpression(p_expression);
    TRawDataType* data_itr = p_expression->begin();

    const IndexType flattened_size = rContainerExpression.GetExpression().GetItemComponentCount();
    const IndexType total_size = number_of_entities * flattened_size;

    IndexPartition<IndexType>(total_size).for_each([pBegin, data_itr](const IndexType Index) {
        data_itr[Index] = pBegin[Index];
    });

    KRATOS_CATCH("");
}

template <class TRawDataType, class TContainerExpressionType>
void MoveFrom(
    TContainerExpressionType& rContainerExpression,
    TRawDataType* pBegin,
    const int NumberOfEntities,
    int const* pShapeBegin,
    const int ShapeSize)
{
    KRATOS_TRY

    const IndexType number_of_entities = NumberOfEntities;

    KRATOS_ERROR_IF_NOT(number_of_entities == rContainerExpression.GetContainer().size())
        << "Number of entities does not match with the local container size. [ "
           "NumberOfEntities = "
        << NumberOfEntities
        << ", local container size = " << rContainerExpression.GetContainer().size() << " ].\n";

    std::vector<IndexType> shape(ShapeSize);
    std::copy(pShapeBegin, pShapeBegin + ShapeSize, shape.begin());

    auto p_expression = LiteralFlatExpression<TRawDataType>::Create(pBegin, number_of_entities, shape);
    rContainerExpression.SetExpression(p_expression);

    KRATOS_CATCH("");
}

} // namespace ContainerExpressionHelperUtilities

template <class TContainerType, class TMeshType>
ContainerExpression<TContainerType, TMeshType>::ContainerExpression(ModelPart& rModelPart)
    : mpExpression(),
      mpModelPart(&rModelPart)
{
}

template <class TContainerType, class TMeshType>
ContainerExpression<TContainerType, TMeshType>::ContainerExpression(
    const ContainerExpression& rOther)
    : mpExpression(rOther.mpExpression),
      mpModelPart(rOther.mpModelPart)
{
}
template <class TContainerType, class TMeshType>
void ContainerExpression<TContainerType, TMeshType>::CopyFrom(
    const ContainerExpression<TContainerType, TMeshType>& rOther)
{
    KRATOS_ERROR_IF(this->GetContainer().size() != rOther.GetContainer().size())
        << "Mismatching model parts found with different number of entities in CopyFrom operation.\n"
        << this->mpModelPart->FullName()
        << ", origin model part name: " << rOther.GetModelPart().FullName() << " ].\n";

    mpExpression = rOther.mpExpression;
}

template <class TContainerType, class TMeshType>
void ContainerExpression<TContainerType, TMeshType>::Read(
    double const* pBegin,
    const int NumberOfEntities,
    int const* pShapeBegin,
    const int ShapeSize)
{
    ContainerExpressionHelperUtilities::Read(*this, pBegin, NumberOfEntities, pShapeBegin, ShapeSize);
}

template <class TContainerType, class TMeshType>
void ContainerExpression<TContainerType, TMeshType>::Read(
    int const* pBegin,
    const int NumberOfEntities,
    int const* pShapeBegin,
    const int ShapeSize)
{
    ContainerExpressionHelperUtilities::Read(*this, pBegin, NumberOfEntities, pShapeBegin, ShapeSize);
}

template <class TContainerType, class TMeshType>
void ContainerExpression<TContainerType, TMeshType>::MoveFrom(
    double* pBegin,
    const int NumberOfEntities,
    int const* pShapeBegin,
    const int ShapeSize)
{
    ContainerExpressionHelperUtilities::MoveFrom(*this, pBegin, NumberOfEntities, pShapeBegin, ShapeSize);
}

template <class TContainerType, class TMeshType>
void ContainerExpression<TContainerType, TMeshType>::MoveFrom(
    int* pBegin,
    const int NumberOfEntities,
    int const* pShapeBegin,
    const int ShapeSize)
{
    ContainerExpressionHelperUtilities::MoveFrom(*this, pBegin, NumberOfEntities, pShapeBegin, ShapeSize);
}

template <class TContainerType, class TMeshType>
void ContainerExpression<TContainerType, TMeshType>::Evaluate(
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

    KRATOS_ERROR_IF_NOT(shape == r_expression.GetItemShape())
        << "Shape mismatch. [ Requested shape  = " << shape
        << ", available shape = " << r_expression.GetItemShape() << " ].\n";

    const IndexType flattened_size = r_expression.GetItemComponentCount();

    IndexPartition<IndexType>(number_of_entities).for_each([pBegin, flattened_size, &r_expression](const IndexType EntityIndex) {
        const IndexType entity_data_begin_index = EntityIndex * flattened_size;
        double* p_input_data_begin = pBegin + entity_data_begin_index;
        for (IndexType i = 0; i < flattened_size; ++i) {
            *(p_input_data_begin+i) = r_expression.Evaluate(EntityIndex, entity_data_begin_index, i);
        }
    });

    KRATOS_CATCH("");
}

template <class TContainerType, class TMeshType>
void ContainerExpression<TContainerType, TMeshType>::SetDataToZero()
{
    mpExpression = LiteralExpression<double>::Create(0.0, GetContainer().size());
}

template <class TContainerType, class TMeshType>
void ContainerExpression<TContainerType, TMeshType>::SetExpression(Expression::Pointer pExpression)
{
    this->mpExpression = pExpression;
}

template <class TContainerType, class TMeshType>
bool ContainerExpression<TContainerType, TMeshType>::HasExpression() const
{
    return mpExpression.has_value();
}

template <class TContainerType, class TMeshType>
const Expression& ContainerExpression<TContainerType, TMeshType>::GetExpression() const
{
    return *(*mpExpression);
}


template <class TContainerType, class TMeshType>
const Expression::Pointer ContainerExpression<TContainerType, TMeshType>::pGetExpression() const
{
    return *mpExpression;
}

template <class TContainerType, class TMeshType>
const std::vector<std::size_t> ContainerExpression<TContainerType, TMeshType>::GetItemShape() const
{
    return this->GetExpression().GetItemShape();
}

template <class TContainerType, class TMeshType>
std::size_t ContainerExpression<TContainerType, TMeshType>::GetItemComponentCount() const
{
    return this->GetExpression().GetItemComponentCount();
}

template <class TContainerType, class TMeshType>
ModelPart& ContainerExpression<TContainerType, TMeshType>::GetModelPart()
{
    return *mpModelPart;
}

template <class TContainerType, class TMeshType>
const ModelPart& ContainerExpression<TContainerType, TMeshType>::GetModelPart() const
{
    return *mpModelPart;
}

template <class TContainerType, class TMeshType>
TContainerType& ContainerExpression<TContainerType, TMeshType>::GetContainer()
{
    return ContainerExpressionHelperUtilities::GetContainer<TContainerType>(ContainerExpressionHelperUtilities::GetMesh<TMeshType>(*mpModelPart));
}

template <class TContainerType, class TMeshType>
const TContainerType& ContainerExpression<TContainerType, TMeshType>::GetContainer() const
{
    return ContainerExpressionHelperUtilities::GetContainer<TContainerType>(ContainerExpressionHelperUtilities::GetMesh<TMeshType>(*mpModelPart));
}

template <class TContainerType, class TMeshType>
std::string ContainerExpression<TContainerType, TMeshType>::Info() const
{
    std::stringstream msg;

    msg << "VariableDataInfo: [ ModelPartName: "
        << this->GetModelPart().FullName()
        << ", Number of entities: " << this->GetContainer().size();

    if (mpExpression.has_value()) {
        msg << ", DataDimension: " << this->GetItemComponentCount();
        msg << ", Expression: " << this->GetExpression();
    } else {
        msg << ", Expression: not initialized";
    }

    msg << " ].\n";

    return msg.str();
}

template <class TContainerType, class TMeshType>
std::string ContainerExpression<TContainerType, TMeshType>::PrintData() const
{
    std::stringstream msg;
    msg << this->Info();
    return msg.str();
}

// template instantiations
template class ContainerExpression<ModelPart::NodesContainerType, MeshType::Local>;
template class ContainerExpression<ModelPart::ConditionsContainerType, MeshType::Local>;
template class ContainerExpression<ModelPart::ElementsContainerType, MeshType::Local>;

template class ContainerExpression<ModelPart::NodesContainerType, MeshType::Ghost>;
template class ContainerExpression<ModelPart::ConditionsContainerType, MeshType::Ghost>;
template class ContainerExpression<ModelPart::ElementsContainerType, MeshType::Ghost>;

template class ContainerExpression<ModelPart::NodesContainerType, MeshType::Interface>;
template class ContainerExpression<ModelPart::ConditionsContainerType, MeshType::Interface>;
template class ContainerExpression<ModelPart::ElementsContainerType, MeshType::Interface>;

} // namespace Kratos