//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes
#include <sstream>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "container_variable_data_holder_base.h"

namespace Kratos {

template <class TContainerType>
ContainerVariableDataHolderBase<TContainerType>::ContainerVariableDataHolderBase(ModelPart& rModelPart)
    : mpModelPart(&rModelPart)
{
}

template <class TContainerType>
ContainerVariableDataHolderBase<TContainerType>::ContainerVariableDataHolderBase(
    const ContainerVariableDataHolderBase& rOther)
    : mpExpression(rOther.mpExpression),
      mpModelPart(rOther.mpModelPart)
{
}
template <class TContainerType>
void ContainerVariableDataHolderBase<TContainerType>::CopyDataFrom(
    const ContainerVariableDataHolderBase<TContainerType>& rOther)
{
    mpExpression = rOther.mpExpression;
}

template <class TContainerType>
void ContainerVariableDataHolderBase<TContainerType>::SetDataToZero()
{
    mpExpression = Kratos::make_intrusive<LiteralDoubleExpression>(0.0);
}

template <class TContainerType>
void ContainerVariableDataHolderBase<TContainerType>::SetExpression(Expression::Pointer pExpression)
{
    this->mpExpression = pExpression;
}

template <class TContainerType>
const Expression& ContainerVariableDataHolderBase<TContainerType>::GetExpression() const
{
    if (mpExpression) {
        return *mpExpression;
    } else {
        KRATOS_ERROR << "Uninitialized expression.\n";
    }
}


template <class TContainerType>
const Expression::Pointer ContainerVariableDataHolderBase<TContainerType>::pGetExpression() const
{
    if (mpExpression) {
        return mpExpression;
    } else {
        KRATOS_ERROR << "Uninitialized expression.\n";
    }
}

template <class TContainerType>
IndexType ContainerVariableDataHolderBase<TContainerType>::GetDataDimension() const
{
    return this->GetExpression().GetDimension();
}

template <class TContainerType>
ModelPart& ContainerVariableDataHolderBase<TContainerType>::GetModelPart()
{
    return *mpModelPart;
}

template <class TContainerType>
const ModelPart& ContainerVariableDataHolderBase<TContainerType>::GetModelPart() const
{
    return *mpModelPart;
}

template <>
ModelPart::NodesContainerType& ContainerVariableDataHolderBase<ModelPart::NodesContainerType>::GetContainer()
{
    return mpModelPart->GetCommunicator().LocalMesh().Nodes();
}

template <>
ModelPart::ConditionsContainerType& ContainerVariableDataHolderBase<ModelPart::ConditionsContainerType>::GetContainer()
{
    return mpModelPart->GetCommunicator().LocalMesh().Conditions();
}

template <>
ModelPart::ElementsContainerType& ContainerVariableDataHolderBase<ModelPart::ElementsContainerType>::GetContainer()
{
    return mpModelPart->GetCommunicator().LocalMesh().Elements();
}

template <>
const ModelPart::NodesContainerType& ContainerVariableDataHolderBase<ModelPart::NodesContainerType>::GetContainer() const
{
    return mpModelPart->GetCommunicator().LocalMesh().Nodes();
}

template <>
const ModelPart::ConditionsContainerType& ContainerVariableDataHolderBase<ModelPart::ConditionsContainerType>::GetContainer() const
{
    return mpModelPart->GetCommunicator().LocalMesh().Conditions();
}

template <>
const ModelPart::ElementsContainerType& ContainerVariableDataHolderBase<ModelPart::ElementsContainerType>::GetContainer() const
{
    return mpModelPart->GetCommunicator().LocalMesh().Elements();
}

template <class TContainerType>
std::string ContainerVariableDataHolderBase<TContainerType>::Info() const
{
    std::stringstream msg;

    msg << "ContainerVariableDataHolderInfo: [ ModelPartName: "
        << this->GetModelPart().FullName()
        << ", Number of entities: " << this->GetContainer().size()
        << ", Data dimension: " << (mpExpression.get() ? this->GetDataDimension() : 0) << " ].\n";

    return msg.str();
}

template <class TContainerType>
std::string ContainerVariableDataHolderBase<TContainerType>::PrintData() const
{
    std::stringstream msg;
    msg << this->Info();
    return msg.str();
}

// template instantiations
template class ContainerVariableDataHolderBase<ModelPart::NodesContainerType>;
template class ContainerVariableDataHolderBase<ModelPart::ConditionsContainerType>;
template class ContainerVariableDataHolderBase<ModelPart::ElementsContainerType>;

} // namespace Kratos