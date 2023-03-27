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
#include "custom_utilities/container_variable_data/expressions/literal/literal_expression.h"

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
void ContainerVariableData<TContainerType>::CopyDataFrom(
    const ContainerVariableData<TContainerType>& rOther)
{
    KRATOS_ERROR_IF(this->mpModelPart != &rOther.GetModelPart())
        << "Model part mismatch. [ Destination model part name: "
        << this->mpModelPart->FullName()
        << ", origin model part name: " << rOther.GetModelPart().FullName() << " ].\n";

    mpExpression = rOther.mpExpression;
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
std::size_t ContainerVariableData<TContainerType>::GetLocalSize() const
{
    return this->GetExpression().GetLocalSize();
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

    msg << "ContainerVariableDataHolderInfo: [ ModelPartName: "
        << this->GetModelPart().FullName()
        << ", Number of entities: " << this->GetContainer().size()
        << ", Data dimension: " << (mpExpression.has_value() ? this->GetLocalSize() : 0) << " ].\n";

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