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

LiteralDoubleExpression::LiteralDoubleExpression(const double Value)
    : mValue(Value)
{
}
double LiteralDoubleExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return mValue;
}

IndexType LiteralDoubleExpression::GetDimension() const
{
    return 1;
}

LiteralArray3Expression::LiteralArray3Expression(
    const array_1d<double, 3>& Value,
    const IndexType Dimension)
    : mValue(Value),
      mDimension(Dimension)
{
}

double LiteralArray3Expression::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return mValue[ComponentIndex];
}

IndexType LiteralArray3Expression::GetDimension() const
{
    return mDimension;
}

LiteralVectorExpression::LiteralVectorExpression(
    std::shared_ptr<Vector> pValue,
    const IndexType Dimension)
    : mpValue(pValue),
      mDimension(Dimension)
{
}

double LiteralVectorExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return (*mpValue)[EntityIndex * mDimension + ComponentIndex];
}

IndexType LiteralVectorExpression::GetDimension() const
{
    return mDimension;
}

BinaryAddExpression::BinaryAddExpression(
    std::shared_ptr<Expression> pLeft,
    std::shared_ptr<Expression> pRight)
    : mpLeft(pLeft),
      mpRight(pRight)
{
    KRATOS_ERROR_IF_NOT(mpLeft->GetDimension() == mpRight->GetDimension() || mpRight->GetDimension() == 1)
        << "Addition should have equal dimensions in left and right side expressions or right hand side should be with dimension = 1. [ Left expresion dimension = "
        << mpLeft->GetDimension() << ", Right expresion dimensions = " << mpRight->GetDimension() << " ].\n";

    KRATOS_ERROR_IF_NOT(mpLeft.get())
        << "Addition is provided with uninitialized left hand side "
           "expression.\n";

    KRATOS_ERROR_IF_NOT(mpRight.get())
        << "Addition is provided with uninitialized right hand side "
           "expression.\n";
}

double BinaryAddExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return mpLeft->Evaluate(EntityIndex, ComponentIndex) + mpRight->Evaluate(EntityIndex, ComponentIndex);
}

IndexType BinaryAddExpression::GetDimension() const
{
    return this->mpLeft->GetDimension();
}

BinarySubstractExpression::BinarySubstractExpression(
    std::shared_ptr<Expression> pLeft,
    std::shared_ptr<Expression> pRight)
    : mpLeft(pLeft),
      mpRight(pRight)
{
    KRATOS_ERROR_IF_NOT(mpLeft->GetDimension() == mpRight->GetDimension() || mpRight->GetDimension() == 1)
        << "Substraction should have equal dimensions in left and right side expressions or right hand side should be with dimension = 1. [ Left expresion dimension = "
        << mpLeft->GetDimension() << ", Right expresion dimensions = " << mpRight->GetDimension() << " ].\n";

    KRATOS_ERROR_IF_NOT(mpLeft.get())
        << "Substraction is provided with uninitialized left hand side "
           "expression.\n";

    KRATOS_ERROR_IF_NOT(mpRight.get())
        << "Substraction is provided with uninitialized right hand side "
           "expression.\n";
}

double BinarySubstractExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return mpLeft->Evaluate(EntityIndex, ComponentIndex) - mpRight->Evaluate(EntityIndex, ComponentIndex);
}

IndexType BinarySubstractExpression::GetDimension() const
{
    return this->mpLeft->GetDimension();
}

BinaryMultiplyExpression::BinaryMultiplyExpression(
    std::shared_ptr<Expression> pLeft,
    std::shared_ptr<Expression> pRight)
    : mpLeft(pLeft),
      mpRight(pRight)
{
    KRATOS_ERROR_IF_NOT(mpLeft->GetDimension() == mpRight->GetDimension() || mpRight->GetDimension() == 1)
        << "Multiplication should have equal dimensions in left and right side expressions or right hand side should be with dimension = 1. [ Left expresion dimension = "
        << mpLeft->GetDimension() << ", Right expresion dimensions = " << mpRight->GetDimension() << " ].\n";

    KRATOS_ERROR_IF_NOT(mpLeft.get())
        << "Multiplication is provided with uninitialized left hand side "
           "expression.\n";

    KRATOS_ERROR_IF_NOT(mpRight.get())
        << "Multiplication is provided with uninitialized right hand side "
           "expression.\n";
}

double BinaryMultiplyExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return mpLeft->Evaluate(EntityIndex, ComponentIndex % this->mpLeft->GetDimension()) * mpRight->Evaluate(EntityIndex, ComponentIndex % this->mpRight->GetDimension());
}

IndexType BinaryMultiplyExpression::GetDimension() const
{
    return this->mpLeft->GetDimension();
}

BinaryDivideExpression::BinaryDivideExpression(
    std::shared_ptr<Expression> pLeft,
    std::shared_ptr<Expression> pRight)
    : mpLeft(pLeft),
      mpRight(pRight)
{
    KRATOS_ERROR_IF_NOT(mpLeft->GetDimension() == mpRight->GetDimension() || mpRight->GetDimension() == 1)
        << "Division should have equal dimensions in left and right side expressions or right hand side should be with dimension = 1. [ Left expresion dimension = "
        << mpLeft->GetDimension() << ", Right expresion dimensions = " << mpRight->GetDimension() << " ].\n";

    KRATOS_ERROR_IF_NOT(mpLeft.get())
        << "Division is provided with uninitialized left hand side "
           "expression.\n";

    KRATOS_ERROR_IF_NOT(mpRight.get())
        << "Division is provided with uninitialized right hand side "
           "expression.\n";
}

double BinaryDivideExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return mpLeft->Evaluate(EntityIndex, ComponentIndex % this->mpLeft->GetDimension()) / mpRight->Evaluate(EntityIndex, ComponentIndex % this->mpRight->GetDimension());
}

IndexType BinaryDivideExpression::GetDimension() const
{
    return this->mpLeft->GetDimension();
}

BinaryPowerExpression::BinaryPowerExpression(
    std::shared_ptr<Expression> pLeft,
    std::shared_ptr<Expression> pRight)
    : mpLeft(pLeft),
      mpRight(pRight)
{
    KRATOS_ERROR_IF_NOT(mpLeft->GetDimension() == mpRight->GetDimension() || mpRight->GetDimension() == 1)
        << "Power should have equal dimensions in left and right side expressions or right hand side should be with dimension = 1. [ Left expresion dimension = "
        << mpLeft->GetDimension() << ", Right expresion dimensions = " << mpRight->GetDimension() << " ].\n";

    KRATOS_ERROR_IF_NOT(mpLeft.get())
        << "Power is provided with uninitialized left hand side "
           "expression.\n";

    KRATOS_ERROR_IF_NOT(mpRight.get())
        << "Power is provided with uninitialized right hand side "
           "expression.\n";
}

double BinaryPowerExpression::Evaluate(
    const IndexType EntityIndex,
    const IndexType ComponentIndex) const
{
    return std::pow(mpLeft->Evaluate(EntityIndex, ComponentIndex), mpRight->Evaluate(EntityIndex, ComponentIndex));
}

IndexType BinaryPowerExpression::GetDimension() const
{
    return this->mpLeft->GetDimension();
}

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
    mpExpression = std::make_shared<LiteralDoubleExpression>(0.0);
}

template <class TContainerType>
void ContainerVariableDataHolderBase<TContainerType>::SetExpression(std::shared_ptr<Expression> pExpression)
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
        << ", Data dimension: " << this->GetDataDimension() << " ].\n";

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