//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes
#include <sstream>
#include <limits>
#include <cmath>

// Project includes
#include "includes/data_communicator.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "custom_utilities/optimization_utils.h"

// Include base h
#include "container_variable_data_holder.h"

namespace Kratos {

namespace ContainerVariableDataHolderHelperUtilities
{

void AssignValueToVector(
    Vector& rOutput,
    const IndexType StartIndex,
    const IndexType ComponentIndex,
    const double& rValue)
{
    rOutput[StartIndex] = rValue;
}

void AssignValueToVector(
    Vector& rOutput,
    const IndexType StartIndex,
    const IndexType ComponentIndex,
    const array_1d<double, 3>& rValue)
{
    rOutput[StartIndex + ComponentIndex] = rValue[ComponentIndex];
}

void AssignValueFromVector(
    double& rOutput,
    const IndexType ValueComponentIndex,
    const IndexType VectoStartingIndex,
    const Vector& rInput)
{
    rOutput = rInput[VectoStartingIndex];
}

void AssignValueFromVector(
    array_1d<double, 3>& rOutput,
    const IndexType ValueComponentIndex,
    const IndexType VectoStartingIndex,
    const Vector& rInput)
{
    rOutput[ValueComponentIndex] = rInput[VectoStartingIndex + ValueComponentIndex];
}

} // namespace ContainerVariableDataHolderHelperUtilities

ContainerVariableDataHolderBase::ContainerVariableDataHolderBase(
    ModelPart& rModelPart,
    const ContainerVariableDataHolderType& rContainerVariableDataHolderType)
    : mrModelPart(rModelPart),
      mDataDimension(0),
      mContainerVariableDataHolderType(rContainerVariableDataHolderType)
{
}

ContainerVariableDataHolderBase::ContainerVariableDataHolderBase(const ContainerVariableDataHolderBase& rOther)
    : mrModelPart(rOther.mrModelPart),
      mContainerVariableDataHolderType(rOther.mContainerVariableDataHolderType)
{
    this->mData.resize(rOther.mData.size(), false);
    this->mDataDimension = rOther.mDataDimension;
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index){
        this->mData[Index] = rOther.mData[Index];
    });
}

IndexType ContainerVariableDataHolderBase::GetDataDimension() const
{
    return mDataDimension;
}

Vector& ContainerVariableDataHolderBase::GetData()
{
    return mData;
}

const Vector& ContainerVariableDataHolderBase::GetData() const
{
    return mData;
}

void ContainerVariableDataHolderBase::CopyDataFrom(const ContainerVariableDataHolderBase& rOther)
{
    KRATOS_ERROR_IF(this->mrModelPart != rOther.mrModelPart)
        << "Modelpart mismatch between origin and destination for copy.\n"
        << "   Destination = " << *this << "\n"
        << "   Origin = " << rOther << "\n.";

    KRATOS_ERROR_IF_NOT(this->IsCompatibleWithContainerVariableDataHolder(rOther))
        << "Copying from unsupported data container types. Followings are supported:\n"
        << "   Copying from same data container type\n"
        << "   Nodal historical         <-> nodal non-historical\n"
        << "   Condition non-historical <-> condition properties\n"
        << "   Element non-historical   <-> element properties\n"
        << "\n"
        << "   Destination data: " << *this << "\n"
        << "   Origin data     : " << rOther << "\n";

    if (this->mData.size() != rOther.mData.size()) {
        this->mData.resize(rOther.mData.size(), false);
    }

    this->mDataDimension = rOther.mDataDimension;

    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index){
        this->mData[Index] = rOther.mData[Index];
    });
}

ModelPart& ContainerVariableDataHolderBase::GetModelPart()
{
    return mrModelPart;
}

const ModelPart& ContainerVariableDataHolderBase::GetModelPart() const
{
    return mrModelPart;
}

bool ContainerVariableDataHolderBase::IsCompatibleWithContainerVariableDataHolder(const ContainerVariableDataHolderBase& rOther) const
{
    const bool valid_same_type_conversion = (this->mContainerVariableDataHolderType == rOther.mContainerVariableDataHolderType);

    const bool valid_nodal_different_conversion =
        (this->mContainerVariableDataHolderType == ContainerVariableDataHolderType::HistoricalContainerVariableDataHolder && rOther.mContainerVariableDataHolderType == ContainerVariableDataHolderType::NodalContainerVariableDataHolder) ||
        (this->mContainerVariableDataHolderType == ContainerVariableDataHolderType::NodalContainerVariableDataHolder && rOther.mContainerVariableDataHolderType == ContainerVariableDataHolderType::HistoricalContainerVariableDataHolder);

    const bool valid_condition_different_conversion =
        (this->mContainerVariableDataHolderType == ContainerVariableDataHolderType::ConditionContainerVariableDataHolder && rOther.mContainerVariableDataHolderType == ContainerVariableDataHolderType::ConditionPropertiesContainerVariableDataHolder) ||
        (this->mContainerVariableDataHolderType == ContainerVariableDataHolderType::ConditionPropertiesContainerVariableDataHolder && rOther.mContainerVariableDataHolderType == ContainerVariableDataHolderType::ConditionContainerVariableDataHolder);

    const bool valid_element_different_conversion =
        (this->mContainerVariableDataHolderType == ContainerVariableDataHolderType::ElementContainerVariableDataHolder && rOther.mContainerVariableDataHolderType == ContainerVariableDataHolderType::ElementPropertiesContainerVariableDataHolder) ||
        (this->mContainerVariableDataHolderType == ContainerVariableDataHolderType::ElementPropertiesContainerVariableDataHolder && rOther.mContainerVariableDataHolderType == ContainerVariableDataHolderType::ElementContainerVariableDataHolder);

    return (this->mrModelPart == rOther.mrModelPart) && (valid_same_type_conversion ||  valid_nodal_different_conversion || valid_condition_different_conversion || valid_element_different_conversion);
}

std::string ContainerVariableDataHolderBase::Info() const
{
    KRATOS_ERROR << "Calling ContainerVariableDataHolderBase::Info. This should "
                    "be implemented in the derrived class.";
    return "";
}

void ContainerVariableDataHolderBase::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

template<class TDataType>
TDataType& HistoricalDataValueContainer::GetValue(
    typename ContainerType::data_type& rEntity,
    const Variable<TDataType>& rVariable)
{
    return rEntity.FastGetSolutionStepValue(rVariable);
}

template<class TDataType>
void HistoricalDataValueContainer::SetValue(
    typename ContainerType::data_type& rEntity,
    const Variable<TDataType>& rVariable,
    const TDataType& rValue)
{
    rEntity.FastGetSolutionStepValue(rVariable) = rValue;
}

template<class TContainerType>
template<class TDataType>
TDataType& NonHistoricalDataValueContainer<TContainerType>::GetValue(
    typename ContainerType::data_type& rEntity,
    const Variable<TDataType>& rVariable)
{
    return rEntity.GetValue(rVariable);
}

template<class TContainerType>
template<class TDataType>
void NonHistoricalDataValueContainer<TContainerType>::SetValue(
    typename ContainerType::data_type& rEntity,
    const Variable<TDataType>& rVariable,
    const TDataType& rValue)
{
    rEntity.SetValue(rVariable, rValue);
}

template<class TContainerType>
template<class TDataType>
TDataType& PropertiesDataValueContainer<TContainerType>::GetValue(
    typename ContainerType::data_type& rEntity,
    const Variable<TDataType>& rVariable)
{
    return rEntity.GetProperties().GetValue(rVariable);
}

template<class TContainerType>
template<class TDataType>
void PropertiesDataValueContainer<TContainerType>::SetValue(
    typename ContainerType::data_type& rEntity,
    const Variable<TDataType>& rVariable,
    const TDataType& rValue)
{
    rEntity.GetProperties().SetValue(rVariable, rValue);
}

template<class TContainerVariableDataHolderType>
ContainerVariableDataHolder<TContainerVariableDataHolderType> ContainerVariableDataHolder<TContainerVariableDataHolderType>::Clone()
{
    return ContainerVariableDataHolder<TContainerVariableDataHolderType>(*this);
}

template<class TContainerVariableDataHolderType>
template<class TDataType>
void ContainerVariableDataHolder<TContainerVariableDataHolderType>::ReadDataFromContainerVariable(const Variable<TDataType>& rVariable)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(this->GetModelPart().GetProcessInfo().Has(DOMAIN_SIZE))
        << "DOMAIN_SIZE variable is not found in process info of model part "
           "with name "
        << this->GetModelPart().FullName() << ".\n";

    const IndexType number_of_entities = this->GetContainer().size();
    mDataDimension = OptimizationUtils::GetLocalSize<TDataType>(this->GetModelPart().GetProcessInfo()[DOMAIN_SIZE]);

    auto& r_data = this->mData;
    auto& r_container = this->GetContainer();

    r_data.resize(number_of_entities * mDataDimension);

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index){
        const auto& values = TContainerVariableDataHolderType::GetValue(*(r_container.begin() + Index), rVariable);
        for (IndexType i = 0; i < mDataDimension; ++i) {
            ContainerVariableDataHolderHelperUtilities::AssignValueToVector(r_data, Index * mDataDimension, i, values);
        }
    });

    KRATOS_CATCH("")
}

template<class TContainerVariableDataHolderType>
template<class TDataType>
void ContainerVariableDataHolder<TContainerVariableDataHolderType>::AssignDataToContainerVariable(const Variable<TDataType>& rVariable)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(this->GetModelPart().GetProcessInfo().Has(DOMAIN_SIZE))
        << "DOMAIN_SIZE variable is not found in process info of model part "
           "with name "
        << this->GetModelPart().FullName() << ".\n";

    const IndexType local_size = OptimizationUtils::GetLocalSize<TDataType>(this->GetModelPart().GetProcessInfo()[DOMAIN_SIZE]);

    KRATOS_ERROR_IF_NOT(local_size == this->GetDataDimension())
        << "Stored data dimension and requested assignment variable data "
           "dimension mismatch. [ Sotred data dimension = "
        << this->GetDataDimension()
        << ", assignment requested variable data diemension = " << local_size << " ].\n";

    const IndexType number_of_entities = this->GetContainer().size();

    const auto& r_values = this->mData;
    auto& r_container = this->GetContainer();

    KRATOS_ERROR_IF_NOT(number_of_entities * local_size == r_values.size())
        << "Container size and values size mismatch. [ Values size = "
        << r_values.size() << ", container size = " << number_of_entities
        << ", local size = " << local_size << " ].\n";

    // initialize the container variables first
    if constexpr(std::is_same_v<TContainerVariableDataHolderType, NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>) {
        // initializes ghost nodes as well for the later synchronization
        // only, the nodal non historical values needs to be set unless
        // they are properly initialized.
        block_for_each(this->mrModelPart.Nodes(), [&](auto& rEntity) {
            TContainerVariableDataHolderType::SetValue(rEntity, rVariable, rVariable.Zero());
        });
    } else {
        // no ghost settings required for elements or conditions nor historical values, hence no special treatment required.
        block_for_each(r_container, [&](auto& rEntity) {
            TContainerVariableDataHolderType::SetValue(rEntity, rVariable, rVariable.Zero());
        });
    }

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index){
        auto& values = TContainerVariableDataHolderType::GetValue(*(r_container.begin() + Index), rVariable);
        for (IndexType i = 0; i < local_size; ++i) {
            ContainerVariableDataHolderHelperUtilities::AssignValueFromVector(values, i, Index * local_size, r_values);
        }
    });

    // synchronize nodal values
    auto& r_communicator = this->mrModelPart.GetCommunicator();

    if constexpr(std::is_same_v<TContainerVariableDataHolderType, HistoricalDataValueContainer>) {
        r_communicator.SynchronizeNodalSolutionStepsData();
    } else if constexpr(std::is_same_v<TContainerVariableDataHolderType, NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>) {
        r_communicator.SynchronizeNonHistoricalVariable(rVariable);
    }

    KRATOS_CATCH("");
}

template<class TContainerVariableDataHolderType>
template<class TDataType>
void ContainerVariableDataHolder<TContainerVariableDataHolderType>::SetDataForContainerVariable(const Variable<TDataType>& rVariable, const TDataType& rValue)
{

    KRATOS_ERROR_IF_NOT(mrModelPart.GetProcessInfo().Has(DOMAIN_SIZE))
        << "DOMAIN_SIZE variable is not found in process info of model part "
           "with name "
        << mrModelPart.FullName() << ".\n";

    mDataDimension = OptimizationUtils::GetLocalSize<TDataType>(mrModelPart.GetProcessInfo()[DOMAIN_SIZE]);

    const IndexType number_of_entities = this->GetContainer().size();

    if (this->mData.size() != mDataDimension * number_of_entities) {
        this->mData.resize(mDataDimension * number_of_entities, false);
    }

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index) {
        for (IndexType i = 0; i < mDataDimension; ++i) {
            ContainerVariableDataHolderHelperUtilities::AssignValueToVector(this->mData, Index * mDataDimension, i, rValue);
        }
    });
}

template<class TContainerVariableDataHolderType>
ContainerVariableDataHolder<TContainerVariableDataHolderType> ContainerVariableDataHolder<TContainerVariableDataHolderType>::operator+(const ContainerVariableDataHolder<TContainerVariableDataHolderType>& rOther) const
{
    KRATOS_ERROR_IF(this->mrModelPart != rOther.mrModelPart)
        << "Mismatching model parts found in addition.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    KRATOS_ERROR_IF_NOT(this->mData.size() == rOther.mData.size())
        << "Data size mismatch in operands for +.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    ContainerVariableDataHolder<TContainerVariableDataHolderType> result(this->mrModelPart);
    result.mDataDimension = this->mDataDimension;
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] + rOther.mData[Index];
    });

    return result;
}

template<class TContainerVariableDataHolderType>
ContainerVariableDataHolder<TContainerVariableDataHolderType>& ContainerVariableDataHolder<TContainerVariableDataHolderType>::operator+=(const ContainerVariableDataHolder<TContainerVariableDataHolderType>& rOther)
{
    KRATOS_ERROR_IF(this->mrModelPart != rOther.mrModelPart)
        << "Mismatching model parts found in addition.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    KRATOS_ERROR_IF_NOT(this->mData.size() == rOther.mData.size())
        << "Data size mismatch in operands for +=.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] += rOther.mData[Index];
    });

    return *this;
}

template<class TContainerVariableDataHolderType>
ContainerVariableDataHolder<TContainerVariableDataHolderType> ContainerVariableDataHolder<TContainerVariableDataHolderType>::operator+(const double Value) const
{
    ContainerVariableDataHolder<TContainerVariableDataHolderType> result(this->mrModelPart);
    result.mDataDimension = this->mDataDimension;
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] + Value;
    });

    return result;
}

template<class TContainerVariableDataHolderType>
ContainerVariableDataHolder<TContainerVariableDataHolderType>& ContainerVariableDataHolder<TContainerVariableDataHolderType>::operator+=(const double Value)
{
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] += Value;
    });

    return *this;
}

template<class TContainerVariableDataHolderType>
ContainerVariableDataHolder<TContainerVariableDataHolderType> ContainerVariableDataHolder<TContainerVariableDataHolderType>::operator-(const ContainerVariableDataHolder<TContainerVariableDataHolderType>& rOther) const
{
    KRATOS_ERROR_IF(this->mrModelPart != rOther.mrModelPart)
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    KRATOS_ERROR_IF_NOT(this->mData.size() == rOther.mData.size())
        << "Data size mismatch in operands for -.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    ContainerVariableDataHolder<TContainerVariableDataHolderType> result(this->mrModelPart);
    result.mDataDimension = this->mDataDimension;
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] - rOther.mData[Index];
    });

    return result;
}

template<class TContainerVariableDataHolderType>
ContainerVariableDataHolder<TContainerVariableDataHolderType>& ContainerVariableDataHolder<TContainerVariableDataHolderType>::operator-=(const ContainerVariableDataHolder<TContainerVariableDataHolderType>& rOther)
{
    KRATOS_ERROR_IF(this->mrModelPart != rOther.mrModelPart)
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    KRATOS_ERROR_IF_NOT(this->mData.size() == rOther.mData.size())
        << "Data size mismatch in operands for -=.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] -= rOther.mData[Index];
    });

    return *this;
}

template<class TContainerVariableDataHolderType>
ContainerVariableDataHolder<TContainerVariableDataHolderType> ContainerVariableDataHolder<TContainerVariableDataHolderType>::operator-(const double Value) const
{
    ContainerVariableDataHolder<TContainerVariableDataHolderType> result(this->mrModelPart);
    result.mDataDimension = this->mDataDimension;
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] - Value;
    });

    return result;
}

template<class TContainerVariableDataHolderType>
ContainerVariableDataHolder<TContainerVariableDataHolderType>& ContainerVariableDataHolder<TContainerVariableDataHolderType>::operator-=(const double Value)
{
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] -= Value;
    });

    return *this;
}

template<class TContainerVariableDataHolderType>
ContainerVariableDataHolder<TContainerVariableDataHolderType> ContainerVariableDataHolder<TContainerVariableDataHolderType>::operator*(const double Value) const
{
    ContainerVariableDataHolder<TContainerVariableDataHolderType> result(this->mrModelPart);
    result.mDataDimension = this->mDataDimension;
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] * Value;
    });

    return result;
}

template<class TContainerVariableDataHolderType>
ContainerVariableDataHolder<TContainerVariableDataHolderType>& ContainerVariableDataHolder<TContainerVariableDataHolderType>::operator*=(const double Value)
{
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] *= Value;
    });

    return *this;
}

template<class TContainerVariableDataHolderType>
ContainerVariableDataHolder<TContainerVariableDataHolderType> ContainerVariableDataHolder<TContainerVariableDataHolderType>::operator/(const double Value) const
{
    KRATOS_ERROR_IF(std::abs(Value) < std::numeric_limits<double>::epsilon())
        << "Division by zero.\n"
        << "      Left operand data : " << *this << "\n"
        << "      divisor           : " << Value << "\n";

    ContainerVariableDataHolder<TContainerVariableDataHolderType> result(this->mrModelPart);
    result.mDataDimension = this->mDataDimension;
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] / Value;
    });

    return result;
}

template<class TContainerVariableDataHolderType>
ContainerVariableDataHolder<TContainerVariableDataHolderType>& ContainerVariableDataHolder<TContainerVariableDataHolderType>::operator/=(const double Value)
{
    KRATOS_ERROR_IF(std::abs(Value) < std::numeric_limits<double>::epsilon())
        << "Division by zero.\n"
        << "      Left operand data : " << *this << "\n"
        << "      divisor           : " << Value << "\n";

    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] /= Value;
    });

    return *this;
}

template<class TContainerVariableDataHolderType>
ContainerVariableDataHolder<TContainerVariableDataHolderType> ContainerVariableDataHolder<TContainerVariableDataHolderType>::operator^(const double Value) const
{
    ContainerVariableDataHolder<TContainerVariableDataHolderType> result(this->mrModelPart);
    result.mDataDimension = this->mDataDimension;
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = std::pow(this->mData[Index], Value);
    });

    return result;
}

template<class TContainerVariableDataHolderType>
ContainerVariableDataHolder<TContainerVariableDataHolderType>& ContainerVariableDataHolder<TContainerVariableDataHolderType>::operator^=(const double Value)
{
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] = std::pow(this->mData[Index], Value);
    });

    return *this;
}

template<class TContainerVariableDataHolderType>
ContainerVariableDataHolder<TContainerVariableDataHolderType>& ContainerVariableDataHolder<TContainerVariableDataHolderType>::operator=(const ContainerVariableDataHolder<TContainerVariableDataHolderType>& rOther)
{
    KRATOS_ERROR_IF(this->mrModelPart != rOther.mrModelPart)
        << "Mismatching model parts found in assignment.\n"
        << "      Assignee data: " << *this << "\n"
        << "      Assignor data: " << rOther << "\n";

    if (this->mData.size() != rOther.mData.size()) {
        this->mData.resize(rOther.mData.size(), false);
    }

    this->mDataDimension = rOther.mDataDimension;

    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] = rOther.mData[Index];
    });

    return *this;
}

template<class TContainerVariableDataHolderType>
typename TContainerVariableDataHolderType::ContainerType& ContainerVariableDataHolder<TContainerVariableDataHolderType>::GetContainer()
{
    using container_type = typename TContainerVariableDataHolderType::ContainerType;

    if constexpr(std::is_same_v<container_type, ModelPart::NodesContainerType>) {
        return this->GetModelPart().GetCommunicator().LocalMesh().Nodes();
    } else if constexpr(std::is_same_v<container_type, ModelPart::ConditionsContainerType>) {
        return this->GetModelPart().GetCommunicator().LocalMesh().Conditions();
    } else if constexpr(std::is_same_v<container_type, ModelPart::ElementsContainerType>) {
        return this->GetModelPart().GetCommunicator().LocalMesh().Elements();
    }
}

template<class TContainerVariableDataHolderType>
const typename TContainerVariableDataHolderType::ContainerType& ContainerVariableDataHolder<TContainerVariableDataHolderType>::GetContainer() const
{
    using container_type = typename TContainerVariableDataHolderType::ContainerType;

    if constexpr(std::is_same_v<container_type, ModelPart::NodesContainerType>) {
        return this->GetModelPart().GetCommunicator().LocalMesh().Nodes();
    } else if constexpr(std::is_same_v<container_type, ModelPart::ConditionsContainerType>) {
        return this->GetModelPart().GetCommunicator().LocalMesh().Conditions();
    } else if constexpr(std::is_same_v<container_type, ModelPart::ElementsContainerType>) {
        return this->GetModelPart().GetCommunicator().LocalMesh().Elements();
    }
}

template<class TContainerVariableDataHolderType>
std::string ContainerVariableDataHolder<TContainerVariableDataHolderType>::Info() const
{
    std::stringstream msg;

    if constexpr(std::is_same_v<TContainerVariableDataHolderType, HistoricalDataValueContainer>) {
        msg << "NodalHistoricalDataContainer";
    } else if constexpr(std::is_same_v<TContainerVariableDataHolderType, NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>) {
        msg << "NodalNonHistoricalDataContainer";
    } else if constexpr(std::is_same_v<TContainerVariableDataHolderType, NonHistoricalDataValueContainer<ModelPart::ConditionsContainerType>>) {
        msg << "ConditionNonHistoricalDataContainer";
    } else if constexpr(std::is_same_v<TContainerVariableDataHolderType, NonHistoricalDataValueContainer<ModelPart::ElementsContainerType>>) {
        msg << "ElementNonHistoricalDataContainer";
    } else if constexpr(std::is_same_v<TContainerVariableDataHolderType, PropertiesDataValueContainer<ModelPart::ConditionsContainerType>>) {
        msg << "ConditionPropertiesDataContainer";
    } else if constexpr(std::is_same_v<TContainerVariableDataHolderType, PropertiesDataValueContainer<ModelPart::ElementsContainerType>>) {
        msg << "ElementPropertiesDataContainer";
    }

    msg << " [ ";
    msg << "ModelPartName = " << this->GetModelPart().FullName() << ", ";
    msg << "Number of entities = " << this->GetContainer().size() << ", ";
    msg << "Data size = " << this->mData.size() << ", ";
    msg << "Data dimension = " << this->GetDataDimension() << " ].\n";

    return msg.str();
}

// template instantiations
template class ContainerVariableDataHolder<HistoricalDataValueContainer>;
template void ContainerVariableDataHolder<HistoricalDataValueContainer>::AssignDataToContainerVariable(const Variable<double>&);
template void ContainerVariableDataHolder<HistoricalDataValueContainer>::AssignDataToContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerVariableDataHolder<HistoricalDataValueContainer>::ReadDataFromContainerVariable(const Variable<double>&);
template void ContainerVariableDataHolder<HistoricalDataValueContainer>::ReadDataFromContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerVariableDataHolder<HistoricalDataValueContainer>::SetDataForContainerVariable(const Variable<double>&, const double&);
template void ContainerVariableDataHolder<HistoricalDataValueContainer>::SetDataForContainerVariable(const Variable<array_1d<double, 3>>&, const array_1d<double, 3>&);

template class ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>;
template void ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>::AssignDataToContainerVariable(const Variable<double>&);
template void ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>::AssignDataToContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>::ReadDataFromContainerVariable(const Variable<double>&);
template void ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>::ReadDataFromContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>::SetDataForContainerVariable(const Variable<double>&, const double&);
template void ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>::SetDataForContainerVariable(const Variable<array_1d<double, 3>>&, const array_1d<double, 3>&);

template class ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::ConditionsContainerType>>;
template void ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::ConditionsContainerType>>::AssignDataToContainerVariable(const Variable<double>&);
template void ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::ConditionsContainerType>>::AssignDataToContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::ConditionsContainerType>>::ReadDataFromContainerVariable(const Variable<double>&);
template void ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::ConditionsContainerType>>::ReadDataFromContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::ConditionsContainerType>>::SetDataForContainerVariable(const Variable<double>&, const double&);
template void ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::ConditionsContainerType>>::SetDataForContainerVariable(const Variable<array_1d<double, 3>>&, const array_1d<double, 3>&);

template class ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::ElementsContainerType>>;
template void ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::ElementsContainerType>>::AssignDataToContainerVariable(const Variable<double>&);
template void ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::ElementsContainerType>>::AssignDataToContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::ElementsContainerType>>::ReadDataFromContainerVariable(const Variable<double>&);
template void ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::ElementsContainerType>>::ReadDataFromContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::ElementsContainerType>>::SetDataForContainerVariable(const Variable<double>&, const double&);
template void ContainerVariableDataHolder<NonHistoricalDataValueContainer<ModelPart::ElementsContainerType>>::SetDataForContainerVariable(const Variable<array_1d<double, 3>>&, const array_1d<double, 3>&);

template class ContainerVariableDataHolder<PropertiesDataValueContainer<ModelPart::ConditionsContainerType>>;
template void ContainerVariableDataHolder<PropertiesDataValueContainer<ModelPart::ConditionsContainerType>>::AssignDataToContainerVariable(const Variable<double>&);
template void ContainerVariableDataHolder<PropertiesDataValueContainer<ModelPart::ConditionsContainerType>>::AssignDataToContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerVariableDataHolder<PropertiesDataValueContainer<ModelPart::ConditionsContainerType>>::ReadDataFromContainerVariable(const Variable<double>&);
template void ContainerVariableDataHolder<PropertiesDataValueContainer<ModelPart::ConditionsContainerType>>::ReadDataFromContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerVariableDataHolder<PropertiesDataValueContainer<ModelPart::ConditionsContainerType>>::SetDataForContainerVariable(const Variable<double>&, const double&);
template void ContainerVariableDataHolder<PropertiesDataValueContainer<ModelPart::ConditionsContainerType>>::SetDataForContainerVariable(const Variable<array_1d<double, 3>>&, const array_1d<double, 3>&);

template class ContainerVariableDataHolder<PropertiesDataValueContainer<ModelPart::ElementsContainerType>>;
template void ContainerVariableDataHolder<PropertiesDataValueContainer<ModelPart::ElementsContainerType>>::AssignDataToContainerVariable(const Variable<double>&);
template void ContainerVariableDataHolder<PropertiesDataValueContainer<ModelPart::ElementsContainerType>>::AssignDataToContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerVariableDataHolder<PropertiesDataValueContainer<ModelPart::ElementsContainerType>>::ReadDataFromContainerVariable(const Variable<double>&);
template void ContainerVariableDataHolder<PropertiesDataValueContainer<ModelPart::ElementsContainerType>>::ReadDataFromContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerVariableDataHolder<PropertiesDataValueContainer<ModelPart::ElementsContainerType>>::SetDataForContainerVariable(const Variable<double>&, const double&);
template void ContainerVariableDataHolder<PropertiesDataValueContainer<ModelPart::ElementsContainerType>>::SetDataForContainerVariable(const Variable<array_1d<double, 3>>&, const array_1d<double, 3>&);

} // namespace Kratos