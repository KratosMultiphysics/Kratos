//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_utilities/optimization_utils.h"

// Include base h
#include "container_variable_data_holder.h"

namespace Kratos {

namespace ContainerVariableDataHolderHelperUtilities
{

using Array3D = array_1d<double, 3>;

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
    const Array3D& rValue)
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
    Array3D& rOutput,
    const IndexType ValueComponentIndex,
    const IndexType VectoStartingIndex,
    const Vector& rInput)
{
    rOutput[ValueComponentIndex] = rInput[VectoStartingIndex + ValueComponentIndex];
}

} // namespace ContainerVariableDataHolderHelperUtilities

template <class TContainerType, class TContainerDataIO>
typename ContainerVariableDataHolder<TContainerType, TContainerDataIO>::Pointer ContainerVariableDataHolder<TContainerType, TContainerDataIO>::Clone() const
{
    return std::make_shared<ContainerVariableDataHolder<TContainerType, TContainerDataIO>>(*this);
}

template <class TContainerType, class TContainerDataIO>
typename ContainerVariableDataHolder<TContainerType, TContainerDataIO>::Pointer ContainerVariableDataHolder<TContainerType, TContainerDataIO>::CloneWithDataInitializedToZero() const
{
    std::shared_ptr<ContainerVariableDataHolder<TContainerType, TContainerDataIO>> result(new ContainerVariableDataHolder<TContainerType, TContainerDataIO>(*(this->mpModelPart)));
    result->mDataDimension = this->mDataDimension;

    auto& r_data = result->mData;
    const IndexType number_of_entities = this->GetContainer().size();

    r_data.resize(number_of_entities * this->mDataDimension);

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index){
        r_data[Index] = 0.0;
    });

    return result;
}

template <class TContainerType, class TContainerDataIO>
template<class TDataType>
void ContainerVariableDataHolder<TContainerType, TContainerDataIO>::ReadDataFromContainerVariable(const Variable<TDataType>& rVariable)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(this->GetModelPart().GetProcessInfo().Has(DOMAIN_SIZE))
        << "DOMAIN_SIZE variable is not found in process info of model part "
           "with name "
        << this->GetModelPart().FullName() << ".\n";

    const IndexType number_of_entities = this->GetContainer().size();
    this->mDataDimension = OptimizationUtils::GetLocalSize<TDataType>(this->GetModelPart().GetProcessInfo()[DOMAIN_SIZE]);

    auto& r_data = this->mData;
    auto& r_container = this->GetContainer();

    r_data.resize(number_of_entities * this->mDataDimension);

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index){
        const auto& values = TContainerDataIO::GetValue(*(r_container.begin() + Index), rVariable);
        for (IndexType i = 0; i < this->mDataDimension; ++i) {
            ContainerVariableDataHolderHelperUtilities::AssignValueToVector(r_data, Index * this->mDataDimension, i, values);
        }
    });

    KRATOS_CATCH("")
}

template <class TContainerType, class TContainerDataIO>
template<class TDataType>
void ContainerVariableDataHolder<TContainerType, TContainerDataIO>::AssignDataToContainerVariable(const Variable<TDataType>& rVariable)
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
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        // initializes ghost nodes as well for the later synchronization
        // only, the nodal non historical values needs to be set unless
        // they are properly initialized.
        block_for_each(this->GetModelPart().Nodes(), [&](auto& rEntity) {
            TContainerDataIO::SetValue(rEntity, rVariable, rVariable.Zero());
        });
    } else {
        // no ghost settings required for elements or conditions nor historical values, hence no special treatment required.
        block_for_each(r_container, [&](auto& rEntity) {
            TContainerDataIO::SetValue(rEntity, rVariable, rVariable.Zero());
        });
    }

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index){
        auto& values = TContainerDataIO::GetValue(*(r_container.begin() + Index), rVariable);
        for (IndexType i = 0; i < local_size; ++i) {
            ContainerVariableDataHolderHelperUtilities::AssignValueFromVector(values, i, Index * local_size, r_values);
        }
    });

    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        // synchronize nodal values
        auto& r_communicator = this->GetModelPart().GetCommunicator();

        if constexpr(std::is_same_v<TContainerDataIO, HistoricalContainerDataIO>) {
            r_communicator.SynchronizeNodalSolutionStepsData();
        } else if constexpr(std::is_same_v<TContainerDataIO, NonHistoricalContainerDataIO>) {
            r_communicator.SynchronizeNonHistoricalVariable(rVariable);
        }
    }

    KRATOS_CATCH("");
}

template <class TContainerType, class TContainerDataIO>
template<class TDataType>
void ContainerVariableDataHolder<TContainerType, TContainerDataIO>::SetDataForContainerVariable(const Variable<TDataType>& rVariable, const TDataType& rValue)
{

    KRATOS_ERROR_IF_NOT(this->GetModelPart().GetProcessInfo().Has(DOMAIN_SIZE))
        << "DOMAIN_SIZE variable is not found in process info of model part "
           "with name "
        << this->GetModelPart().FullName() << ".\n";

    this->mDataDimension = OptimizationUtils::GetLocalSize<TDataType>(this->GetModelPart().GetProcessInfo()[DOMAIN_SIZE]);

    const IndexType number_of_entities = this->GetContainer().size();

    if (this->mData.size() != this->mDataDimension * number_of_entities) {
        this->mData.resize(this->mDataDimension * number_of_entities, false);
    }

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index) {
        for (IndexType i = 0; i < this->mDataDimension; ++i) {
            ContainerVariableDataHolderHelperUtilities::AssignValueToVector(this->mData, Index * this->mDataDimension, i, rValue);
        }
    });
}

template <class TContainerType, class TContainerDataIO>
template<class TDataType>
void ContainerVariableDataHolder<TContainerType, TContainerDataIO>::SetDataForContainerVariableToZero(const Variable<TDataType>& rVariable)
{
    this->SetDataForContainerVariable(rVariable, rVariable.Zero());
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO> ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator+(const ContainerVariableDataHolder<TContainerType, TContainerDataIO>& rOther) const
{
    KRATOS_ERROR_IF(this->GetModelPart() != rOther.GetModelPart())
        << "Mismatching model parts found in addition.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    KRATOS_ERROR_IF_NOT(this->mData.size() == rOther.mData.size())
        << "Data size mismatch in operands for +.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    ContainerVariableDataHolder<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mDataDimension = this->mDataDimension;
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] + rOther.mData[Index];
    });

    return result;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator+=(const ContainerVariableDataHolder<TContainerType, TContainerDataIO>& rOther)
{
    KRATOS_ERROR_IF(this->GetModelPart() != rOther.GetModelPart())
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

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO> ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator+(const double Value) const
{
    ContainerVariableDataHolder<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mDataDimension = this->mDataDimension;
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] + Value;
    });

    return result;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator+=(const double Value)
{
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] += Value;
    });

    return *this;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO> ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator-(const ContainerVariableDataHolder<TContainerType, TContainerDataIO>& rOther) const
{
    KRATOS_ERROR_IF(this->GetModelPart() != rOther.GetModelPart())
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    KRATOS_ERROR_IF_NOT(this->mData.size() == rOther.mData.size())
        << "Data size mismatch in operands for -.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    ContainerVariableDataHolder<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mDataDimension = this->mDataDimension;
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] - rOther.mData[Index];
    });

    return result;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator-=(const ContainerVariableDataHolder<TContainerType, TContainerDataIO>& rOther)
{
    KRATOS_ERROR_IF(this->GetModelPart() != rOther.GetModelPart())
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

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO> ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator-(const double Value) const
{
    ContainerVariableDataHolder<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mDataDimension = this->mDataDimension;
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] - Value;
    });

    return result;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator-=(const double Value)
{
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] -= Value;
    });

    return *this;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO> ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator*(const double Value) const
{
    ContainerVariableDataHolder<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mDataDimension = this->mDataDimension;
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] * Value;
    });

    return result;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator*=(const double Value)
{
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] *= Value;
    });

    return *this;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO> ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator/(const double Value) const
{
    KRATOS_ERROR_IF(std::abs(Value) < std::numeric_limits<double>::epsilon())
        << "Division by zero.\n"
        << "      Left operand data : " << *this << "\n"
        << "      divisor           : " << Value << "\n";

    ContainerVariableDataHolder<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mDataDimension = this->mDataDimension;
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] / Value;
    });

    return result;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator/=(const double Value)
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

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO> ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator^(const double Value) const
{
    ContainerVariableDataHolder<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mDataDimension = this->mDataDimension;
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = std::pow(this->mData[Index], Value);
    });

    return result;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator^=(const double Value)
{
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] = std::pow(this->mData[Index], Value);
    });

    return *this;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator=(const ContainerVariableDataHolder<TContainerType, TContainerDataIO>& rOther)
{
    KRATOS_ERROR_IF(this->GetModelPart() != rOther.GetModelPart())
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

template <class TContainerType, class TContainerDataIO>
std::string ContainerVariableDataHolder<TContainerType, TContainerDataIO>::Info() const
{
    std::stringstream msg;

    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        msg << "Nodal";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        msg << "Condition";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        msg << "Element";
    }

    if constexpr(std::is_same_v<TContainerDataIO, HistoricalContainerDataIO>) {
        msg << "Historical";
    } else if constexpr(std::is_same_v<TContainerDataIO, NonHistoricalContainerDataIO>) {
        msg << "NonHistorical";
    } else if constexpr(std::is_same_v<TContainerDataIO, PropertiesContainerDataIO>) {
        msg << "Properties";
    }

    msg << BaseType::Info();

    return msg.str();
}

//template instantiations
#define INSTANTIATIE_CONTAINER_DATA_METHODS(ContainerType, ContainerDataIOType, DataType)                                                                   \
    template void ContainerVariableDataHolder<ContainerType, ContainerDataIOType>::ReadDataFromContainerVariable(const Variable<DataType>&);                \
    template void ContainerVariableDataHolder<ContainerType, ContainerDataIOType>::SetDataForContainerVariable(const Variable<DataType>&, const DataType&); \
    template void ContainerVariableDataHolder<ContainerType, ContainerDataIOType>::SetDataForContainerVariableToZero(const Variable<DataType>&);            \
    template void ContainerVariableDataHolder<ContainerType, ContainerDataIOType>::AssignDataToContainerVariable(const Variable<DataType>&);

#define INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ContainerType, ContainerDataIOType)                                          \
    template class ContainerVariableDataHolder<ContainerType, ContainerDataIOType>;                                             \
    INSTANTIATIE_CONTAINER_DATA_METHODS(ContainerType, ContainerDataIOType, double)                                             \
    INSTANTIATIE_CONTAINER_DATA_METHODS(ContainerType, ContainerDataIOType, ContainerVariableDataHolderHelperUtilities::Array3D)

INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ModelPart::NodesContainerType, HistoricalContainerDataIO)
INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ModelPart::NodesContainerType, NonHistoricalContainerDataIO)
INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ModelPart::ConditionsContainerType, NonHistoricalContainerDataIO)
INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ModelPart::ElementsContainerType, NonHistoricalContainerDataIO)
INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ModelPart::ConditionsContainerType, PropertiesContainerDataIO)
INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ModelPart::ElementsContainerType, PropertiesContainerDataIO)

#undef INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER
#undef INSTANTIATIE_CONTAINER_DATA_METHODS

}