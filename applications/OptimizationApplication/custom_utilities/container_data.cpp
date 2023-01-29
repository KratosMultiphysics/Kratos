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

// Project includes
#include "includes/data_communicator.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "custom_utilities/optimization_utils.h"

// Include base h
#include "container_data.h"

namespace Kratos {

namespace ContainerDataHelperUtilities
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

} // namespace ContainerDataHelperUtilities

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

ContainerDataBase::ContainerDataBase(
    ModelPart& rModelPart)
    : mrModelPart(rModelPart),
      mDataDimension(0)
{
}

ContainerDataBase::ContainerDataBase(const ContainerDataBase& rOther)
    : mrModelPart(rOther.mrModelPart)
{
    this->mData.resize(rOther.mData.size(), false);
    this->mDataDimension = rOther.mDataDimension;
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index){
        this->mData[Index] = rOther.mData[Index];
    });
}

double ContainerDataBase::NormInf() const
{
    return mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(IndexPartition<IndexType>(this->mData.size()).for_each<MaxReduction<double>>([&](const IndexType Index) {
        return this->mData[Index];
    }));
}

IndexType ContainerDataBase::GetDataDimension() const
{
    return mDataDimension;
}

Vector& ContainerDataBase::GetData()
{
    return mData;
}

const Vector& ContainerDataBase::GetData() const
{
    return mData;
}

void ContainerDataBase::CopyDataFrom(const ContainerDataBase& rOther)
{
    KRATOS_ERROR_IF(this->mrModelPart != rOther.mrModelPart)
        << "Modelpart mismatch between origin and destination for copy.\n"
        << "   Destination = " << *this << "\n"
        << "   Origin = " << rOther << "\n.";

    if (this->mData.size() != rOther.mData.size()) {
        this->mData.resize(rOther.mData.size(), false);
    }

    this->mDataDimension = rOther.mDataDimension;

    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index){
        this->mData[Index] = rOther.mData[Index];
    });
}

ModelPart& ContainerDataBase::GetModelPart()
{
    return mrModelPart;
}

const ModelPart& ContainerDataBase::GetModelPart() const
{
    return mrModelPart;
}

std::string ContainerDataBase::Info() const
{
    KRATOS_ERROR << "Calling ContainerDataBase::Info. This should "
                    "be implemented in the derrived class.";
    return "";
}

void ContainerDataBase::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

template<class TContainerDataType>
ContainerData<TContainerDataType> ContainerData<TContainerDataType>::Clone()
{
    return ContainerData<TContainerDataType>(*this);
}

template<class TContainerDataType>
template<class TDataType>
void ContainerData<TContainerDataType>::ReadDataFromContainerVariable(const Variable<TDataType>& rVariable)
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
        const auto& values = TContainerDataType::GetValue(*(r_container.begin() + Index), rVariable);
        for (IndexType i = 0; i < mDataDimension; ++i) {
            ContainerDataHelperUtilities::AssignValueToVector(r_data, Index * mDataDimension, i, values);
        }
    });

    KRATOS_CATCH("")
}

template<class TContainerDataType>
template<class TDataType>
void ContainerData<TContainerDataType>::AssignDataToContainerVariable(const Variable<TDataType>& rVariable)
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
    block_for_each(r_container, [&](auto& rEntity) {
        TContainerDataType::SetValue(rEntity, rVariable, rVariable.Zero());
    });

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index){
        auto& values = TContainerDataType::GetValue(*(r_container.begin() + Index), rVariable);
        for (IndexType i = 0; i < local_size; ++i) {
            ContainerDataHelperUtilities::AssignValueFromVector(values, i, Index * local_size, r_values);
        }
    });

    KRATOS_CATCH("");
}

template<class TContainerDataType>
template<class TDataType>
void ContainerData<TContainerDataType>::SetDataForContainerVariable(const Variable<TDataType>& rVariable, const TDataType& rValue)
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
            ContainerDataHelperUtilities::AssignValueToVector(this->mData, Index * mDataDimension, i, rValue);
        }
    });
}

template<class TContainerDataType>
double ContainerData<TContainerDataType>::InnerProduct(const ContainerData<TContainerDataType>& rOther) const
{
    KRATOS_ERROR_IF_NOT(this->mData.size() == rOther.mData.size())
        << "Data size mismatch in operands for +.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    // TODO: Make it MPI compatible with SumAll, but need to take care of ghost nodes in the case nodal historical/non historical
    //       containers are used
    return IndexPartition<IndexType>(this->mData.size()).for_each<SumReduction<double>>([&](const IndexType Index) {
        return this->mData[Index] * rOther.mData[Index];
    });
}

template<class TContainerDataType>
ContainerData<TContainerDataType> ContainerData<TContainerDataType>::operator+(const ContainerData<TContainerDataType>& rOther) const
{
    KRATOS_ERROR_IF(this->mrModelPart != rOther.mrModelPart)
        << "Mismatching model parts found in addition.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    KRATOS_ERROR_IF_NOT(this->mData.size() == rOther.mData.size())
        << "Data size mismatch in operands for +.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    ContainerData<TContainerDataType> result(this->mrModelPart);
    result.mDataDimension = this->mDataDimension;
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] + rOther.mData[Index];
    });


    return result;
}

template<class TContainerDataType>
ContainerData<TContainerDataType>& ContainerData<TContainerDataType>::operator+=(const ContainerData<TContainerDataType>& rOther)
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

template<class TContainerDataType>
ContainerData<TContainerDataType> ContainerData<TContainerDataType>::operator-(const ContainerData<TContainerDataType>& rOther) const
{
    KRATOS_ERROR_IF(this->mrModelPart != rOther.mrModelPart)
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    KRATOS_ERROR_IF_NOT(this->mData.size() == rOther.mData.size())
        << "Data size mismatch in operands for -.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    ContainerData<TContainerDataType> result(this->mrModelPart);
    result.mDataDimension = this->mDataDimension;
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] - rOther.mData[Index];
    });

    return result;
}

template<class TContainerDataType>
ContainerData<TContainerDataType>& ContainerData<TContainerDataType>::operator-=(const ContainerData<TContainerDataType>& rOther)
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

template<class TContainerDataType>
ContainerData<TContainerDataType> ContainerData<TContainerDataType>::operator*(const double Value) const
{
    ContainerData<TContainerDataType> result(this->mrModelPart);
    result.mDataDimension = this->mDataDimension;
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] * Value;
    });

    return result;
}

template<class TContainerDataType>
ContainerData<TContainerDataType>& ContainerData<TContainerDataType>::operator*=(const double Value)
{
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] *= Value;
    });

    return *this;
}

template<class TContainerDataType>
ContainerData<TContainerDataType> ContainerData<TContainerDataType>::operator/(const double Value) const
{
    KRATOS_ERROR_IF(std::abs(Value) < std::numeric_limits<double>::epsilon())
        << "Division by zero.\n"
        << "      Left operand data : " << *this << "\n"
        << "      divisor           : " << Value << "\n";

    ContainerData<TContainerDataType> result(this->mrModelPart);
    result.mDataDimension = this->mDataDimension;
    result.mData.resize(this->mData.size(), false);
    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        result.mData[Index] = this->mData[Index] / Value;
    });

    return result;
}

template<class TContainerDataType>
ContainerData<TContainerDataType>& ContainerData<TContainerDataType>::operator/=(const double Value)
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

template<class TContainerDataType>
ContainerData<TContainerDataType>& ContainerData<TContainerDataType>::operator=(const ContainerData<TContainerDataType>& rOther)
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

template<class TContainerDataType>
typename TContainerDataType::ContainerType& ContainerData<TContainerDataType>::GetContainer()
{
    using container_type = typename TContainerDataType::ContainerType;

    if constexpr(std::is_same_v<container_type, ModelPart::NodesContainerType>) {
        return this->GetModelPart().Nodes();
    } else if constexpr(std::is_same_v<container_type, ModelPart::ConditionsContainerType>) {
        return this->GetModelPart().Conditions();
    } else if constexpr(std::is_same_v<container_type, ModelPart::ElementsContainerType>) {
        return this->GetModelPart().Elements();
    }
}

template<class TContainerDataType>
const typename TContainerDataType::ContainerType& ContainerData<TContainerDataType>::GetContainer() const
{
    using container_type = typename TContainerDataType::ContainerType;

    if constexpr(std::is_same_v<container_type, ModelPart::NodesContainerType>) {
        return this->GetModelPart().Nodes();
    } else if constexpr(std::is_same_v<container_type, ModelPart::ConditionsContainerType>) {
        return this->GetModelPart().Conditions();
    } else if constexpr(std::is_same_v<container_type, ModelPart::ElementsContainerType>) {
        return this->GetModelPart().Elements();
    }
}

template<class TContainerDataType>
std::string ContainerData<TContainerDataType>::Info() const
{
    std::stringstream msg;

    if constexpr(std::is_same_v<TContainerDataType, HistoricalDataValueContainer>) {
        msg << "NodalHistoricalDataContainer";
    } else if constexpr(std::is_same_v<TContainerDataType, NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>) {
        msg << "NodalNonHistoricalDataContainer";
    } else if constexpr(std::is_same_v<TContainerDataType, NonHistoricalDataValueContainer<ModelPart::ConditionsContainerType>>) {
        msg << "ConditionNonHistoricalDataContainer";
    } else if constexpr(std::is_same_v<TContainerDataType, NonHistoricalDataValueContainer<ModelPart::ElementsContainerType>>) {
        msg << "ElementNonHistoricalDataContainer";
    } else if constexpr(std::is_same_v<TContainerDataType, PropertiesDataValueContainer<ModelPart::ConditionsContainerType>>) {
        msg << "ConditionPropertiesDataContainer";
    } else if constexpr(std::is_same_v<TContainerDataType, PropertiesDataValueContainer<ModelPart::ElementsContainerType>>) {
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
template class ContainerData<HistoricalDataValueContainer>;
template void ContainerData<HistoricalDataValueContainer>::AssignDataToContainerVariable(const Variable<double>&);
template void ContainerData<HistoricalDataValueContainer>::AssignDataToContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerData<HistoricalDataValueContainer>::ReadDataFromContainerVariable(const Variable<double>&);
template void ContainerData<HistoricalDataValueContainer>::ReadDataFromContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerData<HistoricalDataValueContainer>::SetDataForContainerVariable(const Variable<double>&, const double&);
template void ContainerData<HistoricalDataValueContainer>::SetDataForContainerVariable(const Variable<array_1d<double, 3>>&, const array_1d<double, 3>&);

template class ContainerData<NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>;
template void ContainerData<NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>::AssignDataToContainerVariable(const Variable<double>&);
template void ContainerData<NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>::AssignDataToContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerData<NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>::ReadDataFromContainerVariable(const Variable<double>&);
template void ContainerData<NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>::ReadDataFromContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerData<NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>::SetDataForContainerVariable(const Variable<double>&, const double&);
template void ContainerData<NonHistoricalDataValueContainer<ModelPart::NodesContainerType>>::SetDataForContainerVariable(const Variable<array_1d<double, 3>>&, const array_1d<double, 3>&);

template class ContainerData<NonHistoricalDataValueContainer<ModelPart::ConditionsContainerType>>;
template void ContainerData<NonHistoricalDataValueContainer<ModelPart::ConditionsContainerType>>::AssignDataToContainerVariable(const Variable<double>&);
template void ContainerData<NonHistoricalDataValueContainer<ModelPart::ConditionsContainerType>>::AssignDataToContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerData<NonHistoricalDataValueContainer<ModelPart::ConditionsContainerType>>::ReadDataFromContainerVariable(const Variable<double>&);
template void ContainerData<NonHistoricalDataValueContainer<ModelPart::ConditionsContainerType>>::ReadDataFromContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerData<NonHistoricalDataValueContainer<ModelPart::ConditionsContainerType>>::SetDataForContainerVariable(const Variable<double>&, const double&);
template void ContainerData<NonHistoricalDataValueContainer<ModelPart::ConditionsContainerType>>::SetDataForContainerVariable(const Variable<array_1d<double, 3>>&, const array_1d<double, 3>&);

template class ContainerData<NonHistoricalDataValueContainer<ModelPart::ElementsContainerType>>;
template void ContainerData<NonHistoricalDataValueContainer<ModelPart::ElementsContainerType>>::AssignDataToContainerVariable(const Variable<double>&);
template void ContainerData<NonHistoricalDataValueContainer<ModelPart::ElementsContainerType>>::AssignDataToContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerData<NonHistoricalDataValueContainer<ModelPart::ElementsContainerType>>::ReadDataFromContainerVariable(const Variable<double>&);
template void ContainerData<NonHistoricalDataValueContainer<ModelPart::ElementsContainerType>>::ReadDataFromContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerData<NonHistoricalDataValueContainer<ModelPart::ElementsContainerType>>::SetDataForContainerVariable(const Variable<double>&, const double&);
template void ContainerData<NonHistoricalDataValueContainer<ModelPart::ElementsContainerType>>::SetDataForContainerVariable(const Variable<array_1d<double, 3>>&, const array_1d<double, 3>&);

template class ContainerData<PropertiesDataValueContainer<ModelPart::ConditionsContainerType>>;
template void ContainerData<PropertiesDataValueContainer<ModelPart::ConditionsContainerType>>::AssignDataToContainerVariable(const Variable<double>&);
template void ContainerData<PropertiesDataValueContainer<ModelPart::ConditionsContainerType>>::AssignDataToContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerData<PropertiesDataValueContainer<ModelPart::ConditionsContainerType>>::ReadDataFromContainerVariable(const Variable<double>&);
template void ContainerData<PropertiesDataValueContainer<ModelPart::ConditionsContainerType>>::ReadDataFromContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerData<PropertiesDataValueContainer<ModelPart::ConditionsContainerType>>::SetDataForContainerVariable(const Variable<double>&, const double&);
template void ContainerData<PropertiesDataValueContainer<ModelPart::ConditionsContainerType>>::SetDataForContainerVariable(const Variable<array_1d<double, 3>>&, const array_1d<double, 3>&);

template class ContainerData<PropertiesDataValueContainer<ModelPart::ElementsContainerType>>;
template void ContainerData<PropertiesDataValueContainer<ModelPart::ElementsContainerType>>::AssignDataToContainerVariable(const Variable<double>&);
template void ContainerData<PropertiesDataValueContainer<ModelPart::ElementsContainerType>>::AssignDataToContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerData<PropertiesDataValueContainer<ModelPart::ElementsContainerType>>::ReadDataFromContainerVariable(const Variable<double>&);
template void ContainerData<PropertiesDataValueContainer<ModelPart::ElementsContainerType>>::ReadDataFromContainerVariable(const Variable<array_1d<double, 3>>&);
template void ContainerData<PropertiesDataValueContainer<ModelPart::ElementsContainerType>>::SetDataForContainerVariable(const Variable<double>&, const double&);
template void ContainerData<PropertiesDataValueContainer<ModelPart::ElementsContainerType>>::SetDataForContainerVariable(const Variable<array_1d<double, 3>>&, const array_1d<double, 3>&);

} // namespace Kratos