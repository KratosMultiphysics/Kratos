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
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "container_variable_data_holder.h"

namespace Kratos {

namespace ContainerVariableDataHolderHelperUtilities
{

using Array3D = array_1d<double, 3>;

template<class TDataType>
IndexType GetLocalSize(const IndexType DomainSize);

template<>
IndexType GetLocalSize<double>(const IndexType DomainSize)
{
    return 1;
}

template<>
IndexType GetLocalSize<array_1d<double, 3>>(const IndexType DomainSize)
{
    return DomainSize;
}

void AssignValueToVector(
    Vector& rOutput,
    const IndexType StartIndex,
    const IndexType ComponentIndex,
    double rValue)
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
    const double rInput)
{
    rOutput = rInput;
}

void AssignValueFromVector(
    Array3D& rOutput,
    const IndexType ValueComponentIndex,
    const double rInput)
{
    rOutput[ValueComponentIndex] = rInput;
}

} // namespace ContainerVariableDataHolderHelperUtilities

template <class TContainerType, class TContainerDataIO>
typename ContainerVariableDataHolder<TContainerType, TContainerDataIO>::Pointer ContainerVariableDataHolder<TContainerType, TContainerDataIO>::Clone() const
{
    return Kratos::make_shared<ContainerVariableDataHolder<TContainerType, TContainerDataIO>>(*this);
}

template <class TContainerType, class TContainerDataIO>
typename ContainerVariableDataHolder<TContainerType, TContainerDataIO>::Pointer ContainerVariableDataHolder<TContainerType, TContainerDataIO>::CloneWithDataInitializedToZero() const
{
    std::shared_ptr<ContainerVariableDataHolder<TContainerType, TContainerDataIO>> result(new ContainerVariableDataHolder<TContainerType, TContainerDataIO>(*(this->mpModelPart)));
    result->SetDataToZero();
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
    const auto dimension = ContainerVariableDataHolderHelperUtilities::GetLocalSize<TDataType>(this->GetModelPart().GetProcessInfo()[DOMAIN_SIZE]);

    auto p_data = Kratos::make_shared<Vector>();
    this->mpExpression = LiteralVectorExpression::Create(p_data, dimension);

    auto& r_data = *p_data;
    auto& r_container = this->GetContainer();

    r_data.resize(number_of_entities * dimension);

    IndexPartition<IndexType>(number_of_entities).for_each([&](const IndexType Index){
        const auto& values = TContainerDataIO::GetValue(*(r_container.begin() + Index), rVariable);
        for (IndexType i = 0; i < dimension; ++i) {
            ContainerVariableDataHolderHelperUtilities::AssignValueToVector(r_data, Index * dimension, i, values);
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

    const IndexType local_size = ContainerVariableDataHolderHelperUtilities::GetLocalSize<TDataType>(this->GetModelPart().GetProcessInfo()[DOMAIN_SIZE]);

    KRATOS_ERROR_IF_NOT(local_size == this->GetDataDimension())
        << "Stored data dimension and requested assignment variable data "
           "dimension mismatch. [ Sotred data dimension = "
        << this->GetDataDimension()
        << ", assignment requested variable data diemension = " << local_size << " ].\n";

    auto& r_container = this->GetContainer();
    const IndexType number_of_entities = r_container.size();

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
            const double evaluated_value = this->mpExpression->Evaluate(Index, i);
            ContainerVariableDataHolderHelperUtilities::AssignValueFromVector(values, i, evaluated_value);
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
    static_assert(
        std::is_same_v<TDataType, double>               ||
        std::is_same_v<TDataType, int>                  ||
        std::is_same_v<TDataType, std::size_t>          ||
        std::is_same_v<TDataType, array_1d<double, 3>>  ||
            "Only allowed types are double, int, size_t, array_1d<double, 3>");

    KRATOS_ERROR_IF_NOT(this->GetModelPart().GetProcessInfo().Has(DOMAIN_SIZE))
        << "DOMAIN_SIZE variable is not found in process info of model part "
           "with name "
        << this->GetModelPart().FullName() << ".\n";

    const auto dimension = ContainerVariableDataHolderHelperUtilities::GetLocalSize<TDataType>(this->GetModelPart().GetProcessInfo()[DOMAIN_SIZE]);
    if constexpr(std::is_same_v<TDataType, double> || std::is_same_v<TDataType, int> || std::is_same_v<TDataType, std::size_t>) {
        this->mpExpression = LiteralDoubleExpression::Create(rValue);
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 3>>) {
        this->mpExpression = LiteralArray3Expression::Create(rValue, dimension);
    }

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
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in addition.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    ContainerVariableDataHolder<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinaryAddExpression::Create(this->mpExpression, rOther.mpExpression);
    return result;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator+=(const ContainerVariableDataHolder<TContainerType, TContainerDataIO>& rOther)
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in addition.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    this->mpExpression = BinaryAddExpression::Create(this->mpExpression, rOther.mpExpression);
    return *this;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO> ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator+(const double Value) const
{

    ContainerVariableDataHolder<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinaryAddExpression::Create(this->mpExpression, LiteralDoubleExpression::Create(Value));
    return result;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator+=(const double Value)
{
    this->mpExpression = BinaryAddExpression::Create(this->mpExpression, LiteralDoubleExpression::Create(Value));
    return *this;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO> ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator-(const ContainerVariableDataHolder<TContainerType, TContainerDataIO>& rOther) const
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    ContainerVariableDataHolder<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinarySubstractExpression::Create(this->mpExpression, rOther.mpExpression);
    return result;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator-=(const ContainerVariableDataHolder<TContainerType, TContainerDataIO>& rOther)
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    this->mpExpression = BinarySubstractExpression::Create(this->mpExpression, rOther.mpExpression);
    return *this;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO> ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator-(const double Value) const
{
    ContainerVariableDataHolder<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinarySubstractExpression::Create(this->mpExpression, LiteralDoubleExpression::Create(Value));
    return result;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator-=(const double Value)
{
    this->mpExpression = BinarySubstractExpression::Create(this->mpExpression, LiteralDoubleExpression::Create(Value));
    return *this;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO> ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator*(const ContainerVariableDataHolder<TContainerType, TContainerDataIO>& rOther) const
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    ContainerVariableDataHolder<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinaryMultiplyExpression::Create(this->mpExpression, rOther.mpExpression);
    return result;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator*=(const ContainerVariableDataHolder<TContainerType, TContainerDataIO>& rOther)
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    this->mpExpression = BinaryMultiplyExpression::Create(this->mpExpression, rOther.mpExpression);
    return *this;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO> ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator*(const double Value) const
{
    ContainerVariableDataHolder<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinaryMultiplyExpression::Create(this->mpExpression, LiteralDoubleExpression::Create(Value));
    return result;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator*=(const double Value)
{
    this->mpExpression = BinaryMultiplyExpression::Create(this->mpExpression, LiteralDoubleExpression::Create(Value));
    return *this;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO> ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator/(const ContainerVariableDataHolder<TContainerType, TContainerDataIO>& rOther) const
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    ContainerVariableDataHolder<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinaryDivideExpression::Create(this->mpExpression, rOther.mpExpression);
    return result;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator/=(const ContainerVariableDataHolder<TContainerType, TContainerDataIO>& rOther)
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    this->mpExpression = BinaryDivideExpression::Create(this->mpExpression, rOther.mpExpression);
    return *this;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO> ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator/(const double Value) const
{
    KRATOS_ERROR_IF(std::abs(Value) < std::numeric_limits<double>::epsilon())
        << "Division by zero.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Divisor           : " << Value << "\n";

    ContainerVariableDataHolder<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinaryMultiplyExpression::Create(this->mpExpression, LiteralDoubleExpression::Create(1.0 / Value));
    return result;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator/=(const double Value)
{
    KRATOS_ERROR_IF(std::abs(Value) < std::numeric_limits<double>::epsilon())
        << "Division by zero.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Divisor           : " << Value << "\n";

    this->mpExpression = BinaryMultiplyExpression::Create(this->mpExpression, LiteralDoubleExpression::Create(1.0 / Value));
    return *this;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO> ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator^(const ContainerVariableDataHolder<TContainerType, TContainerDataIO>& rOther) const
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    ContainerVariableDataHolder<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinaryPowerExpression::Create(this->mpExpression, rOther.mpExpression);
    return result;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator^=(const ContainerVariableDataHolder<TContainerType, TContainerDataIO>& rOther)
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    this->mpExpression = BinaryPowerExpression::Create(this->mpExpression, rOther.mpExpression);
    return *this;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO> ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator^(const double Value) const
{
    ContainerVariableDataHolder<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinaryPowerExpression::Create(this->mpExpression, LiteralDoubleExpression::Create(Value));
    return result;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator^=(const double Value)
{
    this->mpExpression = BinaryPowerExpression::Create(this->mpExpression, LiteralDoubleExpression::Create(Value));
    return *this;
}

template <class TContainerType, class TContainerDataIO>
ContainerVariableDataHolder<TContainerType, TContainerDataIO>& ContainerVariableDataHolder<TContainerType, TContainerDataIO>::operator=(const ContainerVariableDataHolder<TContainerType, TContainerDataIO>& rOther)
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in assignment.\n"
        << "      Assignee data: " << *this << "\n"
        << "      Assignor data: " << rOther << "\n";

    this->mpExpression = rOther.mpExpression;
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