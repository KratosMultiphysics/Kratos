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
#include <tuple>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/container_variable_data/container_data_io.h"
#include "custom_utilities/container_variable_data/variable_expression_data_io.h"
#include "custom_utilities/container_variable_data/expressions/literal/literal_expression.h"
#include "custom_utilities/container_variable_data/expressions/literal/literal_flat_expression.h"
#include "custom_utilities/container_variable_data/expressions/binary/binary_expression.h"

// Include base h
#include "specialized_container_variable_data.h"

namespace Kratos {

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO>& SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator=(const SpecializedContainerVariableData& rOther)
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in assignment.\n"
        << "      Assignee data: " << *this << "\n"
        << "      Assignor data: " << rOther << "\n";

    this->mpExpression = rOther.mpExpression;
    return *this;
}

template <class TContainerType, class TContainerDataIO>
typename SpecializedContainerVariableData<TContainerType, TContainerDataIO>::Pointer SpecializedContainerVariableData<TContainerType, TContainerDataIO>::Clone() const
{
    return Kratos::make_shared<SpecializedContainerVariableData<TContainerType, TContainerDataIO>>(*this);
}

template <class TContainerType, class TContainerDataIO>
template<class TDataType>
void SpecializedContainerVariableData<TContainerType, TContainerDataIO>::ReadData(const Variable<TDataType>& rVariable)
{
    KRATOS_TRY

    const auto& r_container = this->GetContainer();
    const IndexType number_of_entities = r_container.size();

    if (number_of_entities != 0) {
        // initialize the shape with the first entity value
        VariableExpressionDataIO<TDataType> variable_flatten_data_io(TContainerDataIO::GetValue(*r_container.begin(), rVariable));

        auto p_expression = LiteralFlatExpression::Create(number_of_entities, variable_flatten_data_io.GetShape());
        auto& r_expression = *p_expression;
        this->mpExpression = p_expression;

        IndexPartition<IndexType>(number_of_entities).for_each([&r_container, &rVariable, &variable_flatten_data_io, &r_expression](const IndexType Index){
            const auto& values = TContainerDataIO::GetValue(*(r_container.begin() + Index), rVariable);
            variable_flatten_data_io.Read(r_expression, Index, values);
        });
    }

    KRATOS_CATCH("")
}

template <class TContainerType, class TContainerDataIO>
template<class TDataType>
void SpecializedContainerVariableData<TContainerType, TContainerDataIO>::AssignData(const Variable<TDataType>& rVariable)
{
    KRATOS_TRY

    auto& r_container = this->GetContainer();
    const IndexType number_of_entities = r_container.size();

    // initialize the container variables first
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType> && std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::NonHistorical>>) {
        // initializes ghost nodes as for the later synchronization
        // only, the nodal non historical values needs to be set unless
        // they are properly initialized. Otherwise, in synchronization, the variables will
        // not be there in the ghost nodes hence seg faults.
        VariableUtils().SetNonHistoricalVariablesToZero(this->GetModelPart().GetCommunicator().GhostMesh().Nodes(), rVariable);
    }

    const auto& r_expression = *this->mpExpression;

    VariableExpressionDataIO<TDataType> variable_flatten_data_io(r_expression.GetShape());

    IndexPartition<IndexType>(number_of_entities).for_each(TDataType{}, [&r_container, &rVariable, &r_expression, &variable_flatten_data_io](const IndexType Index, TDataType& rValue){
        variable_flatten_data_io.Assign(rValue, r_expression, Index);
        TContainerDataIO::SetValue(*(r_container.begin() + Index), rVariable, rValue);
    });

    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        // synchronize nodal values
        auto& r_communicator = this->GetModelPart().GetCommunicator();

        if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::Historical>>) {
            r_communicator.SynchronizeVariable(rVariable);
        } else if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::NonHistorical>>) {
            r_communicator.SynchronizeNonHistoricalVariable(rVariable);
        }
    }

    KRATOS_CATCH("");
}

template <class TContainerType, class TContainerDataIO>
template<class TDataType>
void SpecializedContainerVariableData<TContainerType, TContainerDataIO>::SetData(const TDataType& rValue)
{
    this->mpExpression = LiteralExpression<TDataType>::Create(rValue);
}

template <class TContainerType, class TContainerDataIO>
template<class TDataType>
void SpecializedContainerVariableData<TContainerType, TContainerDataIO>::SetZero(const Variable<TDataType>& rVariable)
{
    this->SetData(rVariable.Zero());
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO> SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator+(const SpecializedContainerVariableData<TContainerType, TContainerDataIO>& rOther) const
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in addition.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    SpecializedContainerVariableData<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinaryExpression<BinaryOperations::Addition>::Create(this->mpExpression, rOther.mpExpression);
    return result;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO>& SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator+=(const SpecializedContainerVariableData<TContainerType, TContainerDataIO>& rOther)
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in addition.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    this->mpExpression = BinaryExpression<BinaryOperations::Addition>::Create(this->mpExpression, rOther.mpExpression);
    return *this;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO> SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator+(const double Value) const
{

    SpecializedContainerVariableData<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinaryExpression<BinaryOperations::Addition>::Create(this->mpExpression, LiteralExpression<double>::Create(Value));
    return result;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO>& SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator+=(const double Value)
{
    this->mpExpression = BinaryExpression<BinaryOperations::Addition>::Create(this->mpExpression, LiteralExpression<double>::Create(Value));
    return *this;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO> SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator-(const SpecializedContainerVariableData<TContainerType, TContainerDataIO>& rOther) const
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    SpecializedContainerVariableData<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinaryExpression<BinaryOperations::Substraction>::Create(this->mpExpression, rOther.mpExpression);
    return result;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO>& SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator-=(const SpecializedContainerVariableData<TContainerType, TContainerDataIO>& rOther)
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    this->mpExpression = BinaryExpression<BinaryOperations::Substraction>::Create(this->mpExpression, rOther.mpExpression);
    return *this;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO> SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator-(const double Value) const
{
    SpecializedContainerVariableData<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinaryExpression<BinaryOperations::Substraction>::Create(this->mpExpression, LiteralExpression<double>::Create(Value));
    return result;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO>& SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator-=(const double Value)
{
    this->mpExpression = BinaryExpression<BinaryOperations::Substraction>::Create(this->mpExpression, LiteralExpression<double>::Create(Value));
    return *this;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO> SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator*(const SpecializedContainerVariableData<TContainerType, TContainerDataIO>& rOther) const
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    SpecializedContainerVariableData<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinaryExpression<BinaryOperations::Multiplication>::Create(this->mpExpression, rOther.mpExpression);
    return result;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO>& SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator*=(const SpecializedContainerVariableData<TContainerType, TContainerDataIO>& rOther)
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    this->mpExpression = BinaryExpression<BinaryOperations::Multiplication>::Create(this->mpExpression, rOther.mpExpression);
    return *this;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO> SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator*(const double Value) const
{
    SpecializedContainerVariableData<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinaryExpression<BinaryOperations::Multiplication>::Create(this->mpExpression, LiteralExpression<double>::Create(Value));
    return result;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO>& SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator*=(const double Value)
{
    this->mpExpression = BinaryExpression<BinaryOperations::Multiplication>::Create(this->mpExpression, LiteralExpression<double>::Create(Value));
    return *this;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO> SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator/(const SpecializedContainerVariableData<TContainerType, TContainerDataIO>& rOther) const
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    SpecializedContainerVariableData<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinaryExpression<BinaryOperations::Division>::Create(this->mpExpression, rOther.mpExpression);
    return result;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO>& SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator/=(const SpecializedContainerVariableData<TContainerType, TContainerDataIO>& rOther)
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    this->mpExpression = BinaryExpression<BinaryOperations::Division>::Create(this->mpExpression, rOther.mpExpression);
    return *this;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO> SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator/(const double Value) const
{
    KRATOS_ERROR_IF(std::abs(Value) < std::numeric_limits<double>::epsilon())
        << "Division by zero.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Divisor           : " << Value << "\n";

    SpecializedContainerVariableData<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinaryExpression<BinaryOperations::Multiplication>::Create(this->mpExpression, LiteralExpression<double>::Create(1.0 / Value));
    return result;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO>& SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator/=(const double Value)
{
    KRATOS_ERROR_IF(std::abs(Value) < std::numeric_limits<double>::epsilon())
        << "Division by zero.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Divisor           : " << Value << "\n";

    this->mpExpression = BinaryExpression<BinaryOperations::Multiplication>::Create(this->mpExpression, LiteralExpression<double>::Create(1.0 / Value));
    return *this;
}

template <class TContainerType, class TContainerDataIO>
std::string SpecializedContainerVariableData<TContainerType, TContainerDataIO>::Info() const
{
    std::stringstream msg;

    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        msg << "Nodal";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        msg << "Condition";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        msg << "Element";
    }

    if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::Historical>>) {
        msg << "Historical";
    } else if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::NonHistorical>>) {
        msg << "NonHistorical";
    } else if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::Properties>>) {
        msg << "Properties";
    }

    msg << BaseType::Info();

    return msg.str();
}

//template instantiations
#define KRATOS_INSTANTIATE_CONTAINER_DATA_METHODS(ContainerType, ContainerDataIOType, ...)                                                                   \
    template void SpecializedContainerVariableData<ContainerType, ContainerDataIOType>::ReadData(const Variable<__VA_ARGS__>&);                \
    template void SpecializedContainerVariableData<ContainerType, ContainerDataIOType>::SetData(const __VA_ARGS__&); \
    template void SpecializedContainerVariableData<ContainerType, ContainerDataIOType>::SetZero(const Variable<__VA_ARGS__>&);            \
    template void SpecializedContainerVariableData<ContainerType, ContainerDataIOType>::AssignData(const Variable<__VA_ARGS__>&);

#define KRATOS_INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ContainerType, ContainerDataIOTag)                                                           \
    template class SpecializedContainerVariableData<ContainerType, ContainerDataIO<ContainerDataIOTag>>;                                             \
    KRATOS_INSTANTIATE_CONTAINER_DATA_METHODS(ContainerType, ContainerDataIO<ContainerDataIOTag>, double)                                             \
    KRATOS_INSTANTIATE_CONTAINER_DATA_METHODS(ContainerType, ContainerDataIO<ContainerDataIOTag>, array_1d<double, 3>)

KRATOS_INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ModelPart::NodesContainerType, ContainerDataIOTags::Historical)
KRATOS_INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ModelPart::NodesContainerType, ContainerDataIOTags::NonHistorical)
KRATOS_INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ModelPart::ConditionsContainerType, ContainerDataIOTags::NonHistorical)
KRATOS_INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ModelPart::ElementsContainerType, ContainerDataIOTags::NonHistorical)
KRATOS_INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ModelPart::ConditionsContainerType, ContainerDataIOTags::Properties)
KRATOS_INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ModelPart::ElementsContainerType, ContainerDataIOTags::Properties)

#undef KRATOS_INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER
#undef KRATOS_INSTANTIATE_CONTAINER_DATA_METHODS

}