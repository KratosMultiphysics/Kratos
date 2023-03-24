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
#include "container_data_io.h"
#include "variable_flatten_data_io.h"

// Include base h
#include "specialized_container_variable_data.h"

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

} // namespace ContainerVariableDataHolderHelperUtilities

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
        KRATOS_ERROR_IF_NOT(this->GetModelPart().GetProcessInfo().Has(DOMAIN_SIZE))
            << "DOMAIN_SIZE variable is not found in process info of model part "
            "with name "
            << this->GetModelPart().FullName() << ".\n";

        const IndexType domain_size = this->GetModelPart().GetProcessInfo()[DOMAIN_SIZE];

        // initialize the shape with the first entity value
        VariableFlattenDataIO<TDataType> variable_flatten_data_io(TContainerDataIO::GetValue(*r_container.begin(), rVariable), domain_size);

        auto p_data = Kratos::make_shared<Vector>();
        this->mpExpression = LiteralVectorExpression::Create(p_data, variable_flatten_data_io.GetShape());
        const IndexType local_size = this->mpExpression->GetLocalSize();

        auto& r_data = *p_data;
        r_data.resize(number_of_entities * local_size);

        IndexPartition<IndexType>(number_of_entities).for_each([&r_container, &rVariable, &r_data, &variable_flatten_data_io, local_size](const IndexType Index){
            const auto& values = TContainerDataIO::GetValue(*(r_container.begin() + Index), rVariable);
            variable_flatten_data_io.Read(&r_data[Index * local_size], values);
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
    const IndexType local_size = r_expression.GetLocalSize();

    struct tls_type
    {
        tls_type(const IndexType LocalSize) { mFlattenedVector.resize(LocalSize); }

        TDataType mValue{};

        Vector mFlattenedVector;
    };

    VariableFlattenDataIO<TDataType> variable_flatten_data_io(r_expression.GetShape());

    IndexPartition<IndexType>(number_of_entities).for_each(tls_type(local_size), [&r_container, &rVariable, &local_size, &r_expression, &variable_flatten_data_io](const IndexType Index, tls_type& rTLS){
        for (IndexType i = 0; i < local_size; ++i) {
            rTLS.mFlattenedVector[i] = r_expression.Evaluate(Index, i);
        }
        variable_flatten_data_io.Assign(rTLS.mValue, &rTLS.mFlattenedVector[0]);
        TContainerDataIO::SetValue(*(r_container.begin() + Index), rVariable, rTLS.mValue);
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

    if constexpr(std::is_same_v<TDataType, double> || std::is_same_v<TDataType, int> || std::is_same_v<TDataType, std::size_t>) {
        this->mpExpression = LiteralDoubleExpression::Create(rValue);
    } else if constexpr(std::is_same_v<TDataType, array_1d<double, 3>>) {
        VariableFlattenDataIO<TDataType> variable_flatten_data_io(rValue, this->GetModelPart().GetProcessInfo()[DOMAIN_SIZE]);
        this->mpExpression = LiteralArray3Expression::Create(rValue, variable_flatten_data_io.GetShape()[0]);
    }

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
    result.mpExpression = BinaryAddExpression::Create(this->mpExpression, rOther.mpExpression);
    return result;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO>& SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator+=(const SpecializedContainerVariableData<TContainerType, TContainerDataIO>& rOther)
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in addition.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    this->mpExpression = BinaryAddExpression::Create(this->mpExpression, rOther.mpExpression);
    return *this;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO> SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator+(const double Value) const
{

    SpecializedContainerVariableData<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinaryAddExpression::Create(this->mpExpression, LiteralDoubleExpression::Create(Value));
    return result;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO>& SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator+=(const double Value)
{
    this->mpExpression = BinaryAddExpression::Create(this->mpExpression, LiteralDoubleExpression::Create(Value));
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
    result.mpExpression = BinarySubstractExpression::Create(this->mpExpression, rOther.mpExpression);
    return result;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO>& SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator-=(const SpecializedContainerVariableData<TContainerType, TContainerDataIO>& rOther)
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    this->mpExpression = BinarySubstractExpression::Create(this->mpExpression, rOther.mpExpression);
    return *this;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO> SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator-(const double Value) const
{
    SpecializedContainerVariableData<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinarySubstractExpression::Create(this->mpExpression, LiteralDoubleExpression::Create(Value));
    return result;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO>& SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator-=(const double Value)
{
    this->mpExpression = BinarySubstractExpression::Create(this->mpExpression, LiteralDoubleExpression::Create(Value));
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
    result.mpExpression = BinaryMultiplyExpression::Create(this->mpExpression, rOther.mpExpression);
    return result;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO>& SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator*=(const SpecializedContainerVariableData<TContainerType, TContainerDataIO>& rOther)
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    this->mpExpression = BinaryMultiplyExpression::Create(this->mpExpression, rOther.mpExpression);
    return *this;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO> SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator*(const double Value) const
{
    SpecializedContainerVariableData<TContainerType, TContainerDataIO> result(*(this->mpModelPart));
    result.mpExpression = BinaryMultiplyExpression::Create(this->mpExpression, LiteralDoubleExpression::Create(Value));
    return result;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO>& SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator*=(const double Value)
{
    this->mpExpression = BinaryMultiplyExpression::Create(this->mpExpression, LiteralDoubleExpression::Create(Value));
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
    result.mpExpression = BinaryDivideExpression::Create(this->mpExpression, rOther.mpExpression);
    return result;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO>& SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator/=(const SpecializedContainerVariableData<TContainerType, TContainerDataIO>& rOther)
{
    KRATOS_ERROR_IF(&this->GetModelPart() != &rOther.GetModelPart())
        << "Mismatching model parts found in substraction.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Right operand data: " << rOther << "\n";

    this->mpExpression = BinaryDivideExpression::Create(this->mpExpression, rOther.mpExpression);
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
    result.mpExpression = BinaryMultiplyExpression::Create(this->mpExpression, LiteralDoubleExpression::Create(1.0 / Value));
    return result;
}

template <class TContainerType, class TContainerDataIO>
SpecializedContainerVariableData<TContainerType, TContainerDataIO>& SpecializedContainerVariableData<TContainerType, TContainerDataIO>::operator/=(const double Value)
{
    KRATOS_ERROR_IF(std::abs(Value) < std::numeric_limits<double>::epsilon())
        << "Division by zero.\n"
        << "      Left operand data : " << *this << "\n"
        << "      Divisor           : " << Value << "\n";

    this->mpExpression = BinaryMultiplyExpression::Create(this->mpExpression, LiteralDoubleExpression::Create(1.0 / Value));
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
#define KRATOS_INSTANTIATE_CONTAINER_DATA_METHODS(ContainerType, ContainerDataIOType, DataType)                                                                   \
    template void SpecializedContainerVariableData<ContainerType, ContainerDataIOType>::ReadData(const Variable<DataType>&);                \
    template void SpecializedContainerVariableData<ContainerType, ContainerDataIOType>::SetData(const DataType&); \
    template void SpecializedContainerVariableData<ContainerType, ContainerDataIOType>::SetZero(const Variable<DataType>&);            \
    template void SpecializedContainerVariableData<ContainerType, ContainerDataIOType>::AssignData(const Variable<DataType>&);

#define KRATOS_INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ContainerType, ContainerDataIOTag)                                                           \
    template class SpecializedContainerVariableData<ContainerType, ContainerDataIO<ContainerDataIOTag>>;                                             \
    KRATOS_INSTANTIATE_CONTAINER_DATA_METHODS(ContainerType, ContainerDataIO<ContainerDataIOTag>, double)                                             \
    KRATOS_INSTANTIATE_CONTAINER_DATA_METHODS(ContainerType, ContainerDataIO<ContainerDataIOTag>, ContainerVariableDataHolderHelperUtilities::Array3D)

KRATOS_INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ModelPart::NodesContainerType, ContainerDataIOTags::Historical)
KRATOS_INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ModelPart::NodesContainerType, ContainerDataIOTags::NonHistorical)
KRATOS_INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ModelPart::ConditionsContainerType, ContainerDataIOTags::NonHistorical)
KRATOS_INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ModelPart::ElementsContainerType, ContainerDataIOTags::NonHistorical)
KRATOS_INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ModelPart::ConditionsContainerType, ContainerDataIOTags::Properties)
KRATOS_INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER(ModelPart::ElementsContainerType, ContainerDataIOTags::Properties)

#undef KRATOS_INSTANTIATE_CONTAINER_VARIABLE_DATA_HOLDER
#undef KRATOS_INSTANTIATE_CONTAINER_DATA_METHODS

}