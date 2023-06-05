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
#include <string>
#include <tuple>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"
#include "containers/container_expression/container_data_io.h"
#include "containers/container_expression/variable_expression_data_io.h"
#include "containers/container_expression/expressions/literal/literal_expression.h"
#include "containers/container_expression/expressions/literal/literal_flat_expression.h"
#include "containers/container_expression/expressions/arithmetic_operators.h"
#include "containers/container_expression/expressions/view_operators.h"
#include "containers/container_expression/expressions/io/variable_expression_io.h"

namespace Kratos {

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator=(const SpecializedContainerExpression& rOther)
{
    KRATOS_ERROR_IF(this->GetContainer().size() != rOther.GetContainer().size())
        << "Mismatching model parts found with different number of entities in assignment operation.\n"
        << "      Assignee data: " << *this << "\n"
        << "      Assignor data: " << rOther << "\n";

    this->mpExpression = rOther.mpExpression;
    return *this;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
typename SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Pointer SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Clone() const
{
    return Kratos::make_shared<SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>>(*this);
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
template<class TDataType>
void SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Read(const Variable<TDataType>& rVariable)
{
    KRATOS_TRY

    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        this->mpExpression = VariableExpressionIO::VariableExpressionInput(*this, rVariable, std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::Historical>>).Execute();
    } else {
        this->mpExpression = VariableExpressionIO::VariableExpressionInput(*this, rVariable).Execute();
    }

    KRATOS_CATCH("")
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
template<class TDataType>
void SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Evaluate(const Variable<TDataType>& rVariable)
{
    KRATOS_TRY

    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        VariableExpressionIO::VariableExpressionOutput(*this, rVariable, std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::Historical>>).Execute(**this->mpExpression);
    } else {
        VariableExpressionIO::VariableExpressionOutput(*this, rVariable).Execute(**this->mpExpression);
    }

    KRATOS_CATCH("");
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
template<class TDataType>
void SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::SetData(const TDataType& rValue)
{
    this->mpExpression = LiteralExpression<TDataType>::Create(rValue, this->GetContainer().size());
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
template<class TDataType>
void SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::SetZero(const Variable<TDataType>& rVariable)
{
    this->SetData(rVariable.Zero());
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Slice(
    const IndexType Offset,
    const IndexType Stride) const
{
    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(*(this->mpModelPart));
    result.mpExpression = Kratos::Slice(*this->mpExpression, Offset, Stride);
    return result;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Reshape(const std::vector<IndexType>& rShape) const
{
    return this->Reshape(rShape.begin(), rShape.end());
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
template<class TIteratorType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Reshape(
    TIteratorType Begin,
    TIteratorType End) const
{
    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(*(this->mpModelPart));
    result.mpExpression = Kratos::Reshape(*this->mpExpression, Begin, End);
    return result;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Comb(const BaseType& rOther) const
{
    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(*(this->mpModelPart));
    std::vector<Expression::Pointer> expressions;
    expressions.push_back(this->pGetExpression());
    expressions.push_back(rOther.pGetExpression());
    result.mpExpression = Kratos::Comb(expressions.begin(), expressions.end());
    return result;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Comb(const std::vector<typename BaseType::Pointer>& rListOfOthers) const
{
    return this->Comb(rListOfOthers.begin(), rListOfOthers.end());
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
template<class TIteratorType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Comb(
    TIteratorType Begin,
    TIteratorType End) const
{
    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(*(this->mpModelPart));
    std::vector<Expression::Pointer> expressions;
    expressions.push_back(this->pGetExpression());
    for (auto itr = Begin; itr != End; ++itr) {
        expressions.push_back((*itr)->pGetExpression());
    }
    result.mpExpression = Kratos::Comb(expressions.begin(), expressions.end());
    return result;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator+(const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rOther) const
{
    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(static_cast<BaseType>(*this) + static_cast<BaseType>(rOther));
    return result;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator+=(const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rOther)
{
    static_cast<BaseType>(*this)+= static_cast<BaseType>(rOther);
    return *this;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator+(const double Value) const
{

    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(static_cast<BaseType>(*this) + Value);
    return result;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator+=(const double Value)
{
    static_cast<BaseType>(*this)+= Value;
    return *this;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator-(const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rOther) const
{
    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(static_cast<BaseType>(*this) - static_cast<BaseType>(rOther));
    return result;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator-=(const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rOther)
{
    static_cast<BaseType>(*this)-= static_cast<BaseType>(rOther);
    return *this;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator-(const double Value) const
{

    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(static_cast<BaseType>(*this) - Value);
    return result;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator-=(const double Value)
{
    static_cast<BaseType>(*this)-= Value;
    return *this;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator*(const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rOther) const
{
    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(static_cast<BaseType>(*this) * static_cast<BaseType>(rOther));
    return result;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator*=(const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rOther)
{
    static_cast<BaseType>(*this)*= static_cast<BaseType>(rOther);
    return *this;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator*(const double Value) const
{

    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(static_cast<BaseType>(*this) * Value);
    return result;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator*=(const double Value)
{
    static_cast<BaseType>(*this)*= Value;
    return *this;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator/(const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rOther) const
{
    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(static_cast<BaseType>(*this) / static_cast<BaseType>(rOther));
    return result;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator/=(const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rOther)
{
    static_cast<BaseType>(*this)/= static_cast<BaseType>(rOther);
    return *this;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator/(const double Value) const
{

    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(static_cast<BaseType>(*this) / Value);
    return result;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator/=(const double Value)
{
    static_cast<BaseType>(*this)/= Value;
    return *this;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Pow(const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rOther) const
{
    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(Kratos::Pow(static_cast<BaseType>(*this), static_cast<BaseType>(rOther)));
    return result;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Pow(const double Value) const
{
    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(Kratos::Pow(static_cast<BaseType>(*this), Value));
    return result;
}

template <class TContainerType, class TContainerDataIO, class TMeshType>
std::string SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Info() const
{
    std::stringstream msg;

    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        msg << "Nodal";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        msg << "Condition";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        msg << "Element";
    }

    msg << TContainerDataIO::mInfo << BaseType::Info();

    return msg.str();
}

}