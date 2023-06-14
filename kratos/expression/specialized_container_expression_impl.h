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
#include "expression/container_data_io.h"
#include "expression/variable_expression_data_io.h"
#include "expression/literal_expression.h"
#include "expression/literal_flat_expression.h"
#include "expression/arithmetic_operators.h"
#include "expression/view_operators.h"
#include "expression/variable_expression_io.h"

namespace Kratos {

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator=(const SpecializedContainerExpression& rOther)
{
    KRATOS_ERROR_IF(this->GetContainer().size() != rOther.GetContainer().size())
        << "Mismatching model parts found with different number of entities in assignment operation.\n"
        << "      Assignee data: " << *this << "\n"
        << "      Assignor data: " << rOther << "\n";

    this->mpExpression = rOther.mpExpression;
    return *this;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
typename SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Pointer SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Clone() const
{
    return Kratos::make_shared<SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>>(*this);
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
template<class TDataType>
void SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Read(const Variable<TDataType>& rVariable)
{
    KRATOS_TRY

        const auto& r_container = this->GetContainer();
    const IndexType number_of_entities = r_container.size();

    using raw_data_type = std::conditional_t<
                            std::is_same_v<TDataType, char>, char,
                            std::conditional_t<
                                std::is_same_v<TDataType, int>, int,
                                double
                            >
                          >;

    if (number_of_entities != 0) {
        // initialize the shape with the first entity value
        VariableExpressionDataIO<TDataType> variable_flatten_data_io(TContainerDataIO::GetValue(*r_container.begin(), rVariable));

        auto p_expression = LiteralFlatExpression<raw_data_type>::Create(number_of_entities, variable_flatten_data_io.GetItemShape());
        auto& r_expression = *p_expression;
        this->mpExpression = p_expression;

        IndexPartition<IndexType>(number_of_entities).for_each([&r_container, &rVariable, &variable_flatten_data_io, &r_expression](const IndexType Index){
            const auto& values = TContainerDataIO::GetValue(*(r_container.begin() + Index), rVariable);
            variable_flatten_data_io.Read(r_expression, Index, values);
        });
    }

    KRATOS_CATCH("")
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
template<class TDataType>
void SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Evaluate(const Variable<TDataType>& rVariable)
{
    KRATOS_TRY

    auto& r_container = this->GetContainer();
    const IndexType number_of_entities = r_container.size();

    if (number_of_entities > 0) {
        const auto& r_expression = this->GetExpression();

        VariableExpressionDataIO<TDataType> variable_flatten_data_io(r_expression.GetItemShape());

        // initialize the container variables first
        if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
            // initializes ghost nodes as for the later synchronization
            // only, the nodal non historical values needs to be set unless
            // they are properly initialized. Otherwise, in synchronization, the variables will
            // not be there in the ghost nodes hence seg faults.

            // the vectors and matrices needs to be initialized in historical and non-historical
            // data containers because they need to be initialized with the correct size for synchronization
            if constexpr(std::is_same_v<TDataType, Vector> || std::is_same_v<TDataType, Matrix>) {
                TDataType dummy_value{};
                variable_flatten_data_io.Assign(dummy_value, r_expression, 0);
                if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::Historical>>) {
                    VariableUtils().SetVariable(rVariable, dummy_value, this->GetModelPart().GetCommunicator().GhostMesh().Nodes());
                } else if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::NonHistorical>>) {
                    VariableUtils().SetNonHistoricalVariable(rVariable, dummy_value, this->GetModelPart().GetCommunicator().GhostMesh().Nodes());
                }
            } else {
                // if it is a static type, then it only needs to be initialized in the non-historical container with zeros.
                // historical container should be initialized with the default values when the container is created.
                if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::NonHistorical>>) {
                    VariableUtils().SetNonHistoricalVariableToZero(rVariable, this->GetModelPart().GetCommunicator().GhostMesh().Nodes());
                }
            }
        }

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
    }

    KRATOS_CATCH("");
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
template<class TDataType>
void SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::SetData(const TDataType& rValue)
{
    this->mpExpression = LiteralExpression<TDataType>::Create(rValue, this->GetContainer().size());
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
template<class TDataType>
void SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::SetZero(const Variable<TDataType>& rVariable)
{
    this->SetData(rVariable.Zero());
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Slice(
    const IndexType Offset,
    const IndexType Stride) const
{
    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(*(this->mpModelPart));
    result.mpExpression = Kratos::Slice(*this->mpExpression, Offset, Stride);
    return result;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Reshape(const std::vector<IndexType>& rShape) const
{
    return this->Reshape(rShape.begin(), rShape.end());
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Comb(const BaseType& rOther) const
{
    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(*(this->mpModelPart));
    std::vector<Expression::ConstPointer> expressions;
    expressions.push_back(this->pGetExpression());
    expressions.push_back(rOther.pGetExpression());
    result.mpExpression = Kratos::Comb(expressions.begin(), expressions.end());
    return result;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Comb(const std::vector<typename BaseType::Pointer>& rListOfOthers) const
{
    return this->Comb(rListOfOthers.begin(), rListOfOthers.end());
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator+(const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rOther) const
{
    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(static_cast<const BaseType&>(*this) + static_cast<const BaseType&>(rOther));
    return result;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator+=(const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rOther)
{
    static_cast<BaseType&>(*this)+= static_cast<const BaseType&>(rOther);
    return *this;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator+(const double Value) const
{

    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(static_cast<const BaseType&>(*this) + Value);
    return result;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator+=(const double Value)
{
    static_cast<BaseType&>(*this)+= Value;
    return *this;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator-(const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rOther) const
{
    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(static_cast<const BaseType&>(*this) - static_cast<const BaseType&>(rOther));
    return result;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator-=(const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rOther)
{
    static_cast<BaseType&>(*this)-= static_cast<const BaseType&>(rOther);
    return *this;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator-(const double Value) const
{

    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(static_cast<const BaseType&>(*this) - Value);
    return result;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator-=(const double Value)
{
    static_cast<BaseType&>(*this)-= Value;
    return *this;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator*(const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rOther) const
{
    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(static_cast<const BaseType&>(*this) * static_cast<const BaseType&>(rOther));
    return result;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator*=(const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rOther)
{
    static_cast<BaseType&>(*this)*= static_cast<const BaseType&>(rOther);
    return *this;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator*(const double Value) const
{

    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(static_cast<const BaseType&>(*this) * Value);
    return result;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator*=(const double Value)
{
    static_cast<BaseType&>(*this)*= Value;
    return *this;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator/(const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rOther) const
{
    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(static_cast<const BaseType&>(*this) / static_cast<const BaseType&>(rOther));
    return result;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator/=(const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rOther)
{
    static_cast<BaseType&>(*this)/= static_cast<const BaseType&>(rOther);
    return *this;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator/(const double Value) const
{

    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(static_cast<const BaseType&>(*this) / Value);
    return result;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::operator/=(const double Value)
{
    static_cast<BaseType&>(*this)/= Value;
    return *this;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Power(const SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>& rOther) const
{
    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(Kratos::Power(static_cast<const BaseType&>(*this), static_cast<const BaseType&>(rOther)));
    return result;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType>::Power(const double Value) const
{
    SpecializedContainerExpression<TContainerType, TContainerDataIO, TMeshType> result(Kratos::Power(static_cast<const BaseType&>(*this), Value));
    return result;
}

template <class TContainerType, class TContainerDataIO, MeshType TMeshType>
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