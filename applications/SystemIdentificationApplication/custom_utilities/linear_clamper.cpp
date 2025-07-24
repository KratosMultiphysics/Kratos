//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <cmath>

// External includes

// Project includes
#include "expression/literal_flat_expression.h"
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "linear_clamper.h"

namespace Kratos {

template<class TContainerType>
LinearClamper<TContainerType>::LinearClamper(
    const double Min,
    const double Max)
    : mMin(Min),
      mMax(Max),
      mDelta(mMax - mMin)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mMin >= mMax) << "Min must be lower than Max";

    KRATOS_CATCH("");
}

template<class TContainerType>
double LinearClamper<TContainerType>::ProjectForward(const double X) const
{
    const double x_tilde = std::clamp(X, 0.0, 1.0);
    return mMin + x_tilde * mDelta;
}

template<class TContainerType>
double LinearClamper<TContainerType>::CalculateForwardProjectionGradient(const double X) const
{
    return mDelta;
}

template<class TContainerType>
double LinearClamper<TContainerType>::ProjectBackward(const double Y) const
{
    return  std::clamp((Y - mMin) / mDelta, 0.0, 1.0);
}

template<class TContainerType>
ContainerExpression<TContainerType> LinearClamper<TContainerType>::ProjectForward(const ContainerExpression<TContainerType>& rInput) const
{
    KRATOS_TRY

    // x*x*(3.0-2.0*x);

    const auto& r_input_exp = rInput.GetExpression();
    const auto number_of_entities = r_input_exp.NumberOfEntities();
    const auto stride = r_input_exp.GetItemComponentCount();

    KRATOS_ERROR_IF_NOT(stride == 1)
        << "LinearClamper only supports scalar expressions. [ Expression stride = " << stride << " ].\n";

    auto p_result_exp = LiteralFlatExpression<double>::Create(number_of_entities, {});

    IndexPartition<IndexType>(number_of_entities).for_each([&](const auto Index) {
        const double x = r_input_exp.Evaluate(Index, Index, 0);
        *(p_result_exp->begin() + Index) = this->ProjectForward(x);
    });

    auto result = rInput;
    result.SetExpression(p_result_exp);
    return result;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> LinearClamper<TContainerType>::CalculateForwardProjectionGradient(const ContainerExpression<TContainerType>& rInput) const
{
    KRATOS_TRY

    // x*x*(3.0-2.0*x);
    // y = 3x^2 - 2x^3
    // dy/dx = 3.2.x - 2.3.x^2

    const auto& r_input_exp = rInput.GetExpression();
    const auto number_of_entities = r_input_exp.NumberOfEntities();
    const auto stride = r_input_exp.GetItemComponentCount();

    KRATOS_ERROR_IF_NOT(stride == 1)
        << "LinearClamper only supports scalar expressions. [ Expression stride = " << stride << " ].\n";

    auto p_result_exp = LiteralFlatExpression<double>::Create(number_of_entities, {});

    IndexPartition<IndexType>(number_of_entities).for_each([&](const auto Index) {
        const double x = r_input_exp.Evaluate(Index, Index, 0);
        *(p_result_exp->begin() + Index) = this->CalculateForwardProjectionGradient(x);
    });

    auto result = rInput;
    result.SetExpression(p_result_exp);
    return result;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> LinearClamper<TContainerType>::ProjectBackward(const ContainerExpression<TContainerType>& rInput) const
{
    KRATOS_TRY

    const auto& r_input_exp = rInput.GetExpression();
    const auto number_of_entities = r_input_exp.NumberOfEntities();
    const auto stride = r_input_exp.GetItemComponentCount();

    KRATOS_ERROR_IF_NOT(stride == 1)
        << "LinearClamper only supports scalar expressions. [ Expression stride = " << stride << " ].\n";

    auto p_result_exp = LiteralFlatExpression<double>::Create(number_of_entities, {});

    IndexPartition<IndexType>(number_of_entities).for_each([&](const auto Index) {
        const double x = r_input_exp.Evaluate(Index, Index, 0);
        *(p_result_exp->begin() + Index) = this->ProjectBackward(x);
    });

    auto result = rInput;
    result.SetExpression(p_result_exp);
    return result;

    KRATOS_CATCH("");
}


// template instantiations
template class LinearClamper<ModelPart::NodesContainerType>;
template class LinearClamper<ModelPart::ConditionsContainerType>;
template class LinearClamper<ModelPart::ElementsContainerType>;

} /* namespace Kratos.*/