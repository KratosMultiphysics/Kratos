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
#include "smooth_clamper.h"

namespace Kratos {

template<class TContainerType>
SmoothClamper<TContainerType>::SmoothClamper(
    const double Min,
    const double Max)
    : mMin(Min),
      mMax(Max)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mMin >= mMax) << "Min must be lower than Max";

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> SmoothClamper<TContainerType>::ProjectForward(const ContainerExpression<TContainerType>& rInput) const
{
    KRATOS_TRY

    // x*x*(3.0-2.0*x);

    const auto min = mMin;
    const auto max = mMax;
    const auto delta = mMax - mMin;
    const auto& r_input_exp = rInput.GetExpression();
    const auto number_of_entities = r_input_exp.NumberOfEntities();
    const auto stride = r_input_exp.GetItemComponentCount();

    KRATOS_ERROR_IF_NOT(stride == 1)
        << "SmoothClamper only supports scalar expressions. [ Expression stride = " << stride << " ].\n";

    auto p_result_exp = LiteralFlatExpression<double>::Create(number_of_entities, {});

    IndexPartition<IndexType>(number_of_entities).for_each([&p_result_exp, &r_input_exp, min, max, delta](const auto Index) {
        const double x = r_input_exp.Evaluate(Index, Index, 0);
        double& value = *(p_result_exp->begin() + Index);
        if (x < 0.0) {
            value = min;
        } else if (x > 1.0) {
            value = max;
        } else {
            const double y = x * x * (3.0 - 2.0 * x);
            value = min + y * delta;
        }
    });

    auto result = rInput;
    result.SetExpression(p_result_exp);
    return result;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> SmoothClamper<TContainerType>::CalculateForwardProjectionGradient(const ContainerExpression<TContainerType>& rInput) const
{
    KRATOS_TRY

    // x*x*(3.0-2.0*x);
    // y = 3x^2 - 2x^3
    // dy/dx = 3.2.x - 2.3.x^2

    const auto delta = mMax - mMin;
    const auto& r_input_exp = rInput.GetExpression();
    const auto number_of_entities = r_input_exp.NumberOfEntities();
    const auto stride = r_input_exp.GetItemComponentCount();

    KRATOS_ERROR_IF_NOT(stride == 1)
        << "SmoothClamper only supports scalar expressions. [ Expression stride = " << stride << " ].\n";

    auto p_result_exp = LiteralFlatExpression<double>::Create(number_of_entities, {});

    IndexPartition<IndexType>(number_of_entities).for_each([&p_result_exp, &r_input_exp, delta](const auto Index) {
        const double x = r_input_exp.Evaluate(Index, Index, 0);
        double& value = *(p_result_exp->begin() + Index);
        if (x < 0.0) {
            value = 0.0;
        } else if (x > 1.0) {
            value = 0.0;
        } else {
            value = (6 * x - 6 * x * x) * delta;
        }
    });

    auto result = rInput;
    result.SetExpression(p_result_exp);
    return result;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> SmoothClamper<TContainerType>::ProjectBackward(const ContainerExpression<TContainerType>& rInput) const
{
    KRATOS_TRY

    const auto min = mMin;
    const auto max = mMax;
    const auto delta = mMax - mMin;
    const auto& r_input_exp = rInput.GetExpression();
    const auto number_of_entities = r_input_exp.NumberOfEntities();
    const auto stride = r_input_exp.GetItemComponentCount();

    KRATOS_ERROR_IF_NOT(stride == 1)
        << "SmoothClamper only supports scalar expressions. [ Expression stride = " << stride << " ].\n";

    auto p_result_exp = LiteralFlatExpression<double>::Create(number_of_entities, {});

    IndexPartition<IndexType>(number_of_entities).for_each([&p_result_exp, &r_input_exp, min, delta, max](const auto Index) {
        const double value = r_input_exp.Evaluate(Index, Index, 0);
        double& x = *(p_result_exp->begin() + Index);
        if (value < min) {
            x = 0.0;
        } else if (value > max) {
            x = 1.0;
        } else {
            const double y = (value - min) / delta;
            x = 0.5 - std::sin(std::asin(1.0 - 2.0 * y) / 3.0);
        }
    });

    auto result = rInput;
    result.SetExpression(p_result_exp);
    return result;

    KRATOS_CATCH("");
}


// template instantiations
template class SmoothClamper<ModelPart::NodesContainerType>;
template class SmoothClamper<ModelPart::ConditionsContainerType>;
template class SmoothClamper<ModelPart::ElementsContainerType>;

} /* namespace Kratos.*/