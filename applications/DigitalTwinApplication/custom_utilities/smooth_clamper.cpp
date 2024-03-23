//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         DigitalTwinApplication/license.txt
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
ContainerExpression<TContainerType> SmoothClamper<TContainerType>::Clamp(const ContainerExpression<TContainerType>& rInput) const
{
    KRATOS_TRY

    // x*x*(3.0-2.0*x);

    const auto min = mMin;
    const auto delta = mMax - mMin;
    const auto& r_input_exp = rInput.GetExpression();
    const auto number_of_entities = r_input_exp.NumberOfEntities();
    const auto stride = r_input_exp.GetItemComponentCount();

    KRATOS_ERROR_IF_NOT(stride == 1)
        << "SmoothClamper only supports scalar expressions. [ Expression stride = " << stride << " ].\n";

    auto p_result_exp = LiteralFlatExpression<double>::Create(number_of_entities, {});

    IndexPartition<IndexType>(number_of_entities).for_each([&p_result_exp, &r_input_exp, min, delta](const auto Index) {
        const double value = r_input_exp.Evaluate(Index, Index, 0);
        const double x = (value - min) / delta;
        double& y = *(p_result_exp->begin() + Index);
        if (x < 0.0) {
            y = 0.0;
        } else if (x > 1.0) {
            y = 1.0;
        } else {
            y = x * x * (3.0 - 2.0 * x);
        }
    });

    auto result = rInput;
    result.SetExpression(p_result_exp);
    return result;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> SmoothClamper<TContainerType>::ClampDerivative(const ContainerExpression<TContainerType>& rInput) const
{
    KRATOS_TRY

    // x*x*(3.0-2.0*x);
    // y = 3x^2 - 2x^3
    // dy/dx = 3.2.x - 2.3.x^2

    const auto min = mMin;
    const auto delta = mMax - mMin;
    const auto& r_input_exp = rInput.GetExpression();
    const auto number_of_entities = r_input_exp.NumberOfEntities();
    const auto stride = r_input_exp.GetItemComponentCount();

    KRATOS_ERROR_IF_NOT(stride == 1)
        << "SmoothClamper only supports scalar expressions. [ Expression stride = " << stride << " ].\n";

    auto p_result_exp = LiteralFlatExpression<double>::Create(number_of_entities, {});

    IndexPartition<IndexType>(number_of_entities).for_each([&p_result_exp, &r_input_exp, min, delta](const auto Index) {
        const double value = r_input_exp.Evaluate(Index, Index, 0);
        const double x = (value - min) / delta;
        double& y = *(p_result_exp->begin() + Index);
        if (x < 0.0) {
            y = 0.0;
        } else if (x > 1.0) {
            y = 0.0;
        } else {
            y = (6 * x - 6 * x * x) / delta;
        }
    });

    auto result = rInput;
    result.SetExpression(p_result_exp);
    return result;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> SmoothClamper<TContainerType>::InverseClamp(const ContainerExpression<TContainerType>& rInput) const
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
        double& value = *(p_result_exp->begin() + Index);
        const double y = r_input_exp.Evaluate(Index, Index, 0);
        if (y < 0.0) {
            value = min;
        } else if (y > 1.0) {
            value = max;
        } else {
            const double x = 0.5 - std::sin(std::asin(1.0 - 2.0 * y) / 3.0);
            value = x * delta + min;
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