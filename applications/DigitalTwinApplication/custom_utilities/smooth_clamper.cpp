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
      mMax(Max),
      mDelta(mMax - mMin)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mMin >= mMax) << "Min must be lower than Max";

    KRATOS_CATCH("");
}

template<class TContainerType>
double  SmoothClamper<TContainerType>::Clamp(const double X) const
{
    double Y;
    if (X < 0.0) {
        Y = mMin;
    } else if (X > 1.0) {
        Y = mMax;
    } else {
        Y = mMin + X * X * (3.0 - 2.0 * X) * mDelta;
    }
    return Y;
}

template<class TContainerType>
double  SmoothClamper<TContainerType>::ClampDerivative(const double X) const
{
    double dY_dX;
    if (X < 0.0) {
        dY_dX = 0.0;
    } else if (X > 1.0) {
        dY_dX = 0.0;
    } else {
        dY_dX = (6 * X - 6 * X * X) * mDelta;
    }
    return dY_dX;
}

template<class TContainerType>
double  SmoothClamper<TContainerType>::InverseClamp(const double Y) const
{
    double x;
    if (Y < mMin) {
        x = 0.0;
    } else if (Y > mMax) {
        x = 1.0;
    } else {
        const double y = (Y - mMin) / mDelta;
        x = 0.5 - std::sin(std::asin(1.0 - 2.0 * y) / 3.0);
    }
    return x;
}

template<class TContainerType>
ContainerExpression<TContainerType> SmoothClamper<TContainerType>::Clamp(const ContainerExpression<TContainerType>& rInput) const
{
    KRATOS_TRY

    // x*x*(3.0-2.0*x);

    const auto& r_input_exp = rInput.GetExpression();
    const auto number_of_entities = r_input_exp.NumberOfEntities();
    const auto stride = r_input_exp.GetItemComponentCount();

    KRATOS_ERROR_IF_NOT(stride == 1)
        << "SmoothClamper only supports scalar expressions. [ Expression stride = " << stride << " ].\n";

    auto p_result_exp = LiteralFlatExpression<double>::Create(number_of_entities, {});

    IndexPartition<IndexType>(number_of_entities).for_each([&](const auto Index) {
        const double x = r_input_exp.Evaluate(Index, Index, 0);
        *(p_result_exp->begin() + Index) = this->Clamp(x);
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

    const auto& r_input_exp = rInput.GetExpression();
    const auto number_of_entities = r_input_exp.NumberOfEntities();
    const auto stride = r_input_exp.GetItemComponentCount();

    KRATOS_ERROR_IF_NOT(stride == 1)
        << "SmoothClamper only supports scalar expressions. [ Expression stride = " << stride << " ].\n";

    auto p_result_exp = LiteralFlatExpression<double>::Create(number_of_entities, {});

    IndexPartition<IndexType>(number_of_entities).for_each([&](const auto Index) {
        const double x = r_input_exp.Evaluate(Index, Index, 0);
        *(p_result_exp->begin() + Index) = this->ClampDerivative(x);
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

    const auto& r_input_exp = rInput.GetExpression();
    const auto number_of_entities = r_input_exp.NumberOfEntities();
    const auto stride = r_input_exp.GetItemComponentCount();

    KRATOS_ERROR_IF_NOT(stride == 1)
        << "SmoothClamper only supports scalar expressions. [ Expression stride = " << stride << " ].\n";

    auto p_result_exp = LiteralFlatExpression<double>::Create(number_of_entities, {});

    IndexPartition<IndexType>(number_of_entities).for_each([&](const auto Index) {
        const double value = r_input_exp.Evaluate(Index, Index, 0);
        *(p_result_exp->begin() + Index) = this->InverseClamp(value);
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