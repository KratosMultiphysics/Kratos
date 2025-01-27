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
#include <limits>

// External includes

// Project includes
#include "expression/literal_expression.h"
#include "expression/literal_flat_expression.h"

// Application includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Include base h
#include "boltzmann_operator.h"

namespace Kratos {

template<class TContainerType>
BoltzmannOperator<TContainerType>::BoltzmannOperator(
    const double Beta,
    const double ScalingFactor)
    : mBeta(Beta),
      mScalingFactor(ScalingFactor)
{
}

template<class TContainerType>
double BoltzmannOperator<TContainerType>::CalculateValue() const
{
    if (mDenominator > std::numeric_limits<double>::epsilon()) {
        return mScalingFactor * mNumerator / mDenominator;
    } else {
        return 0.0;
    }
}

template<class TContainerType>
typename ContainerExpression<TContainerType>::Pointer BoltzmannOperator<TContainerType>::CalculateGradient() const
{
    KRATOS_TRY

    const auto& r_expression = mpContainerExpression->GetExpression();
    const IndexType number_of_elements = r_expression.NumberOfEntities();
    const IndexType number_of_components =  r_expression.GetItemComponentCount();

    auto result = mpContainerExpression->Clone();

    if (mDenominator > std::numeric_limits<double>::epsilon()) {
        auto p_result_exp = LiteralFlatExpression<double>::Create(number_of_elements, r_expression.GetItemShape());

        IndexPartition<IndexType>(number_of_elements).for_each([&](const auto Index) {
            const IndexType entity_data_begin_index = (Index * number_of_components);

            for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
                const double value = r_expression.Evaluate(Index, entity_data_begin_index, i_comp) / mScalingFactor;

                const double exp_beta_x_k = std::exp(mBeta * value);
                const double numerator_gradient = exp_beta_x_k * (1 + mBeta * value);
                const double denominator_gradient = mBeta * exp_beta_x_k;

                *(p_result_exp->begin() + entity_data_begin_index + i_comp) =
                    (mDenominator * numerator_gradient - mNumerator * denominator_gradient) / (mDenominator * mDenominator);
            }

        });
        result->SetExpression(p_result_exp);

    } else {
        result->SetExpression(LiteralExpression<double>::Create(0.0, number_of_elements));
    }

    return result;

    KRATOS_CATCH("");
}

template<class TContainerType>
void BoltzmannOperator<TContainerType>::Update(typename ContainerExpression<TContainerType>::Pointer pContainerExpression)
{
    KRATOS_TRY

    mpContainerExpression = pContainerExpression->Clone();

    const auto& r_expression = mpContainerExpression->GetExpression();
    const IndexType number_of_elements = r_expression.NumberOfEntities();
    const IndexType number_of_components =  r_expression.GetItemComponentCount();

    std::tie(mNumerator, mDenominator) = IndexPartition<IndexType>(number_of_elements).for_each<CombinedReduction<SumReduction<double>, SumReduction<double>>>([&](const auto Index) {
        const IndexType entity_data_begin_index = (Index * number_of_components);

        double numerator{0.0}, denominator{0.0};
        for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
            const double value = r_expression.Evaluate(Index, entity_data_begin_index, i_comp) / mScalingFactor;
            numerator += value * std::exp(mBeta * value);
            denominator += std::exp(mBeta * value);
        }

        return std::make_tuple(numerator, denominator);
    });

    KRATOS_CATCH("");
}

// template instantiations
template class BoltzmannOperator<ModelPart::NodesContainerType>;
template class BoltzmannOperator<ModelPart::ConditionsContainerType>;
template class BoltzmannOperator<ModelPart::ElementsContainerType>;

} // namespace Kratos