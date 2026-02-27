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
#include <ranges>

// External includes

// Project includes

// Application includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Include base h
#include "boltzmann_operator.h"

namespace Kratos {

BoltzmannOperator::BoltzmannOperator(const double Beta)
    : mBeta(Beta)
{
}

double BoltzmannOperator::CalculateValue() const
{
    if (mDenominator > std::numeric_limits<double>::epsilon()) {
        return mNumerator / mDenominator;
    } else {
        return 0.0;
    }
}

TensorAdaptor<double>::Pointer BoltzmannOperator::CalculateGradient() const
{
    KRATOS_TRY

    auto p_result = mpTensorAdaptor->Clone();
    auto output_data_view = p_result->ViewData();
    const auto number_of_entities = p_result->Size();

    if (mDenominator > std::numeric_limits<double>::epsilon()) {
        const double denominator_square = mDenominator * mDenominator;

        IndexPartition<IndexType>(number_of_entities).for_each([this, &output_data_view, denominator_square](const auto Index) {
            const double value = output_data_view[Index];

            const double exp_beta_x_k = std::exp(this->mBeta * value - this->mExtremeValue);
            const double numerator_gradient = exp_beta_x_k * (1 + this->mBeta * value);
            const double denominator_gradient = this->mBeta * exp_beta_x_k;

            output_data_view[Index] = (this->mDenominator * numerator_gradient - this->mNumerator * denominator_gradient) / denominator_square;

        });
    } else {
        std::fill(output_data_view.begin(), output_data_view.end(), 0.0);
    }

    return p_result;

    KRATOS_CATCH("");
}

void BoltzmannOperator::Update(const TensorAdaptor<double>& rTensorAdaptor)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rTensorAdaptor.Shape().size() == 1)
        << "Boltzmann operator can only be used with scalar tensor adaptors [ tensor adaptor = "
        << rTensorAdaptor << " ].\n";

    const auto number_of_entities = rTensorAdaptor.Size();
    const auto data_view = rTensorAdaptor.ViewData();

    // this is identified to avoid std::exp giving inf for intermediate values of the summation.
    mExtremeValue = (mBeta > 0 ? std::ranges::max(data_view) : std::ranges::min(data_view)) * this->mBeta;

    std::tie(mNumerator, mDenominator) = IndexPartition<IndexType>(number_of_entities).for_each<CombinedReduction<SumReduction<double>, SumReduction<double>>>([this, &data_view](const auto Index) {
        const auto value = data_view[Index];
        const double temp = std::exp(this->mBeta * value - this->mExtremeValue);
        return std::make_tuple(value * temp, temp);
    });

    mpTensorAdaptor = rTensorAdaptor.Clone();

    KRATOS_CATCH("");
}

} // namespace Kratos