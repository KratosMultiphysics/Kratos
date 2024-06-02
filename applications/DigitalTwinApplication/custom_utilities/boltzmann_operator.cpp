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
#include <limits>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes

// Include base h
#include "boltzmann_operator.h"

namespace Kratos {

BoltzmannOperator::BoltzmannOperator(const double Beta)
    : mBeta(Beta)
{
}

void BoltzmannOperator::CalculateCoefficients()
{
    KRATOS_TRY

    const IndexType n = mValues.size();

    std::tie(mNumerator, mDenominator) = IndexPartition<IndexType>(n).for_each<CombinedReduction<SumReduction<double>, SumReduction<double>>>([&](const auto Index) {
        const double value = mValues[Index];
        const double exp_value = std::exp(mBeta * value);
        return std::make_tuple(value * exp_value, exp_value);
    });

    KRATOS_CATCH("");
}

Vector& BoltzmannOperator::GetData()
{
    return mValues;
}

const Vector& BoltzmannOperator::GetData() const
{
    return mValues;
}

double BoltzmannOperator::GetValue() const
{
    if (mDenominator > std::numeric_limits<double>::epsilon()) {
        return mNumerator / mDenominator;
    } else {
        return 0.0;
    }
}

double BoltzmannOperator::GetGradient(const Vector& rGradients) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rGradients.size() == mValues.size())
        << "Gradient and values size mismatch [ "
        << "number of values = " << mValues.size()
        << ", number of gradients = " << rGradients.size() << " ].";

    const IndexType n = rGradients.size();

    double numerator_gradient{0.0}, denominator_gradient{0.0};

    for (IndexType i = 0; i < n; ++i) {
        const double value = mValues[i];
        const double exp_value = std::exp(mBeta * value);
        const double value_gradient = rGradients[i];

        numerator_gradient += exp_value * value_gradient * (1.0 + mBeta * value);
        denominator_gradient += mBeta * value_gradient * exp_value;
    }

    if (mDenominator > std::numeric_limits<double>::epsilon()) {
        return numerator_gradient / mDenominator - mNumerator * denominator_gradient / (mDenominator * mDenominator);
    } else {
        return 0.0;
    }

    KRATOS_CATCH("");
}

///@} // Kratos Classes

} /* namespace Kratos.*/