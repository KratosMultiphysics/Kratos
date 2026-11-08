//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl
//                   Suneth Warnakulasuriya
//

// System includes
#include <algorithm>

// Project includes
#include "includes/define.h"
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "sigmoidal_projection_utils.h"

namespace SigmoidalValueProjectionUtils{

using IndexType = std::size_t;

bool HasVectorDuplicates(const std::vector<double>& rValues)
{
    std::set<double> unique_values(rValues.begin(), rValues.end());
    return unique_values.size() != rValues.size();
};

void CheckXYVectors(
    const std::vector<double>& rXValues,
    const std::vector<double>& rYValues)
{
    KRATOS_ERROR_IF(rXValues.size() != rYValues.size())
        << "SigmoidalProjectionUtils: rXLimits and rYLimits should have the same size.\n";

    KRATOS_ERROR_IF(rXValues.size() < 2)
        << "SigmoidalProjectionUtils: rXLimits and rYLimits should have at least two entries.\n";

    KRATOS_ERROR_IF(!std::is_sorted(rXValues.begin(), rXValues.end()))
        << "SigmoidalProjectionUtils: rXLimits should be sorted ascending.\n";

    KRATOS_ERROR_IF(!std::is_sorted(rYValues.begin(), rYValues.end()))
        << "SigmoidalProjectionUtils: rYValues should be sorted ascending.\n";

    KRATOS_ERROR_IF(HasVectorDuplicates(rXValues))
        << "SigmoidalProjectionUtils: rXValues have duplications.\n";

    KRATOS_ERROR_IF(HasVectorDuplicates(rYValues))
        << "SigmoidalProjectionUtils: rYValues have duplications.\n";
}

IndexType GetUpperValueRangeIndex(
    const double Value,
    const std::vector<double>& rRanges)
{
    const  IndexType size = rRanges.size();
    IndexType upper_index;
    for (upper_index = 0; upper_index < size; ++upper_index) {
        if (Value < rRanges[upper_index]) {
            break;
        }
    }

    return std::clamp<IndexType>(upper_index, 1UL, size - 1);
}

double ProjectValueForward(
    const double xValue,
    const std::vector<double>& rXLimits,
    const std::vector<double>& rYLimits,
    const double Beta,
    const int PenaltyFactor)
{
    const IndexType upper_index = GetUpperValueRangeIndex(xValue, rXLimits);

    const double x1 = rXLimits[upper_index - 1];
    const double x2 = rXLimits[upper_index];
    const double y1 = rYLimits[upper_index - 1];
    const double y2 = rYLimits[upper_index];

    const double limit = std::log1p(std::numeric_limits<double>::max());
    double pow_val = -2.0 * Beta * (xValue - (x1 + x2) / 2);
    pow_val = std::clamp(pow_val, -limit, limit);
    return (y2 - y1) / std::pow((1.0 + std::exp(pow_val)), PenaltyFactor) + y1;
}

double ProjectValueBackward(
    const double yValue,
    const std::vector<double>& rXLimits,
    const std::vector<double>& rYLimits,
    const double Beta,
    const int PenaltyFactor)
{
    IndexType size = rXLimits.size();
    KRATOS_ERROR_IF(yValue > rYLimits[size - 1] || yValue < rYLimits[0])
        << "SigmoidalProjectionUtils::ProjectValueBackward: yValue "
        << yValue << " is out of the given range " << rYLimits << "\n";

    const IndexType upper_index = GetUpperValueRangeIndex(yValue, rYLimits);

    const double x1 = rXLimits[upper_index - 1];
    const double x2 = rXLimits[upper_index];
    const double y1 = rYLimits[upper_index - 1];
    const double y2 = rYLimits[upper_index];

    if (std::abs(yValue - y1) < std::numeric_limits<double>::epsilon()) {
        return x1;
    } else if (std::abs(yValue - y2) < std::numeric_limits<double>::epsilon()) {
        return x2;
    } else {
        return ((x2 + x1) / 2.0) +
               (1.0 / (-2.0 * Beta)) *
                   std::log(std::pow((y2 - y1) / (yValue - y1), 1.0 / PenaltyFactor) - 1);
    }
}

double ComputeFirstDerivativeAtValue(
    const double xValue,
    const std::vector<double>& rXLimits,
    const std::vector<double>& rYLimits,
    const double Beta,
    const int PenaltyFactor)
{
    const IndexType upper_index = GetUpperValueRangeIndex(xValue, rXLimits);

    const double x1 = rXLimits[upper_index - 1];
    const double x2 = rXLimits[upper_index];
    const double y1 = rYLimits[upper_index - 1];
    const double y2 = rYLimits[upper_index];

    const double limit = std::log1p(std::numeric_limits<double>::max());
    double pow_val = -2.0 * Beta * (xValue - (x1 + x2) / 2);
    pow_val = std::clamp(pow_val, -limit, limit);
    return (y2 - y1) * (1.0 / std::pow(1 + std::exp(pow_val), PenaltyFactor + 1)) *
           PenaltyFactor * 2.0 * Beta * std::exp(pow_val);
}
} // namespace SigmoidalValueProjectionUtils

namespace Kratos
{

TensorAdaptor<double>::Pointer SigmoidalProjectionUtils::ProjectForward(
    const TensorAdaptor<double>& rInputTensorAdaptor,
    const std::vector<double>& rXValues,
    const std::vector<double>& rYValues,
    const double Beta,
    const int PenaltyFactor)
{
    KRATOS_TRY

    SigmoidalValueProjectionUtils::CheckXYVectors(rXValues,rYValues);

    const IndexType size = rInputTensorAdaptor.Size();
    auto input_data_view = rInputTensorAdaptor.ViewData();

    auto p_nd_data = Kratos::make_shared<NDData<double>>(rInputTensorAdaptor.Shape());
    auto p_result_tensor_adaptor = Kratos::make_shared<TensorAdaptor<double>>(rInputTensorAdaptor.GetContainer(), p_nd_data, false);
    auto result_data_view = p_result_tensor_adaptor->ViewData();

    IndexPartition<IndexType>(size).for_each([&input_data_view, &result_data_view, &rXValues, &rYValues, Beta, PenaltyFactor](const IndexType Index) {
        result_data_view[Index] = SigmoidalValueProjectionUtils::ProjectValueForward(input_data_view[Index], rXValues, rYValues, Beta, PenaltyFactor);
    });

    return p_result_tensor_adaptor;

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer SigmoidalProjectionUtils::ProjectBackward(
    const TensorAdaptor<double>& rInputTensorAdaptor,
    const std::vector<double>& rXValues,
    const std::vector<double>& rYValues,
    const double Beta,
    const int PenaltyFactor)
{
    KRATOS_TRY

    SigmoidalValueProjectionUtils::CheckXYVectors(rXValues,rYValues);

    const IndexType size = rInputTensorAdaptor.Size();
    auto input_data_view = rInputTensorAdaptor.ViewData();

    auto p_nd_data = Kratos::make_shared<NDData<double>>(rInputTensorAdaptor.Shape());
    auto p_result_tensor_adaptor = Kratos::make_shared<TensorAdaptor<double>>(rInputTensorAdaptor.GetContainer(), p_nd_data, false);
    auto result_data_view = p_result_tensor_adaptor->ViewData();

    IndexPartition<IndexType>(size).for_each([&input_data_view, &result_data_view, &rXValues, &rYValues, Beta, PenaltyFactor](const IndexType Index) {
        result_data_view[Index] = SigmoidalValueProjectionUtils::ProjectValueBackward(input_data_view[Index], rXValues, rYValues, Beta, PenaltyFactor);
    });

    return p_result_tensor_adaptor;

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer SigmoidalProjectionUtils::CalculateForwardProjectionGradient(
    const TensorAdaptor<double>& rInputTensorAdaptor,
    const std::vector<double>& rXValues,
    const std::vector<double>& rYValues,
    const double Beta,
    const int PenaltyFactor)
{
    KRATOS_TRY

    SigmoidalValueProjectionUtils::CheckXYVectors(rXValues,rYValues);

    const IndexType size = rInputTensorAdaptor.Size();
    auto input_data_view = rInputTensorAdaptor.ViewData();

    auto p_nd_data = Kratos::make_shared<NDData<double>>(rInputTensorAdaptor.Shape());
    auto p_result_tensor_adaptor = Kratos::make_shared<TensorAdaptor<double>>(rInputTensorAdaptor.GetContainer(), p_nd_data, false);
    auto result_data_view = p_result_tensor_adaptor->ViewData();

    IndexPartition<IndexType>(size).for_each([&input_data_view, &result_data_view, &rXValues, &rYValues, Beta, PenaltyFactor](const IndexType Index) {
        result_data_view[Index] = SigmoidalValueProjectionUtils::ComputeFirstDerivativeAtValue(input_data_view[Index], rXValues, rYValues, Beta, PenaltyFactor);
    });

    return p_result_tensor_adaptor;

    KRATOS_CATCH("");
}

}