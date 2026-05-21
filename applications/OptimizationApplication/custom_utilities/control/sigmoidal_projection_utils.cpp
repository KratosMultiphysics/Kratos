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
#include "expression/literal_flat_expression.h"
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

template<class TContainerType>
ContainerExpression<TContainerType> SigmoidalProjectionUtils::ProjectForward(
    const ContainerExpression<TContainerType>& rInputExpression,
    const std::vector<double>& rXValues,
    const std::vector<double>& rYValues,
    const double Beta,
    const int PenaltyFactor)
{
    KRATOS_TRY

    SigmoidalValueProjectionUtils::CheckXYVectors(rXValues,rYValues);

    const auto& r_input_expression = rInputExpression.GetExpression();
    const IndexType local_size = rInputExpression.GetItemComponentCount();
    const IndexType number_of_entities = rInputExpression.GetContainer().size();

    ContainerExpression<TContainerType> output_container(*rInputExpression.pGetModelPart());
    auto p_flat_data_expression = LiteralFlatExpression<double>::Create(number_of_entities, rInputExpression.GetItemShape());
    output_container.SetExpression(p_flat_data_expression);
    auto& r_output_expression = *p_flat_data_expression;

    IndexPartition<IndexType>(number_of_entities).for_each([&r_input_expression, &r_output_expression, &rXValues, &rYValues, Beta, PenaltyFactor, local_size](const IndexType EntityIndex) {
            const IndexType local_data_begin_index = EntityIndex * local_size;
            for (IndexType i = 0; i < local_size; ++i) {
                const double input_value = r_input_expression.Evaluate(EntityIndex, local_data_begin_index, i);
                const double projected_value = SigmoidalValueProjectionUtils::ProjectValueForward(input_value, rXValues, rYValues, Beta, PenaltyFactor);
                r_output_expression.SetData(local_data_begin_index, i, projected_value);
            }
        });

    return output_container;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> SigmoidalProjectionUtils::ProjectBackward(
    const ContainerExpression<TContainerType>& rInputExpression,
    const std::vector<double>& rXValues,
    const std::vector<double>& rYValues,
    const double Beta,
    const int PenaltyFactor)
{
    KRATOS_TRY

    SigmoidalValueProjectionUtils::CheckXYVectors(rXValues,rYValues);

    const auto& r_input_expression = rInputExpression.GetExpression();
    const IndexType local_size = rInputExpression.GetItemComponentCount();
    const IndexType number_of_entities = rInputExpression.GetContainer().size();

    ContainerExpression<TContainerType> output_container(*rInputExpression.pGetModelPart());
    auto p_flat_data_expression = LiteralFlatExpression<double>::Create(number_of_entities, rInputExpression.GetItemShape());
    output_container.SetExpression(p_flat_data_expression);
    auto& r_output_expression = *p_flat_data_expression;

    IndexPartition<IndexType>(number_of_entities).for_each([&r_input_expression, &r_output_expression, &rXValues, &rYValues, Beta, PenaltyFactor, local_size](const IndexType EntityIndex) {
            const IndexType local_data_begin_index = EntityIndex * local_size;
            for (IndexType i = 0; i < local_size; ++i) {
                const double input_value = r_input_expression.Evaluate(EntityIndex, local_data_begin_index, i);
                const double projected_value = SigmoidalValueProjectionUtils::ProjectValueBackward(input_value, rXValues, rYValues, Beta, PenaltyFactor);
                r_output_expression.SetData(local_data_begin_index, i, projected_value);
            }
        });

    return output_container;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> SigmoidalProjectionUtils::CalculateForwardProjectionGradient(
    const ContainerExpression<TContainerType>& rInputExpression,
    const std::vector<double>& rXValues,
    const std::vector<double>& rYValues,
    const double Beta,
    const int PenaltyFactor)
{
    KRATOS_TRY

    SigmoidalValueProjectionUtils::CheckXYVectors(rXValues,rYValues);

    const auto& r_input_expression = rInputExpression.GetExpression();
    const IndexType local_size = rInputExpression.GetItemComponentCount();
    const IndexType number_of_entities = rInputExpression.GetContainer().size();

    ContainerExpression<TContainerType> output_container(*rInputExpression.pGetModelPart());
    auto p_flat_data_expression = LiteralFlatExpression<double>::Create(number_of_entities, rInputExpression.GetItemShape());
    output_container.SetExpression(p_flat_data_expression);
    auto& r_output_expression = *p_flat_data_expression;

    IndexPartition<IndexType>(number_of_entities).for_each([&r_input_expression, &r_output_expression, &rXValues, &rYValues, Beta, PenaltyFactor, local_size](const IndexType EntityIndex) {
            const IndexType local_data_begin_index = EntityIndex * local_size;
            for (IndexType i = 0; i < local_size; ++i) {
                const double input_value = r_input_expression.Evaluate(EntityIndex, local_data_begin_index, i);
                const double derivative_value = SigmoidalValueProjectionUtils::ComputeFirstDerivativeAtValue(input_value, rXValues, rYValues, Beta, PenaltyFactor);
                r_output_expression.SetData(local_data_begin_index, i, derivative_value);
            }
        });

    return output_container;

    KRATOS_CATCH("");
}

// template instantiations
#define KRATOS_INSTANTIATE_SIGMOIDAL_PROJECTION_UTIL_METHODS(CONTAINER_TYPE)                                                            \
    template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<CONTAINER_TYPE> SigmoidalProjectionUtils::ProjectForward(         \
        const ContainerExpression<CONTAINER_TYPE>&, const std::vector<double>&,                                                         \
        const std::vector<double>&, const double, const int);                                                                           \
    template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<CONTAINER_TYPE> SigmoidalProjectionUtils::ProjectBackward(        \
        const ContainerExpression<CONTAINER_TYPE>&, const std::vector<double>&,                                                         \
        const std::vector<double>&, const double, const int);                                                                           \
    template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<CONTAINER_TYPE> SigmoidalProjectionUtils::CalculateForwardProjectionGradient( \
        const ContainerExpression<CONTAINER_TYPE>&, const std::vector<double>&,                                                         \
        const std::vector<double>&, const double, const int);

KRATOS_INSTANTIATE_SIGMOIDAL_PROJECTION_UTIL_METHODS(ModelPart::NodesContainerType)
KRATOS_INSTANTIATE_SIGMOIDAL_PROJECTION_UTIL_METHODS(ModelPart::ConditionsContainerType)
KRATOS_INSTANTIATE_SIGMOIDAL_PROJECTION_UTIL_METHODS(ModelPart::ElementsContainerType)

#undef KRATOS_INSTANTIATE_SIGMOIDAL_PROJECTION_UTIL_METHODS
///@}
}