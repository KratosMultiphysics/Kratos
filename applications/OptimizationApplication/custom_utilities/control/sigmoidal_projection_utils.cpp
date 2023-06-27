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
#include <sstream>
#include <type_traits>

// Project includes
#include "expression/variable_expression_io.h"
#include "expression/variable_expression_data_io.h"
#include "expression/container_data_io.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "optimization_application_variables.h"

// Include base h
#include "sigmoidal_projection_utils.h"

namespace SigmoidalValueProjectionUtils{

bool HasVectorDuplicates(
    std::vector<double> values){
        std::sort(values.begin(), values.end());
        return std::adjacent_find(values.begin(), values.end()) != values.end();
};

void CheckXYVectors(std::vector<double> rXValues, std::vector<double> rYValues){
    if (rXValues.size() != rYValues.size())
        KRATOS_ERROR << "SigmoidalProjectionUtils: xLimits and yLimits should have the same size.\n";

    if (rXValues.size() < 2)
        KRATOS_ERROR << "SigmoidalProjectionUtils: xLimits and yLimits should have at least two entries.\n";

    if (!std::is_sorted(rXValues.begin(), rXValues.end()))
        KRATOS_ERROR << "SigmoidalProjectionUtils: xLimits should be sorted ascending.\n";

    if (!std::is_sorted(rYValues.begin(), rYValues.end()))
        KRATOS_ERROR << "SigmoidalProjectionUtils: rYValues should be sorted ascending.\n";

    if (HasVectorDuplicates(rXValues))
        KRATOS_ERROR << "SigmoidalProjectionUtils: rXValues have duplications.\n";

    if (HasVectorDuplicates(rYValues))
        KRATOS_ERROR << "SigmoidalProjectionUtils: rYValues have duplications.\n";
}

double ProjectValueForward(
        double xValue,
        std::vector<double> xLimits,
        std::vector<double> yLimits,
        double Beta,
        int PenaltyFactor)
{

    size_t size = xLimits.size();
    double x1 = xLimits[0], x2 = xLimits[1], y1 = yLimits[0], y2 = yLimits[1];

    for (size_t i = 0; i < size - 1; ++i) {
        if (xValue >= xLimits[i] && xValue <= xLimits[i + 1]) {
            x1 = xLimits[i];
            x2 = xLimits[i + 1];
            y1 = yLimits[i];
            y2 = yLimits[i + 1];
            break;
        }
    }

    if (xValue > xLimits[size - 1]) {
        x1 = xLimits[size - 2];
        x2 = xLimits[size - 1];
        y1 = yLimits[size - 2];
        y2 = yLimits[size - 1];
    }

    double pow_val = -2.0 * Beta * (xValue - (x1 + x2) / 2)  * PenaltyFactor;
    double limit = std::log1p(std::numeric_limits<double>::max());
    pow_val = std::clamp(pow_val, -limit, limit);
    return (y2 - y1) / std::pow((1.0 + std::exp(pow_val)),PenaltyFactor) + y1;

}

double ProjectValueBackward(
        double yValue,
        std::vector<double> xLimits,
        std::vector<double> yLimits,
        double Beta,
        int PenaltyFactor)
{

    size_t size = xLimits.size();
    if (yValue > yLimits[size - 1] || yValue < yLimits[0])
        KRATOS_ERROR << "SigmoidalProjectionUtils::ProjectValueBackward: yValue "
                     <<yValue<<" is out of the given range "<< yLimits <<"\n";

    if (std::abs(yValue)>std::numeric_limits<double>::max())
        KRATOS_ERROR << "SigmoidalProjectionUtils::ProjectValueBackward: yValue "
                     <<yValue<<" is too big for backward projection !!! "<< yLimits <<"\n";

    double x1 = xLimits[0], x2 = xLimits[1], y1 = yLimits[0], y2 = yLimits[1];
    for (size_t i = 0; i < size - 1; ++i) {
        if (yValue >= yLimits[i] && yValue <= yLimits[i + 1]) {
            x1 = xLimits[i];
            x2 = xLimits[i + 1];
            y1 = yLimits[i];
            y2 = yLimits[i + 1];
            break;
        }
    }

    if (std::abs(yValue-y1)<std::numeric_limits<double>::epsilon())
        return x1;
    else if (std::abs(yValue-y2)<std::numeric_limits<double>::epsilon())
        return x2;
    else
        return ((x2+x1)/2.0) + (1.0/(-2.0*Beta)) * std::log(std::pow((y2-y1)/(yValue-y1),1.0/PenaltyFactor)-1);
}

double ComputeFirstDerivativeAtValue(
        double xValue,
        std::vector<double> xLimits,
        std::vector<double> yLimits,
        double Beta,
        int PenaltyFactor)
{

    size_t size = xLimits.size();
    double x1 = xLimits[0], x2 = xLimits[1], y1 = yLimits[0], y2 = yLimits[1];

    for (size_t i = 0; i < size - 1; ++i) {
        if (xValue >= xLimits[i] && xValue <= xLimits[i + 1]) {
            x1 = xLimits[i];
            x2 = xLimits[i + 1];
            y1 = yLimits[i];
            y2 = yLimits[i + 1];
            break;
        }
    }

    if (xValue > xLimits[size - 1]) {
        x1 = xLimits[size - 2];
        x2 = xLimits[size - 1];
        y1 = yLimits[size - 2];
        y2 = yLimits[size - 1];
    }

    double pow_val = -2.0 * Beta * (xValue - (x1 + x2) / 2)  * PenaltyFactor;
    double limit = std::log1p(std::numeric_limits<double>::max());
    pow_val = std::clamp(pow_val, -limit, limit);
    return (y2-y1) * (1.0/std::pow(1+std::exp(pow_val),PenaltyFactor+1)) * PenaltyFactor * 2.0 * Beta * std::exp(pow_val);

}
}

namespace Kratos
{

///@name Kratos Classes
///@{



template<class TContainerType>
ContainerExpression<TContainerType> SigmoidalProjectionUtils::ProjectForward(
    ContainerExpression<TContainerType>& rInputExpression,
    const std::vector<double>& rXValues,
    const std::vector<double>& rYValues,
    const double Beta,
    const int PenaltyFactor)
{
    KRATOS_TRY

    SigmoidalValueProjectionUtils::CheckXYVectors(rXValues,rYValues);

    const auto& r_input_expression = rInputExpression.GetExpression();
    const IndexType local_size_1 = rInputExpression.GetItemComponentCount();
    const IndexType number_of_entities_1 = rInputExpression.GetContainer().size();

    ContainerExpression<TContainerType> output_container(rInputExpression.GetModelPart());
    auto p_flat_data_expression = LiteralFlatExpression<double>::Create(number_of_entities_1, rInputExpression.GetItemShape());
    output_container.SetExpression(p_flat_data_expression);
    auto& r_output_expression = *p_flat_data_expression;


    IndexPartition<IndexType>(number_of_entities_1).for_each([&r_input_expression, &r_output_expression, &rXValues, &rYValues, &Beta, &PenaltyFactor, local_size_1](const IndexType EntityIndex) {
            const IndexType local_data_begin_index = EntityIndex * local_size_1;
            for (IndexType i = 0; i < local_size_1; ++i) {
                double input_value = r_input_expression.Evaluate(EntityIndex, local_data_begin_index, i);
                double projected_value = SigmoidalValueProjectionUtils::ProjectValueForward(input_value,rXValues,rYValues,Beta,PenaltyFactor);
                r_output_expression.SetData(local_data_begin_index, i, projected_value);
            }
        });

    return output_container;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> SigmoidalProjectionUtils::ProjectBackward(
    ContainerExpression<TContainerType>& rInputExpression,
    const std::vector<double>& rXValues,
    const std::vector<double>& rYValues,
    const double Beta,
    const int PenaltyFactor)
{
    KRATOS_TRY

    SigmoidalValueProjectionUtils::CheckXYVectors(rXValues,rYValues);

    const auto& r_input_expression = rInputExpression.GetExpression();
    const IndexType local_size_1 = rInputExpression.GetItemComponentCount();
    const IndexType number_of_entities_1 = rInputExpression.GetContainer().size();

    ContainerExpression<TContainerType> output_container(rInputExpression.GetModelPart());
    auto p_flat_data_expression = LiteralFlatExpression<double>::Create(number_of_entities_1, rInputExpression.GetItemShape());
    output_container.SetExpression(p_flat_data_expression);
    auto& r_output_expression = *p_flat_data_expression;


    IndexPartition<IndexType>(number_of_entities_1).for_each([&r_input_expression, &r_output_expression, &rXValues, &rYValues, &Beta, &PenaltyFactor, local_size_1](const IndexType EntityIndex) {
            const IndexType local_data_begin_index = EntityIndex * local_size_1;
            for (IndexType i = 0; i < local_size_1; ++i) {
                double input_value = r_input_expression.Evaluate(EntityIndex, local_data_begin_index, i);
                double projected_value = SigmoidalValueProjectionUtils::ProjectValueBackward(input_value,rXValues,rYValues,Beta,PenaltyFactor);
                r_output_expression.SetData(local_data_begin_index, i, projected_value);
            }
        });

    return output_container;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> SigmoidalProjectionUtils::ComputeFirstDerivative(
    ContainerExpression<TContainerType>& rInputExpression,
    const std::vector<double>& rXValues,
    const std::vector<double>& rYValues,
    const double Beta,
    const int PenaltyFactor)
{
    KRATOS_TRY

    SigmoidalValueProjectionUtils::CheckXYVectors(rXValues,rYValues);

    const auto& r_input_expression = rInputExpression.GetExpression();
    const IndexType local_size_1 = rInputExpression.GetItemComponentCount();
    const IndexType number_of_entities_1 = rInputExpression.GetContainer().size();

    ContainerExpression<TContainerType> output_container(rInputExpression.GetModelPart());
    auto p_flat_data_expression = LiteralFlatExpression<double>::Create(number_of_entities_1, rInputExpression.GetItemShape());
    output_container.SetExpression(p_flat_data_expression);
    auto& r_output_expression = *p_flat_data_expression;


    IndexPartition<IndexType>(number_of_entities_1).for_each([&r_input_expression, &r_output_expression, &rXValues, &rYValues, &Beta, &PenaltyFactor, local_size_1](const IndexType EntityIndex) {
            const IndexType local_data_begin_index = EntityIndex * local_size_1;
            for (IndexType i = 0; i < local_size_1; ++i) {
                double input_value = r_input_expression.Evaluate(EntityIndex, local_data_begin_index, i);
                double derivative_value = SigmoidalValueProjectionUtils::ComputeFirstDerivativeAtValue(input_value,rXValues,rYValues,Beta,PenaltyFactor);
                r_output_expression.SetData(local_data_begin_index, i, derivative_value);
            }
        });

    return output_container;

    KRATOS_CATCH("");
}

template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::NodesContainerType> SigmoidalProjectionUtils::ProjectForward(ContainerExpression<ModelPart::NodesContainerType>&, const std::vector<double>&, const std::vector<double>&, const double, const int);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::ConditionsContainerType> SigmoidalProjectionUtils::ProjectForward(ContainerExpression<ModelPart::ConditionsContainerType>&, const std::vector<double>&, const std::vector<double>&, const double, const int);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::ElementsContainerType> SigmoidalProjectionUtils::ProjectForward(ContainerExpression<ModelPart::ElementsContainerType>&, const std::vector<double>&, const std::vector<double>&, const double, const int);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::NodesContainerType> SigmoidalProjectionUtils::ProjectBackward(ContainerExpression<ModelPart::NodesContainerType>&, const std::vector<double>&, const std::vector<double>&, const double, const int);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::ConditionsContainerType> SigmoidalProjectionUtils::ProjectBackward(ContainerExpression<ModelPart::ConditionsContainerType>&, const std::vector<double>&, const std::vector<double>&, const double, const int);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::ElementsContainerType> SigmoidalProjectionUtils::ProjectBackward(ContainerExpression<ModelPart::ElementsContainerType>&, const std::vector<double>&, const std::vector<double>&, const double, const int);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::NodesContainerType> SigmoidalProjectionUtils::ComputeFirstDerivative(ContainerExpression<ModelPart::NodesContainerType>&, const std::vector<double>&, const std::vector<double>&, const double, const int);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::ConditionsContainerType> SigmoidalProjectionUtils::ComputeFirstDerivative(ContainerExpression<ModelPart::ConditionsContainerType>&, const std::vector<double>&, const std::vector<double>&, const double, const int);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::ElementsContainerType> SigmoidalProjectionUtils::ComputeFirstDerivative(ContainerExpression<ModelPart::ElementsContainerType>&, const std::vector<double>&, const std::vector<double>&, const double, const int);
///@}
}