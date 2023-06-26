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

namespace Kratos
{

///@name Kratos Classes
///@{

bool HasVectorDuplicates(
    std::vector<double> values){
        std::sort(values.begin(), values.end());
        return std::adjacent_find(values.begin(), values.end()) != values.end();
};

void CheckXYVectors(std::vector<double> xValues, std::vector<double> yValues){
    if (xValues.size() != yValues.size())
        KRATOS_ERROR << "SigmoidalProjectionUtils: xLimits and yLimits should have the same size.\n";

    if (xValues.size() < 2)
        KRATOS_ERROR << "SigmoidalProjectionUtils: xLimits and yLimits should have at least two entries.\n";

    if (!std::is_sorted(xValues.begin(), xValues.end()))
        KRATOS_ERROR << "SigmoidalProjectionUtils: xLimits should be sorted ascending.\n";

    if (!std::is_sorted(yValues.begin(), yValues.end()))
        KRATOS_ERROR << "SigmoidalProjectionUtils: yValues should be sorted ascending.\n";

    if (HasVectorDuplicates(xValues))
        KRATOS_ERROR << "SigmoidalProjectionUtils: xValues have duplications.\n";

    if (HasVectorDuplicates(yValues))
        KRATOS_ERROR << "SigmoidalProjectionUtils: yValues have duplications.\n";
}

double ProjectValueForward(
        double xValue,
        std::vector<double> xLimits,
        std::vector<double> yLimits,
        double beta,
        int penalFac)
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

    double pow_val = -2.0 * beta * (xValue - (x1 + x2) / 2)  * penalFac;
    double limit = std::log1p(std::numeric_limits<double>::max());
    pow_val = std::clamp(pow_val, -limit, limit);
    return (y2 - y1) / std::pow((1.0 + std::exp(pow_val)),penalFac) + y1;

}

double ProjectValueBackward(
        double yValue,
        std::vector<double> xLimits,
        std::vector<double> yLimits,
        double beta,
        int penalFac)
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
        return ((x2+x1)/2.0) + (1.0/(-2.0*beta)) * std::log(std::pow((y2-y1)/(yValue-y1),1.0/penalFac)-1);
}

double ComputeFirstDerivativeAtValue(
        double xValue,
        std::vector<double> xLimits,
        std::vector<double> yLimits,
        double beta,
        int penalFac)
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

    double pow_val = -2.0 * beta * (xValue - (x1 + x2) / 2)  * penalFac;
    double limit = std::log1p(std::numeric_limits<double>::max());
    pow_val = std::clamp(pow_val, -limit, limit);
    return (y2-y1) * (1.0/std::pow(1+std::exp(pow_val),penalFac+1)) * penalFac * 2.0 * beta * std::exp(pow_val);

}

template<class TContainerType>
ContainerExpression<TContainerType> SigmoidalProjectionUtils::ProjectForward(
    ContainerExpression<TContainerType>& rInputExpression,
    std::vector<double> xValues,
    std::vector<double> yValues,
    double beta,
    int penalFac)
{
    KRATOS_TRY

    CheckXYVectors(xValues,yValues);

    const auto& r_input_expression = rInputExpression.GetExpression();
    const IndexType local_size_1 = rInputExpression.GetItemComponentCount();
    const IndexType number_of_entities_1 = rInputExpression.GetContainer().size();

    ContainerExpression<TContainerType> output_container(rInputExpression.GetModelPart());
    auto p_flat_data_expression = LiteralFlatExpression<double>::Create(number_of_entities_1, rInputExpression.GetItemShape());
    output_container.SetExpression(p_flat_data_expression);
    auto& r_output_expression = *p_flat_data_expression;


    IndexPartition<IndexType>(number_of_entities_1).for_each([&r_input_expression, &r_output_expression, &xValues, &yValues, &beta, &penalFac, local_size_1](const IndexType EntityIndex) {
            const IndexType local_data_begin_index = EntityIndex * local_size_1;
            for (IndexType i = 0; i < local_size_1; ++i) {
                double input_value = r_input_expression.Evaluate(EntityIndex, local_data_begin_index, i);
                double projected_value = ProjectValueForward(input_value,xValues,yValues,beta,penalFac);
                r_output_expression.SetData(local_data_begin_index, i, projected_value);
            }
        });

    return output_container;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> SigmoidalProjectionUtils::ProjectBackward(
    ContainerExpression<TContainerType>& rInputExpression,
    std::vector<double> xValues,
    std::vector<double> yValues,
    double beta,
    int penalFac)
{
    KRATOS_TRY

    CheckXYVectors(xValues,yValues);

    const auto& r_input_expression = rInputExpression.GetExpression();
    const IndexType local_size_1 = rInputExpression.GetItemComponentCount();
    const IndexType number_of_entities_1 = rInputExpression.GetContainer().size();

    ContainerExpression<TContainerType> output_container(rInputExpression.GetModelPart());
    auto p_flat_data_expression = LiteralFlatExpression<double>::Create(number_of_entities_1, rInputExpression.GetItemShape());
    output_container.SetExpression(p_flat_data_expression);
    auto& r_output_expression = *p_flat_data_expression;


    IndexPartition<IndexType>(number_of_entities_1).for_each([&r_input_expression, &r_output_expression, &xValues, &yValues, &beta, &penalFac, local_size_1](const IndexType EntityIndex) {
            const IndexType local_data_begin_index = EntityIndex * local_size_1;
            for (IndexType i = 0; i < local_size_1; ++i) {
                double input_value = r_input_expression.Evaluate(EntityIndex, local_data_begin_index, i);
                double projected_value = ProjectValueBackward(input_value,xValues,yValues,beta,penalFac);
                r_output_expression.SetData(local_data_begin_index, i, projected_value);
            }
        });

    return output_container;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> SigmoidalProjectionUtils::ComputeFirstDerivative(
    ContainerExpression<TContainerType>& rInputExpression,
    std::vector<double> xValues,
    std::vector<double> yValues,
    double beta,
    int penalFac)
{
    KRATOS_TRY

    CheckXYVectors(xValues,yValues);

    const auto& r_input_expression = rInputExpression.GetExpression();
    const IndexType local_size_1 = rInputExpression.GetItemComponentCount();
    const IndexType number_of_entities_1 = rInputExpression.GetContainer().size();

    ContainerExpression<TContainerType> output_container(rInputExpression.GetModelPart());
    auto p_flat_data_expression = LiteralFlatExpression<double>::Create(number_of_entities_1, rInputExpression.GetItemShape());
    output_container.SetExpression(p_flat_data_expression);
    auto& r_output_expression = *p_flat_data_expression;


    IndexPartition<IndexType>(number_of_entities_1).for_each([&r_input_expression, &r_output_expression, &xValues, &yValues, &beta, &penalFac, local_size_1](const IndexType EntityIndex) {
            const IndexType local_data_begin_index = EntityIndex * local_size_1;
            for (IndexType i = 0; i < local_size_1; ++i) {
                double input_value = r_input_expression.Evaluate(EntityIndex, local_data_begin_index, i);
                double derivative_value = ComputeFirstDerivativeAtValue(input_value,xValues,yValues,beta,penalFac);
                r_output_expression.SetData(local_data_begin_index, i, derivative_value);
            }
        });

    return output_container;

    KRATOS_CATCH("");
}

template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::NodesContainerType> SigmoidalProjectionUtils::ProjectForward(ContainerExpression<ModelPart::NodesContainerType>&, std::vector<double>, std::vector<double>, double, int);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::ConditionsContainerType> SigmoidalProjectionUtils::ProjectForward(ContainerExpression<ModelPart::ConditionsContainerType>&, std::vector<double>, std::vector<double>, double, int);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::ElementsContainerType> SigmoidalProjectionUtils::ProjectForward(ContainerExpression<ModelPart::ElementsContainerType>&, std::vector<double>, std::vector<double>, double, int);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::NodesContainerType> SigmoidalProjectionUtils::ProjectBackward(ContainerExpression<ModelPart::NodesContainerType>&, std::vector<double>, std::vector<double>, double, int);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::ConditionsContainerType> SigmoidalProjectionUtils::ProjectBackward(ContainerExpression<ModelPart::ConditionsContainerType>&, std::vector<double>, std::vector<double>, double, int);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::ElementsContainerType> SigmoidalProjectionUtils::ProjectBackward(ContainerExpression<ModelPart::ElementsContainerType>&, std::vector<double>, std::vector<double>, double, int);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::NodesContainerType> SigmoidalProjectionUtils::ComputeFirstDerivative(ContainerExpression<ModelPart::NodesContainerType>&, std::vector<double>, std::vector<double>, double, int);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::ConditionsContainerType> SigmoidalProjectionUtils::ComputeFirstDerivative(ContainerExpression<ModelPart::ConditionsContainerType>&, std::vector<double>, std::vector<double>, double, int);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::ElementsContainerType> SigmoidalProjectionUtils::ComputeFirstDerivative(ContainerExpression<ModelPart::ElementsContainerType>&, std::vector<double>, std::vector<double>, double, int);
///@}
}