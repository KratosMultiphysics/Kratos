//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes
#include <cmath>

// Project includes
#include "expression/binary_expression.h"
#include "expression/literal_expression.h"
#include "expression/literal_flat_expression.h"
#include "expression/unary_abs_expression.h"
#include "expression/unary_combine_expression.h"
#include "expression/unary_statistics_expression.h"
#include "expression/unary_slice_expression.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Include base h
#include "expression_utils.h"

namespace Kratos
{

Expression::ConstPointer ExpressionUtils::Collapse(const Expression::ConstPointer& rpExpression)
{   KRATOS_TRY

    if (rpExpression->GetMaxDepth() == 1) {
        // no need to collapse.
        return rpExpression;
    } else {
        auto p_expression = LiteralFlatExpression<double>::Create(rpExpression->NumberOfEntities(), rpExpression->GetItemShape());
        const auto number_of_components = rpExpression->GetItemComponentCount();

        IndexPartition<IndexType>(rpExpression->NumberOfEntities()).for_each([&](const auto Index) {
            const auto data_begin_index = Index * number_of_components;
            for (IndexType i = 0; i < number_of_components; ++i) {
                p_expression->SetData(data_begin_index, i, rpExpression->Evaluate(Index, data_begin_index, i));
            }
        });

        return p_expression;
    }

    KRATOS_CATCH("");
}

Expression::ConstPointer ExpressionUtils::Slice(
    const Expression::ConstPointer& rpExpression,
    std::size_t Offset,
    std::size_t Stride)
{
    return UnarySliceExpression::Create(rpExpression, Offset, Stride);
}

Expression::ConstPointer ExpressionUtils::Reshape(
    const Expression::ConstPointer& rpExpression,
    const std::vector<IndexType>& rNewShape)
{
    return Reshape(rpExpression, rNewShape.begin(), rNewShape.end());
}

Expression::ConstPointer ExpressionUtils::Comb(const std::vector<Expression::ConstPointer>& rExpressions)
{
    return Comb(rExpressions.begin(), rExpressions.end());
}

Expression::ConstPointer ExpressionUtils::Abs(const Expression::ConstPointer& rpExpression)
{
    return UnaryAbsExpression::Create(rpExpression);
}

Expression::ConstPointer ExpressionUtils::Pow(
    const Expression::ConstPointer& rpExpression,
    const double Power)
{
    if (Power == 1) {
        return rpExpression;
    } else if (Power == 0) {
        return LiteralExpression<double>::Create(1, rpExpression->NumberOfEntities());
    } else {
        return BinaryExpression<BinaryOperations::Power>::Create(rpExpression, LiteralExpression<double>::Create(Power, rpExpression->NumberOfEntities()));
    }
}

Expression::ConstPointer ExpressionUtils::Pow(
    const Expression::ConstPointer& rpExpression,
    const Expression::ConstPointer& rpPowerExpression)
{
    return BinaryExpression<BinaryOperations::Power>::Create(rpExpression, rpPowerExpression);
}

Expression::ConstPointer ExpressionUtils::Scale(
    const Expression::ConstPointer& rpExpression,
    const double Scale)
{
    if (Scale == 1.0) {
        return rpExpression;
    } else if (Scale == 0.0) {
        return LiteralExpression<double>::Create(0, rpExpression->NumberOfEntities());
    } else {
        return BinaryExpression<BinaryOperations::Multiplication>::Create(rpExpression, LiteralExpression<double>::Create(Scale, rpExpression->NumberOfEntities()));
    }
}

Expression::ConstPointer ExpressionUtils::Scale(
    const Expression::ConstPointer& rpExpression,
    const Expression::ConstPointer& rpScaleExpression)
{
    return BinaryExpression<BinaryOperations::Multiplication>::Create(rpExpression, rpScaleExpression);
}

Expression::ConstPointer ExpressionUtils::EntityMin(const Expression::ConstPointer& rpExpression)
{
    if (rpExpression->GetItemShape().empty()) {
        return rpExpression;
    } else {
        return UnaryStatisticsExpression<UnaryStatisticOperations::Min>::Create(rpExpression);
    }
}

Expression::ConstPointer ExpressionUtils::EntityMax(const Expression::ConstPointer& rpExpression)
{
    if (rpExpression->GetItemShape().empty()) {
        return rpExpression;
    } else {
        return UnaryStatisticsExpression<UnaryStatisticOperations::Max>::Create(rpExpression);
    }
}

Expression::ConstPointer ExpressionUtils::EntitySum(const Expression::ConstPointer& rpExpression)
{
    if (rpExpression->GetItemShape().empty()) {
        return rpExpression;
    } else {
        return UnaryStatisticsExpression<UnaryStatisticOperations::Sum>::Create(rpExpression);
    }
}

double ExpressionUtils::Sum(
    const Expression::ConstPointer& rpExpression,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    const auto number_of_components = rpExpression->GetItemComponentCount();

    return rDataCommunicator.SumAll(IndexPartition<IndexType>(rpExpression->NumberOfEntities()).for_each<SumReduction<double>>([&rpExpression, number_of_components](const auto Index) {
        const auto data_begin_index = number_of_components * Index;
        double value = 0.0;
        for (IndexType i = 0; i < number_of_components; ++i) {
            value += rpExpression->Evaluate(Index, data_begin_index, i);
        }
        return value;
    }));

    KRATOS_CATCH("");
}

double ExpressionUtils::NormInf(
    const Expression::ConstPointer& rpExpression,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    const auto number_of_components = rpExpression->GetItemComponentCount();

    return rDataCommunicator.MaxAll(IndexPartition<IndexType>(rpExpression->NumberOfEntities()).for_each<MaxReduction<double>>([&rpExpression, number_of_components](const auto Index) {
        const auto data_begin_index = number_of_components * Index;
        double value = 0.0;
        for (IndexType i = 0; i < number_of_components; ++i) {
            value = std::max(value, std::abs(rpExpression->Evaluate(Index, data_begin_index, i)));
        }
        return value;
    }));

    KRATOS_CATCH("");
}

double ExpressionUtils::NormL2(
    const Expression::ConstPointer& rpExpression,
    const DataCommunicator& rDataCommunicator)
{
    return NormP(rpExpression, 2, rDataCommunicator);
}

double ExpressionUtils::NormP(
    const Expression::ConstPointer& rpExpression,
    const double P,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    const auto number_of_components = rpExpression->GetItemComponentCount();

    const auto local_sum = rDataCommunicator.SumAll(IndexPartition<IndexType>(rpExpression->NumberOfEntities()).for_each<SumReduction<double>>([&rpExpression, number_of_components, P](const auto Index) {
        const auto data_begin_index = number_of_components * Index;
        double value = 0.0;
        for (IndexType i = 0; i < number_of_components; ++i) {
            value += std::pow(std::abs(rpExpression->Evaluate(Index, data_begin_index, i)), P);
        }
        return value;
    }));

    return std::pow(local_sum, 1 / P);

    KRATOS_CATCH("");
}

double ExpressionUtils::InnerProduct(
    const Expression::ConstPointer& rpExpression1,
    const Expression::ConstPointer& rpExpression2,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    const auto number_of_entities = rpExpression1->NumberOfEntities();
    const auto number_of_components = rpExpression1->GetItemComponentCount();

    KRATOS_ERROR_IF_NOT(number_of_components == rpExpression2->GetItemComponentCount())
        << "Flattened shape size mismatch [ rpExpression1 shape = "
        << rpExpression1->GetItemShape() << ", rpExpression1 flattend size = "
        << rpExpression1->GetItemComponentCount() << ", "
        << "rpExpression2 shape = " << rpExpression2->GetItemShape()
        << ", rpExpression2 flattend size = " << rpExpression2->GetItemComponentCount() << " ].\n";

    KRATOS_ERROR_IF_NOT(number_of_entities == rpExpression2->NumberOfEntities())
        << "Number of entities mismatch [ rpExpression1 # = " << rpExpression1->NumberOfEntities()
        << ", rpExpression2 # = " << rpExpression2->NumberOfEntities() << " ].\n";

    return rDataCommunicator.SumAll(IndexPartition<IndexType>(number_of_entities).for_each<SumReduction<double>>([&rpExpression1, &rpExpression2, number_of_components](const auto Index) {
        const auto data_begin_index = Index * number_of_components;
        double value = 0.0;
        for (IndexType i = 0; i < number_of_components; ++i) {
            value += rpExpression1->Evaluate(Index, data_begin_index, i) * rpExpression2->Evaluate(Index, data_begin_index, i);
        }
        return value;
    }));

    KRATOS_CATCH("");
}

} // namespace Kratos
