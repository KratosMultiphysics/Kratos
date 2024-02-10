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
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"
#include "expression/literal_flat_expression.h"

// Application includes

// Include base h
#include "optimization_utils.h"

namespace Kratos
{

template<class TContainerType>
GeometryData::KratosGeometryType OptimizationUtils::GetContainerEntityGeometryType(
    const TContainerType& rContainer,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    if (rContainer.size() > 0) {
        const auto first_geometry_type = rContainer.begin()->GetGeometry().GetGeometryType();
        const bool local_value = block_for_each<MinReduction<bool>>(rContainer, [&](const auto& rEntity) -> IndexType {
            return first_geometry_type == rEntity.GetGeometry().GetGeometryType();
        });

        if (rDataCommunicator.AndReduceAll(local_value)) {
            return first_geometry_type;
        } else {
            return GeometryData::KratosGeometryType::Kratos_generic_type;
        }
    } else {
        return GeometryData::KratosGeometryType::Kratos_generic_type;
    }

    KRATOS_CATCH("")
}

template<class TContainerType, class TDataType>
bool OptimizationUtils::IsVariableExistsInAllContainerProperties(
    const TContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    const bool local_value = block_for_each<MinReduction<bool>>(rContainer, [&](const auto& rEntity){
        return rEntity.GetProperties().Has(rVariable);
    });

    return rDataCommunicator.AndReduceAll(local_value);

    KRATOS_CATCH("");
}

template<class TContainerType, class TDataType>
bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(
    const TContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    const bool local_value = block_for_each<MaxReduction<bool>>(rContainer, [&](const auto& rEntity){
        return rEntity.GetProperties().Has(rVariable);
    });

    return rDataCommunicator.OrReduceAll(local_value);

    KRATOS_CATCH("");
}

template<class TContainerType>
void OptimizationUtils::CreateEntitySpecificPropertiesForContainer(
    ModelPart& rModelPart,
    TContainerType& rContainer)
{
    KRATOS_TRY

    auto element_properties_id = block_for_each<MaxReduction<IndexType>>(rContainer, [](const auto& rEntity) {
        return rEntity.GetProperties().Id();
    });

    auto properties_id = block_for_each<MaxReduction<IndexType>>(rModelPart.PropertiesArray(), [](const auto pProperties) {
        return pProperties->Id();
    });

    properties_id = std::max(element_properties_id, properties_id);

    // creation of properties is done in serial
    for (auto& r_entity : rContainer) {
        auto p_properties = rModelPart.CreateNewProperties(++properties_id);
        const auto& element_properties = r_entity.GetProperties();
        *p_properties = element_properties;
        p_properties->SetId(properties_id);
        r_entity.SetProperties(p_properties);
    }

    KRATOS_CATCH("");
}

template<>
KRATOS_API(OPTIMIZATION_APPLICATION) IndexType OptimizationUtils::GetVariableDimension(
    const Variable<double>& rVariable,
    const IndexType DomainSize)
{
    return 1;
}

template<>
KRATOS_API(OPTIMIZATION_APPLICATION) IndexType OptimizationUtils::GetVariableDimension(
    const Variable<array_1d<double, 3>>& rVariable,
    const IndexType DomainSize)
{
    return DomainSize;
}

void OptimizationUtils::CopySolutionStepVariablesList(
    ModelPart& rDestinationModelPart,
    const ModelPart& rOriginModelPart)
{
    rDestinationModelPart.GetNodalSolutionStepVariablesList() = rOriginModelPart.GetNodalSolutionStepVariablesList();
}


IndexType OptimizationUtils::Factorial(const IndexType N)
{
    return (N == 1 || N == 0) ? 1 : Factorial(N - 1) * N;
}

IndexType OptimizationUtils::NChooseK(
    const IndexType N,
    const IndexType K)
{
    return Factorial(N) / (Factorial(K) * Factorial(N - K));
}

template<class TContainerType>
ContainerExpression<TContainerType> OptimizationUtils::SmoothClamp(
    const ContainerExpression<TContainerType>& rInputExpression,
    const double MinValue,
    const double MaxValue,
    const IndexType NumberOfContinuousDerivatives)
{
    KRATOS_TRY

    const auto& r_input_expression = rInputExpression.GetExpression();
    const auto number_of_entities = r_input_expression.NumberOfEntities();
    const auto number_of_components = r_input_expression.GetItemComponentCount();

    std::vector<IndexType> coeffs(NumberOfContinuousDerivatives + 1);
    for (IndexType i = 0; i < NumberOfContinuousDerivatives + 1; ++i) {
        coeffs[i] = NChooseK(NumberOfContinuousDerivatives + i, i) * NChooseK(2 * NumberOfContinuousDerivatives + 1, NumberOfContinuousDerivatives - i);
    }

    auto p_expression = LiteralFlatExpression<double>::Create(number_of_entities, r_input_expression.GetItemShape());

    IndexPartition<IndexType>(number_of_entities).for_each([&p_expression, &coeffs, &r_input_expression, MinValue, MaxValue, number_of_components](const auto Index) {
        const auto data_begin_index = Index * number_of_components;
        for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
            const auto x = r_input_expression.Evaluate(Index, data_begin_index, i_comp);
            const auto clamped_x = std::clamp((x - MinValue) / (MaxValue - MinValue), 0.0, 1.0);

            auto& result = *(p_expression->begin() + data_begin_index + i_comp);
            result = 0.0;
            for (IndexType j = 0; j < coeffs.size(); ++j) {
                result += coeffs[j] * std::pow(-clamped_x, j);
            }
            result *= std::pow(clamped_x, coeffs.size());
        }
    });

    auto copy = rInputExpression;
    copy.SetExpression(p_expression);
    return copy;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> OptimizationUtils::SmoothClampGradient(
    const ContainerExpression<TContainerType>& rInputExpression,
    const double MinValue,
    const double MaxValue,
    const IndexType NumberOfContinuousDerivatives)
{
    KRATOS_TRY

    const auto& r_input_expression = rInputExpression.GetExpression();
    const auto number_of_entities = r_input_expression.NumberOfEntities();
    const auto number_of_components = r_input_expression.GetItemComponentCount();

    std::vector<IndexType> coeffs(NumberOfContinuousDerivatives + 1);
    for (IndexType i = 0; i < NumberOfContinuousDerivatives + 1; ++i) {
        coeffs[i] = NChooseK(NumberOfContinuousDerivatives + i, i) * NChooseK(2 * NumberOfContinuousDerivatives + 1, NumberOfContinuousDerivatives - i);
    }

    auto p_expression = LiteralFlatExpression<double>::Create(number_of_entities, r_input_expression.GetItemShape());

    IndexPartition<IndexType>(number_of_entities).for_each([&p_expression, &coeffs, &r_input_expression, MinValue, MaxValue, number_of_components](const auto Index) {
        const auto data_begin_index = Index * number_of_components;
        for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
            const auto x = r_input_expression.Evaluate(Index, data_begin_index, i_comp);
            auto& result = *(p_expression->begin() + data_begin_index + i_comp);

            result = 0.0;
            if (x > MinValue && x < MaxValue) {
                double v1 = 0.0;
                for (IndexType j = 0; j < coeffs.size(); ++j) {
                    v1 += coeffs[j] * std::pow(-x, j);
                }

                double v2 = 0.0;
                for (IndexType j = 1; j < coeffs.size(); ++j) {
                    v2 -= coeffs[j] * j * std::pow(-x, j-1);
                }

                result += v1 * coeffs.size() * std::pow(x, coeffs.size() - 1);
                result += v2 * std::pow(x, coeffs.size());
            }
        }
    });

    auto copy = rInputExpression;
    copy.SetExpression(p_expression);
    return copy;

    KRATOS_CATCH("");
}

// template instantiations
template KRATOS_API(OPTIMIZATION_APPLICATION) GeometryData::KratosGeometryType OptimizationUtils::GetContainerEntityGeometryType(const ModelPart::ConditionsContainerType&, const DataCommunicator&);
template KRATOS_API(OPTIMIZATION_APPLICATION) GeometryData::KratosGeometryType OptimizationUtils::GetContainerEntityGeometryType(const ModelPart::ElementsContainerType&, const DataCommunicator&);

template KRATOS_API(OPTIMIZATION_APPLICATION) void OptimizationUtils::CreateEntitySpecificPropertiesForContainer(ModelPart&, ModelPart::ConditionsContainerType&);
template KRATOS_API(OPTIMIZATION_APPLICATION) void OptimizationUtils::CreateEntitySpecificPropertiesForContainer(ModelPart&, ModelPart::ElementsContainerType&);

template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ElementsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);
template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);

template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ElementsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);
template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);

template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::NodesContainerType> OptimizationUtils::SmoothClamp(const ContainerExpression<ModelPart::NodesContainerType>&, const double, const double, const IndexType);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::ConditionsContainerType> OptimizationUtils::SmoothClamp(const ContainerExpression<ModelPart::ConditionsContainerType>&, const double, const double, const IndexType);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::ElementsContainerType> OptimizationUtils::SmoothClamp(const ContainerExpression<ModelPart::ElementsContainerType>&, const double, const double, const IndexType);

template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::NodesContainerType> OptimizationUtils::SmoothClampGradient(const ContainerExpression<ModelPart::NodesContainerType>&, const double, const double, const IndexType);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::ConditionsContainerType> OptimizationUtils::SmoothClampGradient(const ContainerExpression<ModelPart::ConditionsContainerType>&, const double, const double, const IndexType);
template KRATOS_API(OPTIMIZATION_APPLICATION) ContainerExpression<ModelPart::ElementsContainerType> OptimizationUtils::SmoothClampGradient(const ContainerExpression<ModelPart::ElementsContainerType>&, const double, const double, const IndexType);

}