//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Wranakulasuriya
//

// System includes
#include <algorithm>
#include <numeric>
#include <cmath>

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/reduction_utilities.h"
#include "utilities/parallel_utilities.h"
#include "expression/literal_flat_expression.h"

// Application includes

// Include base h
#include "sensor_utils.h"

namespace Kratos {

template<class TContainerType>
ModelPart& SensorUtils::GetSensorViewsModelPart(const std::vector<typename SensorView<TContainerType>::Pointer>& rSensorViews)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rSensorViews.size() == 0) << "Empty sensor views list provided.";

    auto& p_first_sensor_view = rSensorViews.front();
    auto p_model_part = p_first_sensor_view->GetContainerExpression()->pGetModelPart();

    const auto total_sum = block_for_each<SumReduction<IndexType>>(rSensorViews, [&p_model_part](const auto& pSensorView) -> IndexType{
        return p_model_part == pSensorView->GetContainerExpression()->pGetModelPart();
    });

    KRATOS_ERROR_IF_NOT(total_sum == rSensorViews.size()) << "The sensor views list should only have container expressions with the same model part.";

    return *p_model_part;

    KRATOS_CATCH("");
}

template<class TContainerType>
void SensorUtils::IdentifyBestSensorViewForEveryEntity(
    std::vector<typename SensorView<TContainerType>::Pointer>& rOutputSensorViews,
    const std::vector<typename SensorView<TContainerType>::Pointer>& rNormalizedSensorViews)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rNormalizedSensorViews.size() == 0) << "Empty normalized sensor views list provided.";

    const auto& p_first_sensor_view = rNormalizedSensorViews.front();
    const auto& r_container = p_first_sensor_view->GetContainerExpression()->GetContainer();

    rOutputSensorViews.resize(r_container.size());

    IndexPartition<IndexType>(r_container.size()).for_each([&rOutputSensorViews, &rNormalizedSensorViews](const auto Index){
        double max_entity_contribution = 0.0;
        typename SensorView<TContainerType>::Pointer p_max_contributing_sensor_view;
        for (const auto& p_normalized_sensor_view :  rNormalizedSensorViews) {
            const double current_entity_contribution = p_normalized_sensor_view->GetContainerExpression()->GetExpression().Evaluate(Index, Index, 0);
            if (current_entity_contribution > max_entity_contribution) {
                max_entity_contribution = current_entity_contribution;
                p_max_contributing_sensor_view = p_normalized_sensor_view;
            }
        }

        rOutputSensorViews[Index] = p_max_contributing_sensor_view;
    });

    KRATOS_CATCH("");
}

template<class TContainerType>
double SensorUtils::GetDomainSize(
    const TContainerType& rContainer,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        return 0.0;
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType> || std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        return rDataCommunicator.SumAll(block_for_each<SumReduction<double>>(rContainer, [](const auto& rEntity) { return rEntity.GetGeometry().DomainSize(); }));
    } else {
        static_assert(std::is_same_v<TContainerType, TContainerType>, "Unsupported type.");
    }

    KRATOS_CATCH("");
}

template<class TContainerType>
void SensorUtils::AssignEntitiesToClustersBasedOnOptimalSensor(
    const std::vector<typename SensorViewCluster<TContainerType>::Pointer>& rSensorViewClusters,
    const std::vector<typename ContainerExpression<TContainerType>::Pointer>& rContainerExpressionsList)
{
    KRATOS_TRY

    const auto number_of_clusters = rSensorViewClusters.size();

    KRATOS_ERROR_IF_NOT(number_of_clusters == rContainerExpressionsList.size())
        << "Number of clusters and number of container expressions mismatch [ "
        << "number of clusters = " << number_of_clusters << ", number of container expressions = "
        << rContainerExpressionsList.size() << " ].\n";

    KRATOS_ERROR_IF(number_of_clusters == 0) << "No clusters provided.";

    auto& r_container = rContainerExpressionsList[0]->GetContainer();
    const auto number_of_entities = r_container.size();
    std::vector<IndexType> indices(number_of_entities);

    IndexPartition<IndexType>(number_of_entities).for_each([&indices, &rContainerExpressionsList, number_of_clusters](const auto Index) {
        double best_value = 0.0;
        IndexType best_index = 0;
        for (IndexType j = 0; j < number_of_clusters; ++j) {
            const double current_value = rContainerExpressionsList[j]->GetExpression().Evaluate(Index, Index, 0);
            if (current_value > best_value) {
                best_value = current_value;
                best_index = j;
            }
        }
        indices[Index] = best_index;
    });

    // clean existing lists in the clusters
    for (auto& p_sensor_view_cluster : rSensorViewClusters) {
        p_sensor_view_cluster->GetEntities().clear();
    }

    // now fill the entities containers in serial
    for (IndexType i = 0; i < number_of_entities; ++i) {
        rSensorViewClusters[indices[i]]->GetEntities().push_back(&*(r_container.begin() + i));
    }

    // now sort and make unique
    for (auto& p_sensor_view_cluster : rSensorViewClusters) {
        p_sensor_view_cluster->GetEntities().Unique();
    }

    KRATOS_CATCH("");
}

template<class TContainerType>
void SensorUtils::GetThresholdSensorViews(
    const double ThresholdValue,
    const std::string& rExpressionName,
    std::vector<typename SensorView<TContainerType>::Pointer>& rOutputSensorViews,
    const std::vector<typename SensorView<TContainerType>::Pointer>& rNormalizedSensorViews)
{
    KRATOS_TRY

    rOutputSensorViews.resize(rNormalizedSensorViews.size());

    for (IndexType i = 0; i < rNormalizedSensorViews.size(); ++i) {
        auto& p_normalized_sensor_view = rNormalizedSensorViews[i];
        auto& r_input_expression = p_normalized_sensor_view->GetContainerExpression()->GetExpression();
        auto p_output_expression = LiteralFlatExpression<int>::Create(r_input_expression.NumberOfEntities(), r_input_expression.GetItemShape());
        const auto number_of_components = r_input_expression.GetItemComponentCount();

        IndexPartition<IndexType>(r_input_expression.NumberOfEntities()).for_each([&p_output_expression, &r_input_expression, ThresholdValue, number_of_components](const auto Index){
            const IndexType data_begin_index = Index * number_of_components;
            for (IndexType i = 0; i < number_of_components; ++i) {
                *(p_output_expression->begin() + i) = (r_input_expression.Evaluate(Index, data_begin_index, i) > ThresholdValue);
            }
        });

        auto p_copy = p_normalized_sensor_view->GetContainerExpression()->Clone();
        p_copy->SetExpression(p_output_expression);

        p_normalized_sensor_view->GetSensor()->template AddContainerExpression<TContainerType>(rExpressionName, p_copy);
        rOutputSensorViews[i] = Kratos::make_shared<SensorView<TContainerType>>(p_normalized_sensor_view->GetSensor(), rExpressionName);
    }

    KRATOS_CATCH("");
}

template<class TContainerType>
std::pair<IndexType, typename ContainerExpression<TContainerType>::Pointer> SensorUtils::GetEntityCoverageMask(
    const SensorView<TContainerType>& rSensorView)
{
    KRATOS_TRY

    const auto& p_input_container_expression = rSensorView.GetContainerExpression();
    const auto& r_input_expression = p_input_container_expression->GetExpression();
    const auto number_of_entities = r_input_expression.NumberOfEntities();

    struct Data
    {
        IndexType mIndex;
        double mValue;
        bool operator<(const Data& rRight) const { return mValue < rRight.mValue; }
    };

    std::vector<Data> index_value_pairs_vector;
    index_value_pairs_vector.resize(number_of_entities);

    IndexPartition<IndexType>(number_of_entities).for_each([&index_value_pairs_vector, &r_input_expression](const auto Index) {
        index_value_pairs_vector[Index].mIndex = Index;
        index_value_pairs_vector[Index].mValue = r_input_expression.Evaluate(Index, Index, 0);
    });

    // now sort expression values
    std::sort(index_value_pairs_vector.begin(), index_value_pairs_vector.end(), [](const auto& rV1, const auto& rV2){
        return rV1.mValue > rV2.mValue;
    });

    // now find the coverage
    const auto& r_data = IndexPartition<IndexType>(number_of_entities).for_each<MaxReduction<Data>>([&index_value_pairs_vector](const auto Index){
        const double current_sum = std::accumulate(index_value_pairs_vector.begin(), index_value_pairs_vector.begin() + Index, 0.0, [](const auto& rValue, const auto& rIndexValuePair) {
                                                return rValue + rIndexValuePair.mValue;
                                            });
        return Data{Index, current_sum / std::sqrt(Index)};
    });

    auto p_expression = LiteralFlatExpression<int>::Create(number_of_entities, {});
    IndexPartition<IndexType>(number_of_entities).for_each([&p_expression](const auto Index) {
        *(p_expression->begin() + Index) = 0;
    });
    IndexPartition<IndexType>(r_data.mIndex + 1).for_each([&p_expression, &index_value_pairs_vector](const auto Index){
        *(p_expression->begin() + index_value_pairs_vector[Index].mIndex) = 1;
    });

    auto output_container_expression = rSensorView.GetContainerExpression()->Clone();
    output_container_expression->SetExpression(p_expression);

    return std::make_pair(r_data.mIndex + 1, output_container_expression);

    KRATOS_CATCH("");
}

template<class TContainerType>
IndexType SensorUtils::CountWithInBounds(
    const ContainerExpression<TContainerType>& rContainer,
    const double LowerBound,
    const double UpperBound)
{
    KRATOS_TRY

    const auto& r_expression = rContainer.GetExpression();
    return IndexPartition<IndexType>(r_expression.NumberOfEntities()).for_each<SumReduction<IndexType>>([&r_expression, LowerBound, UpperBound](const auto Index){
        const double value = r_expression.Evaluate(Index, Index, 0);
        return (value > LowerBound && value <= UpperBound);
    });

    KRATOS_CATCH("");
}

template<class TContainerType>
double SensorUtils::Min(
    const ContainerExpression<TContainerType>& rContainer)
{
    KRATOS_TRY

    const auto& r_expression = rContainer.GetExpression();
    return IndexPartition<IndexType>(r_expression.NumberOfEntities()).for_each<MinReduction<double>>([&r_expression](const auto Index){
        return r_expression.Evaluate(Index, Index, 0);
    });

    KRATOS_CATCH("");
}

template<class TContainerType>
double SensorUtils::Max(
    const ContainerExpression<TContainerType>& rContainer)
{
    KRATOS_TRY

    const auto& r_expression = rContainer.GetExpression();
    return IndexPartition<IndexType>(r_expression.NumberOfEntities()).for_each<MaxReduction<double>>([&r_expression](const auto Index){
        return r_expression.Evaluate(Index, Index, 0);
    });

    KRATOS_CATCH("");
}

template<class TContainerType>
double SensorUtils::Sum(
    const ContainerExpression<TContainerType>& rContainer)
{
    KRATOS_TRY

    const auto& r_expression = rContainer.GetExpression();
    return IndexPartition<IndexType>(r_expression.NumberOfEntities()).for_each<SumReduction<double>>([&r_expression](const auto Index){
        return r_expression.Evaluate(Index, Index, 0);
    });

    KRATOS_CATCH("");
}

// template instantiations
#ifndef KRATOS_SENSOR_UTILS
#define KRATOS_SENSOR_UTILS(CONTAINER_TYPE)                                                  \
    template void SensorUtils::IdentifyBestSensorViewForEveryEntity<CONTAINER_TYPE>(         \
        std::vector<SensorView<CONTAINER_TYPE>::Pointer>&,                                   \
        const std::vector<SensorView<CONTAINER_TYPE>::Pointer>&);                            \
    template ModelPart& SensorUtils::GetSensorViewsModelPart<CONTAINER_TYPE>(                \
        const std::vector<SensorView<CONTAINER_TYPE>::Pointer>&);                            \
    template double SensorUtils::GetDomainSize<CONTAINER_TYPE>(                              \
        const CONTAINER_TYPE&, const DataCommunicator&);                                     \
    template void SensorUtils::AssignEntitiesToClustersBasedOnOptimalSensor<CONTAINER_TYPE>( \
        const std::vector<typename SensorViewCluster<CONTAINER_TYPE>::Pointer>&,             \
        const std::vector<typename ContainerExpression<CONTAINER_TYPE>::Pointer>&);          \
    template void SensorUtils::GetThresholdSensorViews<CONTAINER_TYPE>(                      \
        const double, const std::string&,                                                    \
        std::vector<typename SensorView<CONTAINER_TYPE>::Pointer>&,                          \
        const std::vector<typename SensorView<CONTAINER_TYPE>::Pointer>&);                   \
    template std::pair<IndexType, typename ContainerExpression<CONTAINER_TYPE>::Pointer>     \
    SensorUtils::GetEntityCoverageMask<CONTAINER_TYPE>(const SensorView<CONTAINER_TYPE>&);   \
    template IndexType SensorUtils::CountWithInBounds<CONTAINER_TYPE>(                       \
        const ContainerExpression<CONTAINER_TYPE>&, const double, const double);             \
    template double SensorUtils::Min<CONTAINER_TYPE>(                                        \
        const ContainerExpression<CONTAINER_TYPE>&);                                         \
    template double SensorUtils::Max<CONTAINER_TYPE>(                                        \
        const ContainerExpression<CONTAINER_TYPE>&);                                         \
    template double SensorUtils::Sum<CONTAINER_TYPE>(                                        \
        const ContainerExpression<CONTAINER_TYPE>&);

#endif

KRATOS_SENSOR_UTILS(ModelPart::NodesContainerType)
KRATOS_SENSOR_UTILS(ModelPart::ConditionsContainerType)
KRATOS_SENSOR_UTILS(ModelPart::ElementsContainerType)

#undef KRATOS_SENSOR_UTILS

} /* namespace Kratos.*/