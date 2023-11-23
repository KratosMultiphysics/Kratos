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

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/reduction_utilities.h"
#include "utilities/parallel_utilities.h"

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

// template instantiations
#ifndef KRATOS_SENSOR_UTILS
#define KRATOS_SENSOR_UTILS(CONTAINER_TYPE)                                          \
    template void SensorUtils::IdentifyBestSensorViewForEveryEntity<CONTAINER_TYPE>( \
        std::vector<SensorView<CONTAINER_TYPE>::Pointer>&,                           \
        const std::vector<SensorView<CONTAINER_TYPE>::Pointer>&);                    \
    template ModelPart& SensorUtils::GetSensorViewsModelPart<CONTAINER_TYPE>(        \
        const std::vector<SensorView<CONTAINER_TYPE>::Pointer>&);                    \
    template double SensorUtils::GetDomainSize<CONTAINER_TYPE>(                      \
        const CONTAINER_TYPE&, const DataCommunicator&);

#endif

KRATOS_SENSOR_UTILS(ModelPart::NodesContainerType)
KRATOS_SENSOR_UTILS(ModelPart::ConditionsContainerType)
KRATOS_SENSOR_UTILS(ModelPart::ElementsContainerType)

#undef KRATOS_SENSOR_UTILS

} /* namespace Kratos.*/