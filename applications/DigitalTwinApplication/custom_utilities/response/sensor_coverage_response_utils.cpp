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
#include <vector>
#include <tuple>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "expression/literal_flat_expression.h"

// Application includes
#include "custom_utilities/smooth_clamper.h"

// Include base h
#include "sensor_coverage_response_utils.h"

namespace Kratos {

template<class TContainerType>
double SensorCoverageResponseUtils::CalculateValue(const SensorMaskStatus<TContainerType>& rSensorMaskStatus)
{
    KRATOS_TRY

    const auto& r_mask_statuses = rSensorMaskStatus.GetMaskStatuses();
    const auto& r_mask_local_container = rSensorMaskStatus.GetMaskLocalContainer();

    SmoothClamper<ModelPart::NodesContainerType> clamper(0.0, 2.0);

    std::vector<double> local_values(2, 0.0);
    std::tie(local_values[0], local_values[1]) = IndexPartition<IndexType>(r_mask_statuses.size1()).for_each<CombinedReduction<SumReduction<double>, SumReduction<double>>>([&r_mask_local_container, &r_mask_statuses, &clamper](const auto iEntity) {
        double value = 0.0;
        const Vector& r_sensor_statuses = row(r_mask_statuses, iEntity);
        for (IndexType i_sensor = 0; i_sensor < r_mask_statuses.size2(); ++i_sensor) {
            value += r_sensor_statuses[i_sensor];
        }
        const double domain_size = (r_mask_local_container.begin() + iEntity)->GetGeometry().DomainSize();
        return std::make_tuple(domain_size * clamper.Clamp(value) * 0.5, domain_size);
    });

    const auto& global_values = rSensorMaskStatus.GetDataCommunicator().SumAll(local_values);

    return global_values[0] / global_values[1];

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<ModelPart::NodesContainerType> SensorCoverageResponseUtils::CalculateGradient(const SensorMaskStatus<TContainerType>& rSensorMaskStatus)
{
    KRATOS_TRY

    const auto& r_mask_statuses = rSensorMaskStatus.GetMaskStatuses();
    const auto& r_mask_statuses_gradient = rSensorMaskStatus.GetMasks();
    const auto& r_mask_local_container = rSensorMaskStatus.GetMaskLocalContainer();
    const auto& r_data_communicator = rSensorMaskStatus.GetDataCommunicator();

    SmoothClamper<ModelPart::NodesContainerType> clamper(0.0, 1.0);

    auto p_expression = LiteralFlatExpression<double>::Create(r_mask_statuses.size2(), {});
    auto& r_expression = *p_expression;

    std::vector<std::pair<double, double>> entity_values(r_mask_statuses.size1());
    const double local_domain_size = IndexPartition<IndexType>(r_mask_statuses.size1()).for_each<SumReduction<double>>([&r_mask_local_container, &r_mask_statuses, &entity_values](const auto iEntity) {
        auto& r_pair = entity_values[iEntity];
        double& value = std::get<0>(r_pair);
        value = 0.0;
        const Vector& r_sensor_statuses = row(r_mask_statuses, iEntity);
        for (IndexType i_sensor = 0; i_sensor < r_mask_statuses.size2(); ++i_sensor) {
            value += r_sensor_statuses[i_sensor];
        }
        const double domain_size = (r_mask_local_container.begin() + iEntity)->GetGeometry().DomainSize();
        std::get<1>(r_pair) = domain_size;
        return domain_size;
    });

    const double global_domain_size = r_data_communicator.SumAll(local_domain_size);

    IndexPartition<IndexType>(r_mask_statuses.size2()).for_each([&r_expression, &r_mask_statuses, &r_mask_statuses_gradient, &r_mask_local_container, &clamper, &entity_values, &r_data_communicator, global_domain_size](const auto iSensor) {
        double& value = *(r_expression.begin() + iSensor);
        value = 0.0;
        for (IndexType i_entity = 0; i_entity < r_mask_local_container.size(); ++i_entity) {
            const auto& r_pair = entity_values[i_entity];
            value += std::get<1>(r_pair) * clamper.ClampDerivative(std::get<0>(r_pair)) * r_mask_statuses_gradient(i_entity, iSensor);
        }

        value = r_data_communicator.SumAll(value) / global_domain_size;
    });

    ContainerExpression<ModelPart::NodesContainerType> result(rSensorMaskStatus.GetSensorModelPart());
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

// template instantiations
template KRATOS_API(DIGITAL_TWIN_APPLICATION) double SensorCoverageResponseUtils::CalculateValue(const SensorMaskStatus<ModelPart::ConditionsContainerType>&);
template KRATOS_API(DIGITAL_TWIN_APPLICATION) ContainerExpression<ModelPart::NodesContainerType> SensorCoverageResponseUtils::CalculateGradient(const SensorMaskStatus<ModelPart::ConditionsContainerType>&);

template KRATOS_API(DIGITAL_TWIN_APPLICATION) double SensorCoverageResponseUtils::CalculateValue(const SensorMaskStatus<ModelPart::ElementsContainerType>&);
template KRATOS_API(DIGITAL_TWIN_APPLICATION) ContainerExpression<ModelPart::NodesContainerType> SensorCoverageResponseUtils::CalculateGradient(const SensorMaskStatus<ModelPart::ElementsContainerType>&);


} /* namespace Kratos.*/