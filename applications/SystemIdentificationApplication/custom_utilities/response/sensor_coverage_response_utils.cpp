//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: SystemIdentificationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
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

double SensorCoverageResponseUtils::CalculateValue(const SensorMaskStatus& rSensorMaskStatus)
{
    KRATOS_TRY

    const auto& r_mask_statuses = rSensorMaskStatus.GetMaskStatuses();

    SmoothClamper<ModelPart::NodesContainerType> clamper(0.0, 1.0);

    std::vector<double> local_values(2, 0.0);
    std::tie(local_values[0], local_values[1]) = std::visit([&r_mask_statuses, &clamper](const auto& pMaskLocalContainer) {
        const auto& r_mask_local_container = *pMaskLocalContainer;

        return IndexPartition<IndexType>(r_mask_statuses.size1()).for_each<CombinedReduction<SumReduction<double>, SumReduction<double>>>([&r_mask_local_container, &r_mask_statuses, &clamper](const auto iEntity) {
                double value = 0.0;

                const Vector& r_sensor_statuses = row(r_mask_statuses, iEntity);
                for (IndexType i_sensor = 0; i_sensor < r_mask_statuses.size2(); ++i_sensor) {
                    value += r_sensor_statuses[i_sensor];
                }

                if constexpr(std::is_same_v<std::remove_cv_t<std::decay_t<decltype(r_mask_local_container)>>, ModelPart::NodesContainerType>) {
                    KRATOS_ERROR << "TODO: Fix for nodal expression still required.";
                    return std::make_tuple(0.0, 0.0);
                } else {
                    const double domain_size = (r_mask_local_container.begin() + iEntity)->GetGeometry().DomainSize();
                    return std::make_tuple(domain_size * clamper.ProjectForward(value), domain_size);
                }
            });
    }, rSensorMaskStatus.pGetMaskContainer());

    const auto& global_values = rSensorMaskStatus.GetSensorModelPart().GetCommunicator().GetDataCommunicator().SumAll(local_values);

    return global_values[0] / global_values[1];

    KRATOS_CATCH("");
}

ContainerExpression<ModelPart::NodesContainerType> SensorCoverageResponseUtils::CalculateGradient(const SensorMaskStatus& rSensorMaskStatus)
{
    KRATOS_TRY

    const auto& r_mask_statuses = rSensorMaskStatus.GetMaskStatuses();
    const auto& r_mask_statuses_gradient = rSensorMaskStatus.GetMasks();
    const auto& r_data_communicator = rSensorMaskStatus.GetSensorModelPart().GetCommunicator().GetDataCommunicator();

    SmoothClamper<ModelPart::NodesContainerType> clamper(0.0, 1.0);

    std::vector<std::pair<double, double>> entity_values(r_mask_statuses.size1());
    const double local_domain_size = std::visit([&r_mask_statuses, &clamper, &entity_values](const auto& pMaskLocalContainer) {
        const auto& r_mask_local_container = *pMaskLocalContainer;

        return IndexPartition<IndexType>(r_mask_statuses.size1()).for_each<SumReduction<double>>([&r_mask_local_container, &r_mask_statuses, &entity_values](const auto iEntity) {
            auto& r_pair = entity_values[iEntity];
            double& value = std::get<0>(r_pair);
            value = 0.0;
            const Vector& r_sensor_statuses = row(r_mask_statuses, iEntity);
            for (IndexType i_sensor = 0; i_sensor < r_mask_statuses.size2(); ++i_sensor) {
                value += r_sensor_statuses[i_sensor];
            }

            if constexpr(std::is_same_v<std::remove_cv_t<std::decay_t<decltype(r_mask_local_container)>>, ModelPart::NodesContainerType>) {
                KRATOS_ERROR << "TODO: Fix for nodal expression still required.";
                return 0.0;
            } else {
                const double domain_size = (r_mask_local_container.begin() + iEntity)->GetGeometry().DomainSize();
                std::get<1>(r_pair) = domain_size;
                return domain_size;
            }
        });

    }, rSensorMaskStatus.pGetMaskContainer());

    const double global_domain_size = r_data_communicator.SumAll(local_domain_size);

    auto p_expression = LiteralFlatExpression<double>::Create(r_mask_statuses.size2(), {});
    auto& r_expression = *p_expression;

    std::visit([&r_expression, &r_mask_statuses, &r_mask_statuses_gradient, &clamper, &entity_values, &r_data_communicator, global_domain_size](const auto& pMaskLocalContainer) {
        const auto& r_mask_local_container = *pMaskLocalContainer;

        IndexPartition<IndexType>(r_mask_statuses.size2()).for_each([&r_expression, &r_mask_statuses, &r_mask_statuses_gradient, &r_mask_local_container, &clamper, &entity_values, &r_data_communicator, global_domain_size](const auto iSensor) {
            double& value = *(r_expression.begin() + iSensor);
            value = 0.0;
            for (IndexType i_entity = 0; i_entity < r_mask_local_container.size(); ++i_entity) {
                const auto& r_pair = entity_values[i_entity];
                value += std::get<1>(r_pair) * clamper.CalculateForwardProjectionGradient(std::get<0>(r_pair)) * r_mask_statuses_gradient(i_entity, iSensor);
            }

            value = r_data_communicator.SumAll(value) / global_domain_size;
        });
    }, rSensorMaskStatus.pGetMaskContainer());

    ContainerExpression<ModelPart::NodesContainerType> result(*rSensorMaskStatus.pGetSensorModelPart());
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/