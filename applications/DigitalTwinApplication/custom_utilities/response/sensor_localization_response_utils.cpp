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
#include <limits>
#include <cmath>

// External includes
#include "flann/flann.hpp"

// Project includes
#include "expression/literal_flat_expression.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/atomic_utilities.h"

// Application includes
#include "custom_utilities/smooth_clamper.h"
#include "digital_twin_application_variables.h"

// Include base h
#include "sensor_localization_response_utils.h"

namespace Kratos {

SensorLocalizationResponseUtils::SensorLocalizationResponseUtils(
    SensorMaskStatusKDTree<ModelPart::ElementsContainerType>::Pointer pSensorMaskKDTree,
    const double P)
    : mpSensorMaskStatusKDTree(pSensorMaskKDTree),
      mP(P)
{
    KRATOS_TRY

    // get the total domain size and domain size ratios
    const auto& r_container = pSensorMaskKDTree->GetSensorMaskStatus()->GetMaskLocalContainer();
    mDomainSizeRatio.resize(r_container.size());
    const double local_domain_size = IndexPartition<IndexType>(r_container.size()).for_each<SumReduction<double>>([&](const auto Index) {
        const double domain_size = (r_container.begin() + Index)->GetGeometry().DomainSize();
        mDomainSizeRatio[Index] = domain_size;
        return domain_size;
    });
    const double total_domain_size = pSensorMaskKDTree->GetSensorMaskStatus()->GetDataCommunicator().SumAll(local_domain_size);
    block_for_each(mDomainSizeRatio, [total_domain_size](auto& rValue) {
        rValue /= total_domain_size;
    });

    mClusterSizes.resize(r_container.size());

    KRATOS_CATCH("");
}

double SensorLocalizationResponseUtils::CalculateValue()
{
    KRATOS_TRY

    SmoothClamper<ModelPart::ElementsContainerType> clamper(0.0, 1.0);

    // possible number of maximum clusters is the number of elements.
    const IndexType number_of_elements = mpSensorMaskStatusKDTree->GetSensorMaskStatus()->GetMaskLocalContainer().size();

    // getting neighbours for all the elements which are within the radius 1.0 ("0.99999999999999999" is used to make sure that
    // we have all the neighbours within the radius = 1.0, but not the neighbours with 1.0). All other elements which has distance >= 1.0
    // are not relevant since the clamper will anyways make those contribution to zero.
    mpSensorMaskStatusKDTree->GetEntitiesWithinRadius(mNeighbourIndices, mNeighbourSquareDistances, mpSensorMaskStatusKDTree->GetSensorMaskStatus()->GetMaskStatuses(), 1.0 - std::numeric_limits<double>::epsilon());

    // TODO: this will calculate some repeated cluster sizes. Try to avoid them in future.
    const auto local_value = IndexPartition<IndexType>(number_of_elements).for_each<SumReduction<double>>([&](const auto iElement) {
        double& cluster_size = mClusterSizes[iElement];
        cluster_size = 0.0;
        const auto& r_neighbour_square_distances = mNeighbourSquareDistances[iElement];
        const auto& r_neighbour_indices = mNeighbourIndices[iElement];
        for (IndexType i_neighbour = 0; i_neighbour < r_neighbour_square_distances.size(); ++i_neighbour) {
            cluster_size += mDomainSizeRatio[r_neighbour_indices[i_neighbour]] * (1 - clamper.Clamp(r_neighbour_square_distances[i_neighbour]));
        }
        return std::pow(cluster_size, mP);
    });

    mValue = mpSensorMaskStatusKDTree->GetSensorMaskStatus()->GetDataCommunicator().SumAll(local_value);

    return std::pow(mValue, 1 / mP);

    KRATOS_CATCH("");
}

ContainerExpression<ModelPart::NodesContainerType> SensorLocalizationResponseUtils::CalculateGradient() const
{
    KRATOS_TRY

    const auto& r_mask_status = *mpSensorMaskStatusKDTree->GetSensorMaskStatus();
    const auto& r_mask_statuses = r_mask_status.GetMaskStatuses();
    const auto& r_mask_statuses_gradient = r_mask_status.GetMasks();
    const auto& r_data_communicator = r_mask_status.GetDataCommunicator();
    const auto number_of_elements = r_mask_statuses.size1();

    SmoothClamper<ModelPart::NodesContainerType> clamper(0.0, 1.0);

    auto p_expression = LiteralFlatExpression<double>::Create(r_mask_statuses.size2(), {});
    auto& r_expression = *p_expression;

    const double coeff = std::pow(mValue, 1 / mP - 1);

    IndexPartition<IndexType>(r_mask_statuses.size2()).for_each([&](const auto iSensor) {
        double& value = *(r_expression.begin() + iSensor);

        value = 0.0;
        for (IndexType i_element = 0; i_element < number_of_elements; ++i_element) {
            const auto& r_neighbour_square_distances = mNeighbourSquareDistances[i_element];
            const auto& r_neighbour_indices = mNeighbourIndices[i_element];
            const bool i_mask_value = static_cast<bool>(r_mask_statuses_gradient(i_element, iSensor));
            const double cluster_size_derivative = std::pow(mClusterSizes[i_element], mP - 1);
            for (IndexType j_neighbour_element = 0; j_neighbour_element < r_neighbour_square_distances.size(); ++j_neighbour_element) {
                const IndexType neighbour_index = r_neighbour_indices[j_neighbour_element];
                const bool j_mask_value = static_cast<bool>(r_mask_statuses_gradient(neighbour_index, iSensor));
                value -= mDomainSizeRatio[neighbour_index] * clamper.ClampDerivative(r_neighbour_square_distances[j_neighbour_element]) * (i_mask_value ^ j_mask_value) * cluster_size_derivative;
            }
        }

        value = r_data_communicator.SumAll(value) * coeff;
    });

    ContainerExpression<ModelPart::NodesContainerType> result(r_mask_status.GetSensorModelPart());
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/