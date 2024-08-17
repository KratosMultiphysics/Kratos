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
#include <limits>
#include <cmath>

// External includes

// Project includes
#include "expression/literal_flat_expression.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/atomic_utilities.h"
#include "input_output/vtu_output.h"

// Application includes
#include "custom_utilities/smooth_clamper.h"
#include "system_identification_application_variables.h"

// Include base h
#include "sensor_localization_response_utils.h"

namespace Kratos {

SensorLocalizationResponseUtils::SensorLocalizationResponseUtils(
    SensorMaskStatusKDTree::Pointer pSensorMaskKDTree,
    const double Beta,
    const double AllowedDissimilarity)
    : mpSensorMaskStatusKDTree(pSensorMaskKDTree),
      mBeta(Beta),
      mAllowedDissimilarity(AllowedDissimilarity)
{
    KRATOS_TRY

    std::visit([&](const auto& pContainer) {
        const auto& r_container = *pContainer;

        mDomainSizeRatio.resize(r_container.size());
        mClusterSizes.resize(r_container.size());
        mNeighbourData.resize(r_container.size());

        const double local_domain_size = IndexPartition<IndexType>(r_container.size()).for_each<SumReduction<double>>([&](const auto Index) {
            const double domain_size = (r_container.begin() + Index)->GetGeometry().DomainSize();
            mDomainSizeRatio[Index] = domain_size;
            return domain_size;
        });
        const double total_domain_size = pSensorMaskKDTree->GetSensorMaskStatus()->pGetMaskModelPart()->GetCommunicator().GetDataCommunicator().SumAll(local_domain_size);
        block_for_each(mDomainSizeRatio, [total_domain_size, AllowedDissimilarity](auto& rValue) {
            rValue /= (total_domain_size * AllowedDissimilarity);
        });
    }, mpSensorMaskStatusKDTree->GetSensorMaskStatus()->pGetMaskContainer());

    KRATOS_CATCH("");
}

double SensorLocalizationResponseUtils::CalculateValue()
{
    KRATOS_TRY

    SmoothClamper<ModelPart::NodesContainerType> clamper(0.0, mAllowedDissimilarity);

    // possible number of maximum clusters is the number of elements.
    const IndexType number_of_elements = mDomainSizeRatio.size();

    // get the sensor mask statuses
    const auto& r_sensor_mask_statuses = mpSensorMaskStatusKDTree->GetSensorMaskStatus()->GetMaskStatuses();

    // TODO: this will calculate some repeated cluster sizes. Try to avoid them in future.
    std::tie(mNumerator, mDenominator) = IndexPartition<IndexType>(number_of_elements).for_each<CombinedReduction<SumReduction<double>, SumReduction<double>>>([&](const auto iElement) {
        double& cluster_size = mClusterSizes[iElement];
        cluster_size = 0.0;

        auto& r_result = mNeighbourData[iElement];

        // getting neighbours for all the elements which are within the radius mAllowedDissimilarity ("0.99999999999999999" is used to make sure that
        // we have all the neighbours within the radius = mAllowedDissimilarity, but not the neighbours with mAllowedDissimilarity). All other elements which has distance >= mAllowedDissimilarity
        // are not relevant since the clamper will anyways make those contribution to zero.
        mpSensorMaskStatusKDTree->RadiusSearch(row(r_sensor_mask_statuses, iElement), mAllowedDissimilarity - std::numeric_limits<double>::epsilon(), r_result);

        for (const auto& r_neighbour_data : r_result) {
            const auto r_neighbour_index = r_neighbour_data.first;
            const auto r_neighbour_squared_distance = r_neighbour_data.second;
            cluster_size += mDomainSizeRatio[r_neighbour_index] * (mAllowedDissimilarity - clamper.ProjectForward(r_neighbour_squared_distance));
        }

        return std::make_tuple(cluster_size * std::exp(mBeta * cluster_size), std::exp(mBeta * cluster_size));
    });

    if (mDenominator > std::numeric_limits<double>::epsilon()) {
        return mNumerator / mDenominator;
    } else {
        return 0.0;
    }

    KRATOS_CATCH("");
}

ContainerExpression<ModelPart::NodesContainerType> SensorLocalizationResponseUtils::CalculateGradient() const
{
    KRATOS_TRY

    const auto& r_mask_status = *mpSensorMaskStatusKDTree->GetSensorMaskStatus();
    const auto& r_mask_statuses = r_mask_status.GetMaskStatuses();
    const auto& r_mask_statuses_gradient = r_mask_status.GetMasks();
    const auto number_of_elements = r_mask_statuses.size1();

    SmoothClamper<ModelPart::NodesContainerType> clamper(0.0, mAllowedDissimilarity);

    auto p_expression = LiteralFlatExpression<double>::Create(r_mask_statuses.size2(), {});
    auto& r_expression = *p_expression;

    if (mDenominator > std::numeric_limits<double>::epsilon()) {
        const double coeff_1 = 1.0 / mDenominator;
        const double coeff_2 = mNumerator / (mDenominator * mDenominator);

        IndexPartition<IndexType>(r_mask_statuses.size2()).for_each([&](const auto iSensor) {
            double numerator_derivative = 0.0;
            double denominator_derivative = 0.0;

            for (IndexType i_element = 0; i_element < number_of_elements; ++i_element) {
                const auto& r_neighbour_data_list = mNeighbourData[i_element];
                const auto i_mask_value = r_mask_statuses_gradient(i_element, iSensor);
                const double cluster_size = mClusterSizes[i_element];
                double cluster_size_derivative = 0.0;

                for (const auto& r_neighbour_data : r_neighbour_data_list) {
                    const auto j_mask_value = r_mask_statuses_gradient(r_neighbour_data.first, iSensor);
                    cluster_size_derivative -= mDomainSizeRatio[r_neighbour_data.first] * clamper.CalculateForwardProjectionGradient(r_neighbour_data.second) * std::abs(i_mask_value - j_mask_value);
                }

                const double temp_1 = std::exp(mBeta * cluster_size);
                numerator_derivative += temp_1 * cluster_size_derivative * (1 + cluster_size * mBeta);
                denominator_derivative += mBeta * cluster_size_derivative * temp_1;
            }

            *(r_expression.begin() + iSensor) = coeff_1 * numerator_derivative - denominator_derivative * coeff_2;
        });
    } else {
        IndexPartition<IndexType>(r_mask_statuses.size2()).for_each([&](const auto Index) {
            *(p_expression->begin() + Index) = 0.0;
        });
    }

    ContainerExpression<ModelPart::NodesContainerType> result(*r_mask_status.pGetSensorModelPart());
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

ContainerExpression<ModelPart::ElementsContainerType> SensorLocalizationResponseUtils::GetClusterSizes() const
{
    KRATOS_TRY

    auto p_expression = LiteralFlatExpression<double>::Create(mClusterSizes.size(), {});
    auto& r_expression = *p_expression;

    IndexPartition<IndexType>(mClusterSizes.size()).for_each([&](const auto Index) {
        *(r_expression.begin() + Index) = std::log(mClusterSizes[Index]) / mBeta;
    });

    ContainerExpression<ModelPart::ElementsContainerType> result(*mpSensorMaskStatusKDTree->GetSensorMaskStatus()->pGetSensorModelPart());
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/