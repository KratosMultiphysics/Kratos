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

// External includes

// Project includes
#include "utilities/data_type_traits.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "containers/nd_data.h"

// Application includes
#include "custom_utilities/smooth_clamper.h"
#include "custom_utilities/control/sigmoidal_projection_utils.h"
#include "system_identification_application_variables.h"

// Include base h
#include "sensor_localization_response_utils.h"

namespace Kratos {

SensorLocalizationResponseUtils::SensorLocalizationResponseUtils(
    SensorMaskStatusKDTree::Pointer pSensorMaskKDTree,
    const double BoltzmannBeta,
    const double Epsilon)
    : mpSensorMaskStatusKDTree(pSensorMaskKDTree),
      mBoltzmannOperator(BoltzmannBeta),
      mEpsilon(Epsilon)
{
    KRATOS_TRY

    std::visit([&](const auto& pContainer) {
        const auto& r_container = *pContainer;

        using container_type = BareType<decltype(r_container)>;

        mDomainSizeRatio.resize(r_container.size());
        mClusterSizeRatios.resize(r_container.size());

        const double local_domain_size = IndexPartition<IndexType>(r_container.size()).for_each<SumReduction<double>>([&](const auto Index) {
            if constexpr(IsInList<container_type, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>) {
                const double domain_size = (r_container.begin() + Index)->GetGeometry().DomainSize();
                mDomainSizeRatio[Index] = domain_size;
                return domain_size;
            } else {
                KRATOS_ERROR << "TODO: Fix for other container types except condition and element are required.";
                return 0.0;
            }
        });

        const double total_domain_size = pSensorMaskKDTree->GetSensorMaskStatus()->pGetSensorModelPart()->GetCommunicator().GetDataCommunicator().SumAll(local_domain_size);

        block_for_each(mDomainSizeRatio, [total_domain_size](auto& rValue) {
            rValue /= total_domain_size;
        });

    }, mpSensorMaskStatusKDTree->GetSensorMaskStatus()->pGetMaskContainer());

    auto p_sensor_container = mpSensorMaskStatusKDTree->GetSensorMaskStatus()->pGetSensorModelPart()->pNodes();
    mpSensitivities = Kratos::make_shared<TensorAdaptor<double>>(p_sensor_container, Kratos::make_shared<NDData<double>>(DenseVector<unsigned int>(1, p_sensor_container->size())), false);

    KRATOS_CATCH("");
}

double SensorLocalizationResponseUtils::CalculateValue(const double AllowedDissimilarity)
{
    KRATOS_TRY

    using tls_type = std::vector<SensorMaskStatusKDTree::ResultType>;

    const double inverse_allowed_dissimilarity = 1.0 / AllowedDissimilarity;
    const double search_radius = -AllowedDissimilarity * std::log(mEpsilon);

    KRATOS_INFO("SensorLocalizationResponseUtils")
        << "The search radius is " << search_radius << " for allowed dissimilarity of "
        << AllowedDissimilarity << std::endl;

    // possible number of maximum clusters is the number of elements.
    const IndexType number_of_elements = mDomainSizeRatio.size();

    const auto p_sensor_mask_status = mpSensorMaskStatusKDTree->GetSensorMaskStatus();

    // get the sensor mask statuses
    const auto& r_sensor_mask_statuses = p_sensor_mask_status->GetMaskStatuses();

    const auto& r_mask_statuses_gradient = p_sensor_mask_status->GetMasks();

    const IndexType number_of_sensors = r_sensor_mask_statuses.size2();

    Matrix cluster_size_derivatives(number_of_elements, number_of_sensors, 0.0);

    // TODO: this will calculate some repeated cluster sizes. Try to avoid them in future.
    return std::visit([&](const auto& pContainer) {
        const auto& r_container = *pContainer;

        using container_type = BareType<decltype(r_container)>;

        IndexPartition<IndexType>(number_of_elements).for_each(tls_type{}, [&](const auto iElement, auto& rResult) {
            double& cluster_size_ratio = mClusterSizeRatios[iElement];
            cluster_size_ratio = 0.0;
            auto i_mask_row = row(r_mask_statuses_gradient, iElement);

            // getting neighbours for all the elements which are within the radius current_dissimilarity ("0.99999999999999999" is used to make sure that
            // we have all the neighbours within the radius = current_dissimilarity, but not the neighbours with current_dissimilarity). All other elements which has distance >= current_dissimilarity
            // are not relevant since the similarity_clamper will anyways make those contribution to zero.
            mpSensorMaskStatusKDTree->RadiusSearch(row(r_sensor_mask_statuses, iElement), search_radius, rResult);

            for (const auto& r_neighbour_data : rResult) {
                const auto r_neighbour_index = r_neighbour_data.first;
                const auto r_neighbour_squared_distance = r_neighbour_data.second;
                const double weighted_coeff = mDomainSizeRatio[r_neighbour_index] * std::exp(-r_neighbour_squared_distance * inverse_allowed_dissimilarity);

                cluster_size_ratio += weighted_coeff;

                const double derivative_weight = weighted_coeff * inverse_allowed_dissimilarity;
                auto j_mask_row = row(r_mask_statuses_gradient, r_neighbour_index);

                for (IndexType i_sensor = 0; i_sensor < number_of_sensors; ++i_sensor) {
                    cluster_size_derivatives(iElement, i_sensor) -= derivative_weight * std::abs(i_mask_row(i_sensor) - j_mask_row(i_sensor));
                }
            }

            if constexpr(IsInList<container_type, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>) {
                (r_container.begin() + iElement)->SetValue(CLUSTER_SIZE_RATIO, cluster_size_ratio);
            }
        });

        auto p_nd_data = Kratos::make_shared<NDData<double>>(&mClusterSizeRatios[0], DenseVector<unsigned int>(1, mClusterSizeRatios.size()), false);

        mBoltzmannOperator.Update(p_nd_data);
        const double value = mBoltzmannOperator.CalculateValue();

        // now compute the actual gradients
        auto sensitivity_data = mpSensitivities->ViewData();

        auto p_boltzmann_operator_gradient = mBoltzmannOperator.CalculateGradient();
        auto boltzmann_operator_gradient_view = p_boltzmann_operator_gradient->ViewData();

        IndexPartition<IndexType>(sensitivity_data.size()).for_each([&](const auto iSensor) {
            double& sensitivity = sensitivity_data[iSensor];
            sensitivity = 0.0;
            for (IndexType i_element = 0; i_element < number_of_elements; ++i_element) {
                sensitivity += cluster_size_derivatives(i_element, iSensor) * boltzmann_operator_gradient_view[i_element];
            }
        });

        return value;

    }, p_sensor_mask_status->pGetMaskContainer());

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer SensorLocalizationResponseUtils::CalculateGradient(const double AllowedDissimilarity) const
{
    KRATOS_TRY

    return mpSensitivities->Clone();

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer SensorLocalizationResponseUtils::GetClusterSizes() const
{
    KRATOS_TRY

    auto p_result = Kratos::make_shared<TensorAdaptor<double>>(mpSensorMaskStatusKDTree->GetSensorMaskStatus()->pGetSensorModelPart()->pNodes(), Kratos::make_shared<NDData<double>>(DenseVector<unsigned int>(1, mpSensorMaskStatusKDTree->GetSensorMaskStatus()->GetSensorModelPart().NumberOfNodes())));
    auto result_data_view = p_result->ViewData();

    IndexPartition<IndexType>(mClusterSizeRatios.size()).for_each([&](const auto Index) {
        result_data_view[Index] = mClusterSizeRatios[Index];
    });

    return p_result;

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/