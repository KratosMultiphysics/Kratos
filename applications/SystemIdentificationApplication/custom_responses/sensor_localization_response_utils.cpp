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

    KRATOS_CATCH("");
}

double SensorLocalizationResponseUtils::CalculateValue(const double AllowedDissimilarity)
{
    KRATOS_TRY

    using tls_type = std::vector<SensorMaskStatusKDTree::ResultType>;

    const double search_radius = -AllowedDissimilarity * std::log(mEpsilon);

    KRATOS_INFO("SensorLocalizationResponseUtils")
        << "The search radius is " << search_radius << " for allowed dissimilarity of "
        << AllowedDissimilarity << std::endl;

    // possible number of maximum clusters is the number of elements.
    const IndexType number_of_elements = mDomainSizeRatio.size();

    // get the sensor mask statuses
    const auto& r_sensor_mask_statuses = mpSensorMaskStatusKDTree->GetSensorMaskStatus()->GetMaskStatuses();

    // TODO: this will calculate some repeated cluster sizes. Try to avoid them in future.
    return std::visit([&](const auto& pContainer) {
        const auto& r_container = *pContainer;

        using container_type = BareType<decltype(r_container)>;

        IndexPartition<IndexType>(number_of_elements).for_each(tls_type{}, [&](const auto iElement, auto& rResult) {
            double& cluster_size_ratio = mClusterSizeRatios[iElement];
            cluster_size_ratio = 0.0;

            // getting neighbours for all the elements which are within the radius current_dissimilarity ("0.99999999999999999" is used to make sure that
            // we have all the neighbours within the radius = current_dissimilarity, but not the neighbours with current_dissimilarity). All other elements which has distance >= current_dissimilarity
            // are not relevant since the similarity_clamper will anyways make those contribution to zero.
            mpSensorMaskStatusKDTree->RadiusSearch(row(r_sensor_mask_statuses, iElement), search_radius, rResult);

            for (const auto& r_neighbour_data : rResult) {
                const auto r_neighbour_index = r_neighbour_data.first;
                const auto r_neighbour_squared_distance = r_neighbour_data.second;
                cluster_size_ratio += mDomainSizeRatio[r_neighbour_index] * std::exp(-r_neighbour_squared_distance / AllowedDissimilarity);
            }

            if constexpr(IsInList<container_type, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>) {
                (r_container.begin() + iElement)->SetValue(CLUSTER_SIZE_RATIO, cluster_size_ratio);
            }
        });

        auto p_nd_data = Kratos::make_shared<NDData<double>>(&mClusterSizeRatios[0], DenseVector<unsigned int>(1, mClusterSizeRatios.size()), false);

        mBoltzmannOperator.Update(p_nd_data);
        return mBoltzmannOperator.CalculateValue();

    }, mpSensorMaskStatusKDTree->GetSensorMaskStatus()->pGetMaskContainer());

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer SensorLocalizationResponseUtils::CalculateGradient(const double AllowedDissimilarity) const
{
    KRATOS_TRY

    using tls_type = std::vector<SensorMaskStatusKDTree::ResultType>;

    const auto& r_mask_status = *mpSensorMaskStatusKDTree->GetSensorMaskStatus();
    const auto& r_mask_statuses = r_mask_status.GetMaskStatuses();
    const auto& r_mask_statuses_gradient = r_mask_status.GetMasks();
    const auto number_of_elements = r_mask_statuses.size1();

    const double search_radius = -AllowedDissimilarity * std::log(mEpsilon);

    auto p_result = Kratos::make_shared<TensorAdaptor<double>>(r_mask_status.pGetSensorModelPart()->pNodes(), Kratos::make_shared<NDData<double>>(DenseVector<unsigned int>(1, r_mask_status.GetSensorModelPart().NumberOfNodes())));
    auto result_data_view = p_result->ViewData();

    IndexPartition<IndexType>(result_data_view.size()).for_each([&result_data_view](const auto Index) {
        result_data_view[Index] = 0.0;
    });

    // get the sensor mask statuses
    const auto& r_sensor_mask_statuses = mpSensorMaskStatusKDTree->GetSensorMaskStatus()->GetMaskStatuses();

    auto p_boltzmann_operator_gradient = mBoltzmannOperator.CalculateGradient();
    auto boltzmann_operator_gradient_view = p_boltzmann_operator_gradient->ViewData();

    IndexPartition<IndexType>(number_of_elements).for_each(tls_type{}, [&](const auto iElement, auto& rResult){

        mpSensorMaskStatusKDTree->RadiusSearch(row(r_sensor_mask_statuses, iElement), search_radius, rResult);

        for (IndexType i_sensor = 0; i_sensor < r_mask_statuses.size2(); ++i_sensor) {
            const auto i_mask_value = r_mask_statuses_gradient(iElement, i_sensor);

            double cluster_size_derivative = 0.0;

            for (const auto& r_neighbour_data : rResult) {
                const auto j_mask_value = r_mask_statuses_gradient(r_neighbour_data.first, i_sensor);
                cluster_size_derivative -= mDomainSizeRatio[r_neighbour_data.first] * std::exp(-r_neighbour_data.second / AllowedDissimilarity) * std::abs(i_mask_value - j_mask_value) / AllowedDissimilarity;
            }

            result_data_view[i_sensor] += cluster_size_derivative * boltzmann_operator_gradient_view[iElement];
        }
    });

    return p_result;

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