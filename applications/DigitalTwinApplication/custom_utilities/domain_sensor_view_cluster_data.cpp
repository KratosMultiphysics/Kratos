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
#include <cmath>
#include <type_traits>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_utilities/sensor_utils.h"

// Include base h
#include "domain_sensor_view_cluster_data.h"

namespace Kratos {

template<class TContainerType>
TContainerType& DomainSensorViewClusterData<TContainerType>::GetContainer() const
{
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        return mpModelPart->GetCommunicator().LocalMesh().Nodes();
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        return mpModelPart->GetCommunicator().LocalMesh().Conditions();
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        return mpModelPart->GetCommunicator().LocalMesh().Elements();
    } else {
        static_assert(std::is_same_v<TContainerType, TContainerType>, "Unsupported type.");
    }
}

template<class TContainerType>
DomainSensorViewClusterData<TContainerType>::DomainSensorViewClusterData(
    const SensorViewVectorType& rSensorViewsList,
    const SensorViewVectorType& rRepresentativeSensorViewsForEntitiesList)
    : mpModelPart(&SensorUtils::GetSensorViewsModelPart<TContainerType>(rSensorViewsList)),
      mSensorViewPointersList(rSensorViewsList),
      mRepresentativeSensorsForEntities(rRepresentativeSensorViewsForEntitiesList)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mRepresentativeSensorsForEntities.size() == GetContainer().size())
        << "The representative sensor views for entities list size and number of entities size mismatch "
        << "[ representative sensor views for entities list size = " << mRepresentativeSensorsForEntities.size()
        << ", number of entities size = " << GetContainer().size() << " ].\n";

    KRATOS_CATCH("");
}

template<class TContainerType>
void DomainSensorViewClusterData<TContainerType>::AddDistances(
    const std::string& rDistancesName,
    const std::vector<double>& rDistances)
{
    KRATOS_TRY

    auto p_itr = mDistancesMap.find(rDistancesName);
    KRATOS_ERROR_IF_NOT(p_itr == mDistancesMap.end()) << "A distance measure with name = \"" << rDistancesName << "\" already exists.";

    const auto m = rDistances.size();
    const auto n = static_cast<IndexType>((1.0 + std::sqrt(8.0*m + 1.0)) / 2);

    KRATOS_ERROR_IF_NOT(n == mSensorViewPointersList.size())
        << "The provided distances with name = \"" << rDistancesName
        << "\" does not contain required number of distances. [ provided distances size = " << rDistances.size()
        << ", required distances size = "
        << static_cast<IndexType>(mSensorViewPointersList.size() * (mSensorViewPointersList.size() - 1) / 2) << " ].\n";

    mDistancesMap[rDistancesName] = rDistances;

    KRATOS_CATCH("");
}

template<class TContainerType>
std::vector<double> DomainSensorViewClusterData<TContainerType>::GetDistancesForIndices(
    const std::vector<IndexType>& rIndices,
    const std::string& rDistancesName) const
{
    KRATOS_TRY

    const auto n = mSensorViewPointersList.size();

    auto p_itr = mDistancesMap.find(rDistancesName);
    KRATOS_ERROR_IF(p_itr == mDistancesMap.end())
        << "The provided distances name = \"" << rDistancesName << "\" is not found.";

    auto& r_distances = p_itr->second;

    const auto number_of_clusters = static_cast<int>(rIndices.size());
    const auto number_of_distances = static_cast<int>((number_of_clusters * (number_of_clusters - 1)) / 2);

    std::vector<double> cluster_distances;
    cluster_distances.resize(number_of_distances);

    IndexPartition<int>(number_of_distances).for_each([&cluster_distances, &r_distances, &rIndices, n, number_of_clusters](const auto Index){
        for (int cluster_i = 0; cluster_i < number_of_clusters; ++cluster_i) {
            const auto cluster_j = Index + static_cast<int>((cluster_i + 2) * (cluster_i + 1) / 2) - number_of_clusters * cluster_i;
            if (cluster_i < cluster_j && cluster_j >= 0 && cluster_j < number_of_clusters) {
                const auto global_i = rIndices[cluster_i];
                const auto global_j = rIndices[cluster_j];
                // found the correct cluster_j
                cluster_distances[Index] = r_distances[n * global_i + global_j - static_cast<IndexType>((global_i + 2) * (global_i + 1) / 2)];
                break;
            }
        }
    });

    return cluster_distances;

    KRATOS_CATCH("");
}

template<class TContainerType>
typename DomainSensorViewClusterData<TContainerType>::SensorViewVectorType DomainSensorViewClusterData<TContainerType>::GetSensorViewsForIndices(const std::vector<IndexType>& rIndices) const
{
    KRATOS_TRY

    SensorViewVectorType output(rIndices.size());

    IndexPartition<IndexType>(rIndices.size()).for_each([&](const auto Index) {
        output[Index] = mSensorViewPointersList[rIndices[Index]];
    });

    return output;

    KRATOS_CATCH("");
}

// template instantiations
template class DomainSensorViewClusterData<ModelPart::NodesContainerType>;
template class DomainSensorViewClusterData<ModelPart::ConditionsContainerType>;
template class DomainSensorViewClusterData<ModelPart::ElementsContainerType>;

} /* namespace Kratos.*/