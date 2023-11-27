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

// External includes

// Project includes
#include "utilities/reduction_utilities.h"
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "sensor_view_cluster.h"

namespace Kratos {

template<class TContainerType>
SensorViewCluster<TContainerType>::SensorViewCluster(
    const IndexType ClusterId,
    typename DomainSensorViewClusterData<TContainerType>::Pointer pDataContainer)
    : IndexedObject(ClusterId),
      mpDataContainer(pDataContainer)
{
}

template<class TContainerType>
typename SensorViewCluster<TContainerType>::Pointer SensorViewCluster<TContainerType>::Clone() const
{
    auto p_cluster = Kratos::make_shared<SensorViewCluster<TContainerType>>(this->Id(), mpDataContainer);
    p_cluster->SetEntities(this->GetEntities());
    p_cluster->SetSensorViews(this->GetSensorViews());
    return p_cluster;
}

template<class TContainerType>
void SensorViewCluster<TContainerType>::Clear()
{
    this->mEntities.clear();
    this->mSensorViewIndices.clear();
}

template<class TContainerType>
void SensorViewCluster<TContainerType>::SetSensorViews(const SensorViewVectorType& rSensorViews)
{
    KRATOS_TRY

    mSensorViewIndices.resize(rSensorViews.size(), std::numeric_limits<IndexType>::max());

    const auto& r_global_sensor_views = mpDataContainer->GetSensorViews();

    // first find the sensor view global index
    IndexPartition<IndexType>(rSensorViews.size()).for_each([&](const auto Index) {
        for (IndexType g_index = 0; g_index < r_global_sensor_views.size(); ++g_index) {
            if (*(r_global_sensor_views.begin() + g_index) == rSensorViews[Index]) {
                mSensorViewIndices[Index] = g_index;
                break;
            }
        }

        KRATOS_ERROR_IF(mSensorViewIndices[Index] == std::numeric_limits<IndexType>::max())
            << "The sensor view \"" << rSensorViews[Index]->GetSensor()->GetName()
            << "\" not found in the domain sensor view cluster data.";
    });

    // sort the sensor indices
    std::sort(mSensorViewIndices.begin(), mSensorViewIndices.end());

    KRATOS_CATCH("");
}

template<class TContainerType>
typename SensorViewCluster<TContainerType>::SensorViewVectorType SensorViewCluster<TContainerType>::GetSensorViews() const
{
    KRATOS_TRY

    return mpDataContainer->GetSensorViewsForIndices(mSensorViewIndices);

    KRATOS_CATCH("");
}

template<class TContainerType>
std::vector<double> SensorViewCluster<TContainerType>::GetDistances(const std::string& rDistancesName) const
{
    return mpDataContainer->GetDistancesForIndices(mSensorViewIndices, rDistancesName);
}

// template instantiations
template class SensorViewCluster<ModelPart::NodesContainerType>;
template class SensorViewCluster<ModelPart::ConditionsContainerType>;
template class SensorViewCluster<ModelPart::ElementsContainerType>;

} /* namespace Kratos.*/