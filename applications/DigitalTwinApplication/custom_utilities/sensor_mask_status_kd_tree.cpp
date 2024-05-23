//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <cmath>
#include <limits>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/builtin_timer.h"

// Application includes
#include "digital_twin_application_variables.h"

// Include base h
#include "sensor_mask_status_kd_tree.h"

namespace Kratos {

template<class TContainerType>
SensorMaskStatusKDTree<TContainerType>::SensorMaskStatusKDTree(
    typename SensorMaskStatus<TContainerType>::Pointer pSensorMaskStatus,
    const IndexType NumberOfParallelTrees,
    const bool ExactSearch,
    const IndexType EchoLevel)
    : mpSensorMaskStatus(pSensorMaskStatus),
      mExactSearch(ExactSearch),
      mEchoLevel(EchoLevel),
      mData(mpSensorMaskStatus->GetMaskStatuses().data().begin(), mpSensorMaskStatus->GetMaskStatuses().size1(), mpSensorMaskStatus->GetMaskStatuses().size2()),
      mKDTree(mData, flann::KDTreeIndexParams(1)) // since we are using an exact search, there is no point in using more than one tree.
{
}

template<class TContainerType>
void SensorMaskStatusKDTree<TContainerType>::GetEntitiesWithinRadius(
    std::vector<std::vector<long unsigned int>>& rIndices,
    std::vector<std::vector<double>>& rDistances,
    Matrix& rQueries,
    const double Radius)
{
    KRATOS_TRY

    BuiltinTimer timer;

    flann::Matrix<double> queries(rQueries.data().begin(), rQueries.size1(), rQueries.size2());

    flann::SearchParams params;
    params.cores = 0; // uses all the cores available.
    if (mExactSearch) {
        params.checks = flann::FLANN_CHECKS_UNLIMITED; // for an exact search within the radius.
    }
    mKDTree.radiusSearch(queries, rIndices, rDistances, Radius, params);

    KRATOS_INFO_IF("SensorMaskStatusKDTree", mEchoLevel > 0) << "Found neighbours in: " << timer.ElapsedSeconds() << " s" << std::endl;

    KRATOS_CATCH("");
}

template<class TContainerType>
void SensorMaskStatusKDTree<TContainerType>::Update()
{
    KRATOS_TRY

    BuiltinTimer timer;

    mKDTree.buildIndex();

    KRATOS_INFO_IF("SensorMaskStatusKDTree", mEchoLevel > 0) << "Updated KD tree in: " << timer.ElapsedSeconds() << " s" << std::endl;

    KRATOS_CATCH("");
}

template<class TContainerType>
typename SensorMaskStatus<TContainerType>::Pointer SensorMaskStatusKDTree<TContainerType>::GetSensorMaskStatus()
{
    return mpSensorMaskStatus;
}

// template instantiations
template class SensorMaskStatusKDTree<ModelPart::NodesContainerType>;
template class SensorMaskStatusKDTree<ModelPart::ConditionsContainerType>;
template class SensorMaskStatusKDTree<ModelPart::ElementsContainerType>;

} // namespace Kratos