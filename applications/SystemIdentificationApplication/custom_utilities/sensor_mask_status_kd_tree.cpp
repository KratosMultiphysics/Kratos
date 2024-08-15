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

// External includes

// Project includes

// Application includes

// Include base h
#include "sensor_mask_status_kd_tree.h"

namespace Kratos {
///@name Kratos Classes
///@{

SensorMaskStatusKDTree::SensorMaskStatusKDTree(
    SensorMaskStatus::Pointer pSensorMaskStatus,
    const IndexType LeafMaxSize)
    : mpSensorMaskStatus(pSensorMaskStatus),
      mLeafMaxSize(LeafMaxSize)
{
}

SensorMaskStatus::Pointer SensorMaskStatusKDTree::GetSensorMaskStatus() const
{
    return mpSensorMaskStatus;
}

void SensorMaskStatusKDTree::RadiusSearch(
    const Vector& rQueryPoint,
    const double Radius,
    std::vector<ResultType> &IndicesDistances) const
{
    KRATOS_TRY

    mpKDTreeIndex->radiusSearch(rQueryPoint.data().begin(), Radius, IndicesDistances, nanoflann::SearchParameters());

    KRATOS_CATCH("");
}

void SensorMaskStatusKDTree::Update()
{
    KRATOS_TRY

    const Matrix& r_sensor_mask_status = mpSensorMaskStatus->GetMaskStatuses();

    mpKratosMatrixKDTreeAdapter = Kratos::make_unique<KratosMatrixKDTreeAdapter>(r_sensor_mask_status);

    mpKDTreeIndex = Kratos::make_unique<KDTreeIndexType>(
        r_sensor_mask_status.size2(), *mpKratosMatrixKDTreeAdapter,
        nanoflann::KDTreeSingleIndexAdaptorParams(mLeafMaxSize, nanoflann::KDTreeSingleIndexAdaptorFlags::None, 0));

    mpKDTreeIndex->buildIndex();

    KRATOS_CATCH("");
}

} // namespace Kratos