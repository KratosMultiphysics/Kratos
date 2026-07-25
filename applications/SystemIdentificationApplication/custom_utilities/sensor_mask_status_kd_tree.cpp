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
    const IndexType LeafMaxSize,
    const IndexType EchoLevel,
    const bool UseKDTree)
    : mpSensorMaskStatus(pSensorMaskStatus),
      mLeafMaxSize(LeafMaxSize),
      mEchoLevel(EchoLevel),
      mUseKDTree(UseKDTree)
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

    if (mUseKDTree) {
        mpKDTreeIndex->radiusSearch(rQueryPoint.data().begin(), Radius, IndicesDistances, nanoflann::SearchParameters());
    } else {
        // do the radius search manually.
        const Matrix& r_sensor_mask_status = mpSensorMaskStatus->GetMaskStatuses();

        IndicesDistances.clear();
        IndicesDistances.reserve(r_sensor_mask_status.size1());

        for (unsigned int i_element = 0; i_element < r_sensor_mask_status.size1(); ++i_element) {
            const Vector& r_element_mask = row(r_sensor_mask_status, i_element);
            const double l1_distance = norm_1(rQueryPoint - r_element_mask);
            if (l1_distance <= Radius) {
                IndicesDistances.emplace_back(i_element, l1_distance);
            }
        }
    }

    KRATOS_CATCH("");
}

void SensorMaskStatusKDTree::Update()
{
    KRATOS_TRY

    const Matrix& r_sensor_mask_status = mpSensorMaskStatus->GetMaskStatuses();

    if (mUseKDTree) {

        mpKratosMatrixKDTreeAdapter = Kratos::make_unique<KratosMatrixKDTreeAdapter>(r_sensor_mask_status);

        mpKDTreeIndex = Kratos::make_unique<KDTreeIndexType>(
            r_sensor_mask_status.size2(), *mpKratosMatrixKDTreeAdapter,
            nanoflann::KDTreeSingleIndexAdaptorParams(mLeafMaxSize, nanoflann::KDTreeSingleIndexAdaptorFlags::None, 0));

        mpKDTreeIndex->buildIndex();

        KRATOS_INFO_IF("SensorMaskStatusKDTree", mEchoLevel > 0)
            << "Updated sensor mask status kd tree in "
            << mpSensorMaskStatus->GetSensorModelPart().FullName() << ".";
    }

    KRATOS_CATCH("");
}

} // namespace Kratos