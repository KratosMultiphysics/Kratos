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
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "system_identification_application_variables.h"

// Include base h
#include "sensor_mask_status.h"

namespace Kratos {

SensorMaskStatus::SensorMaskStatus(
    const ModelPart& rSensorModelPart,
    const MasksListType& rMaskPointersList)
    : mpSensorModelPart(&rSensorModelPart),
      mMaskPointersList(rMaskPointersList)
{
}

const Matrix& SensorMaskStatus::GetMaskStatuses() const
{
    return mSensorMaskStatuses;
}

const Matrix& SensorMaskStatus::GetMasks() const
{
    return mSensorMasks;
}

const ModelPart& SensorMaskStatus::GetSensorModelPart() const
{
    return *mpSensorModelPart;
}

void SensorMaskStatus::Update()
{
    KRATOS_TRY

    KRATOS_ERROR_IF(std::visit([&](const auto& rMasksPointersList)
                               { return rMasksPointersList.empty(); }, mMaskPointersList))
        << "Please provide non-empty masks list.";

    const auto number_of_masks = std::visit([](const auto& rMasksPointersList)
                                            { return rMasksPointersList.size(); }, mMaskPointersList);

    KRATOS_ERROR_IF_NOT(number_of_masks == mpSensorModelPart->NumberOfNodes())
        << "Number of nodes and number of masks mismatch [ number of nodes = "
        << mpSensorModelPart->NumberOfNodes() << ", number of masks = "
        << number_of_masks << " ].";

    auto& r_mask_statuses = mSensorMasks;

    const auto number_of_entities = std::visit([&r_mask_statuses](const auto& rMasksPointersList) {
        const IndexType number_of_entities = rMasksPointersList.front()->GetContainer().size();
        const IndexType number_of_masks = rMasksPointersList.size();

        if (r_mask_statuses.size1() != number_of_entities || r_mask_statuses.size2() != number_of_masks) {
            r_mask_statuses.resize(number_of_entities, number_of_masks, false);
        }

        IndexPartition<IndexType>(number_of_entities).for_each([&rMasksPointersList, &r_mask_statuses, number_of_masks](const auto iEntity) {
            for (IndexType i_mask = 0; i_mask < number_of_masks; ++i_mask) {
                r_mask_statuses(iEntity, i_mask) = rMasksPointersList[i_mask]->GetExpression().Evaluate(iEntity, iEntity, 0);
            }
        });
        return number_of_entities;

    }, mMaskPointersList);

    Vector sensor_status(mSensorMasks.size2());

    IndexPartition<IndexType>(sensor_status.size()).for_each([&](const auto iSensor) {
        sensor_status[iSensor] = (mpSensorModelPart->NodesBegin() + iSensor)->GetValue(SENSOR_STATUS);
    });

    if (mSensorMaskStatuses.size1() != number_of_entities || mSensorMaskStatuses.size2() != number_of_masks) {
        mSensorMaskStatuses.resize(number_of_entities, number_of_masks, false);
    }

    IndexPartition<IndexType>(mSensorMaskStatuses.size1()).for_each([&](const auto iEntity) {
        const Vector& r_values = row(mSensorMasks, iEntity);
        for (IndexType i_sensor = 0; i_sensor < mSensorMaskStatuses.size2(); ++i_sensor) {
            mSensorMaskStatuses(iEntity, i_sensor) = r_values[i_sensor] * sensor_status[i_sensor];
        }
    });

    KRATOS_CATCH("");
}

} // namespace Kratos