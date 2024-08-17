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
    ModelPart& rSensorModelPart,
    const MasksListType& rMaskPointersList,
    const IndexType EchoLevel)
    : mpSensorModelPart(&rSensorModelPart),
      mMaskPointersList(rMaskPointersList),
      mEchoLevel(EchoLevel)
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

ModelPart * SensorMaskStatus::pGetSensorModelPart() const
{
    return mpSensorModelPart;
}

ModelPart * SensorMaskStatus::pGetMaskModelPart() const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(std::visit([](const auto& rMasksPointersList)
                               { return rMasksPointersList.empty(); }, mMaskPointersList))
        << "Please provide non-empty masks list.";

    return std::visit([](const auto& rMasksPointersList) {
        return rMasksPointersList.front()->pGetModelPart();
    }, mMaskPointersList);

    KRATOS_CATCH("");
}


SensorMaskStatus::MaskContainerPointerType SensorMaskStatus::pGetMaskContainer() const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(std::visit([](const auto& rMasksPointersList)
                               { return rMasksPointersList.empty(); }, mMaskPointersList))
        << "Please provide non-empty masks list.";

    return std::visit([](const auto& rMasksPointersList) -> MaskContainerPointerType {
        auto& r_container = rMasksPointersList.front()->GetContainer();
        using container_type = std::remove_reference_t<decltype(r_container)>;
        if constexpr(std::is_same_v<container_type, ModelPart::ConditionsContainerType>) {
            return rMasksPointersList.front()->pGetModelPart()->pConditions();
        } else if constexpr(std::is_same_v<container_type, ModelPart::ElementsContainerType>) {
            return rMasksPointersList.front()->pGetModelPart()->pElements();
        } else {
            static_assert(!std::is_same_v<container_type, container_type>)
                << "Unsupported container type.";
            return rMasksPointersList.front()->pGetModelPart()->pElements();
        }
    }, mMaskPointersList);

    KRATOS_CATCH("");
}

void SensorMaskStatus::Update()
{
    KRATOS_TRY

    KRATOS_ERROR_IF(std::visit([](const auto& rMasksPointersList)
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

    KRATOS_INFO_IF("SensorMaskStatus", mEchoLevel > 0)
        << "Updated sensor mask status in "
        << mpSensorModelPart->FullName() << ".";

    KRATOS_CATCH("");
}

} // namespace Kratos