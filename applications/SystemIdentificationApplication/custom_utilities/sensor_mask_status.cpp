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
#include "utilities/data_type_traits.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "system_identification_application_variables.h"

// Include base h
#include "sensor_mask_status.h"

namespace Kratos {

SensorMaskStatus::SensorMaskStatus(
    ModelPart& rSensorModelPart,
    const std::vector<TensorAdaptor<double>::Pointer>& rMaskPointersList,
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

TensorAdaptor<double>::ContainerPointerType SensorMaskStatus::pGetMaskContainer() const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mMaskPointersList.empty()) << "Please provide non-empty masks list.";

    return mMaskPointersList.front()->GetContainer();

    KRATOS_CATCH("");
}

void SensorMaskStatus::Update()
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mMaskPointersList.empty()) << "Please provide non-empty masks list.";

    const auto number_of_masks = mMaskPointersList.size();
    const auto number_of_entities = mMaskPointersList.front()->Size();

    KRATOS_ERROR_IF_NOT(number_of_masks == mpSensorModelPart->NumberOfNodes())
        << "Number of nodes and number of masks mismatch [ number of nodes = "
        << mpSensorModelPart->NumberOfNodes() << ", number of masks = "
        << number_of_masks << " ].";

    mSensorMasks.resize(number_of_entities, number_of_masks, false);
    mSensorMaskStatuses.resize(number_of_entities, number_of_masks, false);

    IndexPartition<IndexType>(number_of_masks).for_each([this, number_of_entities](const auto iMask) {
        const auto& r_mask = *(this->mMaskPointersList[iMask]);

        KRATOS_ERROR_IF_NOT(r_mask.Shape().size() == 1)
            << "Masks should be scalar tensor adaptors [ mask = " << r_mask << " ].\n";

        KRATOS_ERROR_IF_NOT(r_mask.Size() == number_of_entities)
            << "Mask number of entities mismatch [ required number of entities = "
            << number_of_entities << ", mask = " << r_mask << " ].\n";

        std::visit([this, &r_mask](auto pReferenceContainer, auto pCurrentContainer){
            using ref_container_type = BareType<decltype(*pReferenceContainer)>;
            using cur_container_type = BareType<decltype(*pCurrentContainer)>;

            if constexpr(std::is_same_v<ref_container_type, cur_container_type>) {
                KRATOS_ERROR_IF_NOT(&*pReferenceContainer == &*pCurrentContainer)
                    << "Container mismatch in tensor adaptors [ reference tensor adaptor = "
                    << *this->mMaskPointersList.front() << ", current tensor adaptor = " << r_mask << " ].\n";
            } else {
                KRATOS_ERROR
                    << "Container mismatch in tensor adaptors [ reference tensor adaptor = "
                    << *this->mMaskPointersList.front() << ", current tensor adaptor = " << r_mask << " ].\n";
            }

        }, this->mMaskPointersList.front()->GetContainer(), r_mask.GetContainer());

        const auto mask_data = r_mask.ViewData();

        for (IndexType i_entity = 0; i_entity < number_of_entities; ++i_entity) {
            this->mSensorMasks(i_entity, iMask) = mask_data[i_entity];
        }
    });

    IndexPartition<IndexType>(mSensorMaskStatuses.size1()).for_each([&](const auto iEntity) {
        const Vector& r_values = row(mSensorMasks, iEntity);
        for (IndexType i_sensor = 0; i_sensor < mSensorMaskStatuses.size2(); ++i_sensor) {
            mSensorMaskStatuses(iEntity, i_sensor) = r_values[i_sensor] * (this->mpSensorModelPart->NodesBegin() + i_sensor)->GetValue(SENSOR_STATUS);
        }
    });

    KRATOS_INFO_IF("SensorMaskStatus", mEchoLevel > 0)
        << "Updated sensor mask status in "
        << mpSensorModelPart->FullName() << ".";

    KRATOS_CATCH("");
}

} // namespace Kratos