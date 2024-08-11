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

template<class TContainerType>
SensorMaskStatus<TContainerType>::SensorMaskStatus(
    ModelPart& rSensorModelPart,
    const std::vector<ContainerExpression<TContainerType>>& rMasksList)
    : mpSensorModelPart(&rSensorModelPart)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rMasksList.empty())
        << "Please provide non-empty masks list.";

    KRATOS_ERROR_IF_NOT(mpSensorModelPart->NumberOfNodes() == rMasksList.size())
        << "Number of sensors and number of masks mismatch [ number of sensors = "
        << mpSensorModelPart->NumberOfNodes() << ", number of masks = "
        << rMasksList.size() << " ].";

    mpMaskModelPart = &*(rMasksList.front().pGetModelPart());
    mpMaskContainer = &(rMasksList.front().GetContainer());
    mpDataCommunicator = &(rMasksList.front().GetModelPart().GetCommunicator().GetDataCommunicator());

    const IndexType number_of_sensors = rMasksList.size();
    const IndexType number_of_entities = mpMaskContainer->size();

    mSensorMasks.resize(number_of_entities, number_of_sensors, false);
    mSensorMaskStatuses.resize(number_of_entities, number_of_sensors, false);

    IndexPartition<IndexType>(number_of_entities).for_each([&](const auto iEntity) {
        for (IndexType i_sensor = 0; i_sensor < number_of_sensors; ++i_sensor) {
            mSensorMasks(iEntity, i_sensor) = rMasksList[i_sensor].GetExpression().Evaluate(iEntity, iEntity, 0);
        }
    });

    KRATOS_CATCH("");
}

template<class TContainerType>
Matrix& SensorMaskStatus<TContainerType>::GetMaskStatuses()
{
    return mSensorMaskStatuses;
}

template<class TContainerType>
const Matrix& SensorMaskStatus<TContainerType>::GetMaskStatuses() const
{
    return mSensorMaskStatuses;
}

template<class TContainerType>
const Matrix& SensorMaskStatus<TContainerType>::GetMasks() const
{
    return mSensorMasks;
}

template<class TContainerType>
ModelPart& SensorMaskStatus<TContainerType>::GetSensorModelPart() const
{
    return *mpSensorModelPart;
}

template<class TContainerType>
const TContainerType& SensorMaskStatus<TContainerType>::GetMaskLocalContainer() const
{
    return *mpMaskContainer;
}

template<class TContainerType>
ModelPart& SensorMaskStatus<TContainerType>::GetMaskModelPart()
{
    return *mpMaskModelPart;
}

template<class TContainerType>
const DataCommunicator& SensorMaskStatus<TContainerType>::GetDataCommunicator() const
{
    return *mpDataCommunicator;
}

template<class TContainerType>
void SensorMaskStatus<TContainerType>::Update()
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mpSensorModelPart->NumberOfNodes() == mSensorMasks.size2())
        << "Number of sensors in the model part and number of stored sensor masks mismatch [ number of sensors in model part = "
        << mpSensorModelPart->NumberOfNodes() << ", number of stored sensor masks = "
        << mSensorMasks.size2() << " ].";

    Vector sensor_status(mSensorMasks.size2());

    IndexPartition<IndexType>(sensor_status.size()).for_each([&](const auto iSensor) {
        sensor_status[iSensor] = (mpSensorModelPart->NodesBegin() + iSensor)->GetValue(SENSOR_STATUS);
    });

    IndexPartition<IndexType>(mSensorMaskStatuses.size1()).for_each([&](const auto iEntity) {
        const Vector& r_values = row(mSensorMasks, iEntity);
        for (IndexType i_sensor = 0; i_sensor < mSensorMaskStatuses.size2(); ++i_sensor) {
            mSensorMaskStatuses(iEntity, i_sensor) = r_values[i_sensor] * sensor_status[i_sensor];
        }
    });

    KRATOS_CATCH("");
}

// template instantiations
template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) class SensorMaskStatus<ModelPart::NodesContainerType>;
template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) class SensorMaskStatus<ModelPart::ConditionsContainerType>;
template KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) class SensorMaskStatus<ModelPart::ElementsContainerType>;

} // namespace Kratos