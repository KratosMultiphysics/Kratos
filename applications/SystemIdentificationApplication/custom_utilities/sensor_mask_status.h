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

#pragma once

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "includes/model_part.h"
#include "includes/data_communicator.h"
#include "expression/container_expression.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

template<class TContainerType>
class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) SensorMaskStatus {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(SensorMaskStatus);

    ///@}
    ///@name Life cycle
    ///@{

    SensorMaskStatus(
        ModelPart& rSensorModelPart,
        const std::vector<ContainerExpression<TContainerType>>& rMasksList);

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Get the mask statuses.
     * @details This method returns the masks multiplied with the corresponding SENSOR_STATUS.
     *          The rows of the return matrix returns number of entities, and columns of the return
     *          matrix represents number of sensors.
     * @return const Matrix&        Mask statuses matrix.
     */
    Matrix& GetMaskStatuses();

    /**
     * @brief Get the mask statuses.
     * @details This method returns the masks multiplied with the corresponding SENSOR_STATUS.
     *          The rows of the return matrix returns number of entities, and columns of the return
     *          matrix represents number of sensors.
     * @return const Matrix&        Mask statuses matrix.
     */
    const Matrix& GetMaskStatuses() const;

    /**
     * @brief Get the masks matrix.
     * @details This method returns the masks. The rows of the return matrix returns number of entities
     *          and columns of it represents the number of sensors.
     * @return const Matrix&        Masks matrix.
     */
    const Matrix& GetMasks() const;

    /**
     * @brief Get the Sensor Model Part
     */
    ModelPart& GetSensorModelPart() const;

    /**
     * @brief Get local container for all the masks
     */
    const TContainerType& GetMaskLocalContainer() const;

    ModelPart& GetMaskModelPart();

    /**
     * @brief Get the Data Communicator used in the mask model part.
     */
    const DataCommunicator& GetDataCommunicator() const;

    /**
     * @brief Updates the masks with the corresponding SENSOR_STATUS.
     */
    void Update();

    ///@}

private:
    ///@name Private member variables
    ///@{

    ModelPart * const mpSensorModelPart;

    ModelPart * mpMaskModelPart;

    TContainerType const * mpMaskContainer;

    DataCommunicator const * mpDataCommunicator;

    Matrix mSensorMaskStatuses;

    Matrix mSensorMasks;

    ///@}
};

} // namespace Kratos