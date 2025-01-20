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

class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) SensorMaskStatus {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using MasksListType = std::variant<
                                std::vector<ContainerExpression<ModelPart::NodesContainerType>::Pointer>,
                                std::vector<ContainerExpression<ModelPart::ConditionsContainerType>::Pointer>,
                                std::vector<ContainerExpression<ModelPart::ElementsContainerType>::Pointer>
                            >;

    using MaskContainerPointerType = std::variant<
                                        ModelPart::NodesContainerType::Pointer,
                                        ModelPart::ConditionsContainerType::Pointer,
                                        ModelPart::ElementsContainerType::Pointer
                                    >;

    KRATOS_CLASS_POINTER_DEFINITION(SensorMaskStatus);

    ///@}
    ///@name Life cycle
    ///@{

    SensorMaskStatus(
        ModelPart& rSensorModelPart,
        const MasksListType& rMaskPointersList,
        const IndexType EchoLevel);

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
    const Matrix& GetMaskStatuses() const;

    /**
     * @brief Get the masks matrix.
     * @details This method returns the masks. The rows of the return matrix returns number of entities
     *          and columns of it represents the number of sensors.
     * @return const Matrix&        Masks matrix.
     */
    const Matrix& GetMasks() const;

    /**
     * @brief Get the Sensor Model Part object
     */
    const ModelPart& GetSensorModelPart() const;

    /**
     * @brief Get the mask model part
     */
    ModelPart* pGetMaskModelPart() const;

    /**
     * @brief Get the Sensor Model Part object
     */
    ModelPart* pGetSensorModelPart() const;

    /**
     * @brief Returns a pointer to the local container of the masks.
     */
    MaskContainerPointerType pGetMaskContainer() const;

    /**
     * @brief Updates the masks with the corresponding SENSOR_STATUS.
     */
    void Update();

    ///@}

private:
    ///@name Private member variables
    ///@{

    ModelPart * const mpSensorModelPart;

    const MasksListType mMaskPointersList;

    const IndexType mEchoLevel;

    Matrix mSensorMaskStatuses;

    Matrix mSensorMasks;

    ///@}
};

} // namespace Kratos