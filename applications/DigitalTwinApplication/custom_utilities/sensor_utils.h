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

#pragma once

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/data_communicator.h"
#include "expression/container_expression.h"

// Application includes
#include "custom_sensors/sensor_view.h"

namespace Kratos {
///@name Kratos Classes
///@{

class SensorUtils
{
public:
    ///@name Public static operations
    ///@{

    template<class TContainerType>
    static ModelPart& GetSensorViewsModelPart(const std::vector<typename SensorView<TContainerType>::Pointer>& rSensorViews);

    template<class TContainerType>
    static void IdentifyBestSensorViewForEveryEntity(
        std::vector<typename SensorView<TContainerType>::Pointer>& rOutputSensorViews,
        const std::vector<typename SensorView<TContainerType>::Pointer>& rNormalizedSensorViews);

    template<class TContainerType>
    static double GetDomainSize(
        const TContainerType& rContainer,
        const DataCommunicator& rDataCommunicator);

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/