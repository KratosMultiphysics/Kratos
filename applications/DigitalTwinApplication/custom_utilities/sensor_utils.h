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
#include "custom_utilities/sensor_view_cluster.h"


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

    template<class TContainerType>
    static void AssignEntitiesToClustersBasedOnOptimalSensor(
        const std::vector<typename SensorViewCluster<TContainerType>::Pointer>& rSensorViewClusers,
        const std::vector<typename ContainerExpression<TContainerType>::Pointer>& rContainerExpressionsList);

    template<class TContainerType>
    static void GetThresholdSensorViews(
        const double ThresholdValue,
        const std::string& rExpressionName,
        std::vector<typename SensorView<TContainerType>::Pointer>& rOutputSensorViews,
        const std::vector<typename SensorView<TContainerType>::Pointer>& rNormalizedSensorViews);
    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/