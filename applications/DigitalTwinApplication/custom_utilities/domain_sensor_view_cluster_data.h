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
#include <string>
#include <unordered_map>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "custom_sensors/sensor_view.h"

namespace Kratos {
///@name Kratos Classes
///@{

template<class TContainerType>
class DomainSensorViewClusterData
{
public:
    ///@name Type definitions
    ///@{

    using SensorViewPointerType = typename SensorView<TContainerType>::Pointer;

    using SensorViewVectorType = std::vector<SensorViewPointerType>;

    KRATOS_CLASS_POINTER_DEFINITION(DomainSensorViewClusterData);

    ///@}
    ///@name Life cycle
    ///@{

    DomainSensorViewClusterData(
        const SensorViewVectorType& rSensorViewsList,
        const SensorViewVectorType& rRepresentativeSensorViewsForEntitiesList);

    ///@}
    ///@name Public operations
    ///@{

    void AddDistances(
        const std::string& rDistancesName,
        const std::vector<double>& rDistances);


    TContainerType& GetContainer() const;

    std::vector<double> GetDistancesForIndices(
        const std::vector<IndexType>& rIndices,
        const std::string& rDistancesName) const;

    SensorViewVectorType GetSensorViewsForIndices(const std::vector<IndexType>& rIndices) const;

    SensorViewVectorType& GetSensorViews() { return mSensorViewPointersList; };

    SensorViewVectorType& GetRepresentativeSensorViewsForEntities() { return mRepresentativeSensorsForEntities; };

    ModelPart& GetModelPart() const { return *mpModelPart; };

    ///@}

private:
    ///@name Private member variables
    ///@{

    ModelPart * const mpModelPart;

    SensorViewVectorType mSensorViewPointersList;

    SensorViewVectorType mRepresentativeSensorsForEntities;

    std::unordered_map<std::string, std::vector<double>> mDistancesMap;

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/