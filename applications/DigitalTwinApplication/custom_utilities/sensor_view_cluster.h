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
#include <set>
#include <string>
#include <unordered_map>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/indexed_object.h"
#include "includes/model_part.h"

// Application includes
#include "custom_utilities/domain_sensor_view_cluster_data.h"

namespace Kratos {
///@name Kratos Classes
///@{

template<class TContainerType>
class SensorViewCluster: public IndexedObject
{
public:
    ///@name Type definitions
    ///@{

    using SensorViewPointerType = typename SensorView<TContainerType>::Pointer;

    using SensorViewVectorType = std::vector<SensorViewPointerType>;

    KRATOS_CLASS_POINTER_DEFINITION(SensorViewCluster);

    ///@}
    ///@name Life cycle
    ///@{

    SensorViewCluster(
        const IndexType ClusterId,
        typename DomainSensorViewClusterData<TContainerType>::Pointer pDataContainer);

    ///@}
    ///@name Public operations
    ///@{

    typename SensorViewCluster<TContainerType>::Pointer Clone() const;

    void SetSensorViews(const SensorViewVectorType& rSensorViews);

    SensorViewVectorType GetSensorViews() const;

    std::vector<double> GetDistances(const std::string& rDistancesName) const;

    void SetEntities(const TContainerType& rContainer) { mEntities = rContainer; }

    TContainerType& GetEntities() { return mEntities; }

    const TContainerType& GetEntities() const { return mEntities; }

    typename DomainSensorViewClusterData<TContainerType>::Pointer GetDataContainer() { return mpDataContainer; }

    void Clear();

    ///@}
private:
    ///@name Private member variables
    ///@{

    typename DomainSensorViewClusterData<TContainerType>::Pointer mpDataContainer;

    std::vector<IndexType> mSensorViewIndices;

    TContainerType mEntities;

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/