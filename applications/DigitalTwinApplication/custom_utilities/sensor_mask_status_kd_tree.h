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
#include "flann/flann.hpp"

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "includes/model_part.h"
#include "includes/data_communicator.h"
#include "expression/container_expression.h"

// Application includes
#include "custom_utilities/sensor_mask_status.h"

namespace Kratos {
///@name Kratos Classes
///@{

template<class TContainerType>
class SensorMaskStatusKDTree {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(SensorMaskStatusKDTree);

    ///@}
    ///@name Life cycle
    ///@{

    SensorMaskStatusKDTree(
        typename SensorMaskStatus<TContainerType>::Pointer pSensorMaskStatus,
        const IndexType NumberOfParallelTrees,
        const bool ExactSearch,
        const IndexType EchoLevel);

    ///@}
    ///@name Public operations
    ///@{

    void GetEntitiesWithinRadius(
        std::vector<std::vector<long unsigned int>>& rIndices,
        std::vector<std::vector<double>>& rDistances,
        Matrix& rQueries,
        const double Radius);

    /**
     * @brief Updates the masks with the corresponding SENSOR_STATUS.
     */
    void Update();

    typename SensorMaskStatus<TContainerType>::Pointer GetSensorMaskStatus();

    ///@}

private:
    ///@name Private member variables
    ///@{

    typename SensorMaskStatus<TContainerType>::Pointer mpSensorMaskStatus;

    const bool mExactSearch;

    const IndexType mEchoLevel;

    flann::Matrix<double> mData;

    flann::KDTreeIndex<flann::L1<double>> mKDTree;

    ///@}
};

} // namespace Kratos