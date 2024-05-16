//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "expression/container_expression.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"

// Application includes
#include "custom_utilities/sensor_mask_status_kd_tree.h"

namespace Kratos {
///@name Kratos Classes
///@{

class SensorLocalizationResponseUtils
{
public:
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SensorLocalizationResponseUtils);

    ///@}
    ///@name Life cycle
    ///@{

    SensorLocalizationResponseUtils(
        SensorMaskStatusKDTree<ModelPart::ElementsContainerType>::Pointer pSensorMaskKDTree,
        const double Beta,
        const double AllowedDissimilarity);

    ///@}
    ///@name Public operations
    ///@{

    double CalculateValue();

    ContainerExpression<ModelPart::NodesContainerType> CalculateGradient() const;

    ContainerExpression<ModelPart::ElementsContainerType> GetClusterSizes() const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    SensorMaskStatusKDTree<ModelPart::ElementsContainerType>::Pointer mpSensorMaskStatusKDTree;

    const double mBeta;

    const double mAllowedDissimilarity;

    double mNumerator;

    double mDenominator;

    double mValue;

    std::vector<double> mDomainSizeRatio;

    std::vector<double> mClusterSizes;

    std::vector<std::vector<long unsigned int>> mNeighbourIndices;

    std::vector<std::vector<double>> mNeighbourSquareDistances;

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/