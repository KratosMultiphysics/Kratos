//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <variant>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_utilities/sensor_mask_status.h"
#include "custom_utilities/filtering/explicit_filter_utils.h"

namespace Kratos {
///@name Kratos Classes
///@{

class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) SensorResolutionMatrixResponseUtils
{
public:
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SensorResolutionMatrixResponseUtils);

    ///@}
    ///@name Life cycle
    ///@{

    SensorResolutionMatrixResponseUtils(
        SensorMaskStatus::Pointer pSensorMaskStatus,
        const double StepSize,
        const double FilterRadius,
        ModelPart& rModelPart,
        const std::string& rKernelFunctionType,
        const IndexType MaxLeafSize = 10,
        const IndexType EchoLevel = 0,
        const bool NodeCloudMesh = false,
        const bool StoreFilterMatrix = true);

    ///@}
    ///@name Public operations
    ///@{

    double CalculateValue();

    TensorAdaptor<double>::Pointer CalculateGradient() const;

    std::variant<
        ExplicitFilterUtils<ModelPart::NodesContainerType>::Pointer,
        ExplicitFilterUtils<ModelPart::ConditionsContainerType>::Pointer,
        ExplicitFilterUtils<ModelPart::ElementsContainerType>::Pointer> GetFilter();

    ///@}

private:
    ///@name Private member variables
    ///@{

    SensorMaskStatus::Pointer mpSensorMaskStatus;

    double mStepSize;

    std::variant<
        ExplicitFilterUtils<ModelPart::NodesContainerType>::Pointer,
        ExplicitFilterUtils<ModelPart::ConditionsContainerType>::Pointer,
        ExplicitFilterUtils<ModelPart::ElementsContainerType>::Pointer> mpFilter;

    Matrix mResolutionMatrix;

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/