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

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_utilities/boltzmann_operator.h"
#include "custom_utilities/sensor_mask_status_kd_tree.h"

namespace Kratos {
///@name Kratos Classes
///@{

class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) SensorLocalizationResponseUtils
{
public:
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SensorLocalizationResponseUtils);

    ///@}
    ///@name Life cycle
    ///@{

    SensorLocalizationResponseUtils(
        SensorMaskStatusKDTree::Pointer pSensorMaskKDTree,
        const double BoltzmannBeta,
        const double SigmoidalBeta,
        const double PenaltyFactor,
        const double InitialDissimilarityMultiplier,
        const double DissimilarityDecayingFactor,
        const IndexType DissimilarityDecayingPeriod,
        const double AllowedDissimilarity);

    ///@}
    ///@name Public operations
    ///@{

    double CalculateValue(const IndexType Step);

    TensorAdaptor<double>::Pointer CalculateGradient(const IndexType Step) const;

    TensorAdaptor<double>::Pointer GetClusterSizes() const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    SensorMaskStatusKDTree::Pointer mpSensorMaskStatusKDTree;

    BoltzmannOperator mBoltzmannOperator;

    const double mSigmoidalBeta;

    const double mPenaltyFactor;

    const double mInitialDissimilarityMultiplier;

    const double mDissimilarityDecayingFactor;

    const IndexType mDissimilarityDecayingPeriod;

    const double mAllowedDissimilarity;

    std::vector<double> mDomainSizeRatio;

    std::vector<double> mClusterSizeRatios;

    std::vector<std::vector<SensorMaskStatusKDTree::ResultType>> mNeighbourData;

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/