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
#include "includes/model_part.h"
#include "containers/nd_data.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) BoltzmannOperator
{
public:
    ///@name Life cycle
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(BoltzmannOperator);

    BoltzmannOperator(const double Beta);

    ///@}
    ///@name Public operations
    ///@{

    double CalculateValue() const;

    NDData<double>::Pointer CalculateGradient() const;

    void Update(NDData<double>::Pointer pNDData);

    ///@}
private:

    ///@name Private member variables
    ///@{

    double mBeta;

    double mNumerator;

    double mDenominator;

    double mExtremeValue;

    NDData<double>::Pointer mpNDData;

    ///@}
};

///@}
} // namespace Kratos