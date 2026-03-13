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
#include "tensor_adaptors/tensor_adaptor.h"

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

    TensorAdaptor<double>::Pointer CalculateGradient() const;

    void Update(const TensorAdaptor<double>& rInputTensorAdaptor);

    ///@}
private:

    ///@name Private member variables
    ///@{

    double mBeta;

    double mNumerator;

    double mDenominator;

    double mExtremeValue;

    TensorAdaptor<double>::Pointer mpTensorAdaptor;

    ///@}
};

///@}
} // namespace Kratos