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

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class BoltzmannOperator
{
public:
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(BoltzmannOperator);

    using IndexType = std::size_t;

    ///@}
    ///@name Life cycle
    ///@{

    BoltzmannOperator(const double Beta);

    ///@}
    ///@name Public operations
    ///@{

    void CalculateCoefficients();

    Vector& GetData();

    const Vector& GetData() const;

    double GetValue() const;

    double GetGradient(const Vector& rGradients) const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const double mBeta;

    double mNumerator;

    Vector mValues;

    double mDenominator;

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/