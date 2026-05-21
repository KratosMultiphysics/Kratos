//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//                   Suneth Warnakulasuriya
//

#pragma once

// System includes

// Project includes
#include "includes/define.h"

// Application includes
#include "filter_function.h"

namespace Kratos {

///@name Kratos Classes
///@{

/// Short class definition.
/**
 * DampingFunction implementations.
 */

class KRATOS_API(OPTIMIZATION_APPLICATION) DampingFunction : public FilterFunction {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DampingFunction
    KRATOS_CLASS_POINTER_DEFINITION(DampingFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DampingFunction(
        const std::string& rKernelFunctionType,
        const double DampingDistanceMultiplier);

    ///@}
    ///@name Operations
    ///@{

    double ComputeWeight(
        const double Radius,
        const double Distance) const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const double mDampingDistanceMultiplier;

    ///@}

}; // Class DampingFunction

///@}

} // namespace Kratos.
