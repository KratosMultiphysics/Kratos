//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel
//                   Aditya Ghantasala
//                   Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <string>
#include <functional>

// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/**
 * FilterFunction implementations.
*/

class KRATOS_API(OPTIMIZATION_APPLICATION) FilterFunction
{
  public:
    ///@name Type Definitions
    ///@{

    using Array3DType = array_1d<double, 3>;

    /// Pointer definition of FilterFunction
    KRATOS_CLASS_POINTER_DEFINITION(FilterFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    FilterFunction(const std::string& rKernelFunctionType);

    ///@}
    ///@name Operations
    ///@{

    double ComputeWeight(
        const Array3DType& ICoord,
        const Array3DType& JCoord,
        const double Radius) const;

    ///@}

  private:
    ///@name Member Variables
    ///@{

    std::function<double (const double, const double)> mFilterFunctional;

    ///@}
    ///@name Private Operations
    ///@{

    double inline GetDistance(
        const Array3DType& ICoord,
        const Array3DType& JCoord) const;

    ///@}

}; // Class FilterFunction

///@}

} // namespace Kratos.

