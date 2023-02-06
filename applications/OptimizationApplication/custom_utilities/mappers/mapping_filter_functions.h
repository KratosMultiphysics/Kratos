//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//                   Aditya Ghantasala, https://github.com/adityaghantasala
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
 * MappingFilterFunction implementations.
*/

class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) MappingFilterFunction
{
  public:
    ///@name Type Definitions
    ///@{

    using Array3DType = array_1d<double, 3>;

    /// Pointer definition of MappingFilterFunction
    KRATOS_CLASS_POINTER_DEFINITION(MappingFilterFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    MappingFilterFunction(const std::string& MappingFilterFunctionType);

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

    std::function<double (const double, const double)> mMappingFilterFunctional;

    ///@}
    ///@name Private Operations
    ///@{

    double inline GetDistance(
        const Array3DType& ICoord,
        const Array3DType& JCoord) const;

    ///@}

}; // Class MappingFilterFunction

///@}

} // namespace Kratos.

