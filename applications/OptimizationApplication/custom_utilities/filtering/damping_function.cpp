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

// System includes

// Project includes
#include "damping_function.h"

namespace Kratos
{

DampingFunction::DampingFunction(const std::string& rKernelFunctionType)
    : FilterFunction(rKernelFunctionType)
{
}

double DampingFunction::ComputeWeight(
    const Array3DType& ICoord,
    const Array3DType& JCoord,
    const double Radius) const
{
    KRATOS_TRY;

    // Compute distance vector
    const double distance = GetDistance(ICoord, JCoord);

    if (distance <= Radius) {
        return 0.0;
    } else {
        return 1.0 - mFilterFunctional(Radius, distance - Radius);
    }

    KRATOS_CATCH("");
}

} // namespace Kratos.
