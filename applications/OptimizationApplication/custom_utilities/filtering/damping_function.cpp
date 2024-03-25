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

DampingFunction::DampingFunction(
    const std::string& rKernelFunctionType,
    const double DampingDistanceMultiplier)
    : FilterFunction(rKernelFunctionType),
      mDampingDistanceMultiplier(DampingDistanceMultiplier)
{
}

double DampingFunction::ComputeWeight(
    const double Radius,
    const double Distance) const
{
    KRATOS_TRY;

    if (Distance <= Radius) {
        return 0.0;
    } else {
        return 1.0 - mFilterFunctional(Radius, (Distance - Radius) * mDampingDistanceMultiplier);
    }

    KRATOS_CATCH("");
}

} // namespace Kratos.
