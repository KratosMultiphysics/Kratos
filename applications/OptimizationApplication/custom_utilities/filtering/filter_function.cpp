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

// System includes
#include <string>

// Project includes
#include "filter_function.h"

namespace Kratos
{

FilterFunction::FilterFunction(const std::string& rKernelFunctionType)
{
    // Set type of weighting function

    if (rKernelFunctionType == "gaussian") {
        // Type 1: Gaussian function
        // at distance = radius, the filter value will be 1e-8
        mFilterFunctional =  [](double radius, double distance) { return std::max(0.0, exp(-(8 * 2.3025850929940455 * distance * distance) / (radius * radius))); };
    } else if (rKernelFunctionType == "linear") {
        // Type 2: Linear function
        mFilterFunctional =  [](double radius, double distance) { return std::max(0.0, (radius - distance) / radius); };
    } else if (rKernelFunctionType == "constant") {
        // Type 3: Constant function
        mFilterFunctional = [](double radius, double distance) { return (distance < radius) ? 1.0: 0.0; };
    } else if (rKernelFunctionType == "cosine") {
        // Type 4: Cosine function
        mFilterFunctional = [](double radius, double distance) { return std::max(0.0, 1 - 0.5 * (1 - std::cos(Globals::Pi / radius * std::min(distance, radius))));};
    } else if (rKernelFunctionType == "quartic") {
        // Type 5: Quartic function
        mFilterFunctional = [](double radius, double distance) { return std::max(0.0, (std::pow(std::min(distance, radius) - radius, 4.0) / std::pow(radius, 4.0))); };
    } else {
        // Throw error message in case of wrong specification
        KRATOS_ERROR << "Specified kernel function of type : "
                     << rKernelFunctionType << " is not recognized. \n \t Options are:"
                     << "\n\tconstant"
                     << "\n\tlinear"
                     << "\n\tgaussian"
                     << "\n\tcosine"
                     << "\n\tquartic.\n";
    }
}

double FilterFunction::GetDistance(
    const Array3DType& ICoord,
    const Array3DType& JCoord) const
{
    const Array3DType dist_vector = ICoord - JCoord;
    return sqrt(dist_vector[0] * dist_vector[0] + dist_vector[1] * dist_vector[1] + dist_vector[2] * dist_vector[2]);
}

double FilterFunction::ComputeWeight(
    const Array3DType& ICoord,
    const Array3DType& JCoord,
    const double Radius) const
{
    KRATOS_TRY;

    // Compute distance vector
    const double distance = GetDistance(ICoord, JCoord);

    // Depending on which weighting function is chosen, compute weight
    return mFilterFunctional(Radius, distance);

    KRATOS_CATCH("");
}

} // namespace Kratos.
