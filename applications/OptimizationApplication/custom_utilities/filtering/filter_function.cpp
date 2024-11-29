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
#include <cmath>

// Project includes
#include "includes/global_variables.h"

// Include base h
#include "filter_function.h"

namespace Kratos
{

FilterFunction::FilterFunction(const std::string& rKernelFunctionType)
{
    // Set type of weighting function

    if (rKernelFunctionType == "gaussian") {
        // Type 1: Gaussian function
        // at distance = radius, the filter value will be 1e-8
        mFilterFunctional =  [](double radius, double distance) { return std::max<double>(0.0, exp(-(8 * 2.3025850929940455 * distance * distance) / (radius * radius))); };
    } else if (rKernelFunctionType == "linear") {
        // Type 2: Linear function
        mFilterFunctional =  [](double radius, double distance) { return std::max<double>(0.0, (radius - distance) / radius); };
    } else if (rKernelFunctionType == "constant") {
        // Type 3: Constant function
        mFilterFunctional = [](double radius, double distance) { return (distance < radius) ? 1.0: 0.0; };
    } else if (rKernelFunctionType == "cosine") {
        // Type 4: Cosine function
        mFilterFunctional = [](double radius, double distance) { return std::max<double>(0.0, 1 - 0.5 * (1 - std::cos(Globals::Pi / radius * std::min(distance, radius))));};
    } else if (rKernelFunctionType == "quartic") {
        // Type 5: Quartic function
        mFilterFunctional = [](double radius, double distance) { return std::max<double>(0.0, (std::pow(std::min(distance, radius) - radius, 4.0) / std::pow(radius, 4.0))); };
    } else if (rKernelFunctionType == "sigmoidal") {
        // Type 6: Sigmoidal function
        mFilterFunctional = [](double radius, double distance) {
                // in the following line, order of multiplication matters. If the max is used in between (not at the
                // last position of multiplication, then pow_val will be "nan", which is not recognised by the
                // following clamp method.)
                const double pow_val = 2.0 * (distance - radius) * std::numeric_limits<double>::max();
                const double exp_val = std::clamp(std::exp(pow_val), 0.0, std::numeric_limits<double>::max());
                return 1.0 / (1.0 + exp_val);
            };
    } else {
        // Throw error message in case of wrong specification
        KRATOS_ERROR << "Specified kernel function of type : "
                     << rKernelFunctionType << " is not recognized. \n \t Options are:"
                     << "\n\tconstant"
                     << "\n\tlinear"
                     << "\n\tgaussian"
                     << "\n\tcosine"
                     << "\n\tquartic"
                     << "\n\tsigmoidal\n";
    }
}

double FilterFunction::ComputeWeight(
    const double Radius,
    const double Distance) const
{
    KRATOS_TRY;

    // Depending on which weighting function is chosen, compute weight
    return mFilterFunctional(Radius, Distance);

    KRATOS_CATCH("");
}

} // namespace Kratos.
