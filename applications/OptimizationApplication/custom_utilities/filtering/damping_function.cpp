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
//

// System includes
#include <string>

// Project includes
#include "damping_function.h"

namespace Kratos
{

DampingFunction::DampingFunction(const std::string& rKernelFunctionType):
mKernelFunctionType(rKernelFunctionType)
{
    // Set type of weighting function

    std::function<double (const double, const double)> filter_function;

    if (mKernelFunctionType == "gaussian") {
        // Type 1: Gaussian function
        filter_function = [](double radius, double distance) {return std::max(0.0, exp(-(distance*distance)));};
    } else if (mKernelFunctionType == "linear") {
        // Type 2: Linear function
        filter_function =  [](double radius, double distance) {return std::max(0.0, (radius - distance) / radius);};
    } else if (mKernelFunctionType == "constant") {
        // Type 3: Constant function
        filter_function = [](double radius, double distance) {return 1.0;};
    } else if (mKernelFunctionType == "cosine") {
        // Type 4: Cosine function
        filter_function = [](double radius, double distance) {return std::max(0.0, 1-0.5*(1-std::cos(Globals::Pi/radius*distance)));};
    } else if (mKernelFunctionType == "quartic") {
        // Type 5: Quartic function
        filter_function = [](double radius, double distance) {return std::max(0.0, (pow(distance-radius,4.0)/pow(radius,4.0)));};
    } else if (mKernelFunctionType == "sigmoidal") {
        // Type 6: Sigmoidal function
        filter_function = [](double radius, double distance) {
            const double limit = std::log1p(std::numeric_limits<double>::max());
            double pow_val = -2.0 * std::numeric_limits<double>::max() * (distance - radius);
            pow_val = std::clamp(pow_val, -limit, limit);
            return 1.0 / (1.0 + std::exp(pow_val));};
    } else {
        // Throw error message in case of wrong specification
        KRATOS_ERROR << "Specified kernel function of type : "
                     << rKernelFunctionType << " is not recognized. \n \t Options are:"
                     << "\n\tconstant"
                     << "\n\tlinear"
                     << "\n\tgaussian"
                     << "\n\tcosine"
                     << "\n\tsigmoidal"
                     << "\n\tquartic.\n";
    }

    if (mKernelFunctionType == "sigmoidal") {
        mFilterFunctional = filter_function;
    } else {
        mFilterFunctional = [filter_function](const double radius, const double distance) {
            if (distance <= radius) {
                return 0.0;
            } else if (distance >= 2 * radius) {
                return 1.0;
            } else {
                return filter_function(radius, distance / 2.0);
            }
        };
    }
}

double DampingFunction::GetDistance(
    const Array3DType& ICoord,
    const Array3DType& JCoord) const
{
    const Array3DType dist_vector = ICoord - JCoord;
    return sqrt(dist_vector[0] * dist_vector[0] + dist_vector[1] * dist_vector[1] + dist_vector[2] * dist_vector[2]);
}

double DampingFunction::ComputeWeight(
    const Array3DType& ICoord,
    const Array3DType& JCoord,
    const double Radius) const
{
    KRATOS_TRY;

    // Compute distance vector
    const double distance = GetDistance(ICoord, JCoord);

    return mFilterFunctional(Radius, distance);

    KRATOS_CATCH("");
}

} // namespace Kratos.
