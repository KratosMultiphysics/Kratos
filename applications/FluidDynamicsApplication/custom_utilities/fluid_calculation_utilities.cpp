//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// External includes

// Project includes

// Include base h
#include "custom_utilities/fluid_calculation_utilities.h"

namespace Kratos
{

template<>
KRATOS_API(FLUID_DYNAMICS_APPLICATION)
void FluidCalculationUtilities::AssignValue(
    const array_1d<double, 3>& rInput,
    array_1d<double, 2>& rOutput)
{
    rOutput[0] = rInput[0];
    rOutput[1] = rInput[1];
}

template<class TOutputDataType, class TInputDataType>
void FluidCalculationUtilities::AssignValue(
    const TInputDataType& rInput,
    TOutputDataType& rOutput)
{
    rOutput = rInput;
}

template<>
KRATOS_API(FLUID_DYNAMICS_APPLICATION)
void FluidCalculationUtilities::UpdateValue(
    const double& rInput,
    double& rOutput)
{
    rOutput += rInput;
}

template<>
KRATOS_API(FLUID_DYNAMICS_APPLICATION)
void FluidCalculationUtilities::UpdateValue(
    const array_1d<double, 3>& rInput,
    array_1d<double, 2>& rOutput)
{
    rOutput[0] += rInput[0];
    rOutput[1] += rInput[1];
}

template<class TOutputDataType, class TInputDataType>
void FluidCalculationUtilities::UpdateValue(
    const TInputDataType& rInput,
    TOutputDataType& rOutput)
{
    noalias(rOutput) += rInput;
}

double FluidCalculationUtilities::CalculateLogarithmicYPlusLimit(
    const double Kappa,
    const double Beta,
    const int MaxIterations,
    const double Tolerance)
{
    double y_plus = 11.06;
    const double inv_kappa = 1.0 / Kappa;
    double dx = 0.0;
    for (int i = 0; i < MaxIterations; ++i) {
        const double value = inv_kappa * std::log(y_plus) + Beta;
        dx = value - y_plus;

        if (std::abs(dx) < Tolerance) {
            return y_plus;
        }

        y_plus = value;
    }

    KRATOS_WARNING("LogarithmicYPlusLimit")
        << "Logarithmic y_plus limit reached max iterations with dx > "
           "Tolerance [ "
        << dx << " > " << Tolerance << ", MaxIterations = " << MaxIterations << " ].\n";
    return y_plus;
}

double FluidCalculationUtilities::CalculateLogarithmicYPlus(
    const double WallVelocityMagnitude,
    const double WallHeight,
    const double KinematicViscosity,
    const double Kappa,
    const double Beta,
    const double YPlusLimit,
    const int MaxIterations,
    const double Tolerance)
{
    // linear region
    double u_tau = std::sqrt(WallVelocityMagnitude * KinematicViscosity / WallHeight);
    double y_plus = u_tau * WallHeight / KinematicViscosity;
    const double inv_kappa = 1.0 / Kappa;

    // log region
    if (y_plus > YPlusLimit) {
        int iter = 0;
        double dx = 1e10;
        double u_plus = inv_kappa * std::log(y_plus) + Beta;

        while (iter < MaxIterations && std::abs(dx) > Tolerance * u_tau) {
            // Newton-Raphson iteration
            double f = u_tau * u_plus - WallVelocityMagnitude;
            double df = u_plus + inv_kappa;
            dx = f / df;

            // Update variables
            u_tau -= dx;
            y_plus = WallHeight * u_tau / KinematicViscosity;
            u_plus = inv_kappa * std::log(y_plus) + Beta;
            ++iter;
        }
        if (iter == MaxIterations) {
            std::cout << "Warning: wall condition Newton-Raphson did not "
                         "converge. Residual is "
                      << dx << std::endl;
        }
    }

    return y_plus;
}

// template instantiations
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidCalculationUtilities::AssignValue<double>(const double&, double&);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidCalculationUtilities::AssignValue<array_1d<double, 3>>(const array_1d<double, 3>&, array_1d<double, 3>&);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidCalculationUtilities::AssignValue<Vector>(const Vector&, Vector&);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidCalculationUtilities::AssignValue<Matrix>(const Matrix&, Matrix&);

template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidCalculationUtilities::UpdateValue<array_1d<double, 3>>(const array_1d<double, 3>&, array_1d<double, 3>&);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidCalculationUtilities::UpdateValue<Vector>(const Vector&, Vector&);
template KRATOS_API(FLUID_DYNAMICS_APPLICATION) void FluidCalculationUtilities::UpdateValue<Matrix>(const Matrix&, Matrix&);

} // namespace Kratos
