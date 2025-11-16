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

// System includes
#include <cmath>

// Project includes

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "element_data_utilities.h"

namespace Kratos
{
namespace KEpsilonElementData
{

double CalculateTurbulentViscosity(
    const Geometry<Node>& rGeometry,
    const Vector& rN,
    const double Cmu)
{
    // we calculate tke and epsilon based on the values
    // stored in non-historical nodal container
    // historical container values are transferred to non-historical
    // container after each coupling step of two equation turbulence model
    // this is done to achieve convergence easily.
    double tke, epsilon;
    FluidCalculationUtilities::EvaluateNonHistoricalInPoint(
        rGeometry, rN,
        std::tie(tke, TURBULENT_KINETIC_ENERGY),
        std::tie(epsilon, TURBULENT_ENERGY_DISSIPATION_RATE));

    if (epsilon > 0.0) {
        return std::max(Cmu * std::pow(tke, 2) / epsilon, 1e-12);
    } else {
        return 1e-12;
    }
}

template <unsigned int TDim>
double CalculateProductionTerm(
    const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
    const double TurbulentKinematicViscosity)
{
    const double velocity_divergence =
        RansCalculationUtilities::CalculateMatrixTrace<TDim>(rVelocityGradient);
    identity_matrix<double> identity(TDim);

    BoundedMatrix<double, TDim, TDim> symmetric_velocity_gradient;
    noalias(symmetric_velocity_gradient) = rVelocityGradient + trans(rVelocityGradient);

    BoundedMatrix<double, TDim, TDim> reynolds_stress_tensor;

    noalias(reynolds_stress_tensor) =
        TurbulentKinematicViscosity *
        (symmetric_velocity_gradient - (2.0 / 3.0) * velocity_divergence * identity);

    double source = 0.0;
    for (unsigned int i = 0; i < TDim; ++i) {
        for (unsigned int j = 0; j < TDim; ++j) {
            source += reynolds_stress_tensor(i, j) * rVelocityGradient(i, j);
        }
    }

    return source;
}

double CalculateGamma(
    const double Cmu,
    const double TurbulentKineticEnergy,
    const double TurbulentKinematicViscosity)
{
    return std::max(Cmu * TurbulentKineticEnergy / TurbulentKinematicViscosity, 0.0);
}

// template instantiations

template double CalculateProductionTerm<2>(const BoundedMatrix<double, 2, 2>&, const double);
template double CalculateProductionTerm<3>(const BoundedMatrix<double, 3, 3>&, const double);

} // namespace KEpsilonElementData

} // namespace Kratos
