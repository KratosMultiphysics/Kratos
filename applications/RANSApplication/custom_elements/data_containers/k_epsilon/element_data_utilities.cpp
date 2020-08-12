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
#include "custom_utilities/rans_calculation_utilities.h"

// Include base h
#include "element_data_utilities.h"

namespace Kratos
{
namespace KEpsilonElementData
{
double CalculateTurbulentViscosity(
    const double Cmu,
    const double TurbulentKineticEnergy,
    const double TurbulentEnergyDissipationRate)
{
    return Cmu * std::pow(TurbulentKineticEnergy, 2) / TurbulentEnergyDissipationRate;
}

template <unsigned int TDim>
double CalculateSourceTerm(
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

template double CalculateSourceTerm<2>(const BoundedMatrix<double, 2, 2>&, const double);
template double CalculateSourceTerm<3>(const BoundedMatrix<double, 3, 3>&, const double);

} // namespace KEpsilonElementData

} // namespace Kratos
