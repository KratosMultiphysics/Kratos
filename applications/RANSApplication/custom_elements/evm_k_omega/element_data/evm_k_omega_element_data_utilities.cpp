//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Dharmin Shah (https://github.com/sdharmin)
//                   Bence Rochlitz (https://github.com/bencerochlitz)
//
//  Supervised by:   Jordi Cotela (https://github.com/jcotela)
//                   Suneth Warnakulasuriya (https://github.com/sunethwarna)

// System includes
#include <cmath>
#include <iostream>
#include <limits>

// Project includes
#include "containers/array_1d.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"

// Include base h
#include "evm_k_omega_element_data_utilities.h"

namespace Kratos
{
namespace EvmKOmegaElementDataUtilities
{
double CalculateTurbulentViscosity(const double TurbulentKineticEnergy,
                                   const double TurbulentSpecificEnergyDissipationRate)
{
    return TurbulentKineticEnergy / TurbulentSpecificEnergyDissipationRate;
}

// we need to calculate our own parameters, which are: Beta and Beta_Star --> parameters of the CDR equations.
template <unsigned int TDim>
double CalculateBeta(const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
                     const double VelocityDivergence,
                     const double TurbulentSpecificDissipationRate,
                     const double BetaZero,
                     const double BetaStar) // omega element
{
    BoundedMatrix<double, TDim, TDim> symmetric_velocity_gradient;
    noalias(symmetric_velocity_gradient) =
        0.5 * (rVelocityGradient + trans(rVelocityGradient));

    BoundedMatrix<double, TDim, TDim> antisymmetric_velocity_gradient;
    noalias(antisymmetric_velocity_gradient) = rVelocityGradient - symmetric_velocity_gradient;

    IdentityMatrix identity(TDim);
    BoundedMatrix<double, TDim, TDim> s_tilde;
    noalias(s_tilde) = symmetric_velocity_gradient - identity * (0.5 * VelocityDivergence);

    BoundedMatrix<double, TDim, TDim> result;
    noalias(result) = prod(antisymmetric_velocity_gradient, antisymmetric_velocity_gradient);

    double x_omega = 0.0;
    for (IndexType i = 0; i < TDim; ++i)
    {
        for (IndexType j = 0; j < TDim; ++j)
        {
            x_omega += result(i, j) * s_tilde(j, i);
        }
    }

    x_omega = std::abs(x_omega / std::pow(TurbulentSpecificDissipationRate * BetaStar, 3));

    const double f_beta = (1 + 85 * x_omega) / (1 + 100 * x_omega);

    return BetaZero * f_beta;
}

double CalculateSigmaD(const array_1d<double, 3>& rTurbulentKineticEnergyGradient,
                       const array_1d<double, 3>& rTurbulentSpecificEnergyDissipationRateGradient)
{
    const double value = inner_prod(rTurbulentKineticEnergyGradient,
                                    rTurbulentSpecificEnergyDissipationRateGradient);
    if (value > 0.0)
    {
        return 1.0 / 8.0;
    }
    else
    {
        return 0.0;
    }
}

template double CalculateBeta<2>(const BoundedMatrix<double, 2, 2>&,
                                 const double,
                                 const double,
                                 const double,
                                 const double);

template double CalculateBeta<3>(const BoundedMatrix<double, 3, 3>&,
                                 const double,
                                 const double,
                                 const double,
                                 const double); // omega element

} // namespace EvmKOmegaElementDataUtilities

///@}

} // namespace Kratos