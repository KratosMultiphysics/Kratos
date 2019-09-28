#include "evm_k_epsilon_utilities.h"
#include <cmath>
#include <iostream>
#include <limits>

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

namespace EvmKepsilonModelUtilities
{
double CalculateTurbulentViscosity(const double C_mu,
                                   const double turbulent_kinetic_energy,
                                   const double turbulent_energy_dissipation_rate,
                                   const double f_mu)
{
    return C_mu * f_mu * std::pow(turbulent_kinetic_energy, 2) / turbulent_energy_dissipation_rate;
}

double CalculateFmu(const double y_plus)
{
    return 1.0 - std::exp(-0.0115 * y_plus);
}

double CalculateF2(const double turbulent_kinetic_energy,
                   const double kinematic_viscosity,
                   const double turbulent_energy_dissipation_rate)
{
    if (turbulent_energy_dissipation_rate == 0.0)
        return 1.0;
    else
    {
        const double Re_t = std::pow(turbulent_kinetic_energy, 2) /
                            (kinematic_viscosity * turbulent_energy_dissipation_rate);
        const double f2 = 1.0 - 0.22 * std::exp(-1.0 * std::pow(Re_t * (1.0 / 6.0), 2));

        return f2;
    }
}

template <unsigned int TDim>
double CalculateSourceTerm(const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
                           const double turbulent_kinematic_viscosity,
                           const double turbulent_kinetic_energy)
{
    const double velocity_divergence =
        RansCalculationUtilities().CalculateMatrixTrace<TDim>(rVelocityGradient);
    identity_matrix<double> identity(TDim);

    BoundedMatrix<double, TDim, TDim> symmetric_velocity_gradient;
    noalias(symmetric_velocity_gradient) = rVelocityGradient + trans(rVelocityGradient);

    BoundedMatrix<double, TDim, TDim> reynolds_stress_tensor;

    noalias(reynolds_stress_tensor) =
        turbulent_kinematic_viscosity *
        (symmetric_velocity_gradient - (2.0 / 3.0) * velocity_divergence * identity);

    double source = 0.0;
    for (unsigned int i = 0; i < TDim; ++i)
        for (unsigned int j = 0; j < TDim; ++j)
            source += reynolds_stress_tensor(i, j) * rVelocityGradient(i, j);

    return source;
}

double CalculateGamma(const double C_mu,
                      const double f_mu,
                      const double turbulent_kinetic_energy,
                      const double turbulent_kinematic_viscosity)
{
    return std::max<double>(
        0.0, C_mu * f_mu * turbulent_kinetic_energy / turbulent_kinematic_viscosity);
}

void CalculateTurbulentValues(double& turbulent_kinetic_energy,
                              double& turbulent_energy_dissipation_rate,
                              const double y_plus,
                              const double kinematic_viscosity,
                              const double wall_distance,
                              const double c_mu,
                              const double von_karman)
{
    const double u_tau = y_plus * kinematic_viscosity / wall_distance;
    turbulent_kinetic_energy = std::pow(u_tau, 2) / std::sqrt(c_mu);
    turbulent_energy_dissipation_rate = std::pow(u_tau, 3) / (von_karman * wall_distance);
}

void CalculateTurbulentValues(double& turbulent_kinetic_energy,
                              double& turbulent_energy_dissipation_rate,
                              const double velocity_mag,
                              const double turbulence_intensity,
                              const double mixing_length,
                              const double c_mu)
{
    turbulent_kinetic_energy = 1.5 * std::pow(velocity_mag * turbulence_intensity, 2);
    turbulent_energy_dissipation_rate =
        c_mu * std::pow(turbulent_kinetic_energy, 1.5) / mixing_length;
}

template double CalculateSourceTerm<2>(const BoundedMatrix<double, 2, 2>&, const double, const double);
template double CalculateSourceTerm<3>(const BoundedMatrix<double, 3, 3>&, const double, const double);

} // namespace EvmKepsilonModelUtilities

///@}

} // namespace Kratos
