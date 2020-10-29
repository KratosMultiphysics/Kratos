//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// Project includes

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"

// Include base h
#include "adjoint_element_data_utilities.h"

namespace Kratos
{
///@name  Functions
///@{

namespace KEpsilonElementData
{

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointUtilities<TDim, TNumNodes>::CalculateProductionVelocitySensitivities(
    BoundedVector<double, TNumNodes * TDim>& rOutput,
    const double TurbulentKinematicViscosity,
    const double ProductionTerm,
    const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
    const Matrix& rShapeDerivatives)
{
    rOutput.clear();

    double velocity_divergence = 0.0;
    velocity_divergence =
        RansCalculationUtilities::CalculateMatrixTrace<TDim>(rVelocityGradient);
    identity_matrix<double> identity(TDim);

    BoundedMatrix<double, TDim, TDim> reynolds_stress_tensor;
    noalias(reynolds_stress_tensor) = rVelocityGradient + trans(rVelocityGradient) -
                                      (2.0 / 3.0) * velocity_divergence * identity;

    for (IndexType c = 0; c < TNumNodes; ++c) {
        for (IndexType k = 0; k < TDim; ++k) {
            double value = 0.0;

            for (IndexType j = 0; j < TDim; ++j) {
                value += rShapeDerivatives(c, j) * rVelocityGradient(k, j);
                value += rShapeDerivatives(c, j) * rVelocityGradient(j, k);
                value -= rShapeDerivatives(c, k) * (2.0 / 3.0) * rVelocityGradient(j, j);
                value += rShapeDerivatives(c, j) * reynolds_stress_tensor(k, j);
            }

            rOutput[c * TDim + k] += TurbulentKinematicViscosity * value;
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointUtilities<TDim, TNumNodes>::CalculateProductionScalarSensitivities(
    BoundedVector<double, TNumNodes>& rOutput,
    const double TurbulentKinematicViscosity,
    const double ProductionTerm,
    const BoundedVector<double, TNumNodes>& rGaussNutScalarDerivatives)
{
    noalias(rOutput) =
        rGaussNutScalarDerivatives * (ProductionTerm / TurbulentKinematicViscosity);
}

template <unsigned int TDim, unsigned int TNumNodes>
double AdjointUtilities<TDim, TNumNodes>::CalculateProductionShapeSensitivities(
    const double TurbulentKinematicViscosity,
    const double TurbulentKinematicViscosityDerivative,
    const double ProductionTerm,
    const BoundedMatrix<double, TNumNodes, TDim>& rNodalVelocity,
    const Matrix& rShapeDerivatives,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_DxDerivatives)
{
    double output = 0.0;

    const double negative_two_thirds = -2.0 / 3.0;

    for (std::size_t a = 0; a < TNumNodes; ++a)
    {
        const Vector& r_dna_dx_derivative = row(rDN_DxDerivatives, a);
        const Vector& r_velocity_a = row(rNodalVelocity, a);
        const Vector& r_dna_dx = row(rShapeDerivatives, a);
        for (std::size_t b = 0; b < TNumNodes; ++b)
        {
            const Vector& r_dnb_dx_derivative = row(rDN_DxDerivatives, b);
            const Vector& r_dnb_dx = row(rShapeDerivatives, b);
            const Vector& r_velocity_b = row(rNodalVelocity, b);

            const double uai_ubi = inner_prod(r_velocity_a, r_velocity_b);

            output += inner_prod(r_dna_dx_derivative, r_dnb_dx) * uai_ubi;
            output += inner_prod(r_dna_dx, r_dnb_dx_derivative) * uai_ubi;
            output += inner_prod(r_velocity_a, r_dnb_dx_derivative) *
                      inner_prod(r_velocity_b, r_dna_dx);
            output += inner_prod(r_velocity_a, r_dnb_dx) *
                      inner_prod(r_velocity_b, r_dna_dx_derivative);
            output += negative_two_thirds * inner_prod(r_velocity_a, r_dna_dx_derivative) *
                      inner_prod(r_velocity_b, r_dnb_dx);
            output += negative_two_thirds * inner_prod(r_velocity_a, r_dna_dx) *
                      inner_prod(r_velocity_b, r_dnb_dx_derivative);
        }
    }

    output *= TurbulentKinematicViscosity;
    output += TurbulentKinematicViscosityDerivative * ProductionTerm / TurbulentKinematicViscosity;

    return output;
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointUtilities<TDim, TNumNodes>::CalculateNodalNutTKESensitivity(
    BoundedVector<double, TNumNodes>& rOutput,
    const double Cmu,
    const BoundedVector<double, TNumNodes>& rNodalTurbulentKineticEnergy,
    const BoundedVector<double, TNumNodes>& rNodalTurbulentEnergyDissipationRate)
{
    for (IndexType i = 0; i < TNumNodes; ++i) {
        if (rNodalTurbulentEnergyDissipationRate[i] > 0.0) {
            rOutput[i] = 2.0 * Cmu * rNodalTurbulentKineticEnergy[i] /
                         rNodalTurbulentEnergyDissipationRate[i];
        } else {
            rOutput[i] = 0.0;
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointUtilities<TDim, TNumNodes>::CalculateNodalNutEpsilonSensitivity(
    BoundedVector<double, TNumNodes>& rOutput,
    const double Cmu,
    const BoundedVector<double, TNumNodes>& rNodalTurbulentKineticEnergy,
    const BoundedVector<double, TNumNodes>& rNodalTurbulentEnergyDissipationRate)
{
    for (IndexType i = 0; i < TNumNodes; ++i) {
        if (rNodalTurbulentEnergyDissipationRate[i] > 0.0) {
            rOutput[i] = -1.0 * Cmu *
                         std::pow(rNodalTurbulentKineticEnergy[i] /
                                      rNodalTurbulentEnergyDissipationRate[i],
                                  2);
        } else {
            rOutput[i] = 0.0;
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointUtilities<TDim, TNumNodes>::CalculateGaussGammaTkeSensitivity(
    BoundedVector<double, TNumNodes>& rOutput,
    const double Cmu,
    const double Gamma,
    const double TurbulentKinematicViscosity,
    const BoundedVector<double, TNumNodes>& rGaussNutTKESensitivities,
    const Vector& rShapeFunctions)
{
    if (Gamma > 0.0) {
        noalias(rOutput) = rShapeFunctions * (Cmu / TurbulentKinematicViscosity);
        noalias(rOutput) -=
            rGaussNutTKESensitivities * (Gamma / TurbulentKinematicViscosity);
    } else {
        rOutput.clear();
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointUtilities<TDim, TNumNodes>::CalculateGaussGammaEpsilonSensitivity(
    BoundedVector<double, TNumNodes>& rOutput,
    const double Gamma,
    const double TurbulentKinematicViscosity,
    const BoundedVector<double, TNumNodes>& rGaussNutEpsilonSensitivities)
{
    if (Gamma > 0.0) {
        noalias(rOutput) = rGaussNutEpsilonSensitivities *
                           (-1.0 * Gamma / TurbulentKinematicViscosity);
    } else {
        rOutput.clear();
    }
}

// template instantiations
template class AdjointUtilities<2, 3>;
template class AdjointUtilities<3, 4>;

}

///@}

} // namespace Kratos