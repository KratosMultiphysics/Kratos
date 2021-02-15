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

// Project includes

// Application includes

// Include base h
#include "element_data_derivative_utilities.h"

namespace Kratos
{

namespace KEpsilonElementData
{

template <unsigned int TDim, unsigned int TNumNodes>
double AdjointUtilities<TDim, TNumNodes>::CalculateProductionVelocityDerivative(
    const IndexType NodeIndex,
    const IndexType DirectionIndex,
    const double ProductionTerm,
    const double TurbulentKinematicViscosity,
    const double TurbulentKinematicViscosityDerivative,
    const MatrixDD& rVelocityGradient,
    const Matrix& rdNdX)
{
    if (TurbulentKinematicViscosity != 0.0) {
        double value = 0.0;

        for (IndexType j = 0; j < TDim; ++j) {
            value += 2.0 * rdNdX(NodeIndex, j) * rVelocityGradient(DirectionIndex, j);
            value += 2.0 * rdNdX(NodeIndex, j) * rVelocityGradient(j, DirectionIndex);
            value -= rdNdX(NodeIndex, DirectionIndex) * (4.0 / 3.0) * rVelocityGradient(j, j);
        }
        return value * TurbulentKinematicViscosity + ProductionTerm *
               TurbulentKinematicViscosityDerivative / TurbulentKinematicViscosity;
    } else {
        return 0.0;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double AdjointUtilities<TDim, TNumNodes>::CalculateProductionScalarDerivative(
    const double TurbulentKinematicViscosity,
    const double ProductionTerm,
    const double GaussNutScalarDerivative)
{
    return GaussNutScalarDerivative * (ProductionTerm / TurbulentKinematicViscosity);
}

template <unsigned int TDim, unsigned int TNumNodes>
double AdjointUtilities<TDim, TNumNodes>::CalculateProductionShapeDerivative(
    const double TurbulentKinematicViscosity,
    const double TurbulentKinematicViscosityDerivative,
    const double ProductionTerm,
    const MatrixND& rNodalVelocity,
    const Matrix& rdNdX,
    const Matrix& rdNdXDerivative)
{
    double output = 0.0;

    const double negative_two_thirds = -2.0 / 3.0;

    for (std::size_t a = 0; a < TNumNodes; ++a)
    {
        const Vector& r_dna_dx_derivative = row(rdNdXDerivative, a);
        const Vector& r_velocity_a = row(rNodalVelocity, a);
        const Vector& r_dna_dx = row(rdNdX, a);
        for (std::size_t b = 0; b < TNumNodes; ++b)
        {
            const Vector& r_dnb_dx_derivative = row(rdNdXDerivative, b);
            const Vector& r_dnb_dx = row(rdNdX, b);
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
double AdjointUtilities<TDim, TNumNodes>::CalculateGammaKDerivative(
    const IndexType NodeIndex,
    const double Cmu,
    const double Gamma,
    const double TurbulentKinematicViscosity,
    const double GaussNutTKEDerivative,
    const Vector& rShapeFunctions)
{
    if (Gamma > 0.0) {
        double value = 0.0;
        value = rShapeFunctions[NodeIndex] * (Cmu / TurbulentKinematicViscosity);
        value -= GaussNutTKEDerivative * (Gamma / TurbulentKinematicViscosity);
        return value;
    } else {
        return 0.0;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double AdjointUtilities<TDim, TNumNodes>::CalculateGammaEpsilonDerivative(
    const double Gamma,
    const double TurbulentKinematicViscosity,
    const double GaussNutEpsilonDerivative)
{
    if (Gamma > 0.0) {
        return GaussNutEpsilonDerivative *
                           (-1.0 * Gamma / TurbulentKinematicViscosity);
    } else {
        return 0.0;
    }
}

// template instantiations
template class AdjointUtilities<2, 3>;
template class AdjointUtilities<2, 4>;

template class AdjointUtilities<3, 4>;
template class AdjointUtilities<3, 8>;

} // namespace KEpsilonElementData

} // namespace Kratos