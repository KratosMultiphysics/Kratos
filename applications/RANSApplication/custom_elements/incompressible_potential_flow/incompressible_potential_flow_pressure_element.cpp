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

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/variables.h"

// Application includes
#include "rans_application_variables.h"
#include "custom_utilities/rans_calculation_utilities.h"

// Include base h
#include "incompressible_potential_flow_pressure_element.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowPressureElement<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Calculate RHS
    this->CalculateRightHandSideVelocityContribution(rRightHandSideVector);

    // Calculate LHS
    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    BoundedVector<double, TNumNodes> values;
    this->GetValuesArray(values);

    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, values);
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowPressureElement<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    this->CalculateRightHandSideVelocityContribution(rRightHandSideVector);

    BoundedMatrix<double, TNumNodes, TNumNodes> lhs;
    this->CalculateBoundedLeftHandSide(lhs, rCurrentProcessInfo);

    BoundedVector<double, TNumNodes> values;
    this->GetValuesArray(values);

    noalias(rRightHandSideVector) -= prod(lhs, values);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
int IncompressiblePotentialFlowPressureElement<TDim, TNumNodes>::Check(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int check = BaseType::Check(rCurrentProcessInfo);

    for (IndexType i_node = 0; i_node < this->GetGeometry().size(); ++i_node) {
        const auto& r_node = this->GetGeometry()[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY_POTENTIAL, r_node);
    }

    return check;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowPressureElement<TDim, TNumNodes>::CalculateRightHandSideVelocityContribution(
    Vector& rRHS) const
{
    if (rRHS.size() != TNumNodes) {
        rRHS.resize(TNumNodes, false);
    }
    noalias(rRHS) = ZeroVector(TNumNodes);

    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
    const IndexType num_gauss_points = gauss_weights.size();

    const auto& r_geometry = this->GetGeometry();

    for (IndexType g = 0; g < num_gauss_points; ++g) {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector& gauss_shape_functions = row(shape_functions, g);

        const double density = RansCalculationUtilities::EvaluateInPoint(
            r_geometry, DENSITY, gauss_shape_functions);

        array_1d<double, 3> kinetic_energy_gradient;
        // VELOCITY_POTENTIAL contains velocity magnitude square
        RansCalculationUtilities::CalculateGradient(
            kinetic_energy_gradient, r_geometry, VELOCITY_POTENTIAL, r_shape_derivatives);
        noalias(kinetic_energy_gradient) =
            kinetic_energy_gradient * (gauss_weights[g] * 0.5 * density);

        for (IndexType a = 0; a < TNumNodes; ++a) {
            double value = 0.0;
            for (IndexType d = 0; d < TDim; ++d) {
                value += r_shape_derivatives(a, d) * kinetic_energy_gradient[d];
            }

            rRHS[a] -= value;
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& IncompressiblePotentialFlowPressureElement<TDim, TNumNodes>::GetVariable() const
{
    return PRESSURE_POTENTIAL;
}

// template instantiations
template class IncompressiblePotentialFlowPressureElement<2, 3>;
template class IncompressiblePotentialFlowPressureElement<3, 4>;

} // namespace Kratos.
