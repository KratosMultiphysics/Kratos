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
#include "custom_utilities/fluid_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "incompressible_potential_flow_velocity_element.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePotentialFlowVelocityElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == VELOCITY) {
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const IndexType num_gauss_points = gauss_weights.size();

        if (rValues.size() != num_gauss_points) {
            rValues.resize(num_gauss_points);
        }

        const auto& r_geometry = this->GetGeometry();

        for (IndexType g = 0; g < num_gauss_points; ++g) {
            const Matrix& r_shape_derivatives = shape_derivatives[g];

            array_1d<double, 3> velocity;
            FluidCalculationUtilities::EvaluateGradientInPoint(
                r_geometry, r_shape_derivatives,
                std::tie(velocity, VELOCITY_POTENTIAL));

            rValues[g] = velocity;
        }
    } else {
        KRATOS_ERROR << "CalculateOnIntegrationPoints for variable "
                     << rVariable.Name() << " not defined for " << this->Info();
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& IncompressiblePotentialFlowVelocityElement<TDim, TNumNodes>::GetVariable() const
{
    return VELOCITY_POTENTIAL;
}

// template instantiations
template class IncompressiblePotentialFlowVelocityElement<2, 3>;
template class IncompressiblePotentialFlowVelocityElement<3, 4>;

} // namespace Kratos.
