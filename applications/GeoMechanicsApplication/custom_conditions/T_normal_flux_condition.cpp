// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Mohamed Nabi
//                   John van Esch
//

// Application includes
#include "custom_conditions/T_normal_flux_condition.h"
#include "utilities/math_utils.h"

namespace Kratos {

template <unsigned int TDim, unsigned int TNumNodes>
GeoTNormalFluxCondition<TDim, TNumNodes>::GeoTNormalFluxCondition()
    : GeoTCondition<TDim, TNumNodes>()
{
}

template <unsigned int TDim, unsigned int TNumNodes>
GeoTNormalFluxCondition<TDim, TNumNodes>::GeoTNormalFluxCondition(IndexType NewId,
                                                            GeometryType::Pointer pGeometry)
    : GeoTCondition<TDim, TNumNodes>(NewId, pGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
GeoTNormalFluxCondition<TDim, TNumNodes>::GeoTNormalFluxCondition(IndexType NewId,
                                                            GeometryType::Pointer pGeometry,
                                                            PropertiesType::Pointer pProperties)
    : GeoTCondition<TDim, TNumNodes>(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
GeoTNormalFluxCondition<TDim, TNumNodes>::~GeoTNormalFluxCondition() = default;

template <unsigned int TDim, unsigned int TNumNodes>
void GeoTNormalFluxCondition<TDim, TNumNodes>::CalculateRHS(VectorType& rRightHandSideVector,
                                                         const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& r_geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points =
        r_geom.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int num_integration_points = r_integration_points.size();
    const unsigned int local_dim = r_geom.LocalSpaceDimension();

    const Matrix& r_N_container =
        r_geom.ShapeFunctionsValues(this->GetIntegrationMethod());
    GeometryType::JacobiansType j_container(num_integration_points);

    for (auto& x : j_container) {
        x.resize(TDim, local_dim, false);
    }

    r_geom.Jacobian(j_container, this->GetIntegrationMethod());

    array_1d<double, TNumNodes> normal_flux_vector;
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        normal_flux_vector[i] = r_geom[i].FastGetSolutionStepValue(NORMAL_HEAT_FLUX);
    }
    NormalFluxVariables variables;

    for (unsigned int integration_point = 0;
         integration_point < num_integration_points; ++integration_point) {
        // Obtain N
        noalias(variables.N) = row(r_N_container, integration_point);

        // Interpolation of nodal normal flux to integration point normal flux.
        variables.NormalFlux = MathUtils<>::Dot(variables.N, normal_flux_vector);

        // Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(
            variables.IntegrationCoefficient, j_container[integration_point],
            r_integration_points[integration_point].Weight());

        // Contributions to the right hand side
        this->CalculateAndAddRHS(rRightHandSideVector, variables);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoTNormalFluxCondition<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                               NormalFluxVariables& rVariables)
{
    noalias(rVariables.FluxVector) =
        rVariables.NormalFlux * rVariables.N * rVariables.IntegrationCoefficient;

    GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(
        rRightHandSideVector, rVariables.FluxVector);
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoTNormalFluxCondition<TDim, TNumNodes>::CalculateIntegrationCoefficient(
    double& rIntegrationCoefficient, const Matrix& rJacobian, double Weight)
{
    Vector normal_vector = ZeroVector(TDim);

    if constexpr (TDim == 2) {
        normal_vector = column(rJacobian, 0);
    }
    else if constexpr (TDim == 3) {
        MathUtils<double>::CrossProduct(normal_vector, column(rJacobian, 0),
                                        column(rJacobian, 1));
    }
    rIntegrationCoefficient = Weight * MathUtils<double>::Norm(normal_vector);
}

template class GeoTNormalFluxCondition<2, 2>;
template class GeoTNormalFluxCondition<2, 3>;
template class GeoTNormalFluxCondition<2, 4>;
template class GeoTNormalFluxCondition<2, 5>;
template class GeoTNormalFluxCondition<3, 3>;
template class GeoTNormalFluxCondition<3, 4>;
template class GeoTNormalFluxCondition<3, 6>;
template class GeoTNormalFluxCondition<3, 8>;
template class GeoTNormalFluxCondition<3, 9>;

} // namespace Kratos
