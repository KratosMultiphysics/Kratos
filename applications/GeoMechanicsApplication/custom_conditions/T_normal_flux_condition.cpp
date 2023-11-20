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
TNormalFluxCondition<TDim, TNumNodes>::TNormalFluxCondition()
    : TCondition<TDim, TNumNodes>()
{
}

template <unsigned int TDim, unsigned int TNumNodes>
TNormalFluxCondition<TDim, TNumNodes>::TNormalFluxCondition(IndexType NewId,
                                                            GeometryType::Pointer pGeometry)
    : TCondition<TDim, TNumNodes>(NewId, pGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
TNormalFluxCondition<TDim, TNumNodes>::TNormalFluxCondition(IndexType NewId,
                                                            GeometryType::Pointer pGeometry,
                                                            PropertiesType::Pointer pProperties)
    : TCondition<TDim, TNumNodes>(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
TNormalFluxCondition<TDim, TNumNodes>::~TNormalFluxCondition() = default;

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer TNormalFluxCondition<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new TNormalFluxCondition(
        NewId, this->GetGeometry().Create(rThisNodes), pProperties));
}

template <unsigned int TDim, unsigned int TNumNodes>
void TNormalFluxCondition<TDim, TNumNodes>::CalculateRHS(VectorType& rRightHandSideVector,
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

        // Compute normal flux
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
void TNormalFluxCondition<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                               NormalFluxVariables& rVariables)
{
    noalias(rVariables.FluxVector) =
        rVariables.NormalFlux * rVariables.N * rVariables.IntegrationCoefficient;

    GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(
        rRightHandSideVector, rVariables.FluxVector);
}

template <unsigned int TDim, unsigned int TNumNodes>
void TNormalFluxCondition<TDim, TNumNodes>::CalculateIntegrationCoefficient(
    double& rIntegrationCoefficient, const Matrix& rJacobian, const double& Weight)
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

template class TNormalFluxCondition<2, 2>;
template class TNormalFluxCondition<2, 3>;
template class TNormalFluxCondition<2, 4>;
template class TNormalFluxCondition<2, 5>;
template class TNormalFluxCondition<3, 3>;
template class TNormalFluxCondition<3, 4>;
template class TNormalFluxCondition<3, 6>;
template class TNormalFluxCondition<3, 8>;
template class TNormalFluxCondition<3, 9>;

} // namespace Kratos
