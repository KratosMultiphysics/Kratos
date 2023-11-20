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
    const GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        rGeom.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int NumGPoints = IntegrationPoints.size();
    const unsigned int LocalDim = rGeom.LocalSpaceDimension();

    const Matrix& NContainer = rGeom.ShapeFunctionsValues(this->GetIntegrationMethod());
    GeometryType::JacobiansType JContainer(NumGPoints);

    for (auto& x : JContainer) {
        x.resize(TDim, LocalDim, false);
    }

    rGeom.Jacobian(JContainer, this->GetIntegrationMethod());

    array_1d<double, TNumNodes> normalFluxVector;
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        normalFluxVector[i] = rGeom[i].FastGetSolutionStepValue(NORMAL_HEAT_FLUX);
    }
    NormalFluxVariables variables;

    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        // Obtain N
        noalias(variables.N) = row(NContainer, GPoint);

        // Compute normal flux
        variables.normalFlux = MathUtils<>::Dot(variables.N, normalFluxVector);

        // Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(variables.IntegrationCoefficient,
                                              JContainer[GPoint],
                                              IntegrationPoints[GPoint].Weight());

        // Contributions to the right hand side
        this->CalculateAndAddRHS(rRightHandSideVector, variables);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void TNormalFluxCondition<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                               NormalFluxVariables& rVariables)
{
    noalias(rVariables.fluxVector) =
        rVariables.normalFlux * rVariables.N * rVariables.IntegrationCoefficient;

    GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(
        rRightHandSideVector, rVariables.fluxVector);
}

template <unsigned int TDim, unsigned int TNumNodes>
void TNormalFluxCondition<TDim, TNumNodes>::CalculateIntegrationCoefficient(
    double& rIntegrationCoefficient, const Matrix& Jacobian, const double& Weight)
{
    Vector normalVector = ZeroVector(TDim);

    if constexpr (TDim == 2) {
        normalVector = column(Jacobian, 0);
    }
    else if constexpr (TDim == 3) {
        MathUtils<double>::CrossProduct(normalVector, column(Jacobian, 0),
                                        column(Jacobian, 1));
    }
    rIntegrationCoefficient = Weight * MathUtils<double>::Norm(normalVector);
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
