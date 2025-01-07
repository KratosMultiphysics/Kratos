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
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi,
//                   Aron Noordam
//

// Application includes
#include "custom_conditions/Pw_normal_flux_condition.hpp"
#include "custom_utilities/condition_utilities.hpp"
#include <custom_utilities/variables_utilities.hpp>

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer PwNormalFluxCondition<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                  NodesArrayType const& ThisNodes,
                                                                  PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new PwNormalFluxCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

template <unsigned int TDim, unsigned int TNumNodes>
void PwNormalFluxCondition<TDim, TNumNodes>::CalculateRHS(VectorType&        rRightHandSideVector,
                                                          const ProcessInfo& CurrentProcessInfo)
{
    // Previous definitions
    const GeometryType&                             Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        Geom.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int NumGPoints = IntegrationPoints.size();
    const unsigned int local_dim  = Geom.LocalSpaceDimension();

    // Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues(this->GetIntegrationMethod());
    GeometryType::JacobiansType j_container(NumGPoints);
    for (auto& j : j_container) {
        j.resize(TDim, local_dim, false);
    }
    Geom.Jacobian(j_container, this->GetIntegrationMethod());

    // Condition variables
    Vector normal_flux_vector(TNumNodes);
    VariablesUtilities::GetNodalValues(Geom, NORMAL_FLUID_FLUX, normal_flux_vector.begin());

    NormalFluxVariables Variables;

    // Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        // Obtain Np
        noalias(Variables.Np) = row(NContainer, GPoint);

        // Interpolation of nodal normal flux to integration point normal flux.
        Variables.NormalFlux = MathUtils<>::Dot(Variables.Np, normal_flux_vector);

        // Compute weighting coefficient for integration
        Variables.IntegrationCoefficient = ConditionUtilities::CalculateIntegrationCoefficient<TDim, TNumNodes>(
            j_container[GPoint], IntegrationPoints[GPoint].Weight());

        // Contributions to the right hand side
        this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void PwNormalFluxCondition<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                                NormalFluxVariables& rVariables)
{
    noalias(rVariables.PVector) = -rVariables.NormalFlux * rVariables.Np * rVariables.IntegrationCoefficient;

    rRightHandSideVector += rVariables.PVector;
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string PwNormalFluxCondition<TDim, TNumNodes>::Info() const
{
    return "PwNormalFluxCondition";
}

template class PwNormalFluxCondition<2, 2>;
template class PwNormalFluxCondition<2, 3>;
template class PwNormalFluxCondition<2, 4>;
template class PwNormalFluxCondition<2, 5>;
template class PwNormalFluxCondition<3, 3>;
template class PwNormalFluxCondition<3, 4>;

} // Namespace Kratos.
