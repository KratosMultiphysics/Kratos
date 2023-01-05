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
//                   
//                   
//


// Application includes
#include "custom_conditions/T_normal_flux_condition.hpp"

namespace Kratos
{

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer TNormalFluxCondition<TDim,TNumNodes>::Create(
    IndexType NewId,NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new TNormalFluxCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TNormalFluxCondition<TDim,TNumNodes>::CalculateRHS(
    VectorType& rRightHandSideVector,
    const ProcessInfo& CurrentProcessInfo)
{        
    //Previous definitions
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints( mThisIntegrationMethod );
    const unsigned int NumGPoints = IntegrationPoints.size();
    const unsigned int LocalDim = Geom.LocalSpaceDimension();
    
    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    GeometryType::JacobiansType JContainer(NumGPoints);
    for(unsigned int i = 0; i < NumGPoints; ++i)
        (JContainer[i]).resize(TDim,LocalDim,false);
    Geom.Jacobian( JContainer, mThisIntegrationMethod );
    
    //Condition variables
    array_1d<double,TNumNodes> NormalFluxVector;
    for(unsigned int i = 0; i < TNumNodes; ++i)
    {
        NormalFluxVector[i] = Geom[i].FastGetSolutionStepValue(NORMAL_HEAT_FLUX);
    }
    NormalFluxVariables Variables;
    
    //Loop over integration points
    for(unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        //Compute normal flux 
        Variables.NormalFlux = 0.0;
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            Variables.NormalFlux += NContainer(GPoint,i) * NormalFluxVector[i];
        }
        
        //Obtain Np
        noalias(Variables.Np) = row(NContainer,GPoint);
                
        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient,
                                              JContainer[GPoint],
                                              IntegrationPoints[GPoint].Weight() );
                
        //Contributions to the right hand side
        this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TNormalFluxCondition<TDim,TNumNodes>::CalculateAndAddRHS(
    VectorType& rRightHandSideVector,
    NormalFluxVariables& rVariables )
{
    noalias(rVariables.TVector) = rVariables.NormalFlux * rVariables.Np * rVariables.IntegrationCoefficient;

    GeoElementUtilities::
        AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector,
                                                rVariables.TVector);
}

// ============================================================================================
// ============================================================================================
template< >
void TNormalFluxCondition<2, 2>::CalculateIntegrationCoefficient(
    double& rIntegrationCoefficient,
    const Matrix& Jacobian,
    const double& Weight)
{
    const double dx_dxi = Jacobian(0, 0);
    const double dy_dxi = Jacobian(1, 0);

    const double ds = std::sqrt(dx_dxi * dx_dxi + dy_dxi * dy_dxi);

    rIntegrationCoefficient = ds * Weight;
}

// ============================================================================================
// ============================================================================================
template< >
void TNormalFluxCondition<2, 3>::CalculateIntegrationCoefficient(
    double& rIntegrationCoefficient,
    const Matrix& Jacobian,
    const double& Weight)
{
    const double dx_dxi = Jacobian(0, 0);
    const double dy_dxi = Jacobian(1, 0);

    const double ds = std::sqrt(dx_dxi * dx_dxi + dy_dxi * dy_dxi);

    rIntegrationCoefficient = ds * Weight;
}

// ============================================================================================
// ============================================================================================
template< >
void TNormalFluxCondition<3, 3>::CalculateIntegrationCoefficient(
    double& rIntegrationCoefficient,
    const Matrix& Jacobian,
    const double& Weight)
{
    double NormalVector[3];

    NormalVector[0] = Jacobian(1, 0) * Jacobian(2, 1) - Jacobian(2, 0) * Jacobian(1, 1);

    NormalVector[1] = Jacobian(2, 0) * Jacobian(0, 1) - Jacobian(0, 0) * Jacobian(2, 1);

    NormalVector[2] = Jacobian(0, 0) * Jacobian(1, 1) - Jacobian(1, 0) * Jacobian(0, 1);

    const double dA = std::sqrt(NormalVector[0] * NormalVector[0]
        + NormalVector[1] * NormalVector[1]
        + NormalVector[2] * NormalVector[2]);

    rIntegrationCoefficient = dA * Weight;
}

// ============================================================================================
// ============================================================================================
template< >
void TNormalFluxCondition<3, 4>::CalculateIntegrationCoefficient(
    double& rIntegrationCoefficient,
    const Matrix& Jacobian,
    const double& Weight)
{
    double NormalVector[3];

    NormalVector[0] = Jacobian(1, 0) * Jacobian(2, 1) - Jacobian(2, 0) * Jacobian(1, 1);

    NormalVector[1] = Jacobian(2, 0) * Jacobian(0, 1) - Jacobian(0, 0) * Jacobian(2, 1);

    NormalVector[2] = Jacobian(0, 0) * Jacobian(1, 1) - Jacobian(1, 0) * Jacobian(0, 1);

    const double dA = std::sqrt(NormalVector[0] * NormalVector[0]
        + NormalVector[1] * NormalVector[1]
        + NormalVector[2] * NormalVector[2]);

    rIntegrationCoefficient = dA * Weight;
}

// ============================================================================================
// ============================================================================================
template class TNormalFluxCondition<2,2>;
template class TNormalFluxCondition<2,3>;
template class TNormalFluxCondition<3,3>;
template class TNormalFluxCondition<3,4>;

} // Namespace Kratos.
