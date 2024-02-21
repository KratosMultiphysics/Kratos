//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//                   Danilo Cavalcanti
//


// Application includes
#include "custom_conditions/two-phase_flow/U_Pl_Pg_normal_liquid_flux_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPlPgNormalLiquidFluxCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPlPgNormalLiquidFluxCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgNormalLiquidFluxCondition<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{        
    //Previous definitions
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
    const unsigned int NumGPoints = integration_points.size();
    const unsigned int LocalDim = Geom.LocalSpaceDimension();
    
    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    GeometryType::JacobiansType JContainer(NumGPoints);
    for(unsigned int i = 0; i<NumGPoints; i++)
        (JContainer[i]).resize(TDim,LocalDim,false);
    Geom.Jacobian( JContainer, mThisIntegrationMethod );
    
    //Condition variables
    array_1d<double,TNumNodes> NormalLiquidFluxVector;
    this->GetNormalFluidFluxVector(NormalLiquidFluxVector,Geom);
    NormalFluidFluxVariables Variables;
    
    //Loop over integration points
    for(unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute normal flux 
        Variables.NormalFluidFlux = 0.0;
        for(unsigned int i=0; i<TNumNodes; i++)
        {
            Variables.NormalFluidFlux += NContainer(GPoint,i)*NormalLiquidFluxVector[i];
        }
        
        //Obtain Np
        noalias(Variables.Np) = row(NContainer,GPoint);
                
        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, JContainer[GPoint], integration_points[GPoint].Weight() );
                
        //Contributions to the right hand side
        this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgNormalLiquidFluxCondition<TDim,TNumNodes>::GetNormalFluidFluxVector(array_1d<double,TNumNodes>& rNormalFluidFluxVector, const GeometryType& Geom)
{
    for(unsigned int i=0; i<TNumNodes; i++)
    {
        rNormalFluidFluxVector[i] = Geom[i].FastGetSolutionStepValue(NORMAL_LIQUID_FLUX);
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgNormalLiquidFluxCondition<TDim,TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector, NormalFluidFluxVariables& rVariables)
{
    noalias(rVariables.PVector) = -rVariables.NormalFluidFlux * rVariables.Np * rVariables.IntegrationCoefficient;
    
    PoroElementUtilities::AssemblePlBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPlPgNormalLiquidFluxCondition<2,2>;
template class UPlPgNormalLiquidFluxCondition<3,3>;
template class UPlPgNormalLiquidFluxCondition<3,4>;

} // Namespace Kratos.
