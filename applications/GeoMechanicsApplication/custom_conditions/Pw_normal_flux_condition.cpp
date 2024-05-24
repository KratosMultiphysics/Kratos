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

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer PwNormalFluxCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new PwNormalFluxCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void PwNormalFluxCondition<TDim,TNumNodes>::
    CalculateRHS(VectorType& rRightHandSideVector,
                 const ProcessInfo& CurrentProcessInfo)
{        
    //Previous definitions
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int NumGPoints = IntegrationPoints.size();
    const unsigned int LocalDim = Geom.LocalSpaceDimension();
    
    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues(this->GetIntegrationMethod());
    GeometryType::JacobiansType JContainer(NumGPoints);
    for(unsigned int i = 0; i<NumGPoints; ++i)
        (JContainer[i]).resize(TDim,LocalDim,false);
    Geom.Jacobian(JContainer, this->GetIntegrationMethod());
    
    //Condition variables
    array_1d<double,TNumNodes> NormalFluxVector;
    for(unsigned int i=0; i<TNumNodes; ++i)
    {
        NormalFluxVector[i] = Geom[i].FastGetSolutionStepValue(NORMAL_FLUID_FLUX);
    }
    NormalFluxVariables Variables;
    
    //Loop over integration points
    for(unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        //Compute normal flux 
        Variables.NormalFlux = 0.0;
        for (unsigned int i=0; i<TNumNodes; ++i) {
            Variables.NormalFlux += NContainer(GPoint,i)*NormalFluxVector[i];
        }
        
        //Obtain Np
        noalias(Variables.Np) = row(NContainer,GPoint);
                
        //Compute weighting coefficient for integration
        Variables.IntegrationCoefficient = 
            ConditionUtilities::CalculateIntegrationCoefficient<TDim, TNumNodes>(
            JContainer[GPoint], IntegrationPoints[GPoint].Weight());
                
        //Contributions to the right hand side
        this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void PwNormalFluxCondition<TDim,TNumNodes>::
    CalculateAndAddRHS( VectorType& rRightHandSideVector,
                        NormalFluxVariables& rVariables )
{
    noalias(rVariables.PVector) = - rVariables.NormalFlux * rVariables.Np * rVariables.IntegrationCoefficient;

    rRightHandSideVector += rVariables.PVector;
}


template class PwNormalFluxCondition<2,2>;
template class PwNormalFluxCondition<2,3>;
template class PwNormalFluxCondition<2,4>;
template class PwNormalFluxCondition<2,5>;
template class PwNormalFluxCondition<3,3>;
template class PwNormalFluxCondition<3,4>;

} // Namespace Kratos.
