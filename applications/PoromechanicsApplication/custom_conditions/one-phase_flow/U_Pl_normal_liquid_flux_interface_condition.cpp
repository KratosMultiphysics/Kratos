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
//


// Application includes
#include "custom_conditions/one-phase_flow/U_Pl_normal_liquid_flux_interface_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPlNormalLiquidFluxInterfaceCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPlNormalLiquidFluxInterfaceCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlNormalLiquidFluxInterfaceCondition<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
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
    array_1d<double,TNumNodes*TDim> DisplacementVector;
    PoroConditionUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);
    array_1d<double,TNumNodes> NormalLiquidFluxVector;
    for(unsigned int i=0; i<TNumNodes; i++)
    {
        // Multiplied by -1.0 to indicate that positive value = inlet
        NormalLiquidFluxVector[i] = -1.0*Geom[i].FastGetSolutionStepValue(NORMAL_LIQUID_FLUX);
    }
    BoundedMatrix<double,TDim,TDim> RotationMatrix;
    const double& InitialJointWidth = this->GetProperties()[INITIAL_JOINT_WIDTH];
    bool ComputeJointWidth;
    double JointWidth;
    this->CheckJointWidth(JointWidth,ComputeJointWidth,RotationMatrix,InitialJointWidth,Geom);
    BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
    array_1d<double,TDim> LocalRelDispVector;
    array_1d<double,TDim> RelDispVector;
    array_1d<double,TNumNodes> Np;
    array_1d<double,TNumNodes> PVector;
    double NormalLiquidFlux;
    double IntegrationCoefficient;
    
    //Loop over integration points
    for(unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute normal flux
        NormalLiquidFlux = 0.0;
        for(unsigned int i=0; i<TNumNodes; i++)
        {
            NormalLiquidFlux += NContainer(GPoint,i)*NormalLiquidFluxVector[i];
        }
        
        //Obtain Np
        noalias(Np) = row(NContainer,GPoint);
        
        if(ComputeJointWidth==true)
        {
            //Compute Nu Matrix and Joint Width
            InterfaceElementUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);
            this->CalculateJointWidth(JointWidth, Nu, DisplacementVector, RelDispVector, RotationMatrix, LocalRelDispVector, InitialJointWidth,GPoint);
        }
        
        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(IntegrationCoefficient, JContainer[GPoint], integration_points[GPoint].Weight(), JointWidth );
                
        //Contributions to the right hand side
        noalias(PVector) = -NormalLiquidFlux * Np * IntegrationCoefficient;
        PoroElementUtilities::AssemblePBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,PVector,TDim,TNumNodes);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPlNormalLiquidFluxInterfaceCondition<2,2>;
template class UPlNormalLiquidFluxInterfaceCondition<3,4>;

} // Namespace Kratos.
