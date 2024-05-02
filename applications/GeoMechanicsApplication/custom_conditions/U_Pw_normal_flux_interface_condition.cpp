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
//                   Vahid Galavi
//


// Application includes
#include "custom_conditions/U_Pw_normal_flux_interface_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPwNormalFluxInterfaceCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPwNormalFluxInterfaceCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwNormalFluxInterfaceCondition<TDim,TNumNodes>::
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
    array_1d<double,TNumNodes*TDim> DisplacementVector;
    ConditionUtilities::GetDisplacementsVector(DisplacementVector,Geom);
    array_1d<double,TNumNodes> NormalFluxVector;
    for(unsigned int i=0; i<TNumNodes; ++i)
    {
        NormalFluxVector[i] = Geom[i].FastGetSolutionStepValue(NORMAL_FLUID_FLUX);
    }
    BoundedMatrix<double,TDim,TDim> RotationMatrix;
    const double& MinimumJointWidth = this->GetProperties()[MINIMUM_JOINT_WIDTH];
    bool ComputeJointWidth;
    double JointWidth;
    this->CheckJointWidth(JointWidth,ComputeJointWidth,RotationMatrix,MinimumJointWidth,Geom);
    BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
    array_1d<double,TDim> LocalRelDispVector;
    array_1d<double,TDim> RelDispVector;
    array_1d<double,TNumNodes> Np;
    array_1d<double,TNumNodes> PVector;
    double NormalFlux;
    
    //Loop over integration points
    for(unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute normal flux
        NormalFlux = 0.0;
        for(unsigned int i=0; i<TNumNodes; ++i)
        {
            NormalFlux += NContainer(GPoint,i)*NormalFluxVector[i];
        }
        
        //Obtain Np
        noalias(Np) = row(NContainer,GPoint);
        
        if(ComputeJointWidth==true)
        {
            //Compute Nu Matrix and Joint Width
            InterfaceElementUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);
            this->CalculateJointWidth(JointWidth, Nu, DisplacementVector, RelDispVector, RotationMatrix, LocalRelDispVector, MinimumJointWidth,GPoint);
        }
        
        //Compute weighting coefficient for integration
        double IntegrationCoefficient = this->CalculateIntegrationCoefficient(JContainer[GPoint], IntegrationPoints[GPoint].Weight(), JointWidth );
                
        //Contributions to the right hand side
        noalias(PVector) = -NormalFlux * Np * IntegrationCoefficient;
        GeoElementUtilities::AssemblePBlockVector(rRightHandSideVector,PVector);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPwNormalFluxInterfaceCondition<2,2>;
template class UPwNormalFluxInterfaceCondition<3,4>;

} // Namespace Kratos.
