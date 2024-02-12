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
#include "custom_conditions/one-phase_flow/U_Pl_face_load_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPlFaceLoadCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPlFaceLoadCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlFaceLoadCondition<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
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
    array_1d<double,TNumNodes*TDim> FaceLoadVector;
    PoroConditionUtilities::GetNodalVariableVector(FaceLoadVector,Geom,FACE_LOAD);
    BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
    array_1d<double,TDim> TractionVector;
    array_1d<double,TNumNodes*TDim> UVector;
    double IntegrationCoefficient;
    
    //Loop over integration points
    for(unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute traction vector 
        PoroConditionUtilities::InterpolateVariableWithComponents(TractionVector,NContainer,FaceLoadVector,GPoint);
        
        //Compute Nu Matrix
        PoroConditionUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);
        
        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(IntegrationCoefficient, JContainer[GPoint], integration_points[GPoint].Weight());
                
        //Contributions to the right hand side
        noalias(UVector) = prod(trans(Nu),TractionVector) * IntegrationCoefficient;
        PoroConditionUtilities::AssembleUBlockVector(rRightHandSideVector,UVector);
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPlFaceLoadCondition<2,2>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& Weight)
{
    double dx_dxi = Jacobian(0,0);
    double dy_dxi = Jacobian(1,0);

    double ds = sqrt(dx_dxi*dx_dxi + dy_dxi*dy_dxi);

    rIntegrationCoefficient = ds * Weight;
}

//----------------------------------------------------------------------------------------

template< >
void UPlFaceLoadCondition<3,3>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& Weight)
{
    double NormalVector[3];

    NormalVector[0] = Jacobian(1,0) * Jacobian(2,1) - Jacobian(2,0) * Jacobian(1,1);

    NormalVector[1] = Jacobian(2,0) * Jacobian(0,1) - Jacobian(0,0) * Jacobian(2,1);

    NormalVector[2] = Jacobian(0,0) * Jacobian(1,1) - Jacobian(1,0) * Jacobian(0,1);

    double dA = sqrt(NormalVector[0]*NormalVector[0] + NormalVector[1]*NormalVector[1] + NormalVector[2]*NormalVector[2]);

    rIntegrationCoefficient = dA * Weight;
}

//----------------------------------------------------------------------------------------

template< >
void UPlFaceLoadCondition<3,4>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& Weight)
{
    double NormalVector[3];

    NormalVector[0] = Jacobian(1,0) * Jacobian(2,1) - Jacobian(2,0) * Jacobian(1,1);

    NormalVector[1] = Jacobian(2,0) * Jacobian(0,1) - Jacobian(0,0) * Jacobian(2,1);

    NormalVector[2] = Jacobian(0,0) * Jacobian(1,1) - Jacobian(1,0) * Jacobian(0,1);

    double dA = sqrt(NormalVector[0]*NormalVector[0] + NormalVector[1]*NormalVector[1] + NormalVector[2]*NormalVector[2]);

    rIntegrationCoefficient = dA * Weight;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPlFaceLoadCondition<2,2>;
template class UPlFaceLoadCondition<3,3>;
template class UPlFaceLoadCondition<3,4>;

} // Namespace Kratos.
