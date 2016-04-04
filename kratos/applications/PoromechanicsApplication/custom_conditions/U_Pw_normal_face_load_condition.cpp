//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_conditions/U_Pw_normal_face_load_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPwNormalFaceLoadCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPwNormalFaceLoadCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwNormalFaceLoadCondition<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
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
    NormalFaceLoadVariables Variables;
    this->InitializeConditionVariables(Variables,Geom);
    array_1d<double,TDim> TractionVector;
    boost::numeric::ublas::bounded_matrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
    double IntegrationCoefficient;
    array_1d<double,TNumNodes*TDim> UVector;
    
    //Loop over integration points
    for(unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute traction vector
        this->CalculateTractionVector(TractionVector,JContainer[GPoint],NContainer,Variables,GPoint);
        
        //Compute Nu Matrix
        ConditionUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);
        
        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(IntegrationCoefficient, integration_points[GPoint].Weight());
                
        //Contributions to the right hand side
        noalias(UVector) = prod(trans(Nu),TractionVector) * IntegrationCoefficient;
        ConditionUtilities::AssembleUBlockVector(rRightHandSideVector,UVector);
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPwNormalFaceLoadCondition<2,2>::InitializeConditionVariables(NormalFaceLoadVariables& rVariables, const GeometryType& Geom)
{
    for(unsigned int i=0; i<2; i++)
    {
        rVariables.NormalStressVector[i] = Geom[i].FastGetSolutionStepValue(NORMAL_CONTACT_STRESS);
        rVariables.TangentialStressVector[i] = Geom[i].FastGetSolutionStepValue(TANGENTIAL_CONTACT_STRESS);
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPwNormalFaceLoadCondition<3,3>::InitializeConditionVariables(NormalFaceLoadVariables& rVariables, const GeometryType& Geom)
{
    for(unsigned int i=0; i<3; i++)
    {
        rVariables.NormalStressVector[i] = Geom[i].FastGetSolutionStepValue(NORMAL_CONTACT_STRESS);
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPwNormalFaceLoadCondition<3,4>::InitializeConditionVariables(NormalFaceLoadVariables& rVariables, const GeometryType& Geom)
{
    for(unsigned int i=0; i<4; i++)
    {
        rVariables.NormalStressVector[i] = Geom[i].FastGetSolutionStepValue(NORMAL_CONTACT_STRESS);
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPwNormalFaceLoadCondition<2,2>::CalculateTractionVector(array_1d<double,2>& rTractionVector,const Matrix& Jacobian,const Matrix& NContainer,
                                                                const NormalFaceLoadVariables& Variables,const unsigned int& GPoint)
{
    double NormalStress = 0.0;
    double TangentialStress = 0.0;
    for(unsigned int i=0; i<2; i++)
    {
        NormalStress += NContainer(GPoint,i)*Variables.NormalStressVector[i];
        TangentialStress += NContainer(GPoint,i)*Variables.TangentialStressVector[i];
    }
    
    double dx_dxi = Jacobian(0,0);
    double dy_dxi = Jacobian(1,0);
    
    rTractionVector[0] = TangentialStress * dx_dxi - NormalStress * dy_dxi;
    rTractionVector[1] = NormalStress * dx_dxi + TangentialStress * dy_dxi;
}

//----------------------------------------------------------------------------------------

template< >
void UPwNormalFaceLoadCondition<3,3>::CalculateTractionVector(array_1d<double,3>& rTractionVector,const Matrix& Jacobian,const Matrix& NContainer,
                                                                const NormalFaceLoadVariables& Variables,const unsigned int& GPoint)
{
    double NormalStress = 0.0;

    for(unsigned int i=0; i<3; i++)
    {
        NormalStress += NContainer(GPoint,i)*Variables.NormalStressVector[i];
    }
        
    double NormalVector[3];

    NormalVector[0] = Jacobian(1,0) * Jacobian(2,1) - Jacobian(2,0) * Jacobian(1,1);

    NormalVector[1] = Jacobian(2,0) * Jacobian(0,1) - Jacobian(0,0) * Jacobian(2,1);

    NormalVector[2] = Jacobian(0,0) * Jacobian(1,1) - Jacobian(1,0) * Jacobian(0,1);
    
    rTractionVector[0] = NormalStress * NormalVector[0];
    rTractionVector[1] = NormalStress * NormalVector[1];
    rTractionVector[2] = NormalStress * NormalVector[2];
}

//----------------------------------------------------------------------------------------

template< >
void UPwNormalFaceLoadCondition<3,4>::CalculateTractionVector(array_1d<double,3>& rTractionVector,const Matrix& Jacobian,const Matrix& NContainer,
                                                                const NormalFaceLoadVariables& Variables,const unsigned int& GPoint)
{
    double NormalStress = 0.0;

    for(unsigned int i=0; i<4; i++)
    {
        NormalStress += NContainer(GPoint,i)*Variables.NormalStressVector[i];
    }
        
    double NormalVector[3];

    NormalVector[0] = Jacobian(1,0) * Jacobian(2,1) - Jacobian(2,0) * Jacobian(1,1);

    NormalVector[1] = Jacobian(2,0) * Jacobian(0,1) - Jacobian(0,0) * Jacobian(2,1);

    NormalVector[2] = Jacobian(0,0) * Jacobian(1,1) - Jacobian(1,0) * Jacobian(0,1);
    
    rTractionVector[0] = NormalStress * NormalVector[0];
    rTractionVector[1] = NormalStress * NormalVector[1];
    rTractionVector[2] = NormalStress * NormalVector[2];
}

//----------------------------------------------------------------------------------------

template< >
void UPwNormalFaceLoadCondition<2,2>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& Weight)
{
    rIntegrationCoefficient = this->GetProperties()[THICKNESS] * Weight;
}

//----------------------------------------------------------------------------------------

template< >
void UPwNormalFaceLoadCondition<3,3>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& Weight)
{
    rIntegrationCoefficient = Weight;
}

//----------------------------------------------------------------------------------------

template< >
void UPwNormalFaceLoadCondition<3,4>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& Weight)
{
    rIntegrationCoefficient = Weight;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPwNormalFaceLoadCondition<2,2>;
template class UPwNormalFaceLoadCondition<3,3>;
template class UPwNormalFaceLoadCondition<3,4>;

} // Namespace Kratos.
