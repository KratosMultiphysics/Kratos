//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Lorenzo Gracia
//


// Application includes
#include "custom_conditions/UP_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------

template< >
void UPCondition<2,2>::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    const unsigned int condition_size = 2 * (2 + 1);
    unsigned int index = 0;
       
    if (rConditionDofList.size() != condition_size)
      rConditionDofList.resize( condition_size );
    
    for (unsigned int i = 0; i < 2; i++)
    {
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_X);
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Y);
        rConditionDofList[index++] = rGeom[i].pGetDof(PRESSURE);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPCondition<3,3>::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 3 * (3 + 1);
    unsigned int index = 0;
    
    if (rConditionDofList.size() != condition_size)
      rConditionDofList.resize( condition_size );
    
    for (unsigned int i = 0; i < 3; i++)
    {
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_X);
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Y);
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Z);
        rConditionDofList[index++] = rGeom[i].pGetDof(PRESSURE);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPCondition<3,4>::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 4 * (3 + 1);
    unsigned int index = 0;
    
    if (rConditionDofList.size() != condition_size)
      rConditionDofList.resize( condition_size );
    
    for (unsigned int i = 0; i < 4; i++)
    {
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_X);
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Y);
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Z);
        rConditionDofList[index++] = rGeom[i].pGetDof(PRESSURE);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPCondition<TDim,TNumNodes>::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int condition_size = TNumNodes * (TDim + 1);
    
    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != condition_size )
        rLeftHandSideMatrix.resize( condition_size, condition_size, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( condition_size, condition_size );
       
    //Resetting the RHS
    if ( rRightHandSideVector.size() != condition_size )
        rRightHandSideVector.resize( condition_size, false );
    noalias( rRightHandSideVector ) = ZeroVector( condition_size );
      
    this->CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
        
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPCondition<TDim,TNumNodes>::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    unsigned int condition_size = TNumNodes * (TDim + 1);
    
    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != condition_size )
        rLeftHandSideMatrix.resize( condition_size, condition_size, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( condition_size, condition_size );
    
    this->CalculateLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
    
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPCondition<TDim,TNumNodes>::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    
    unsigned int condition_size = TNumNodes * (TDim + 1);
        
    //Resetting the RHS
    if ( rRightHandSideVector.size() != condition_size )
        rRightHandSideVector.resize( condition_size, false );
    noalias( rRightHandSideVector ) = ZeroVector( condition_size );
    
    this->CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPCondition<2,2>::EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 2 * (2 + 1);
    unsigned int index = 0;
    
    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < 2; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPCondition<3,3>::EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 3 * (3 + 1);
    unsigned int index = 0;
    
    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < 3; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index++] = rGeom[i].GetDof(PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPCondition<3,4>::EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 4 * (3 + 1);
    unsigned int index = 0;
    
    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < 4; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index++] = rGeom[i].GetDof(PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPCondition<TDim,TNumNodes>::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{   
	KRATOS_TRY
	
	this->CalculateLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
	 
    this->CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPCondition<TDim,TNumNodes>::CalculateLHS( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo )
{    
    KRATOS_TRY
    
		const PropertiesType& Prop = this->GetProperties();
		const GeometryType& Geom = this->GetGeometry();
		const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
		const unsigned int NumGPoints = integration_points.size();
		const unsigned int LocalDim = Geom.LocalSpaceDimension();
			
		// Components of the Jacobian for computing the normal vector and shape functions container
		const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
		GeometryType::JacobiansType JContainer(NumGPoints);
		for(unsigned int i = 0; i<NumGPoints; i++)
			(JContainer[i]).resize(TDim,LocalDim,false);
		Geom.Jacobian( JContainer, mThisIntegrationMethod );
				
		//Condition variables
		ConditionVariables Variables; 
		this->InitializeConditionVariables(Variables, Geom, Prop, rCurrentProcessInfo);
		             
        for(unsigned int igauss = 0; igauss < NumGPoints; igauss++ )
        {	
			// Compute Np Vector
			noalias(Variables.Np) = row(NContainer,igauss);
			
			//Computing Nu Matrix
			PoroConditionUtilities::CalculateNuMatrix(Variables.Nu,NContainer,igauss);
								
			//Calculating weighting coefficient for integration
			this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, JContainer[igauss], integration_points[igauss].Weight() );
			
			//Calculating the unit normal vector
			this->CalculateNormalVector(Variables.NormalVector, JContainer[igauss]);
			
			//Calculating UP and PU Contributions
			this-> CalculateLHSContribution(rLeftHandSideMatrix, Variables);
			
		}

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPCondition<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{    
    KRATOS_TRY
    
		const PropertiesType& Prop = this->GetProperties();
		const GeometryType& Geom = this->GetGeometry();
		const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
		const unsigned int NumGPoints = integration_points.size();
		const unsigned int LocalDim = Geom.LocalSpaceDimension();
			
		// Components of the Jacobian for computing the normal vector and shape functions container
		const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
		GeometryType::JacobiansType JContainer(NumGPoints);
		for(unsigned int i = 0; i<NumGPoints; i++)
			(JContainer[i]).resize(TDim,LocalDim,false);
		Geom.Jacobian( JContainer, mThisIntegrationMethod );
		
		//Condition variables
		ConditionVariables Variables; 
		this->InitializeConditionVariables(Variables, Geom, Prop, rCurrentProcessInfo);
		             
        for(unsigned int igauss = 0; igauss < NumGPoints; igauss++ )
        {	
			// Compute Np Vector
			noalias(Variables.Np) = row(NContainer,igauss);
			
			//Computing Nu Matrix
			PoroConditionUtilities::CalculateNuMatrix(Variables.Nu,NContainer,igauss);
								
			//Calculating weighting coefficient for integration
			this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, JContainer[igauss], integration_points[igauss].Weight() );
			
			//Calculating the unit normal vector
			this->CalculateNormalVector(Variables.NormalVector, JContainer[igauss]);
			
			//Calculating UP and PU Contributions
			this-> CalculateRHSContribution(rRightHandSideVector, Variables);
			
		}

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPCondition<2,2>::CalculateNormalVector(VectorType& rNormalVector,  const Matrix& Jacobian )
{
	rNormalVector.resize(2,false);
	
	rNormalVector[0] = Jacobian(0,0);
    rNormalVector[1] = Jacobian(1,0);
    
    //Norm of the normal vector
    double norm = norm_2(rNormalVector);

    //Compute the unit normal vector
    if(norm>0)
    {
      rNormalVector   /= norm;
    }	
        
}

//----------------------------------------------------------------------------------------

template< >
void UPCondition<3,3>::CalculateNormalVector(VectorType& rNormalVector,  const Matrix& Jacobian )
{
	rNormalVector.resize(3,false);
	
	rNormalVector[0] = Jacobian(1,0) * Jacobian(2,1) - Jacobian(2,0) * Jacobian(1,1);
    rNormalVector[1] = Jacobian(2,0) * Jacobian(0,1) - Jacobian(0,0) * Jacobian(2,1);
    rNormalVector[2] = Jacobian(0,0) * Jacobian(1,1) - Jacobian(1,0) * Jacobian(0,1);
	
	//Norm of the normal vector
    double norm = norm_2(rNormalVector);

    //Compute the unit normal vector
    if(norm>0)
    {
      rNormalVector   /= norm;
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPCondition<3,4>::CalculateNormalVector(VectorType& rNormalVector,  const Matrix& Jacobian )
{
	rNormalVector.resize(3,false);
	
	rNormalVector[0] = Jacobian(1,0) * Jacobian(2,1) - Jacobian(2,0) * Jacobian(1,1);
    rNormalVector[1] = Jacobian(2,0) * Jacobian(0,1) - Jacobian(0,0) * Jacobian(2,1);
    rNormalVector[2] = Jacobian(0,0) * Jacobian(1,1) - Jacobian(1,0) * Jacobian(0,1);
	
	//Norm of the normal vector
    double norm = norm_2(rNormalVector);

    //Compute the unit normal vector
    if(norm>0)
    {
      rNormalVector   /= norm;
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPCondition<2,2>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& weight)
{
    double NormalVector[2];
	
	NormalVector[0] = Jacobian(0,0);
    NormalVector[1] = Jacobian(1,0);
		
	double detJ = sqrt(NormalVector[0]*NormalVector[0] + NormalVector[1]*NormalVector[1]);
	
    rIntegrationCoefficient = weight * detJ;
    
}

//----------------------------------------------------------------------------------------

template< >
void UPCondition<3,3>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& weight)
{
    double NormalVector[3];
	
	NormalVector[0] = Jacobian(1,0) * Jacobian(2,1) - Jacobian(2,0) * Jacobian(1,1);
    NormalVector[1] = Jacobian(2,0) * Jacobian(0,1) - Jacobian(0,0) * Jacobian(2,1);
    NormalVector[2] = Jacobian(0,0) * Jacobian(1,1) - Jacobian(1,0) * Jacobian(0,1);
	
	double detJ = sqrt(NormalVector[0]*NormalVector[0] + NormalVector[1]*NormalVector[1] + NormalVector[2]*NormalVector[2]);
	
    rIntegrationCoefficient = weight * detJ;
}

//----------------------------------------------------------------------------------------

template< >
void UPCondition<3,4>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& weight)
{
    double NormalVector[3];
	
	NormalVector[0] = Jacobian(1,0) * Jacobian(2,1) - Jacobian(2,0) * Jacobian(1,1);
    NormalVector[1] = Jacobian(2,0) * Jacobian(0,1) - Jacobian(0,0) * Jacobian(2,1);
    NormalVector[2] = Jacobian(0,0) * Jacobian(1,1) - Jacobian(1,0) * Jacobian(0,1);
	
	double detJ = sqrt(NormalVector[0]*NormalVector[0] + NormalVector[1]*NormalVector[1] + NormalVector[2]*NormalVector[2]);
	
    rIntegrationCoefficient = weight * detJ;
}


//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPCondition<TDim,TNumNodes>::InitializeConditionVariables(ConditionVariables& rVariables, const GeometryType& Geom, const PropertiesType& Prop, const ProcessInfo& rCurrentProcessInfo)
{   
    KRATOS_TRY
     
		//Properties variables    
		rVariables.Density = 1000.0;
	 
		//ProcessInfo variables
		// We can use the factor of the pressure since we are using the same temporal scheme for pressure and displacements
		rVariables.AcelerationCoefficient = rCurrentProcessInfo[ACCELERATION_PRESSURE_COEFFICIENT];
		
		//Variables computed at each GP
		rVariables.Np.resize(TNumNodes,false);
		rVariables.NormalVector.resize(TDim,false);
		
		//Nodal Variables
		for(unsigned int i=0; i<TNumNodes; i++)
		{
			rVariables.PressureVector[i] = Geom[i].FastGetSolutionStepValue(PRESSURE);
		}

		this->GetAccelerationVector(rVariables.AccelerationVector,0);

    
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPCondition<TDim,TNumNodes>::GetAccelerationVector( Vector& rValues, int Step )
{
    unsigned int       element_size    = TNumNodes * TDim;
    if ( rValues.size() != element_size ) rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < TNumNodes; i++ )
    {
        unsigned int index = i * TDim;
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

        if ( TDim == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
    }

}
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPCondition<TDim,TNumNodes>::CalculateLHSContribution(MatrixType& rLeftHandSideMatrix, ConditionVariables& rVariables)
{
	KRATOS_TRY
    
    noalias(rVariables.UPMatrix) = -1.0*prod(trans(rVariables.Nu),Matrix(outer_prod(rVariables.NormalVector,rVariables.Np)))*rVariables.IntegrationCoefficient;
    
    PoroConditionUtilities::AssembleUPMatrix(rLeftHandSideMatrix,rVariables.UPMatrix);
        
    noalias(rVariables.PUMatrix) = -1.0*rVariables.AcelerationCoefficient*rVariables.Density*trans(rVariables.UPMatrix);
    
    PoroConditionUtilities::AssemblePUMatrix(rLeftHandSideMatrix,rVariables.PUMatrix);
    
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPCondition<TDim,TNumNodes>::CalculateRHSContribution(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
	KRATOS_TRY
	
    noalias(rVariables.UPMatrix) = -1.0*prod(trans(rVariables.Nu),Matrix(outer_prod(rVariables.NormalVector,rVariables.Np)))*rVariables.IntegrationCoefficient;
    
    noalias(rVariables.UVector) = -1.0*prod(rVariables.UPMatrix,rVariables.PressureVector);
    
    PoroConditionUtilities::AssembleUBlockVector(rRightHandSideVector,rVariables.UVector);
    
    noalias(rVariables.PUMatrix) = -1.0*rVariables.Density*trans(rVariables.UPMatrix);
    
    noalias(rVariables.PVector) = -1.0*prod(rVariables.PUMatrix,rVariables.AccelerationVector);
        
    PoroConditionUtilities::AssemblePBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPCondition<2,2>;
template class UPCondition<3,3>;
template class UPCondition<3,4>;

} // Namespace Kratos.
