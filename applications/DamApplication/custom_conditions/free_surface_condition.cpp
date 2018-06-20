//   
//   Project Name:        			KratosDamApplication $
//   Last Modified by:    $Author:    	  Lorenzo Gracia $
//   Date:                $Date:           	January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_conditions/free_surface_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer FreeSurfaceCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new FreeSurfaceCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void FreeSurfaceCondition<TDim,TNumNodes>::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
       
    GeometryType& rGeom = GetGeometry();
    const unsigned int condition_size = TNumNodes;
    unsigned int index = 0;
    
    if (rConditionDofList.size() != condition_size)
      rConditionDofList.resize( condition_size );
    
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        rConditionDofList[index++] = rGeom[i].pGetDof(PRESSURE);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void FreeSurfaceCondition<TDim,TNumNodes>::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int condition_size = TNumNodes;
    
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
void FreeSurfaceCondition<TDim,TNumNodes>::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    
    unsigned int condition_size = TNumNodes;
        
    //Resetting the RHS
    if ( rLeftHandSideMatrix.size1() != condition_size )
        rLeftHandSideMatrix.resize( condition_size, condition_size, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( condition_size, condition_size );
    
    this->CalculateLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void FreeSurfaceCondition<TDim,TNumNodes>::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    
    unsigned int condition_size = TNumNodes;
        
    //Resetting the RHS
    if ( rRightHandSideVector.size() != condition_size )
        rRightHandSideVector.resize( condition_size, false );
    noalias( rRightHandSideVector ) = ZeroVector( condition_size );
    
    this->CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void FreeSurfaceCondition<TDim,TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = TNumNodes;
    unsigned int index = 0;
    
    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        rResult[index++] = rGeom[i].GetDof(PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void FreeSurfaceCondition<TDim,TNumNodes>::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
{   
	this->CalculateLHS(rLeftHandSideMatrix,CurrentProcessInfo);
	 
    this->CalculateRHS(rRightHandSideVector,CurrentProcessInfo);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void FreeSurfaceCondition<TDim,TNumNodes>::CalculateLHS( MatrixType& rLeftHandSideMatrix, const ProcessInfo& CurrentProcessInfo )
{    
		KRATOS_TRY

		const GeometryType& Geom = this->GetGeometry();
		const unsigned int element_size = TNumNodes;
		const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
		const unsigned int NumGPoints = integration_points.size();
        const unsigned int LocalDim = Geom.LocalSpaceDimension();

		//Resetting the LHS
		if ( rLeftHandSideMatrix.size1() != element_size )
			rLeftHandSideMatrix.resize( element_size, element_size, false );
		noalias( rLeftHandSideMatrix ) = ZeroMatrix( element_size, element_size );
		
		 //Defining the shape functions, the jacobian and the shape functions local gradients Containers
		array_1d<double,TNumNodes> Np;
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        GeometryType::JacobiansType JContainer(NumGPoints);
		for(unsigned int i = 0; i<NumGPoints; i++)
			(JContainer[i]).resize(TDim,LocalDim,false);
		Geom.Jacobian( JContainer, mThisIntegrationMethod );
        
        
		double IntegrationCoefficient;
		double inv_gravity = 1.0/9.81;
               
        for ( unsigned int igauss = 0; igauss < NumGPoints; igauss++ )
        {	
			noalias(Np) = row(NContainer,igauss);
								
			//calculating weighting coefficient for integration
            this->CalculateIntegrationCoefficient(IntegrationCoefficient, JContainer[igauss], integration_points[igauss].Weight() );

			// Mass matrix contribution
			noalias(rLeftHandSideMatrix) += CurrentProcessInfo[ACCELERATION_PRESSURE_COEFFICIENT]*(inv_gravity)*outer_prod(Np,Np)*IntegrationCoefficient;
            
		}
				
		KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void FreeSurfaceCondition<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
{    
		KRATOS_TRY

		const GeometryType& Geom = this->GetGeometry();
		const unsigned int element_size = TNumNodes;
		const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
		const unsigned int NumGPoints = integration_points.size();
        const unsigned int LocalDim = Geom.LocalSpaceDimension();

        
		//Resetting the RHS
		if ( rRightHandSideVector.size() != element_size )
			rRightHandSideVector.resize( element_size, false );
		noalias( rRightHandSideVector ) = ZeroVector( element_size );
		

		BoundedMatrix<double,TNumNodes,TNumNodes> MassMatrix;
		
		 //Defining the shape functions, the jacobian and the shape functions local gradients Containers
		array_1d<double,TNumNodes> Np;
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
		GeometryType::JacobiansType JContainer(NumGPoints);
		for(unsigned int i = 0; i<NumGPoints; i++)
			(JContainer[i]).resize(TDim,LocalDim,false);
		Geom.Jacobian( JContainer, mThisIntegrationMethod );
		double IntegrationCoefficient;
		double inv_gravity = 1.0/9.81;
		
		//Nodal Variables
        array_1d<double,TNumNodes> Dt2PressureVector;
		        
		for(unsigned int i=0; i<TNumNodes; i++)
		{
			Dt2PressureVector[i] = Geom[i].FastGetSolutionStepValue(Dt2_PRESSURE);
		}
               
        for ( unsigned int igauss = 0; igauss < NumGPoints; igauss++ )
        {	
			noalias(Np) = row(NContainer,igauss);
								
			//calculating weighting coefficient for integration
            this->CalculateIntegrationCoefficient(IntegrationCoefficient, JContainer[igauss], integration_points[igauss].Weight() );

			// Mass matrix contribution
			noalias(MassMatrix) = (inv_gravity)*outer_prod(Np,Np)*IntegrationCoefficient;
			noalias(rRightHandSideVector) += -1.0*prod(MassMatrix,Dt2PressureVector); 
            
			
		}
				
		KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


template< >
void FreeSurfaceCondition<2,2>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& weight)
{
    double NormalVector[2];
	
	NormalVector[0] = Jacobian(0,0);
    NormalVector[1] = Jacobian(1,0);
		
	double detJ = sqrt(NormalVector[0]*NormalVector[0] + NormalVector[1]*NormalVector[1]);
	
    rIntegrationCoefficient = weight * detJ;
    
}

//----------------------------------------------------------------------------------------

template< >
void FreeSurfaceCondition<3,3>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& weight)
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
void FreeSurfaceCondition<3,4>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& weight)
{
    double NormalVector[3];
	
	NormalVector[0] = Jacobian(1,0) * Jacobian(2,1) - Jacobian(2,0) * Jacobian(1,1);
    NormalVector[1] = Jacobian(2,0) * Jacobian(0,1) - Jacobian(0,0) * Jacobian(2,1);
    NormalVector[2] = Jacobian(0,0) * Jacobian(1,1) - Jacobian(1,0) * Jacobian(0,1);
	
	double detJ = sqrt(NormalVector[0]*NormalVector[0] + NormalVector[1]*NormalVector[1] + NormalVector[2]*NormalVector[2]);
	
    rIntegrationCoefficient = weight * detJ;
}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class FreeSurfaceCondition<2,2>;
template class FreeSurfaceCondition<3,3>;
template class FreeSurfaceCondition<3,4>;

} // Namespace Kratos.
