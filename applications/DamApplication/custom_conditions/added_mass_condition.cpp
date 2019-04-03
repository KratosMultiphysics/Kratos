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
#include "custom_conditions/added_mass_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer AddedMassCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new AddedMassCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void AddedMassCondition<TDim,TNumNodes>::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    const unsigned int condition_size = TDim * TNumNodes;
    unsigned int index = 0;
       
    if (rConditionDofList.size() != condition_size)
      rConditionDofList.resize( condition_size );
    
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_X);
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Y);
        if( TDim > 2)
        {
            rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Z);
        }
    }

    KRATOS_CATCH( "" )
}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void AddedMassCondition<TDim,TNumNodes>::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int condition_size = TNumNodes * TDim;
    
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
void AddedMassCondition<TDim,TNumNodes>::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    unsigned int condition_size = TNumNodes * TDim;
    
    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != condition_size )
        rLeftHandSideMatrix.resize( condition_size, condition_size, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( condition_size, condition_size );
    
    this->CalculateLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
    
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void AddedMassCondition<TDim,TNumNodes>::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    
    unsigned int condition_size = TNumNodes * TDim;
        
    //Resetting the RHS
    if ( rRightHandSideVector.size() != condition_size )
        rRightHandSideVector.resize( condition_size, false );
    noalias( rRightHandSideVector ) = ZeroVector( condition_size );
    
    this->CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void AddedMassCondition<TDim,TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = TDim * TNumNodes;
    unsigned int index = 0;
    
    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        if( TDim > 2)
        {
            rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();
        }
        
    }

    KRATOS_CATCH( "" )
}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void AddedMassCondition<TDim,TNumNodes>::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{   
	KRATOS_TRY
	
	this->CalculateLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
	 
    this->CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void AddedMassCondition<TDim,TNumNodes>::CalculateLHS( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo )
{    
    KRATOS_TRY
    
		const GeometryType& Geom = this->GetGeometry();
  		const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
		const unsigned int NumGPoints = integration_points.size();
		const unsigned int LocalDim = Geom.LocalSpaceDimension();
			
		// Components of the Jacobian for computing the normal vector and shape functions container
		BoundedMatrix<double,TDim, TNumNodes*TDim> Nu;
        //const Vector& ShapeFunctionsValues = Geom.ShapeFunctionsValues();
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
		GeometryType::JacobiansType JContainer(NumGPoints);
		for(unsigned int i = 0; i<NumGPoints; i++)
			(JContainer[i]).resize(TDim,LocalDim,false);
		Geom.Jacobian( JContainer, mThisIntegrationMethod );

        double IntegrationCoefficient;
        Vector ShapeFunctionsValues;
        ShapeFunctionsValues.resize(TNumNodes,false);
        
        for(unsigned int igauss = 0; igauss < NumGPoints; igauss++ )
        {
            double mass_contribution =0.0;

            // Distributed mass contribution according the surface element
            noalias(ShapeFunctionsValues) = row(NContainer,igauss);
            for ( unsigned int j = 0; j < TNumNodes; j++ )
            {
                 mass_contribution += ShapeFunctionsValues[j] * Geom[j].GetSolutionStepValue(ADDED_MASS);
            }

			//Computing Nu Matrix
			PoroConditionUtilities::CalculateNuMatrix(Nu,NContainer,igauss);
								
			//Calculating weighting coefficient for integration
			this->CalculateIntegrationCoefficient(IntegrationCoefficient, JContainer[igauss], integration_points[igauss].Weight() );
			
            // Mass matrix contribution
			noalias(rLeftHandSideMatrix) += rCurrentProcessInfo[ACCELERATION_PRESSURE_COEFFICIENT]*mass_contribution*prod(trans(Nu),Nu)*IntegrationCoefficient;
			
		}

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void AddedMassCondition<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{    
    KRATOS_TRY
    
     	const GeometryType& Geom = this->GetGeometry();
		const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
		const unsigned int NumGPoints = integration_points.size();
		const unsigned int LocalDim = Geom.LocalSpaceDimension();
			
		// Components of the Jacobian for computing the normal vector and shape functions container
		BoundedMatrix<double,TDim, TNumNodes*TDim> Nu;
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
		GeometryType::JacobiansType JContainer(NumGPoints);
		for(unsigned int i = 0; i<NumGPoints; i++)
			(JContainer[i]).resize(TDim,LocalDim,false);
		Geom.Jacobian( JContainer, mThisIntegrationMethod );
        
        double IntegrationCoefficient;
        BoundedMatrix<double,TNumNodes,TNumNodes> MassMatrix;
        Vector ShapeFunctionsValues;
        ShapeFunctionsValues.resize(TNumNodes,false);
        Vector AccelerationVector;
        this->GetAccelerationVector(AccelerationVector,0);

        for(unsigned int igauss = 0; igauss < NumGPoints; igauss++ )
        {	
            double mass_contribution =0.0;

            // Distributed mass contribution according the surface element
            noalias(ShapeFunctionsValues) = row(NContainer,igauss);
            for ( unsigned int j = 0; j < TNumNodes; j++ )
            {
                 mass_contribution += ShapeFunctionsValues[j] * Geom[j].GetSolutionStepValue(ADDED_MASS);
            }

			//Computing Nu Matrix
			PoroConditionUtilities::CalculateNuMatrix(Nu,NContainer,igauss);
								
			//Calculating weighting coefficient for integration
			this->CalculateIntegrationCoefficient(IntegrationCoefficient, JContainer[igauss], integration_points[igauss].Weight() );
			
    		// Mass matrix contribution
			noalias(MassMatrix) = mass_contribution*prod(trans(Nu),Nu)*IntegrationCoefficient;
            noalias(rRightHandSideVector) += -1.0*prod( MassMatrix, AccelerationVector );

		}

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void AddedMassCondition<2,2>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& weight)
{
    double NormalVector[2];
	
	NormalVector[0] = Jacobian(0,0);
    NormalVector[1] = Jacobian(1,0);
		
	double detJ = sqrt(NormalVector[0]*NormalVector[0] + NormalVector[1]*NormalVector[1]);
	
    rIntegrationCoefficient = weight * detJ;
    
}

//----------------------------------------------------------------------------------------

template< >
void AddedMassCondition<3,3>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& weight)
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
void AddedMassCondition<3,4>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& weight)
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
void AddedMassCondition<TDim,TNumNodes>::GetAccelerationVector( Vector& rValues, int Step )
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

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class AddedMassCondition<2,2>;
template class AddedMassCondition<3,3>;
template class AddedMassCondition<3,4>;

} // Namespace Kratos.