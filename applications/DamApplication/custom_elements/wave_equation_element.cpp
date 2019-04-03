//
//   Project Name:        KratosDamApplication  $
//   Last modified by:    $Author: lgracia		$
//   Date:                $Date: December 2016  $
//   Revision:            $Revision: 1.0        $
//
//
#include "custom_elements/wave_equation_element.hpp"

#include "utilities/math_utils.h"

namespace Kratos
{
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
	Element::Pointer WaveEquationElement<TDim,TNumNodes>::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
	{    
		return Element::Pointer( new WaveEquationElement( NewId, this->GetGeometry().Create( ThisNodes ), pProperties ) );
	}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
	Element::Pointer WaveEquationElement<TDim,TNumNodes>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
	{    
		return Element::Pointer( new WaveEquationElement( NewId, pGeom, pProperties ) );
	}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
int WaveEquationElement<TDim,TNumNodes>::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

	const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();

    // verify nodal variables and dofs
	if ( PRESSURE.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "PRESSURE has Key zero at element", this->Id() )

	for ( unsigned int i = 0; i < TNumNodes; i++ )
    {
		if ( Geom[i].SolutionStepsDataHas( PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable PRESSURE on node ", Geom[i].Id() )
       	if ( Geom[i].SolutionStepsDataHas( Dt_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable Dt_PRESSURE on node ", Geom[i].Id() )  
       	if ( Geom[i].SolutionStepsDataHas( Dt2_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable Dt2_PRESSURE on node ", Geom[i].Id() )
            
        if ( Geom[i].HasDofFor( PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing the dof for the variable PRESSURE on node ", Geom[i].Id() )  
	}
	
	    // Verify ProcessInfo variables
    if ( VELOCITY_PRESSURE_COEFFICIENT.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"VELOCITY_PRESSURE_COEFFICIENT has Key zero at element", this->Id() )
    if ( ACCELERATION_PRESSURE_COEFFICIENT.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"ACCELERATION_PRESSURE_COEFFICIENT has Key zero at element", this->Id() )
        
       // Verify properties
    if ( BULK_MODULUS_FLUID.Key() == 0 || Prop.Has( BULK_MODULUS_FLUID ) == false || Prop[BULK_MODULUS_FLUID] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"BULK_MODULUS_FLUID has Key zero, is not defined or has an invalid value at element", this->Id() )
    if ( DENSITY_WATER.Key() == 0 || Prop.Has( DENSITY_WATER ) == false || Prop[DENSITY_WATER] < 0.0 )
        KRATOS_THROW_ERROR( std::invalid_argument,"DENSITY_WATER has Key zero, is not defined or has an invalid value at element", this->Id() )     
   
	
	return 0;

    KRATOS_CATCH( "" );
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
	void WaveEquationElement<TDim,TNumNodes>::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
	{
		KRATOS_TRY
		
		GeometryType& rGeom = this->GetGeometry();
		const unsigned int element_size = TNumNodes;
		unsigned int index = 0;
		
		if (rElementalDofList.size() != element_size)
		  rElementalDofList.resize( element_size );
		
		for (unsigned int i = 0; i < TNumNodes; i++)
		{
			rElementalDofList[index++] = rGeom[i].pGetDof(PRESSURE);
		}

		KRATOS_CATCH( "" )
	}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
	void WaveEquationElement<TDim,TNumNodes>::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
	{
		KRATOS_TRY

		GeometryType& rGeom = this->GetGeometry();
		const unsigned int element_size = TNumNodes;
		unsigned int index = 0;
		
		if (rResult.size() != element_size)
		  rResult.resize( element_size, false );

		for (unsigned int i = 0; i < TNumNodes; i++)
		{
			rResult[index++] = rGeom[i].GetDof(PRESSURE).EquationId();
		}

		KRATOS_CATCH( "" )
	}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
	void WaveEquationElement<TDim,TNumNodes>::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
	{
		KRATOS_TRY

		const unsigned int element_size = TNumNodes;
		
		//Resetting the LHS
		if ( rLeftHandSideMatrix.size1() != element_size )
			rLeftHandSideMatrix.resize( element_size, element_size, false );
		noalias( rLeftHandSideMatrix ) = ZeroMatrix( element_size, element_size );
		
				//Resetting the RHS
		if ( rRightHandSideVector.size() != element_size )
			rRightHandSideVector.resize( element_size, false );
		noalias( rRightHandSideVector ) = ZeroVector( element_size );
		
		this->CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
					
		KRATOS_CATCH( "" )
	}


//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
	void WaveEquationElement<TDim,TNumNodes>::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
	{
    KRATOS_TRY;
    
		const unsigned int element_size = TNumNodes;
		
		//Resetting the LHS
		if ( rLeftHandSideMatrix.size1() != element_size )
			rLeftHandSideMatrix.resize( element_size, element_size, false );
		noalias( rLeftHandSideMatrix ) = ZeroMatrix( element_size, element_size );
		
		this->CalculateLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
    
    KRATOS_CATCH("");
	}
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
	void WaveEquationElement<TDim,TNumNodes>::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
	{
		KRATOS_TRY
		
		const unsigned int element_size = TNumNodes;
				
		 //Resetting the RHS
		if ( rRightHandSideVector.size() != element_size )
			rRightHandSideVector.resize( element_size, false );
		noalias( rRightHandSideVector ) = ZeroVector( element_size );
		
		this->CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);
		
		KRATOS_CATCH( "" )
         
	}
	
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void WaveEquationElement<TDim,TNumNodes>::GetValuesVector( Vector& rValues, int Step )
{
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int element_size = TNumNodes;
    unsigned int index = 0;

    if ( rValues.size() != element_size )
        rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < TNumNodes; i++ )
    {
        rValues[index++] = Geom[i].FastGetSolutionStepValue( PRESSURE, Step );
    }
    
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void WaveEquationElement<TDim,TNumNodes>::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int element_size = TNumNodes;
    unsigned int index = 0;

    if ( rValues.size() != element_size )
        rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < TNumNodes; i++ )
    {
        rValues[index++] = Geom[i].FastGetSolutionStepValue( Dt_PRESSURE, Step );
    }
    
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void WaveEquationElement<TDim,TNumNodes>::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    const GeometryType& Geom = this->GetGeometry();
    const unsigned int element_size = TNumNodes;
    unsigned int index = 0;

    if ( rValues.size() != element_size )
        rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < TNumNodes; i++ )
    {
        rValues[index++] = Geom[i].FastGetSolutionStepValue( Dt2_PRESSURE, Step );
    }
    
}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void WaveEquationElement<TDim,TNumNodes>::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{    
	//Contributions to the left hand side
    this->CalculateLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
    
    //Contributions to the right hand side
    this->CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);
	
}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void WaveEquationElement<TDim,TNumNodes>::CalculateLHS( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
		KRATOS_TRY

		const PropertiesType& Prop = this->GetProperties();
		const GeometryType& Geom = this->GetGeometry();
		const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
		const unsigned int NumGPoints = integration_points.size();
		
		 //Defining the shape functions, the jacobian and the shape functions local gradients Containers
		array_1d<double,TNumNodes> Np;
		BoundedMatrix<double,TNumNodes,TDim> GradNpT; 
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
		Vector detJContainer(NumGPoints);
		Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,detJContainer,mThisIntegrationMethod);
		double IntegrationCoefficient;

        // Definition of the speed in the fluid
        const double BulkModulus = Prop[BULK_MODULUS_FLUID];
        const double Water_density = Prop[DENSITY_WATER];
        const double inv_c_speed = 1.0 /sqrt(BulkModulus/Water_density);
                
        for ( unsigned int igauss = 0; igauss < NumGPoints; igauss++ )
        {	
			
			noalias(Np) = row(NContainer,igauss);
			
			noalias(GradNpT) = DN_DXContainer[igauss];
            
			//calculating weighting coefficient for integration
			this->CalculateIntegrationCoefficient( IntegrationCoefficient, detJContainer[igauss], integration_points[igauss].Weight() );

			// Mass matrix contribution
			noalias(rLeftHandSideMatrix) += rCurrentProcessInfo[ACCELERATION_PRESSURE_COEFFICIENT]*(inv_c_speed*inv_c_speed)*outer_prod(Np,Np)*IntegrationCoefficient;
            
			// Stiffness matrix contribution
			noalias(rLeftHandSideMatrix) += prod(GradNpT,trans(GradNpT))*IntegrationCoefficient;
            
		}
				
		KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void WaveEquationElement<TDim,TNumNodes>::CalculateRHS(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
		KRATOS_TRY
	
		const PropertiesType& Prop = this->GetProperties();
		const GeometryType& Geom = this->GetGeometry();
		const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
		const unsigned int NumGPoints = integration_points.size();
		
		BoundedMatrix<double,TNumNodes,TNumNodes> MassMatrix;
		BoundedMatrix<double,TNumNodes,TNumNodes> SiffnessMatrix;

		//Defining the shape functions, the jacobian and the shape functions local gradients Containers
		array_1d<double,TNumNodes> Np;
		BoundedMatrix<double,TNumNodes,TDim> GradNpT; 
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer(NumGPoints);
		Vector detJContainer(NumGPoints);
		Geom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer,detJContainer,mThisIntegrationMethod);
		double IntegrationCoefficient;

        // Definition of the speed in the fluid
        const double BulkModulus = Prop[BULK_MODULUS_FLUID];
        const double Water_density = Prop[DENSITY_WATER];
        const double inv_c_speed = 1.0 /sqrt(BulkModulus/Water_density);
        
         //Nodal Variables
        Vector PressureVector;
        Vector Dt2PressureVector;
        this->GetValuesVector(PressureVector, 0);
        this->GetSecondDerivativesVector(Dt2PressureVector, 0);
         
        for ( unsigned int igauss = 0; igauss < NumGPoints; igauss++ )
        {	
			
			noalias(Np) = row(NContainer,igauss);
			
			noalias(GradNpT) = DN_DXContainer[igauss];
						
			//calculating weighting coefficient for integration
			this->CalculateIntegrationCoefficient( IntegrationCoefficient, detJContainer[igauss], integration_points[igauss].Weight() );

			// Mass matrix contribution
			noalias(MassMatrix) = (inv_c_speed*inv_c_speed)*outer_prod(Np,Np)*IntegrationCoefficient;
			noalias(rRightHandSideVector) += -1.0*prod(MassMatrix,Dt2PressureVector); 
			
			// Stiffness matrix contribution
			noalias(SiffnessMatrix) = prod(GradNpT,trans(GradNpT))*IntegrationCoefficient;
			noalias(rRightHandSideVector) += -1.0*prod(SiffnessMatrix,PressureVector);
            
		}
				
		KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void WaveEquationElement<TDim,TNumNodes>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const double& detJ, const double& weight)
{
    rIntegrationCoefficient = weight * detJ;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class WaveEquationElement<2,3>;
template class WaveEquationElement<2,4>;
template class WaveEquationElement<3,4>;
template class WaveEquationElement<3,8>;

} // Namespace Kratos
