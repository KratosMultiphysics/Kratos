//   
//   Project Name:        			KratosDamApplication $
//   Last Modified by:    $Author:    	  Lorenzo Gracia $
//   Date:                $Date:           	January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_conditions/infinite_domain_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer InfiniteDomainCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new InfiniteDomainCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void InfiniteDomainCondition<TDim,TNumNodes>::CalculateLHS( MatrixType& rLeftHandSideMatrix, const ProcessInfo& CurrentProcessInfo )
{    
        
        KRATOS_TRY
		
		//const PropertiesType& Prop = this->GetProperties();
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
		
		// Definition of the speed in the fluid
        //~ const double BulkModulus = Prop[BULK_MODULUS_FLUID];
        //~ const double Water_density = Prop[DENSITY_WATER];
        const double BulkModulus = 2.21e9;
        const double Water_density = 1000.0;
        const double inv_c_speed = 1.0 /sqrt(BulkModulus/Water_density);
               
        for ( unsigned int igauss = 0; igauss < NumGPoints; igauss++ )
        {	
			noalias(Np) = row(NContainer,igauss);
								
			//calculating weighting coefficient for integration
			this->CalculateIntegrationCoefficient( IntegrationCoefficient, JContainer[igauss], integration_points[igauss].Weight() );

			// Mass matrix contribution
			noalias(rLeftHandSideMatrix) += CurrentProcessInfo[VELOCITY_PRESSURE_COEFFICIENT]*(inv_c_speed)*outer_prod(Np,Np)*IntegrationCoefficient;
		}
        
				
		KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void InfiniteDomainCondition<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
{    
		KRATOS_TRY

		//const PropertiesType& Prop = this->GetProperties();
		const GeometryType& Geom = this->GetGeometry();
		const unsigned int element_size = TNumNodes;
		const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
		const unsigned int NumGPoints = integration_points.size();
        const unsigned int LocalDim = Geom.LocalSpaceDimension();

        
		//Resetting the RHS
		if ( rRightHandSideVector.size() != element_size )
			rRightHandSideVector.resize( element_size, false );
		noalias( rRightHandSideVector ) = ZeroVector( element_size );
		

		BoundedMatrix<double,TNumNodes,TNumNodes> DampingMatrix;
		
		 //Defining the shape functions, the jacobian and the shape functions local gradients Containers
		array_1d<double,TNumNodes> Np;
        const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
        GeometryType::JacobiansType JContainer(NumGPoints);
		for(unsigned int i = 0; i<NumGPoints; i++)
			(JContainer[i]).resize(TDim,LocalDim,false);
		Geom.Jacobian( JContainer, mThisIntegrationMethod );
		double IntegrationCoefficient;
		
		// Definition of the speed in the fluid
        //~ const double BulkModulus = Prop[BULK_MODULUS_FLUID];
        //~ const double Water_density = Prop[DENSITY_WATER];
        const double BulkModulus = 2.21e9;
        const double Water_density = 1000.0;
        const double inv_c_speed = 1.0 /sqrt(BulkModulus/Water_density);
        	
		//Nodal Variables
        array_1d<double,TNumNodes> DtPressureVector;
		        
		for(unsigned int i=0; i<TNumNodes; i++)
		{
			DtPressureVector[i] = Geom[i].FastGetSolutionStepValue(Dt_PRESSURE);
		}
               
        for ( unsigned int igauss = 0; igauss < NumGPoints; igauss++ )
        {	
			noalias(Np) = row(NContainer,igauss);
								
			//calculating weighting coefficient for integration
			this->CalculateIntegrationCoefficient( IntegrationCoefficient, JContainer[igauss], integration_points[igauss].Weight() );

			// Mass matrix contribution
			noalias(DampingMatrix) = (inv_c_speed)*outer_prod(Np,Np)*IntegrationCoefficient;
			noalias(rRightHandSideVector) += -1.0*prod(DampingMatrix,DtPressureVector); 
			
		}
        
		KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class InfiniteDomainCondition<2,2>;
template class InfiniteDomainCondition<3,3>;
template class InfiniteDomainCondition<3,4>;

} // Namespace Kratos.
