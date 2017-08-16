//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_conditions/elastic_condition.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//***********************************************************************************
//***********************************************************************************
ElasticCondition::ElasticCondition()
    : BoundaryCondition()
{
  //DO NOT CALL IT: only needed for Register and Serialization!!!
}


//***********************************************************************************
//***********************************************************************************
ElasticCondition::ElasticCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : BoundaryCondition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//***********************************************************************************
//***********************************************************************************
ElasticCondition::ElasticCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : BoundaryCondition(NewId, pGeometry, pProperties)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
ElasticCondition::ElasticCondition( ElasticCondition const& rOther )
    : BoundaryCondition(rOther)

{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer ElasticCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new ElasticCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}


//************************************CLONE*******************************************
//************************************************************************************

Condition::Pointer ElasticCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
  
  ElasticCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

  NewCondition.SetData(this->GetData());
  NewCondition.SetFlags(this->GetFlags());

  
  return Condition::Pointer( new ElasticCondition(NewCondition) );
}


//***********************************************************************************
//***********************************************************************************
ElasticCondition::~ElasticCondition()
{
}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//***********************************************************************************
//***********************************************************************************

void ElasticCondition::CalculateExternalStiffness(ConditionVariables& rVariables)
{
    KRATOS_TRY

    //const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    
    if( rVariables.ExternalVectorValue.size() != dimension )
      rVariables.ExternalVectorValue.resize(dimension,false);

    noalias(rVariables.ExternalVectorValue) = ZeroVector(dimension);

    rVariables.ExternalScalarValue = 0.0;

    // the elasticity modulus from the elastic constraint must be read.
    
    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void ElasticCondition::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
						     ConditionVariables& rVariables,
						     double& rIntegrationWeight)

{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    Vector CurrentValueVector(dimension);
    noalias(CurrentValueVector) = ZeroVector(dimension); 

    int index = 0;
    
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        index = dimension * i;

	//current displacements
	CurrentValueVector = GetNodalCurrentValue( DISPLACEMENT, CurrentValueVector, i );

        for ( unsigned int j = 0; j < dimension; j++ )
        {
	  rRightHandSideVector[index + j] += CurrentValueVector[j] * rVariables.ExternalVectorValue[j] * rIntegrationWeight;
        }

    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

double& ElasticCondition::CalculateAndAddExternalEnergy(double& rEnergy,
							ConditionVariables& rVariables,
							double& rIntegrationWeight)

{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Energy Calculation:
    Vector CurrentValueVector(dimension);
    noalias(CurrentValueVector) = ZeroVector(dimension); 
    Vector Displacements(dimension);
    noalias(Displacements) = ZeroVector(dimension);
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	//current displacements to compute energy
	CurrentValueVector = GetNodalCurrentValue( DISPLACEMENT, CurrentValueVector, i );

	Displacements += rVariables.N[i] * CurrentValueVector;
      }
    //------
    
    Vector ForceVector(dimension);
    noalias(ForceVector) = ZeroVector(dimension);
    
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
 	//current displacements
	CurrentValueVector = GetNodalCurrentValue( DISPLACEMENT, CurrentValueVector, i );
	
        for ( unsigned int j = 0; j < dimension; j++ )
        {
	  ForceVector[j] += CurrentValueVector[j] * rVariables.ExternalVectorValue[j] * rIntegrationWeight;
        }
    }
    
    rEnergy += inner_prod( ForceVector, Displacements );
    
    return rEnergy;

    KRATOS_CATCH( "" )
}

  
//***********************************************************************************
//***********************************************************************************

void ElasticCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BoundaryCondition )
}

void ElasticCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BoundaryCondition )
}


} // Namespace Kratos.
