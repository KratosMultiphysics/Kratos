//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_conditions/elastic_conditions/point_elastic_condition.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{


  //***********************************************************************************
  //***********************************************************************************
  PointElasticCondition::PointElasticCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : ElasticCondition(NewId, pGeometry)
  {

  }

  //***********************************************************************************
  //***********************************************************************************
  PointElasticCondition::PointElasticCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : ElasticCondition(NewId, pGeometry, pProperties)
  {
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
  }

  //************************************************************************************
  //************************************************************************************
  PointElasticCondition::PointElasticCondition( PointElasticCondition const& rOther )
    : ElasticCondition(rOther)     
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  Condition::Pointer PointElasticCondition::Create(
						  IndexType NewId,
						  NodesArrayType const& ThisNodes,
						  PropertiesType::Pointer pProperties) const
  {
    return Condition::Pointer(new PointElasticCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
  }


  //************************************CLONE*******************************************
  //************************************************************************************

  Condition::Pointer PointElasticCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {
    PointElasticCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewCondition.SetData(this->GetData());
    NewCondition.SetFlags(this->GetFlags());

    //-----------//      
    return Condition::Pointer( new PointElasticCondition(NewCondition) );
  }


  //***********************************************************************************
  //***********************************************************************************
  PointElasticCondition::~PointElasticCondition()
  {
  }


  //************************************************************************************
  //************************************************************************************

  void PointElasticCondition::InitializeConditionVariables(ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY
      
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int local_dimension = GetGeometry().LocalSpaceDimension();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rVariables.Initialize(dimension, local_dimension, number_of_nodes);
   
    //Only one node:
    rVariables.N[0] = 1.0;

    KRATOS_CATCH( "" )

  }


  //***********************************************************************************
  //***********************************************************************************

  void PointElasticCondition::CalculateExternalStiffness(ConditionVariables& rVariables)
  {
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if( rVariables.ExternalVectorValue.size() != dimension )
      rVariables.ExternalVectorValue.resize(dimension,false);

    noalias(rVariables.ExternalVectorValue) = ZeroVector(dimension);
   
    //STIFFNESS CONDITION:
    
    //defined on condition
    if( this->Has( POINT_STIFFNESS ) ){
      array_1d<double, 3 > & PointStiffness = this->GetValue( POINT_STIFFNESS );
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
	{
	  for( unsigned int k = 0; k < dimension; k++ )
	    rVariables.ExternalVectorValue[k] += rVariables.N[i] * fabs(PointStiffness[k]);
	}
    }
    
    //defined on condition nodes
    if( this->Has( POINT_STIFFNESS_VECTOR ) ){
      Vector& PointStiffnesss = this->GetValue( POINT_STIFFNESS_VECTOR );
      unsigned int counter = 0;
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
	{
	  counter = i*3;
	  for( unsigned int k = 0; k < dimension; k++ )
	    {
	      rVariables.ExternalVectorValue[k] += rVariables.N[i] * fabs(PointStiffnesss[counter+k]);
	    }
	  
	}
    }
    
    //defined on condition nodes      
    for (unsigned int i = 0; i < number_of_nodes; i++)
      {
	if( GetGeometry()[i].SolutionStepsDataHas( POINT_STIFFNESS ) ){
	  array_1d<double, 3 > & PointStiffness = GetGeometry()[i].FastGetSolutionStepValue( POINT_STIFFNESS );
	  for( unsigned int k = 0; k < dimension; k++ )
	    rVariables.ExternalVectorValue[k] += rVariables.N[i] * fabs(PointStiffness[k]);
 
	}
      }

    KRATOS_CATCH( "" )
  }


  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************

  void PointElasticCondition::CalculateKinematics(ConditionVariables& rVariables,
						 const double& rPointNumber)
  {
    KRATOS_TRY

    rVariables.Jacobian = 1.0;

    this->CalculateExternalStiffness(rVariables);

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void PointElasticCondition::CalculateConditionSystem(LocalSystemComponents& rLocalSystem,
						      const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    //create and initialize condition variables:
    ConditionVariables Variables;
    this->InitializeConditionVariables(Variables,rCurrentProcessInfo);

    //reading integration points
    unsigned int PointNumber = 0;

    //compute element kinematics B, F, DN_DX ...
    this->CalculateKinematics(Variables,PointNumber);
    
    //calculating weights for integration on the "reference configuration"
    double IntegrationWeight = 1;
    
    if ( rLocalSystem.CalculationFlags.Is(BoundaryCondition::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
      {
	//contributions to stiffness matrix calculated on the reference config
	this->CalculateAndAddLHS ( rLocalSystem, Variables, IntegrationWeight );
      }
    
    if ( rLocalSystem.CalculationFlags.Is(BoundaryCondition::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
      {
	
	this->CalculateAndAddRHS ( rLocalSystem, Variables, IntegrationWeight );
      }
    
  
    KRATOS_CATCH( "" )
  }

 
  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************

  int PointElasticCondition::Check( const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    // Perform base condition checks
    int ErrorCode = 0;
    ErrorCode = ElasticCondition::Check(rCurrentProcessInfo);
    
    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(POINT_STIFFNESS);
    KRATOS_CHECK_VARIABLE_KEY(POINT_STIFFNESS_VECTOR);

    return ErrorCode;
    
    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void PointElasticCondition::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ElasticCondition )
  }

  void PointElasticCondition::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ElasticCondition )
  }



} // Namespace Kratos.
