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
#include "custom_conditions/moment_conditions/point_moment_condition.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{


  //***********************************************************************************
  //***********************************************************************************
  PointMomentCondition::PointMomentCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : MomentCondition(NewId, pGeometry)
  {

  }

  //***********************************************************************************
  //***********************************************************************************
  PointMomentCondition::PointMomentCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : MomentCondition(NewId, pGeometry, pProperties)
  {
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
  }

  //************************************************************************************
  //************************************************************************************
  PointMomentCondition::PointMomentCondition( PointMomentCondition const& rOther )
    : MomentCondition(rOther)     
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  Condition::Pointer PointMomentCondition::Create(
						  IndexType NewId,
						  NodesArrayType const& ThisNodes,
						  PropertiesType::Pointer pProperties) const
  {
    return Condition::Pointer(new PointMomentCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
  }


  //************************************CLONE*******************************************
  //************************************************************************************

  Condition::Pointer PointMomentCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {
    PointMomentCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewCondition.SetData(this->GetData());
    NewCondition.SetFlags(this->GetFlags());

    //-----------//      
    return Condition::Pointer( new PointMomentCondition(NewCondition) );
  }


  //***********************************************************************************
  //***********************************************************************************
  PointMomentCondition::~PointMomentCondition()
  {
  }


  //************************************************************************************
  //************************************************************************************

  void PointMomentCondition::InitializeConditionVariables(ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
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

  void PointMomentCondition::CalculateExternalMoment(ConditionVariables& rVariables)
  {
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if( rVariables.ExternalVectorValue.size() != dimension )
      rVariables.ExternalVectorValue.resize(dimension,false);

    noalias(rVariables.ExternalVectorValue) = ZeroVector(dimension);
   
    //MOMENT CONDITION:

    if( dimension == 2 ){
    
      //defined on condition
      if( this->Has( PLANE_POINT_MOMENT ) ){
	double& PointMoment = this->GetValue( PLANE_POINT_MOMENT );
	for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
	      rVariables.ExternalScalarValue += rVariables.N[i] * PointMoment;
	  }
      }
    
      //defined on condition nodes
      if( this->Has( PLANE_POINT_MOMENT_VECTOR ) ){
	Vector& PointMoments = this->GetValue( PLANE_POINT_MOMENT_VECTOR );
	for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
	    rVariables.ExternalScalarValue += rVariables.N[i] * PointMoments[i];	  
	  }
      }
    
      //defined on condition nodes      
      for (unsigned int i = 0; i < number_of_nodes; i++)
	{
	  if( GetGeometry()[i].SolutionStepsDataHas( PLANE_POINT_MOMENT ) ){
	    double& PointMoment = GetGeometry()[i].FastGetSolutionStepValue( PLANE_POINT_MOMENT );
	    rVariables.ExternalScalarValue += rVariables.N[i] * PointMoment;
 	  }
	}

    }
    else{

      //defined on condition
      if( this->Has( POINT_MOMENT ) ){
	array_1d<double, 3 > & PointMoment = this->GetValue( POINT_MOMENT );
	for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
	    for( unsigned int k = 0; k < dimension; k++ )
	      rVariables.ExternalVectorValue[k] += rVariables.N[i] * PointMoment[k];
	  }
      }
    
      //defined on condition nodes
      if( this->Has( POINT_MOMENT_VECTOR ) ){
	Vector& PointMoments = this->GetValue( POINT_MOMENT_VECTOR );
	unsigned int counter = 0;
	for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
	    counter = i*3;
	    for( unsigned int k = 0; k < dimension; k++ )
	      {
		rVariables.ExternalVectorValue[k] += rVariables.N[i] * PointMoments[counter+k];
	      }
	  
	  }
      }
    
      //defined on condition nodes      
      for (unsigned int i = 0; i < number_of_nodes; i++)
	{
	  if( GetGeometry()[i].SolutionStepsDataHas( POINT_MOMENT ) ){
	    array_1d<double, 3 > & PointMoment = GetGeometry()[i].FastGetSolutionStepValue( POINT_MOMENT );
	    for( unsigned int k = 0; k < dimension; k++ )
	      rVariables.ExternalVectorValue[k] += rVariables.N[i] * PointMoment[k];
 
	  }
	}

    }
    
    KRATOS_CATCH( "" )
  }


  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************

  void PointMomentCondition::CalculateKinematics(ConditionVariables& rVariables,
						 const double& rPointNumber)
  {
    KRATOS_TRY

    rVariables.Jacobian = 1.0;

    this->CalculateExternalMoment(rVariables);

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void PointMomentCondition::CalculateConditionSystem(LocalSystemComponents& rLocalSystem,
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

  int PointMomentCondition::Check( const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    // Perform base condition checks
    int ErrorCode = 0;
    ErrorCode = MomentCondition::Check(rCurrentProcessInfo);
    
    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(POINT_MOMENT);
    KRATOS_CHECK_VARIABLE_KEY(POINT_MOMENT_VECTOR);
    KRATOS_CHECK_VARIABLE_KEY(PLANE_POINT_MOMENT);
    KRATOS_CHECK_VARIABLE_KEY(PLANE_POINT_MOMENT_VECTOR);

    return ErrorCode;
    
    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void PointMomentCondition::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MomentCondition )
  }

  void PointMomentCondition::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MomentCondition )
  }



} // Namespace Kratos.
