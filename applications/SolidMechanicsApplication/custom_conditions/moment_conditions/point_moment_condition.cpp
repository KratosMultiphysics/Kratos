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
    return Kratos::make_intrusive<PointMomentCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
  }


  //************************************CLONE*******************************************
  //************************************************************************************

  Condition::Pointer PointMomentCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {
    PointMomentCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewCondition.SetData(this->GetData());
    NewCondition.SetFlags(this->GetFlags());

    //-----------//
    return Kratos::make_intrusive<PointMomentCondition>(NewCondition);
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

    const SizeType number_of_nodes = GetGeometry().size();
    const unsigned int local_dimension = GetGeometry().LocalSpaceDimension();
    const SizeType& dimension       = GetGeometry().WorkingSpaceDimension();

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

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType& dimension       = GetGeometry().WorkingSpaceDimension();

    if( rVariables.ExternalVectorValue.size() != dimension )
      rVariables.ExternalVectorValue.resize(dimension,false);

    noalias(rVariables.ExternalVectorValue) = ZeroVector(dimension);

    //MOMENT CONDITION:

    if( dimension == 2 ){

      //defined on condition
      if( this->Has( PLANE_MOMENT_LOAD ) ){
	double& PointMoment = this->GetValue( PLANE_MOMENT_LOAD );
	for ( SizeType i = 0; i < number_of_nodes; i++ )
	  {
	      rVariables.ExternalScalarValue += rVariables.N[i] * PointMoment;
	  }
      }

      //defined on condition nodes
      if( this->Has( PLANE_MOMENT_LOAD_VECTOR ) ){
	Vector& PointMoments = this->GetValue( PLANE_MOMENT_LOAD_VECTOR );
	for ( SizeType i = 0; i < number_of_nodes; i++ )
	  {
	    rVariables.ExternalScalarValue += rVariables.N[i] * PointMoments[i];
	  }
      }

      //defined on condition nodes
      for (SizeType i = 0; i < number_of_nodes; i++)
	{
	  if( GetGeometry()[i].SolutionStepsDataHas( PLANE_MOMENT_LOAD ) ){
	    double& PointMoment = GetGeometry()[i].FastGetSolutionStepValue( PLANE_MOMENT_LOAD );
	    rVariables.ExternalScalarValue += rVariables.N[i] * PointMoment;
 	  }
	}

    }
    else{

      //defined on condition
      if( this->Has( MOMENT_LOAD ) ){
	array_1d<double, 3 > & PointMoment = this->GetValue( MOMENT_LOAD );
	for ( SizeType i = 0; i < number_of_nodes; i++ )
	  {
	    for( SizeType k = 0; k < dimension; k++ )
	      rVariables.ExternalVectorValue[k] += rVariables.N[i] * PointMoment[k];
	  }
      }

      //defined on condition nodes
      if( this->Has( MOMENT_LOAD_VECTOR ) ){
	Vector& PointMoments = this->GetValue( MOMENT_LOAD_VECTOR );
	unsigned int counter = 0;
	for ( SizeType i = 0; i < number_of_nodes; i++ )
	  {
	    counter = i*3;
	    for( SizeType k = 0; k < dimension; k++ )
	      {
		rVariables.ExternalVectorValue[k] += rVariables.N[i] * PointMoments[counter+k];
	      }

	  }
      }

      //defined on condition nodes
      for (SizeType i = 0; i < number_of_nodes; i++)
	{
	  if( GetGeometry()[i].SolutionStepsDataHas( MOMENT_LOAD ) ){
	    array_1d<double, 3 > & PointMoment = GetGeometry()[i].FastGetSolutionStepValue( MOMENT_LOAD );
	    for( SizeType k = 0; k < dimension; k++ )
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
