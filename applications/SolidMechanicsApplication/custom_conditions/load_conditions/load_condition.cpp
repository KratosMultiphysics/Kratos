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
#include "custom_conditions/load_conditions/load_condition.hpp"
#include "utilities/beam_math_utilities.hpp"

namespace Kratos
{

  //***********************************************************************************
  //***********************************************************************************
  LoadCondition::LoadCondition()
    : BoundaryCondition()
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  LoadCondition::LoadCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : BoundaryCondition(NewId, pGeometry)
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  LoadCondition::LoadCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : BoundaryCondition(NewId, pGeometry, pProperties)
  {
    this->Set(SOLID);
    this->Set(STRUCTURE);
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
  }

  //************************************************************************************
  //************************************************************************************
  LoadCondition::LoadCondition( LoadCondition const& rOther )
    : BoundaryCondition(rOther)
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  Condition::Pointer LoadCondition::Create(IndexType NewId,
					   NodesArrayType const& ThisNodes,
					   PropertiesType::Pointer pProperties) const
  {
    return Kratos::make_shared<LoadCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
  }


  //************************************CLONE*******************************************
  //************************************************************************************

  Condition::Pointer LoadCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {
    std::cout<<" Call base class LOAD CONDITION Clone "<<std::endl;

    LoadCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewCondition.SetData(this->GetData());
    NewCondition.SetFlags(this->GetFlags());


    return Kratos::make_shared<LoadCondition>(NewCondition);
  }


  //***********************************************************************************
  //***********************************************************************************
  LoadCondition::~LoadCondition()
  {
  }


  //************************************************************************************
  //************************************************************************************

  void LoadCondition::InitializeConditionVariables(ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    BoundaryCondition::InitializeConditionVariables(rVariables, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void LoadCondition::CalculateExternalLoad(ConditionVariables& rVariables)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class CalculateExternalLoad method for a load condition... " << std::endl;

    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void LoadCondition::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
						    ConditionVariables& rVariables,
						    double& rIntegrationWeight)

  {
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType& dimension = GetGeometry().WorkingSpaceDimension();

    unsigned int index = 0;
    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
        index = dimension * i;

        for ( SizeType j = 0; j < dimension; j++ )
	  {
	    rRightHandSideVector[index + j] += rVariables.N[i] * rVariables.ExternalVectorValue[j] * rIntegrationWeight;
	  }

      }

    if( this->HasVariableDof(ROTATION) ){

      // external load
      Vector ExternalLoad(3);
      noalias(ExternalLoad) = ZeroVector(3);

      Vector ExternalCouple(3);
      noalias(ExternalCouple) = ZeroVector(3);

      Matrix SkewSymMatrix(3,3);
      noalias(SkewSymMatrix) = ZeroMatrix(3,3);

      Vector IntegrationPointPosition(3);
      noalias(IntegrationPointPosition) = ZeroVector(3);

      Vector CurrentValueVector(3);
      noalias(CurrentValueVector) = ZeroVector(3);

      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  CurrentValueVector = GetGeometry()[i].Coordinates();
	  IntegrationPointPosition += rVariables.N[i] * CurrentValueVector;
	}

      unsigned int RowIndex = 0;
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  RowIndex = i * ( (dimension-1) * 3 );

	  ExternalLoad  = rIntegrationWeight * rVariables.N[i] * rVariables.ExternalVectorValue * rIntegrationWeight;

	  //integration for the moment force vector
	  BeamMathUtils<double>::VectorToSkewSymmetricTensor(ExternalLoad,SkewSymMatrix); // m = f x r = skewF · r
	  CurrentValueVector = GetGeometry()[i].Coordinates();
	  CurrentValueVector -= IntegrationPointPosition;
	  ExternalCouple = prod(SkewSymMatrix,CurrentValueVector);

	  if( dimension == 2 ){
	    ExternalLoad[2] = ExternalCouple[2];
	    BeamMathUtils<double>::AddVector(ExternalCouple,  rRightHandSideVector, RowIndex);
	  }
	  else{
	    BeamMathUtils<double>::AddVector(ExternalLoad,  rRightHandSideVector, RowIndex);
	    BeamMathUtils<double>::AddVector(ExternalCouple, rRightHandSideVector, RowIndex+dimension);
	  }
	}


    }

    //std::cout<<" ExternalForces ["<<this->Id()<<"]"<<rRightHandSideVector<<std::endl;

    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  double& LoadCondition::CalculateAndAddExternalEnergy(double& rEnergy,
						       ConditionVariables& rVariables,
						       double& rIntegrationWeight,
						       const ProcessInfo& rCurrentProcessInfo)

  {
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType& dimension = GetGeometry().WorkingSpaceDimension();

    // Energy Calculation:
    Vector CurrentValueVector(dimension);
    noalias(CurrentValueVector) = ZeroVector(dimension);
    Vector Displacements(dimension);
    noalias(Displacements) = ZeroVector(dimension);
    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	//current displacements to compute energy
	CurrentValueVector = GetNodalCurrentValue( DISPLACEMENT, CurrentValueVector, i );

	Displacements += rVariables.N[i] * CurrentValueVector;
      }
    //------

    Vector ForceVector(dimension);
    noalias(ForceVector) = ZeroVector(dimension);
    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	ForceVector += rVariables.N[i] * rVariables.ExternalVectorValue * rIntegrationWeight;
      }

    rEnergy += inner_prod( ForceVector, Displacements );

    return rEnergy;

    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************


  int LoadCondition::Check( const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    // Perform base condition checks
    int ErrorCode = 0;
    ErrorCode = BoundaryCondition::Check(rCurrentProcessInfo);

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);

    return ErrorCode;

    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void LoadCondition::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BoundaryCondition )
  }

  void LoadCondition::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BoundaryCondition )
  }


} // Namespace Kratos.
