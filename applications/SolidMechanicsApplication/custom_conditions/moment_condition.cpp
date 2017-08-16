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
#include "custom_conditions/moment_condition.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

  //************************************************************************************
  //************************************************************************************
  MomentCondition::MomentCondition()
    : BoundaryCondition()
  {
  }
  
  //************************************************************************************
  //************************************************************************************
  MomentCondition::MomentCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : BoundaryCondition(NewId, pGeometry)
  {
  }

  //************************************************************************************
  //************************************************************************************
  MomentCondition::MomentCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : BoundaryCondition(NewId, pGeometry, pProperties)
  {
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
  }

  //************************************************************************************
  //************************************************************************************
  MomentCondition::MomentCondition( MomentCondition const& rOther )
    : BoundaryCondition(rOther)
  {
  }

  //************************************************************************************
  //************************************************************************************
  Condition::Pointer MomentCondition::Create(IndexType NewId,
					     NodesArrayType const& ThisNodes,
					     PropertiesType::Pointer pProperties) const
  {
    return Condition::Pointer(new MomentCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
  }

  
  //************************************CLONE*******************************************
  //************************************************************************************
  
  Condition::Pointer MomentCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {
    std::cout<<" Call base class MOMENT CONDITION Clone "<<std::endl;
  
    MomentCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewCondition.SetData(this->GetData());
    NewCondition.SetFlags(this->GetFlags());
  
    return Condition::Pointer( new MomentCondition(NewCondition) );
  }
  
  //************************************************************************************
  //************************************************************************************
  MomentCondition::~MomentCondition()
  {
  }


  //************************************************************************************
  //************************************************************************************
  unsigned int MomentCondition::GetDofsSize()
  {
    KRATOS_TRY
     
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();    

    if(dimension == 2)
      {
	return number_of_nodes;
      }
    else
      {
	unsigned int size = number_of_nodes * dimension;
	return size;
      }
    
    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void MomentCondition::GetDofList(DofsVectorType& rConditionDofList,
				     ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    rConditionDofList.resize(0);
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    for (unsigned int i = 0; i < number_of_nodes; i++)
      {
	if( dimension == 2 ){
	  rConditionDofList.push_back(GetGeometry()[i].pGetDof(PLANE_ROTATION));
	}
	else{
	  rConditionDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
	  rConditionDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
	  rConditionDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));
	}
      }


    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void MomentCondition::EquationIdVector(EquationIdVectorType& rResult,
					  ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int condition_size        = number_of_nodes * dimension;

    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    unsigned int index = 0;
    for (unsigned int i = 0; i < number_of_nodes; i++)
      {
        if( dimension == 2 ){
	  rResult[i]     = GetGeometry()[i].GetDof(PLANE_ROTATION).EquationId();
	}
	else{
	  index = i * dimension;
	  rResult[index]     = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
	  rResult[index + 1] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
	  rResult[index + 2] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
	}
      }

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void MomentCondition::GetValuesVector(Vector& rValues, int Step)
  {
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       condition_size  = number_of_nodes * dimension;

    if ( rValues.size() != condition_size ) 
      rValues.resize( condition_size, false );

    unsigned int index = 0;
    for (unsigned int i = 0; i < number_of_nodes; i++)
      {
	if ( dimension == 2 ){
	  rValues[i]     = GetGeometry()[i].GetSolutionStepValue( PLANE_ROTATION, Step );
	}
	else{
	  index = i * dimension;
	  rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ROTATION_X, Step );
	  rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Y, Step );
	  rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Z, Step );
	}
      }
  }

  //***********************************************************************************
  //***********************************************************************************

  void MomentCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
  {
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       condition_size  = number_of_nodes * dimension;

    if ( rValues.size() != condition_size )
      rValues.resize( condition_size, false );
    
    unsigned int index = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	if ( dimension == 2 ){
	  rValues[i]     = GetGeometry()[i].GetSolutionStepValue( PLANE_ANGULAR_VELOCITY, Step );
	}
	else{
	  index = i * dimension;
	  rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_X, Step );
	  rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_Y, Step );
	  rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_Z, Step );
	}
      }
  }
  
  //***********************************************************************************
  //***********************************************************************************

  void MomentCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
  {
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       condition_size    = number_of_nodes * dimension;

    if ( rValues.size() != condition_size )
      rValues.resize( condition_size, false );
      
    unsigned int index = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	if ( dimension == 2 ){
	  rValues[i]     = GetGeometry()[i].GetSolutionStepValue( PLANE_ANGULAR_VELOCITY, Step );
	}
	else{
	  index = i * dimension;
	  rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_X, Step );
	  rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_Y, Step );
	  rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_Z, Step );
	}
	  
      }
  }


  //************************************************************************************
  //************************************************************************************
  void MomentCondition::InitializeExplicitContributions()
  {
    KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	if( GetGeometry()[i].SolutionStepsDataHas(EXTERNAL_MOMENT) && GetGeometry()[i].SolutionStepsDataHas(MOMENT_RESIDUAL) ){
	  
	  array_1d<double, 3 > & ExternalMoment = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_MOMENT);
	  array_1d<double, 3 > & ResidualMoment = GetGeometry()[i].FastGetSolutionStepValue(MOMENT_RESIDUAL);
  
	  GetGeometry()[i].SetLock();
	  ExternalMoment.clear();
	  ResidualMoment.clear();
	  GetGeometry()[i].UnSetLock();

	}

      }

    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void MomentCondition::AddExplicitContribution(const VectorType& rRHS, 
						  const Variable<VectorType>& rRHSVariable,
						  Variable<array_1d<double,3> >& rDestinationVariable, 
						  const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if( rRHSVariable == EXTERNAL_FORCES_VECTOR && rDestinationVariable == EXTERNAL_MOMENT )
      {

	for(unsigned int i=0; i< number_of_nodes; i++)
	  {
	    int index = dimension * i;

	    GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
	    for(unsigned int j=0; j<dimension; j++)
	      {
		ExternalForce[j] += rRHS[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }

    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == MOMENT_RESIDUAL )
      {

	for(unsigned int i=0; i< number_of_nodes; i++)
	  {
	    int index = dimension * i;

	    GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &MomentResidual = GetGeometry()[i].FastGetSolutionStepValue(MOMENT_RESIDUAL);
	    for(unsigned int j=0; j<dimension; j++)
	      {
		MomentResidual[j] += rRHS[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }

    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void MomentCondition::CalculateExternalMoment(ConditionVariables& rVariables)
  {
    KRATOS_TRY
      
    KRATOS_ERROR << "calling the base class CalculateExternalMoment method for a moment condition... " << std::endl;
    
    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void MomentCondition::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
						    ConditionVariables& rVariables,
						    double& rIntegrationWeight)

  {
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    int index = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        index = dimension * i;
	
        for ( unsigned int j = 0; j < dimension; j++ )
	  {
	    rRightHandSideVector[index + j] += rVariables.N[i] * rVariables.ExternalVectorValue[j] * rIntegrationWeight;
	  }

      }

    std::cout<<" ExternalForces ["<<this->Id()<<"]"<<rRightHandSideVector<<std::endl;

    KRATOS_CATCH( "" )
  }  


  //***********************************************************************************
  //***********************************************************************************

  double& MomentCondition::CalculateAndAddExternalEnergy(double& rEnergy,
							 ConditionVariables& rVariables,
							 double& rIntegrationWeight,
							 const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Energy Calculation:
    Vector CurrentValueVector(dimension);
    noalias(CurrentValueVector) = ZeroVector(dimension); 
    Vector PreviousValueVector(dimension);
    noalias(PreviousValueVector) = ZeroVector(dimension); 
    Vector Rotations(dimension);
    noalias(Rotations) = ZeroVector(dimension);
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	//current rotation to compute energy
	if( GetGeometry()[i].SolutionStepsDataHas(ANGULAR_VELOCITY) ){
	  CurrentValueVector  = GetNodalCurrentValue( ANGULAR_VELOCITY, CurrentValueVector, i );
	  PreviousValueVector = GetNodalPreviousValue( ANGULAR_VELOCITY, PreviousValueVector, i );
	  CurrentValueVector += PreviousValueVector;
	  CurrentValueVector *= 0.5 * rCurrentProcessInfo[DELTA_TIME];
	}
	else{
	  CurrentValueVector = GetNodalCurrentValue( ROTATION, CurrentValueVector, i );
	}
	
	Rotations += rVariables.N[i] * CurrentValueVector;
      }
    //------

    Vector MomentVector(dimension);
    noalias(MomentVector) = ZeroVector(dimension);
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	MomentVector += rVariables.N[i] * rVariables.ExternalVectorValue * rIntegrationWeight;
      }

    rEnergy += inner_prod( MomentVector, Rotations );

    return rEnergy;

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void MomentCondition::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BoundaryCondition )
  }
  
  void MomentCondition::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BoundaryCondition )
  }


} // Namespace Kratos


