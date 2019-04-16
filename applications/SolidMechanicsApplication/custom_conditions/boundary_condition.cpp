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
#include "custom_conditions/boundary_condition.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

  /**
   * Flags related to the condition computation
   */
  KRATOS_CREATE_LOCAL_FLAG( BoundaryCondition, COMPUTE_RHS_VECTOR,                 0 );
  KRATOS_CREATE_LOCAL_FLAG( BoundaryCondition, COMPUTE_LHS_MATRIX,                 1 );

  //***********************************************************************************
  //***********************************************************************************
  BoundaryCondition::BoundaryCondition()
    : Condition()
  {
  }


  //***********************************************************************************
  //***********************************************************************************
  BoundaryCondition::BoundaryCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  BoundaryCondition::BoundaryCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
  {
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
  }

  //************************************************************************************
  //************************************************************************************
  BoundaryCondition::BoundaryCondition( BoundaryCondition const& rOther )
    : Condition(rOther)
    ,mThisIntegrationMethod(rOther.mThisIntegrationMethod)
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  Condition::Pointer BoundaryCondition::Create(IndexType NewId,
					       NodesArrayType const& ThisNodes,
					       PropertiesType::Pointer pProperties) const
  {
    return Kratos::make_shared<BoundaryCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
  }


  //************************************CLONE*******************************************
  //************************************************************************************

  Condition::Pointer BoundaryCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {
    std::cout<<" Call base class BOUNDARY CONDITION Clone "<<std::endl;

    BoundaryCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewCondition.SetData(this->GetData());
    NewCondition.SetFlags(this->GetFlags());


    return Kratos::make_shared<BoundaryCondition>(NewCondition);
  }


  //***********************************************************************************
  //***********************************************************************************
  BoundaryCondition::~BoundaryCondition()
  {
  }


  //************************************************************************************
  //************************************************************************************

  bool BoundaryCondition::HasVariableDof(VariableVectorType& rVariable)
  {
    KRATOS_TRY

    typedef VectorComponentAdaptor<array_1d<double,3> >  VectorComponentType;
    const VariableComponent<VectorComponentType>& var_x  = KratosComponents<VariableComponent<VectorComponentType> >::Get(rVariable.Name()+"_X");

    //usually if the dofs do not exist condition adds them, standard conditions do not work like this
    if( GetGeometry()[0].HasDofFor(var_x) == true )
      return true;

    return false;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  bool BoundaryCondition::HasVariableDof(VariableScalarType& rVariable)
  {
    KRATOS_TRY

    //usually if the dofs do not exist condition adds them, standard conditions do not work like this
    if( GetGeometry()[0].HasDofFor(rVariable) == true )
      return true;

    return false;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  unsigned int BoundaryCondition::GetDofsSize()
  {
    KRATOS_TRY

    const SizeType& dimension       = GetGeometry().WorkingSpaceDimension();
    const SizeType number_of_nodes = GetGeometry().PointsNumber();

    unsigned int size = number_of_nodes * dimension;

    if( HasVariableDof(ROTATION) ){
      if(dimension == 2){
	size += number_of_nodes;
      }
      else{
	size += number_of_nodes * dimension;
      }
    }

    return size;

    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void BoundaryCondition::GetDofList(DofsVectorType& rConditionDofList,
				     ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    rConditionDofList.resize(0);
    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType& dimension = GetGeometry().WorkingSpaceDimension();

    for (SizeType i = 0; i < number_of_nodes; i++)
      {
        if( HasVariableDof(DISPLACEMENT) ){
          rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
          rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
          if( dimension == 3 )
            rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }
        else if( HasVariableDof(VELOCITY) ){
          rConditionDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
          rConditionDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
          if( dimension == 3 )
            rConditionDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Z));
        }

	if( HasVariableDof(ROTATION) ){
	  if( dimension == 2 ){
	    rConditionDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));
	  }
	  else{
	    rConditionDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
	    rConditionDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
	    rConditionDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));
	  }
	}

      }


    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void BoundaryCondition::EquationIdVector(EquationIdVectorType& rResult,
					   ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType& dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = this->GetDofsSize();

    if ( rResult.size() != dofs_size )
      rResult.resize( dofs_size, false );

    unsigned int index = 0;

    if( HasVariableDof(ROTATION) && HasVariableDof(DISPLACEMENT) ){
      if( dimension == 2 ){
	for ( SizeType i = 0; i < number_of_nodes; i++ )
	  {
	    index = i * ( (dimension-1) * 3 );
	    rResult[index]   = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
	    rResult[index+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
	    rResult[index+2] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
	  }
      }
      else{
	for ( SizeType i = 0; i < number_of_nodes; i++ )
	  {
	    index = i * ( (dimension-1) * 3 );
	    rResult[index]   = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
	    rResult[index+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
	    rResult[index+2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();

	    rResult[index+3] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
	    rResult[index+4] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
	    rResult[index+5] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
	  }
      }
    }
    else if( HasVariableDof(DISPLACEMENT) ){

      for (SizeType i = 0; i < number_of_nodes; i++)
	{
	  index = i * dimension;
	  rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
	  rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
	  if( dimension == 3)
	    rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
	}
    }
    else if( HasVariableDof(VELOCITY) ){

      for (SizeType i = 0; i < number_of_nodes; i++)
	{
	  index = i * dimension;
	  rResult[index]     = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
	  rResult[index + 1] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
	  if( dimension == 3)
	    rResult[index + 2] = GetGeometry()[i].GetDof(VELOCITY_Z).EquationId();
	}
    }

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void BoundaryCondition::GetValuesVector(Vector& rValues, int Step)
  {
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType& dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = this->GetDofsSize();

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    unsigned int index = 0;

    if( HasVariableDof(ROTATION) ){

      if( dimension == 2 ){
	for ( SizeType i = 0; i < number_of_nodes; i++ )
	  {
	    index = i * ( (dimension-1) * 3 );
	    rValues[index]   = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
	    rValues[index+1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );
	    rValues[index+2] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Z, Step );
	  }
      }
      else{
	for ( SizeType i = 0; i < number_of_nodes; i++ )
	  {
	    index = i * ( (dimension-1) * 3 );
	    rValues[index]   = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
	    rValues[index+1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );
	    rValues[index+2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );

	    rValues[index+3] = GetGeometry()[i].GetSolutionStepValue( ROTATION_X, Step );
	    rValues[index+4] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Y, Step );
	    rValues[index+5] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Z, Step );
	  }
      }

    }
    else{

      for (SizeType i = 0; i < number_of_nodes; i++)
	{
	  index = i * dimension;
	  rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
	  rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

	  if ( dimension == 3 )
	    rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
	}

    }

    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void BoundaryCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
  {
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType& dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = this->GetDofsSize();

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    unsigned int index = 0;

    if( HasVariableDof(ROTATION) ){

      if( dimension == 2 ){
	for ( SizeType i = 0; i < number_of_nodes; i++ )
	  {
	    index = i * ( (dimension-1) * 3 );
	    rValues[index]   = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
	    rValues[index+1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );
	    rValues[index+2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
	  }
      }
      else{
	for ( SizeType i = 0; i < number_of_nodes; i++ )
	  {
	    index = i * ( (dimension-1) * 3 );
	    rValues[index]   = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
	    rValues[index+1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );
	    rValues[index+2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );

	    rValues[index+3] = GetGeometry()[i].GetSolutionStepValue( ROTATION_X, Step );
	    rValues[index+4] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Y, Step );
	    rValues[index+5] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Z, Step );
	  }
      }

    }
    else{

      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  index = i * dimension;
	  rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
	  rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

	  if ( dimension == 3 )
	    rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
	}

    }

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void BoundaryCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
  {
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType& dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = this->GetDofsSize();

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    unsigned int index = 0;

    if( HasVariableDof(ROTATION) ){

      if( dimension == 2 ){

	for ( SizeType i = 0; i < number_of_nodes; i++ )
	  {
	    index = i * ( (dimension-1) * 3 );
	    rValues[index]   = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
	    rValues[index+1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );
	    rValues[index+2] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_Z, Step );
	  }
      }
      else{
	for ( SizeType i = 0; i < number_of_nodes; i++ )
	  {
	    index = i * ( (dimension-1) * 3 );
	    rValues[index]   = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
	    rValues[index+1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );
	    rValues[index+2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );

	    rValues[index+3] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_X, Step );
	    rValues[index+4] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_Y, Step );
	    rValues[index+5] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_Z, Step );
	  }
      }

    }
    else{

      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  index = i * dimension;
	  rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
	  rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

	  if ( dimension == 3 )
	    rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
	}
    }

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void BoundaryCondition::InitializeExplicitContributions()
  {
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	if( GetGeometry()[i].SolutionStepsDataHas(EXTERNAL_FORCE) && GetGeometry()[i].SolutionStepsDataHas(FORCE_RESIDUAL) ){

	  array_1d<double, 3 > & ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
	  array_1d<double, 3 > & ResidualForce = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

	  GetGeometry()[i].SetLock();
	  ExternalForce.clear();
	  ResidualForce.clear();
	  GetGeometry()[i].UnSetLock();
	}

	if( HasVariableDof(ROTATION) ){

	  if( GetGeometry()[i].SolutionStepsDataHas(EXTERNAL_MOMENT) && GetGeometry()[i].SolutionStepsDataHas(MOMENT_RESIDUAL) ){

	    array_1d<double, 3 > & ExternalMoment = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_MOMENT);
	    array_1d<double, 3 > & ResidualMoment = GetGeometry()[i].FastGetSolutionStepValue(MOMENT_RESIDUAL);

	    GetGeometry()[i].SetLock();
	    ExternalMoment.clear();
	    ResidualMoment.clear();
	    GetGeometry()[i].UnSetLock();

	  }

	}

      }

    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void BoundaryCondition::AddExplicitContribution(const VectorType& rRHS,
						  const Variable<VectorType>& rRHSVariable,
						  Variable<array_1d<double,3> >& rDestinationVariable,
						  const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType& dimension       = GetGeometry().WorkingSpaceDimension();

    if( rRHSVariable == EXTERNAL_FORCES_VECTOR && rDestinationVariable == EXTERNAL_FORCE )
      {

	for(SizeType i=0; i< number_of_nodes; i++)
	  {
	    int index = dimension * i;

	    GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
	    for(SizeType j=0; j<dimension; j++)
	      {
		ExternalForce[j] += rRHS[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }

    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL )
      {

	for(SizeType i=0; i< number_of_nodes; i++)
	  {
	    int index = dimension * i;

	    GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &ForceResidual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
	    for(SizeType j=0; j<dimension; j++)
	      {
		ForceResidual[j] += rRHS[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }


    if( HasVariableDof(ROTATION) ){

      if( rRHSVariable == EXTERNAL_FORCES_VECTOR && rDestinationVariable == EXTERNAL_MOMENT )
	{

	  for(SizeType i=0; i< number_of_nodes; i++)
	    {
	      int index = dimension * i;

	      GetGeometry()[i].SetLock();

	      array_1d<double, 3 > &ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
	      for(SizeType j=0; j<dimension; j++)
		{
		  ExternalForce[j] += rRHS[index + j];
		}

	      GetGeometry()[i].UnSetLock();
	    }
	}

      if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == MOMENT_RESIDUAL )
	{

	  for(SizeType i=0; i< number_of_nodes; i++)
	    {
	      int index = dimension * i;

	      GetGeometry()[i].SetLock();

	      array_1d<double, 3 > &MomentResidual = GetGeometry()[i].FastGetSolutionStepValue(MOMENT_RESIDUAL);
	      for(SizeType j=0; j<dimension; j++)
		{
		  MomentResidual[j] += rRHS[index + j];
		}

	      GetGeometry()[i].UnSetLock();
	    }
	}
    }


    KRATOS_CATCH( "" )
  }

  //************* STARTING - ENDING  METHODS
  //***********************************************************************************
  //***********************************************************************************

  void BoundaryCondition::Initialize()
  {
    KRATOS_TRY

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void BoundaryCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    InitializeExplicitContributions();

    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void BoundaryCondition::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void BoundaryCondition::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
						   VectorType& rRightHandSideVector,
						   Flags& rCalculationFlags)

  {

    //resizing as needed the LHS
    unsigned int MatSize = this->GetDofsSize();

    if ( rCalculationFlags.Is(BoundaryCondition::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
      {
        if ( rLeftHandSideMatrix.size1() != MatSize )
	  rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
      }


    //resizing as needed the RHS
    if ( rCalculationFlags.Is(BoundaryCondition::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
      {
        if ( rRightHandSideVector.size() != MatSize )
	  rRightHandSideVector.resize( MatSize, false );

	noalias(rRightHandSideVector) = ZeroVector( MatSize ); //resetting RHS

      }
  }


  //************************************************************************************
  //************************************************************************************

  void BoundaryCondition::InitializeConditionVariables(ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().PointsNumber();
    const SizeType& local_dimension = GetGeometry().LocalSpaceDimension();
    const SizeType& dimension       = GetGeometry().WorkingSpaceDimension();

    rVariables.Initialize(dimension, local_dimension, number_of_nodes);

    //set variables including all integration points values

    //reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    //reading shape functions local gradients
    rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

    KRATOS_CATCH( "" )

  }

  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************

  void BoundaryCondition::CalculateKinematics(ConditionVariables& rVariables,
					      const double& rPointNumber)
  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class CalculateKinematics method for a boundary condition... " << std::endl;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void BoundaryCondition::CalculateConditionSystem(LocalSystemComponents& rLocalSystem,
						   const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    //create and initialize condition variables:
    ConditionVariables Variables;
    this->InitializeConditionVariables(Variables,rCurrentProcessInfo);

    //reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    double IntegrationWeight = 1;

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {
        //compute element kinematics B, F, DN_DX ...
        this->CalculateKinematics(Variables,PointNumber);

        //calculating weights for integration on the "reference configuration" Jacobian respect to the reference configuration
	//take in account in a linear element (Jacobian=2*Area) this is the relation

        IntegrationWeight = Variables.Jacobian * integration_points[PointNumber].Weight();

	// std::cout<<" Variables.Jacobian "<<Variables.Jacobian<<" Weight "<<integration_points[PointNumber].Weight()<<" / "<<std::endl;

        if ( rLocalSystem.CalculationFlags.Is(BoundaryCondition::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
	  {
            //contributions to stiffness matrix calculated on the reference config
	    this->CalculateAndAddLHS ( rLocalSystem, Variables, IntegrationWeight );
	  }

        if ( rLocalSystem.CalculationFlags.Is(BoundaryCondition::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
	  {
            //contribution to external forces
	    this->CalculateAndAddRHS ( rLocalSystem, Variables, IntegrationWeight );
	  }

      }

    //VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
    //std::cout<<" rRightHandSideVector "<<rRightHandSideVector<<std::endl;

    KRATOS_CATCH( "" )
  }


  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************

  void BoundaryCondition::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ConditionVariables& rVariables, double& rIntegrationWeight)
  {

    //contributions of the stiffness matrix calculated on the reference configuration
    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

    // operation performed: add Kg to the rLefsHandSideMatrix
    this->CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    //KRATOS_WATCH( rLeftHandSideMatrix )

  }


  //************************************************************************************
  //************************************************************************************

  void BoundaryCondition::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ConditionVariables& rVariables, double& rIntegrationWeight)
  {
    //contribution of the internal and external forces

    VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();

    // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
    this->CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rIntegrationWeight );

    //std::cout<<" rRightHandSideVectorPart "<<rRightHandSideVector<<std::endl;

  }


  //************************************************************************************
  //************************************************************************************

  void BoundaryCondition::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
  {
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(BoundaryCondition::COMPUTE_LHS_MATRIX);

    VectorType RightHandSideVector = Vector();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(RightHandSideVector);

    //Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );

    //KRATOS_WATCH( rRightHandSideVector )

  }

  //************************************************************************************
  //************************************************************************************

  void BoundaryCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
  {
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(BoundaryCondition::COMPUTE_RHS_VECTOR);

    MatrixType LeftHandSideMatrix = Matrix();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    //Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );

  }


  //************************************************************************************
  //************************************************************************************

  void BoundaryCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
  {
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(BoundaryCondition::COMPUTE_LHS_MATRIX);
    LocalSystem.CalculationFlags.Set(BoundaryCondition::COMPUTE_RHS_VECTOR);

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    //Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );

  }


  //***********************************************************************************
  //***********************************************************************************

  void BoundaryCondition::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    rMassMatrix.resize(0, 0, false);

    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void BoundaryCondition::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    rDampingMatrix.resize(0, 0, false);

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void BoundaryCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
					      ConditionVariables& rVariables,
					      double& rIntegrationWeight)

  {
    KRATOS_TRY

    unsigned int MatSize = this->GetDofsSize();
    if(rLeftHandSideMatrix.size1() != MatSize)
      rLeftHandSideMatrix.resize(MatSize,MatSize,false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize,MatSize);

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void BoundaryCondition::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
							ConditionVariables& rVariables,
							double& rIntegrationWeight)

  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class CalculateAndAddExternalForces method for a boundary condition... " << std::endl;

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  double& BoundaryCondition::CalculateAndAddExternalEnergy(double& rEnergy,
							   ConditionVariables& rVariables,
							   double& rIntegrationWeight,
							   const ProcessInfo& rCurrentProcessInfo)

  {
    KRATOS_TRY

    KRATOS_ERROR << "calling the base class CalculateAndAddExternalEnergy method for a boundary condition... " << std::endl;

    return rEnergy;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void BoundaryCondition::GetNodalDeltaMovements(Vector& rValues, const int& rNode)
  {
    KRATOS_TRY

    const SizeType& dimension = GetGeometry().WorkingSpaceDimension();

    if( rValues.size() != dimension )
      rValues.resize(dimension);

    noalias(rValues) = ZeroVector(dimension);

    Vector CurrentValueVector(3);
    noalias(CurrentValueVector) = ZeroVector(3);
    CurrentValueVector = GetNodalCurrentValue( DISPLACEMENT, CurrentValueVector, rNode );

    Vector PreviousValueVector(3);
    noalias(PreviousValueVector)= ZeroVector(3);
    CurrentValueVector = GetNodalPreviousValue( DISPLACEMENT, CurrentValueVector, rNode );


    rValues[0] = CurrentValueVector[0] - PreviousValueVector[0];
    rValues[1] = CurrentValueVector[1] - PreviousValueVector[1];

    if( dimension == 3 )
      rValues[2] = CurrentValueVector[2] - PreviousValueVector[2];

    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************

  Vector& BoundaryCondition::GetNodalCurrentValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const unsigned int& rNode)
  {
    KRATOS_TRY

    const SizeType& dimension       = GetGeometry().WorkingSpaceDimension();

    if( rValue.size() != dimension )
      rValue.resize(dimension, false);

    rValue = GetGeometry()[rNode].FastGetSolutionStepValue( rVariable );

    return rValue;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  Vector& BoundaryCondition::GetNodalPreviousValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const unsigned int& rNode)
  {
    KRATOS_TRY

    const SizeType& dimension       = GetGeometry().WorkingSpaceDimension();

    if( rValue.size() != dimension )
      rValue.resize(dimension, false);

    rValue = GetGeometry()[rNode].FastGetSolutionStepValue( rVariable, 1 );

    return rValue;

    KRATOS_CATCH( "" )
  }


  //*********************************GET DOUBLE VALUE***********************************
  //************************************************************************************

  void BoundaryCondition::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
						       std::vector<double>& rValues,
						       const ProcessInfo& rCurrentProcessInfo )
  {
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
  }

  //************************************************************************************
  //************************************************************************************

  void BoundaryCondition::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo )
  {

    KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    unsigned int integration_points = 0;
    if( integration_points_number == 0 )
      integration_points = 1;

    if ( rOutput.size() != integration_points )
      rOutput.resize( integration_points, false );


    if ( rVariable == EXTERNAL_ENERGY )
      {

	//create and initialize condition variables:
	ConditionVariables Variables;
	this->InitializeConditionVariables(Variables,rCurrentProcessInfo);

	//reading integration points
	const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

	double Energy = 0;
	double IntegrationWeight = 0;

	for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
	  {
	    //compute element kinematics B, F, DN_DX ...
	    this->CalculateKinematics(Variables,PointNumber);

	    IntegrationWeight = Variables.Jacobian * integration_points[PointNumber].Weight();

	    Energy = 0;

	    Energy = this->CalculateAndAddExternalEnergy( Energy, Variables, IntegrationWeight, rCurrentProcessInfo);
	    rOutput[PointNumber] = Energy;
	  }

      }

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************


  int BoundaryCondition::Check( const ProcessInfo& rCurrentProcessInfo )
  {

    // Perform base element checks
    int ErrorCode = 0;
    ErrorCode = Condition::Check(rCurrentProcessInfo);

    return ErrorCode;
  }

  //***********************************************************************************
  //***********************************************************************************

  void BoundaryCondition::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
    int IntMethod = int(mThisIntegrationMethod);
    rSerializer.save("IntegrationMethod",IntMethod);

  }

  void BoundaryCondition::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);

  }


} // Namespace Kratos.
