//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/rigid_body_element.hpp"
#include "contact_mechanics_application_variables.h"

namespace Kratos
{

/**
 * Flags related to the element computation
 */
KRATOS_CREATE_LOCAL_FLAG( RigidBodyElement, COMPUTE_RHS_VECTOR, 0 );
KRATOS_CREATE_LOCAL_FLAG( RigidBodyElement, COMPUTE_LHS_MATRIX, 1 );

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

RigidBodyElement::RigidBodyElement(IndexType NewId,GeometryType::Pointer pGeometry)
    :Element(NewId, pGeometry)
{
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

RigidBodyElement::RigidBodyElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    :Element(NewId, pGeometry, pProperties)
{
    KRATOS_TRY

    this->Set(RIGID);

    KRATOS_CATCH("")
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************


RigidBodyElement::RigidBodyElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, NodesContainerType::Pointer pNodes)
    :Element(NewId, pGeometry, pProperties)
{
    KRATOS_TRY

    mpNodes = pNodes;

    KRATOS_CATCH("")
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

RigidBodyElement::RigidBodyElement(RigidBodyElement const& rOther)
    :Element(rOther)
    ,mInitialLocalQuaternion(rOther.mInitialLocalQuaternion)
    ,mpNodes(rOther.mpNodes)
{
}

//*********************************CREATE*********************************************
//************************************************************************************

Element::Pointer RigidBodyElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
  return Kratos::make_shared<RigidBodyElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}


//*********************************CLONE**********************************************
//************************************************************************************

Element::Pointer RigidBodyElement::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
{
  RigidBodyElement NewElement( NewId, GetGeometry().Create(ThisNodes), pGetProperties(), mpNodes );

  NewElement.mInitialLocalQuaternion = this->mInitialLocalQuaternion;
  NewElement.SetData(this->GetData());
  NewElement.SetFlags(this->GetFlags());

  return Kratos::make_shared<RigidBodyElement>(NewElement);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

RigidBodyElement::~RigidBodyElement()
{
}

//************************************************************************************
//************************************************************************************

void RigidBodyElement::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
{

    ElementalDofList.resize(0);

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension       = GetGeometry().WorkingSpaceDimension();

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
	ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
	if( dimension == 2 ){
	  ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));
	}
	else{
	  ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
	  ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
	  ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
	  ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));
	}
      }

}

//************************************************************************************
//************************************************************************************

void RigidBodyElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    SizeType dofs_size = this->GetDofsSize();

    if ( rResult.size() != dofs_size )
        rResult.resize(dofs_size, false);

    dofs_size = dimension * (dimension + 1) * 0.5;
    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
      SizeType index = i * (dofs_size);
      rResult[index]   = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
      rResult[index+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
      if( dimension == 2 ){
	rResult[index+2] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
      }
      else{
	rResult[index+2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
	rResult[index+3] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
	rResult[index+4] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
	rResult[index+5] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
      }
    }

}


//*********************************DISPLACEMENT***************************************
//************************************************************************************

void RigidBodyElement::GetValuesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    SizeType dofs_size = this->GetDofsSize();

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    dofs_size = dimension * (dimension + 1) * 0.5;
    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
      SizeType index = i * (dofs_size);
      rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
      rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );
      if( dimension == 2 ){
	rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Z, Step );
      }
      else{
	rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
	rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( ROTATION_X, Step );
	rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Y, Step );
	rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Z, Step );
      }
    }

    KRATOS_CATCH("")
}

//************************************VELOCITY****************************************
//************************************************************************************

void RigidBodyElement::GetFirstDerivativesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    SizeType dofs_size = this->GetDofsSize();

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    dofs_size = dimension * (dimension + 1) * 0.5;
    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
      SizeType index = i * (dofs_size);
      rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
      rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );
      if( dimension == 2 ){
	rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_Z, Step );
      }
      else{
	rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
	rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_X, Step );
	rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_Y, Step );
	rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_Z, Step );
      }
    }


    KRATOS_CATCH("")
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void RigidBodyElement::GetSecondDerivativesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    SizeType dofs_size = this->GetDofsSize();

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    dofs_size = dimension * (dimension + 1) * 0.5;
    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
      SizeType index = i * (dofs_size);
      rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
      rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );
      if( dimension == 2 ){
	rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_Z, Step );
      }
      else{
	rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
	rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_X, Step );
	rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_Y, Step );
	rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_Z, Step );
      }

    }


    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyElement::Initialize()
{
    KRATOS_TRY

    Matrix& LocalAxesMatrix = GetProperties()[LOCAL_AXES_MATRIX];

    mInitialLocalQuaternion = QuaternionType::FromRotationMatrix( LocalAxesMatrix );

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyElement::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyElement::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
     KRATOS_TRY

     KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyElement::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
     KRATOS_TRY

     this->UpdateRigidBodyNodes(rCurrentProcessInfo);

     KRATOS_CATCH("")
}



//************************************************************************************
//************************************************************************************

void RigidBodyElement::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
     KRATOS_TRY

     KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
						VectorType& rRightHandSideVector,
						Flags& rCalculationFlags)

{
    //resizing as needed the LHS
    SizeType MatSize               = this->GetDofsSize();

    if ( rCalculationFlags.Is(RigidBodyElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( rCalculationFlags.Is(RigidBodyElement::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
	    rRightHandSideVector.resize( MatSize, false );

	noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
    }
}

//************************************************************************************
//************************************************************************************

Vector& RigidBodyElement::GetNodalCurrentValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const unsigned int& rNode)
{
    KRATOS_TRY

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    array_1d<double,3> ArrayValue;
    ArrayValue = GetGeometry()[rNode].FastGetSolutionStepValue( rVariable );

    if( rValue.size() != dimension )
      rValue.resize(dimension, false);

    for( SizeType i=0; i<dimension; i++ )
    {
      rValue[i] = ArrayValue[i];
    }

    return rValue;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

Vector& RigidBodyElement::GetNodalPreviousValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const unsigned int& rNode)
{
    KRATOS_TRY

    const SizeType dimension       = GetGeometry().WorkingSpaceDimension();

    array_1d<double,3> ArrayValue;
    ArrayValue = GetGeometry()[rNode].FastGetSolutionStepValue( rVariable, 1 );

    if( rValue.size() != dimension )
      rValue.resize(dimension, false);

    for( SizeType i=0; i<dimension; i++ )
    {
      rValue[i] = ArrayValue[i];
    }

    return rValue;

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyElement::CalculateRigidBodyProperties(RigidBodyProperties & rRigidBody)
{
    KRATOS_TRY


    if( GetProperties().Has(NODAL_MASS) ){
        rRigidBody.Mass = GetProperties()[NODAL_MASS];
    }


    if( GetProperties().Has(LOCAL_INERTIA_TENSOR) )
    {
        Matrix& inertia = GetProperties()[LOCAL_INERTIA_TENSOR];

        rRigidBody.InertiaTensor = inertia;
    }

    KRATOS_CATCH("")

}

//************************************************************************************
//************************************************************************************

void RigidBodyElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
					      ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(RigidBodyElement::COMPUTE_RHS_VECTOR);

    MatrixType LeftHandSideMatrix = Matrix();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices(LeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags);

    //Calculate elemental system
    this->CalculateElementalSystem(LocalSystem, rCurrentProcessInfo);

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                             ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(RigidBodyElement::COMPUTE_LHS_MATRIX);

    VectorType RightHandSideVector = Vector();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices(rLeftHandSideMatrix, RightHandSideVector,  LocalSystem.CalculationFlags);

    //Calculate elemental system
    this->CalculateElementalSystem(LocalSystem, rCurrentProcessInfo);

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(RigidBodyElement::COMPUTE_RHS_VECTOR);
    LocalSystem.CalculationFlags.Set(RigidBodyElement::COMPUTE_LHS_MATRIX);

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices(rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags);

    //Calculate elemental system
    this->CalculateElementalSystem(LocalSystem, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyElement::CalculateElementalSystem(LocalSystemComponents& rLocalSystem,
                                                ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
        
    KRATOS_ERROR << " calling the default method CalculateElementalSystem for a rigid body element " << std::endl;
    
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyElement::CalculateDynamicSystem(LocalSystemComponents& rLocalSystem,
                                              ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    //create and initialize element variables:
    ElementVariables Variables;
    Variables.Initialize(dimension,rCurrentProcessInfo);

    std::cout<<" RigidBodyElement "<<this->Id()<<std::endl;
    std::cout<<" [Displacement "<<GetGeometry()[0].FastGetSolutionStepValue( DISPLACEMENT )<<"]"<<std::endl;
    std::cout<<" [Rotation     "<<GetGeometry()[0].FastGetSolutionStepValue( ROTATION )<<"]"<<std::endl;

    //Compute Rigid Body Properties:
    this->CalculateRigidBodyProperties(Variables.RigidBody);

    // initialize variables short version;

    if ( rLocalSystem.CalculationFlags.Is(RigidBodyElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
      {
	MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

	this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);
      }

    if ( rLocalSystem.CalculationFlags.Is(RigidBodyElement::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
      {
	VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();

	this->CalculateAndAddRHS(rRightHandSideVector, Variables);
      }


    // Note:
    // That means that the standard rotation K = Q·K'·QT and F = Q·F' is the correct transformation

    Matrix InitialLocalMatrix = ZeroMatrix(3,3);
    mInitialLocalQuaternion.ToRotationMatrix(InitialLocalMatrix);

    // Transform Local to Global LHSMatrix:
    if ( rLocalSystem.CalculationFlags.Is(RigidBodyElement::COMPUTE_LHS_MATRIX) ){

      MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

      // works for 2D and 3D case
      BeamMathUtilsType::MapLocalToGlobal3D(InitialLocalMatrix, rLeftHandSideMatrix);

      //std::cout<<"["<<this->Id()<<"] RB RotatedDynamic rLeftHandSideMatrix "<<rLeftHandSideMatrix<<std::endl;
    }

    // Transform Local to Global RHSVector:
    if ( rLocalSystem.CalculationFlags.Is(RigidBodyElement::COMPUTE_RHS_VECTOR) ){

      VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();

      //std::cout<<"["<<this->Id()<<"] RB Dynamic rRightHandSideVector "<<rRightHandSideVector<<std::endl;

      // works for 2D and 3D case
      BeamMathUtilsType::MapLocalToGlobal3D(InitialLocalMatrix, rRightHandSideVector);

      //std::cout<<"["<<this->Id()<<"] RB RotatedDynamic rRightHandSideVector "<<rRightHandSideVector<<std::endl;
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

Vector& RigidBodyElement::MapToInitialLocalFrame(Vector& rVariable)
{
    KRATOS_TRY

    BeamMathUtilsType::MapToCurrentLocalFrame(mInitialLocalQuaternion, rVariable);

    return rVariable;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyElement::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
						     ElementVariables& rVariables)
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size = dimension * (dimension + 1) * 0.5;

    rVariables.VolumeForce = CalculateVolumeForce(rVariables.VolumeForce);

    double DomainSize = rVariables.RigidBody.Mass;

    //gravity load
    Vector GravityLoad(dimension);
    noalias(GravityLoad) = ZeroVector(dimension);

    SizeType RowIndex = 0;
    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
      RowIndex = i * (dofs_size);

      for ( SizeType j = 0; j < dimension; j++ )
      {
        GravityLoad[j] = rVariables.VolumeForce[j] * DomainSize;
      }

      //substract because is added as a component of the InertiaRHS and is substracted again later in the scheme
      BeamMathUtilsType::SubstractVector( GravityLoad, rRightHandSideVector, RowIndex );

    }

    std::cout<<" Rigid Element Gravity "<<GravityLoad<<std::endl;
    std::cout<<" rRightHandSideVector "<<rRightHandSideVector<<std::endl;

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyElement::AddExplicitContribution(const VectorType& rRHSVector,
					       const Variable<VectorType>& rRHSVariable,
					       Variable<array_1d<double,3> >& rDestinationVariable,
					       const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension       = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size       = dimension * (dimension + 1) * 0.5;

    if( (rRHSVariable == RESIDUAL_VECTOR) ){

      if ( rDestinationVariable == FORCE_RESIDUAL )
      {

        for(SizeType i=0; i< number_of_nodes; i++)
        {
          int index = (dofs_size) * i;

          GetGeometry()[i].SetLock();

          array_1d<double, 3 > &ForceResidual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

          for(SizeType j=0; j<dimension; j++)
          {
            ForceResidual[j] += rRHSVector[index + j];
          }

          GetGeometry()[i].UnSetLock();
        }
      }
      else if( rDestinationVariable == MOMENT_RESIDUAL )
      {

        for(SizeType i=0; i< number_of_nodes; i++)
        {
          int index = dimension + (dofs_size) * i;

          GetGeometry()[i].SetLock();

          array_1d<double, 3 > &MomentResidual = GetGeometry()[i].FastGetSolutionStepValue(MOMENT_RESIDUAL);

          if( dimension == 2 ){
            MomentResidual[2] += rRHSVector[index];
          }
          else{
            for(SizeType j=0; j<dimension; j++)
            {
              MomentResidual[j] += rRHSVector[index + j];
            }
          }
          GetGeometry()[i].UnSetLock();
        }

      }


    }

    KRATOS_CATCH("")
}

//************************************CALCULATE VOLUME ACCELERATION*******************
//************************************************************************************

Vector&  RigidBodyElement::CalculateVolumeForce( Vector& rVolumeForce)
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension       = GetGeometry().WorkingSpaceDimension();

    if( rVolumeForce.size() != dimension )
      rVolumeForce.resize(dimension,false);

    noalias(rVolumeForce) = ZeroVector(dimension);

    for ( SizeType j = 0; j < number_of_nodes; j++ )
    {
        if( GetGeometry()[j].SolutionStepsDataHas(VOLUME_ACCELERATION) ) //temporary, will be checked once at the beginning only
            rVolumeForce += GetGeometry()[j].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }

    //rVolumeForce *= GetProperties()[DENSITY];

    //Current Frame is the local frame
    rVolumeForce = MapToInitialLocalFrame( rVolumeForce );

    return rVolumeForce;

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyElement::CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
                                                               VectorType& rRightHandSideVector,
                                                               ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY


    bool ComputeDynamicTangent = false;
    if( rCurrentProcessInfo.Has(COMPUTE_DYNAMIC_TANGENT) ){

      if(rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true){
	ComputeDynamicTangent = true;
      }
    }

    if( ComputeDynamicTangent == true ){

      //create local system components
      LocalSystemComponents LocalSystem;

      //calculation flags
      LocalSystem.CalculationFlags.Set(RigidBodyElement::COMPUTE_RHS_VECTOR);
      LocalSystem.CalculationFlags.Set(RigidBodyElement::COMPUTE_LHS_MATRIX);

      //Initialize sizes for the system components:
      this->InitializeSystemMatrices(rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags);

      //Set Variables to Local system components
      LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
      LocalSystem.SetRightHandSideVector(rRightHandSideVector);

      //Calculate elemental system
      CalculateDynamicSystem( LocalSystem, rCurrentProcessInfo );

    }
    else{

      //1.-Calculate Tangent Inertia Matrix:
      this->CalculateMassMatrix( rLeftHandSideMatrix, rCurrentProcessInfo );

      double MatSize = rLeftHandSideMatrix.size1();

      //2.-Calculate Inertial Forces:
      if ( rRightHandSideVector.size() != MatSize )
	rRightHandSideVector.resize( MatSize, false );

      noalias(rRightHandSideVector) = ZeroVector( MatSize ); //resetting RHS
    }


    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyElement::CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                                     ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    bool ComputeDynamicTangent = false;
    if( rCurrentProcessInfo.Has(COMPUTE_DYNAMIC_TANGENT) )
      if(rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true)
	ComputeDynamicTangent = true;

    if( ComputeDynamicTangent == true ){

      //create local system components
      LocalSystemComponents LocalSystem;

      //calculation flags
      LocalSystem.CalculationFlags.Set(RigidBodyElement::COMPUTE_LHS_MATRIX);

      VectorType RightHandSideVector = Vector();

      //Initialize sizes for the system components:
      this->InitializeSystemMatrices(rLeftHandSideMatrix, RightHandSideVector,  LocalSystem.CalculationFlags);

      //Set Variables to Local system components
      LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
      LocalSystem.SetRightHandSideVector(RightHandSideVector);

      //Calculate elemental system
      CalculateDynamicSystem( LocalSystem, rCurrentProcessInfo );

    }
    else{

      //1.-Calculate Tangent Inertia Matrix:
      this->CalculateMassMatrix( rLeftHandSideMatrix, rCurrentProcessInfo );

    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyElement::CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
                                                     ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    bool ComputeDynamicTangent = false;
    if( rCurrentProcessInfo.Has(COMPUTE_DYNAMIC_TANGENT) )
      if(rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true)
	ComputeDynamicTangent = true;

    if( ComputeDynamicTangent == true ){

      //create local system components
      LocalSystemComponents LocalSystem;

      //calculation flags
      LocalSystem.CalculationFlags.Set(RigidBodyElement::COMPUTE_RHS_VECTOR);

      MatrixType LeftHandSideMatrix = Matrix();

      //Initialize sizes for the system components:
      this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

      //Set Variables to Local system components
      LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
      LocalSystem.SetRightHandSideVector(rRightHandSideVector);

      //Calculate elemental system
      CalculateDynamicSystem( LocalSystem, rCurrentProcessInfo );

    }
    else{

      //create local system components
      LocalSystemComponents LocalSystem;

      //calculation flags
      LocalSystem.CalculationFlags.Set(RigidBodyElement::COMPUTE_RHS_VECTOR);

      MatrixType LeftHandSideMatrix = Matrix();

      //Initialize sizes for the system components:
      this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

      //RHS reseted to zero
    }


    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyElement::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                                          ElementVariables& rVariables)
{
  // add inertia RHS
  this->CalculateAndAddInertiaLHS(rLeftHandSideMatrix, rVariables);
}

//************************************************************************************
//************************************************************************************

void RigidBodyElement::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                          ElementVariables& rVariables)
{
  // calculate and add external forces
  this->CalculateAndAddExternalForces(rRightHandSideVector, rVariables);

  // add inertia RHS
  this->CalculateAndAddInertiaRHS(rRightHandSideVector, rVariables);
}

//************************************************************************************
//************************************************************************************

//Inertia in the SPATIAL configuration
void RigidBodyElement::CalculateAndAddInertiaLHS(MatrixType& rLeftHandSideMatrix,
                                                 ElementVariables& rVariables)
{
    KRATOS_TRY

    const ProcessInfo& rCurrentProcessInfo = rVariables.GetProcessInfo();

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension       = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size       = dimension * (dimension + 1) * 0.5;
    SizeType MatSize               = number_of_nodes * (dofs_size);

    if(rLeftHandSideMatrix.size1() != MatSize)
      rLeftHandSideMatrix.resize (MatSize, MatSize, false);

    rLeftHandSideMatrix = ZeroMatrix( MatSize, MatSize );

    //rCurrentProcessInfo must give it:
    double DeltaTime = rCurrentProcessInfo[DELTA_TIME];

    double Newmark1 = (1.0/ ( DeltaTime * DeltaTime * rCurrentProcessInfo[NEWMARK_BETA] ));
    double Newmark2 = ( DeltaTime * rCurrentProcessInfo[NEWMARK_GAMMA] );

    //block m(1,1) of the mass matrix

    MatrixType m11 = ZeroMatrix(3,3);

    double TotalMass = 0;
    TotalMass = rVariables.RigidBody.Mass;

    //block m(2,2) of the mass matrix

    MatrixType m22 = ZeroMatrix(3,3);

    Vector CurrentCompoundRotationVector(3);
    noalias(CurrentCompoundRotationVector) = ZeroVector(3);
    Vector PreviousCompoundRotationVector(3);
    noalias(PreviousCompoundRotationVector) = ZeroVector(3);

    Vector CurrentStepRotationVector(3);
    noalias(CurrentStepRotationVector) = ZeroVector(3);
    Vector AngularVelocityVector(3);
    noalias(AngularVelocityVector) = ZeroVector(3);
    Vector AngularAccelerationVector(3);
    noalias(AngularAccelerationVector) = ZeroVector(3);
    Vector CurrentAngularAccelerationVector(3);
    noalias(CurrentAngularAccelerationVector) = ZeroVector(3);
    Vector PreviousAngularAccelerationVector(3);
    noalias(PreviousAngularAccelerationVector) = ZeroVector(3);

    Vector CurrentValueVector(3);
    noalias(CurrentValueVector) = ZeroVector(3);
    Vector PreviousValueVector(3);
    noalias(PreviousValueVector) = ZeroVector(3);

    std::vector<QuaternionType> PreviousNodeQuaternions;
    QuaternionType QuaternionValue;

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {

	//Current Compound Rotation Vector
	CurrentValueVector = GetNodalCurrentValue( ROTATION, CurrentValueVector, i );
	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector );

    	CurrentCompoundRotationVector = CurrentValueVector;

    	//Previous Compound Rotation Vector
	PreviousValueVector = GetNodalPreviousValue( ROTATION, PreviousValueVector, i );
	PreviousValueVector = MapToInitialLocalFrame( PreviousValueVector );

	PreviousCompoundRotationVector = PreviousValueVector;

	//Current Step Rotation Vector
	CurrentValueVector = GetNodalCurrentValue( STEP_ROTATION, CurrentValueVector, i );
	CurrentStepRotationVector              = CurrentValueVector;

	//Angular Velocity Vector
	CurrentValueVector = GetNodalCurrentValue( ANGULAR_VELOCITY, CurrentValueVector, i );
	AngularVelocityVector                  = CurrentValueVector;

	//CurrentAngular Acceleration Vector
	CurrentValueVector = GetNodalCurrentValue( ANGULAR_ACCELERATION, CurrentValueVector, i );
	CurrentAngularAccelerationVector       = CurrentValueVector;

        //PreviousAngular Acceleration Vector
	CurrentValueVector = GetNodalPreviousValue( ANGULAR_ACCELERATION, CurrentValueVector, i );
	PreviousAngularAccelerationVector       = CurrentValueVector;

      }

    //Set step variables to local frame (current Frame is the local frame)
    CurrentStepRotationVector         = MapToInitialLocalFrame( CurrentStepRotationVector );
    AngularVelocityVector             = MapToInitialLocalFrame( AngularVelocityVector );
    CurrentAngularAccelerationVector  = MapToInitialLocalFrame( CurrentAngularAccelerationVector );
    PreviousAngularAccelerationVector = MapToInitialLocalFrame( PreviousAngularAccelerationVector );

    double AlphaM = rCurrentProcessInfo[BOSSAK_ALPHA];

    AngularAccelerationVector = (1.0-AlphaM)*CurrentAngularAccelerationVector + AlphaM*(PreviousAngularAccelerationVector);

    //Set step variables to local frame (current Frame is the local frame)
    Matrix CurrentRotationMatrix  = ZeroMatrix(3,3);
    Matrix PreviousRotationMatrix = ZeroMatrix(3,3);


    //1.-Get rotation matrix
    QuaternionType TotalQuaternion;

    //1.1.-Current Rotation Matrix
    TotalQuaternion = QuaternionType::FromRotationVector(CurrentCompoundRotationVector);

    TotalQuaternion.ToRotationMatrix( CurrentRotationMatrix );

    //1.2.-Previous Rotation Matrix

    //option 1:
    TotalQuaternion = QuaternionType::FromRotationVector(PreviousCompoundRotationVector);

    TotalQuaternion.ToRotationMatrix( PreviousRotationMatrix );

    //2.-Get inertia dyadic
    Matrix InertiaDyadic = ZeroMatrix(3,3);
    InertiaDyadic = rVariables.RigidBody.InertiaTensor;
    //Inertia dyadic expressed in the spatial frame
    InertiaDyadic = prod(CurrentRotationMatrix,InertiaDyadic);
    InertiaDyadic = prod(InertiaDyadic,trans(CurrentRotationMatrix));


    // SIMO -----------------------------
    //2.- Compute Term 1:

    Matrix MassTerm1 = ZeroMatrix(3,3);

    Vector InertiaxAngularVelocity     = prod( InertiaDyadic, AngularVelocityVector );
    Vector InertiaxAngularAcceleration = prod( InertiaDyadic, AngularAccelerationVector );

    // previous implementation
    // Vector VectorTerm1 = MathUtils<double>::CrossProduct( AngularVelocityVector, InertiaxAngularVelocity);
    // VectorTerm1 += InertiaxAngularAcceleration;

    Matrix TensorAngularVelocity = ZeroMatrix(3,3);
    BeamMathUtilsType::VectorToSkewSymmetricTensor( AngularVelocityVector, TensorAngularVelocity );

    Vector VectorTerm1(3);
    noalias(VectorTerm1) = ZeroVector(3);

    VectorTerm1  = prod( TensorAngularVelocity, InertiaxAngularVelocity );

    VectorTerm1 += InertiaxAngularAcceleration;

    BeamMathUtilsType::VectorToSkewSymmetricTensor(VectorTerm1, MassTerm1);

    //3.- Compute Term 2:

    Matrix MassTerm2 = ZeroMatrix(3,3);

    // Matrix TensorAngularVelocity = ZeroMatrix(3,3);
    // BeamMathUtilsType::VectorToSkewSymmetricTensor( AngularVelocityVector, TensorAngularVelocity );

    Matrix InertiaxAngularVelocityTensor = ZeroMatrix(3,3);
    BeamMathUtilsType::VectorToSkewSymmetricTensor( InertiaxAngularVelocity, InertiaxAngularVelocityTensor );

    Matrix TensorAngularVelocityxInertia = ZeroMatrix(3,3);
    noalias(TensorAngularVelocityxInertia) = prod( TensorAngularVelocity, InertiaDyadic );

    Matrix InertiaxTensorAngularVelocity = ZeroMatrix(3,3);
    noalias(InertiaxTensorAngularVelocity) = prod( InertiaDyadic, TensorAngularVelocity );

    //MassTerm2 = InertiaDyadic - Newmark2 * InertiaxAngularVelocityTensor + Newmark2 * TensorAngularVelocityxInertia + Newmark2 * InertiaxTensorAngularVelocity; //Cardona

    MassTerm2 = (1.0-AlphaM) * InertiaDyadic - Newmark2 * InertiaxAngularVelocityTensor + Newmark2 * TensorAngularVelocityxInertia; //Simo

    MassTerm2 *= Newmark1;

    MatrixType MassMatrixBlock2 = ZeroMatrix(3,3);

    MassMatrixBlock2 = (-1) * MassTerm1 + MassTerm2;

    // Compute Linear Part of the Step Rotation
    Matrix LinearPartRotationTensor = ZeroMatrix(3,3);

    double NormCurrentStepRotation =  norm_2(CurrentStepRotationVector);
    if( NormCurrentStepRotation != 0 ){

      this->CalculateRotationLinearPartTensor( CurrentStepRotationVector , LinearPartRotationTensor );

    }
    else{

      //std::cout<<" Attention... problem sure "<<std::endl;
      LinearPartRotationTensor = IdentityMatrix(3);
    }

    //this->CalculateRotationLinearPartTensor( CurrentStepRotationVector , LinearPartRotationTensor );
    // std::cout<<" MassMatrixBlock 2 SPA "<<MassMatrixBlock2<<std::endl;
    // std::cout<<" LinearPartRotationTensor SPA "<<LinearPartRotationTensor<<std::endl;

    MassMatrixBlock2 = prod( MassMatrixBlock2, LinearPartRotationTensor );

    SizeType RowIndex = 0;
    SizeType ColIndex = 0;

    Matrix DiagonalMatrix = IdentityMatrix(3);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
    	m11 = ZeroMatrix(3,3);
    	m22 = ZeroMatrix(3,3);

    	RowIndex = i * (dofs_size);

    	for ( SizeType j = 0; j < number_of_nodes; j++ )
    	  {
    	    ColIndex = j * (dofs_size);

    	    m11 = (1.0-AlphaM) * Newmark1 * TotalMass * DiagonalMatrix;

    	    m22 = MassMatrixBlock2;


    	    //Building the Local Tangent Inertia Matrix
    	    BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, m11, RowIndex, ColIndex );
	    if(dimension == 2)
	      rLeftHandSideMatrix(RowIndex+3,ColIndex+3) = m22(2,2);
	    else
	      BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, m22, RowIndex+3, ColIndex+3 );
    	  }
      }

    //std::cout<<" rLeftHandSideMatrix "<<rLeftHandSideMatrix<<std::endl;

    // SIMO -----------------------------


    // GERADIN --------------------------
    // //2.- Compute Damping term:

    // Matrix DampingTerm = ZeroMatrix(3,3);
    // Matrix TensorAngularVelocity = ZeroMatrix(3,3);
    // BeamMathUtilsType::VectorToSkewSymmetricTensor( AngularVelocityVector, TensorAngularVelocity );
    // TensorAngularVelocity *= (-1);
    // Matrix TensorAngularVelocityxInertia = ZeroMatrix(3,3);
    // noalias(TensorAngularVelocityxInertia) = prod( TensorAngularVelocity, InertiaDyadic );


    // Vector InertiaxAngularVelocity     = prod( InertiaDyadic, AngularVelocityVector );
    // Matrix TensorInertiaxAngularVelocity = ZeroMatrix(3,3);
    // BeamMathUtilsType::VectorToSkewSymmetricTensor( InertiaxAngularVelocity, TensorInertiaxAngularVelocity );
    // TensorInertiaxAngularVelocity *= (-1);

    // DampingTerm = TensorAngularVelocityxInertia - TensorInertiaxAngularVelocity;


    // //3.- Compute Stiffness term:

    // Matrix StiffnessTerm = ZeroMatrix(3,3);

    // Matrix TensorAngularAcceleration = ZeroMatrix(3,3);
    // BeamMathUtilsType::VectorToSkewSymmetricTensor( AngularAccelerationVector, TensorAngularAcceleration );
    // TensorAngularAcceleration *= (-1);

    // Matrix InertiaxTensorAngularAcceleration = prod( InertiaDyadic, TensorAngularAcceleration );

    // Matrix VelocityxInertiaxVelocity = prod( InertiaDyadic, TensorAngularVelocity );
    // VelocityxInertiaxVelocity = prod( TensorAngularVelocity, InertiaDyadic );

    // Matrix TensorInertiaxAngularVelocityxVelocity = prod( TensorInertiaxAngularVelocity, TensorAngularVelocity );

    // Matrix InertiaxVelocityxVelocity = prod( InertiaDyadic, TensorAngularVelocity );
    // InertiaxVelocityxVelocity = prod( InertiaxVelocityxVelocity, TensorAngularVelocity );

    // StiffnessTerm = 0.5 * ( VelocityxInertiaxVelocity + InertiaxTensorAngularAcceleration - TensorInertiaxAngularVelocityxVelocity ) - (1.0/3.0) * InertiaxVelocityxVelocity;


    // Matrix MassMatrixBlock2 = Newmark1 * (InertiaDyadic + Newmark2*DampingTerm) + StiffnessTerm;

    // MassMatrixBlock2 = prod( CurrentRotationMatrix, MassMatrixBlock2 );

    // MassMatrixBlock2 = prod( MassMatrixBlock2, trans(PreviousRotationMatrix));

    // // Compute Linear Part of the Step Rotation
    // Matrix LinearPartRotationTensor = ZeroMatrix(3,3);

    // double NormCurrentStepRotation =  norm_2(CurrentStepRotationVector);
    // if( NormCurrentStepRotation != 0 ){

    //   this->CalculateRotationLinearPartTensor( CurrentStepRotationVector , LinearPartRotationTensor );

    // }
    // else{

    //   LinearPartRotationTensor = IdentityMatrix(3);
    //   MassMatrixBlock2 = InertiaDyadic;

    // }

    // //std::cout<<" LinearPartRotationTensor "<<LinearPartRotationTensor<<std::endl;

    // MassMatrixBlock2 = prod( MassMatrixBlock2, LinearPartRotationTensor );



    // SizeType RowIndex = 0;
    // SizeType ColIndex = 0;

    // Matrix DiagonalMatrix = IdentityMatrix(3);

    // for ( SizeType i = 0; i < number_of_nodes; i++ )
    //   {
    // 	m11 = ZeroMatrix(3,3);
    // 	m22 = ZeroMatrix(3,3);

    // 	RowIndex = i * (dofs_size);

    // 	for ( SizeType j = 0; j < number_of_nodes; j++ )
    // 	  {

    // 	    ColIndex = j * (dofs_size);

    // 	    m11 = TotalMass * DiagonalMatrix;

    // 	    m22 = MassMatrixBlock2;


    // 	    //Building the Local Tangent Inertia Matrix
    // 	    BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, m11, RowIndex, ColIndex );
    // 	    BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, m22, RowIndex+3, ColIndex+3 );

    // 	  }
    //   }

    // GERADIN --------------------------

    KRATOS_CATCH("")

}


//************************************************************************************
//************************************************************************************

//Inertia in the SPATIAL configuration
void RigidBodyElement::CalculateAndAddInertiaRHS(VectorType& rRightHandSideVector,
                                                 ElementVariables& rVariables)
{
    KRATOS_TRY

    const ProcessInfo& rCurrentProcessInfo = rVariables.GetProcessInfo();

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension       = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size       = dimension * (dimension + 1) * 0.5;
    SizeType MatSize               = number_of_nodes * (dofs_size);

    if(rRightHandSideVector.size() != MatSize)
      rRightHandSideVector.resize(MatSize, false);

    noalias(rRightHandSideVector) = ZeroVector( MatSize );

    double TotalMass = 0;
    TotalMass = rVariables.RigidBody.Mass;


    //displacements and rotations vector
    Matrix TotalRotationMatrix(3,3);
    noalias(TotalRotationMatrix) = ZeroMatrix(3,3);
    Vector TotalRotationVector(3);
    noalias(TotalRotationVector) = ZeroVector(3);

    Vector LinearAccelerationVector(3);
    noalias(LinearAccelerationVector) = ZeroVector(3);
    Vector CurrentLinearAccelerationVector(3);
    noalias(CurrentLinearAccelerationVector) = ZeroVector(3);
    Vector PreviousLinearAccelerationVector(3);
    noalias(PreviousLinearAccelerationVector) = ZeroVector(3);

    Vector AngularVelocityVector(3);
    noalias(AngularVelocityVector) = ZeroVector(3);
    Vector AngularAccelerationVector(3);
    noalias(AngularAccelerationVector) = ZeroVector(3);
    Vector CurrentAngularAccelerationVector(3);
    noalias(CurrentAngularAccelerationVector) = ZeroVector(3);
    Vector PreviousAngularAccelerationVector(3);
    noalias(PreviousAngularAccelerationVector) = ZeroVector(3);


    Vector CurrentValueVector(3);
    noalias(CurrentValueVector) = ZeroVector(3);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	//Current Compound Rotation Vector
	CurrentValueVector = GetNodalCurrentValue( ROTATION, CurrentValueVector, i );
	TotalRotationVector             =  CurrentValueVector;

	//Current Linear Acceleration Vector
	CurrentValueVector = GetNodalCurrentValue( ACCELERATION, CurrentValueVector, i );
	CurrentLinearAccelerationVector  = CurrentValueVector;

	//Previous Linear Acceleration Vector
	CurrentValueVector = GetNodalPreviousValue( ACCELERATION, CurrentValueVector, i );
	PreviousLinearAccelerationVector = CurrentValueVector;

	//Angular Velocity Vector
	CurrentValueVector = GetNodalCurrentValue( ANGULAR_VELOCITY, CurrentValueVector, i );
	AngularVelocityVector           =  CurrentValueVector;

	//Current Angular Acceleration Vector
	CurrentValueVector = GetNodalCurrentValue( ANGULAR_ACCELERATION, CurrentValueVector, i );
	CurrentAngularAccelerationVector       =  CurrentValueVector;

        //Previous Angular Acceleration Vector
	CurrentValueVector = GetNodalPreviousValue( ANGULAR_ACCELERATION, CurrentValueVector, i );
	PreviousAngularAccelerationVector       =  CurrentValueVector;

      }

    //Set step variables to local frame (current Frame is the local frame)
    TotalRotationVector               = MapToInitialLocalFrame( TotalRotationVector );
    CurrentLinearAccelerationVector   = MapToInitialLocalFrame( CurrentLinearAccelerationVector );
    PreviousLinearAccelerationVector  = MapToInitialLocalFrame( PreviousLinearAccelerationVector );
    AngularVelocityVector             = MapToInitialLocalFrame( AngularVelocityVector );
    CurrentAngularAccelerationVector  = MapToInitialLocalFrame( CurrentAngularAccelerationVector );
    PreviousAngularAccelerationVector = MapToInitialLocalFrame( PreviousAngularAccelerationVector );

    double AlphaM = rCurrentProcessInfo[BOSSAK_ALPHA];
    LinearAccelerationVector  = (1.0-AlphaM) * CurrentLinearAccelerationVector + AlphaM * (PreviousLinearAccelerationVector);
    AngularAccelerationVector = (1.0-AlphaM) * CurrentAngularAccelerationVector + AlphaM * (PreviousAngularAccelerationVector);

    QuaternionType TotalQuaternion = QuaternionType::FromRotationVector(TotalRotationVector);
    Matrix CurrentRotationMatrix   = ZeroMatrix(3,3);
    TotalQuaternion.ToRotationMatrix( CurrentRotationMatrix );

    //-----------------
    //block m(1) of the inertial force vector

    //Compute Linear Term:

    Vector LinearInertialForceVector(3);
    noalias(LinearInertialForceVector) = ZeroVector(3);

    //this transformation is wrong
    //LinearInertialForceVector  = MapToMaterialFrame( TotalQuaternion, LinearInertialForceVector );

    LinearInertialForceVector = TotalMass * LinearAccelerationVector;

    //-----------------
    //block m(2,2) of the inertial force vector (rotations part::to be defined)

    //Get inertia dyadic
    Matrix InertiaDyadic = ZeroMatrix(3,3);
    InertiaDyadic = rVariables.RigidBody.InertiaTensor;
    //Inertia dyadic expressed in the initial frame
    InertiaDyadic = prod(CurrentRotationMatrix,InertiaDyadic);
    InertiaDyadic = prod(InertiaDyadic,trans(CurrentRotationMatrix));

    //Compute Angular Term:

    Vector InertiaxAngularVelocity     = prod( InertiaDyadic, AngularVelocityVector );
    Vector InertiaxAngularAcceleration = prod( InertiaDyadic, AngularAccelerationVector );

    Matrix TensorAngularVelocity = ZeroMatrix(3,3);
    BeamMathUtilsType::VectorToSkewSymmetricTensor( AngularVelocityVector, TensorAngularVelocity );

    Vector AngularInertialForceVector(3);
    noalias(AngularInertialForceVector) = ZeroVector(3);

    AngularInertialForceVector  = prod( TensorAngularVelocity, InertiaxAngularVelocity );

    // CROSS PRODUCT of AxB = prod( skewA, B )  where (skewA) = [Ax] = [A]^v =>  (skewA)^T = hat(A) (nomenclature)
    AngularInertialForceVector += InertiaxAngularAcceleration;

    //compose total acceleration integral function:

    Vector TotalInertialForceVector(3);
    noalias(TotalInertialForceVector) = ZeroVector(dofs_size);

    BeamMathUtilsType::AddVector(LinearInertialForceVector, TotalInertialForceVector, 0);

    if(dimension == 2)
      TotalInertialForceVector[dofs_size-1] = AngularInertialForceVector[2];
    else
      BeamMathUtilsType::AddVector(AngularInertialForceVector, TotalInertialForceVector, 3);

    //Initialize Local Matrices
    VectorType Fi(dofs_size);
    noalias(Fi) = ZeroVector(dofs_size);
    SizeType RowIndex = 0;

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {

    	RowIndex = i * (dofs_size);

    	noalias(Fi) = ZeroVector(dofs_size);

    	//nodal force vector
    	Fi  = TotalInertialForceVector;

	BeamMathUtilsType::AddVector(Fi, rRightHandSideVector, RowIndex);

    	//std::cout<<" Fi "<<Fi<<std::endl;

      }

    //std::cout<<" Rigid Body: rRightHandSideVector "<<rRightHandSideVector<<std::endl;

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyElement::CalculateMassMatrix(MatrixType& rMassMatrix,
                                           ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension       = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size       = dimension * (dimension + 1) * 0.5;
    SizeType MatSize               = number_of_nodes * (dofs_size);

    if(rMassMatrix.size1() != MatSize)
        rMassMatrix.resize (MatSize, MatSize, false);

    rMassMatrix = ZeroMatrix( MatSize, MatSize );

    // Rigid Body Properties
    RigidBodyProperties RigidBody;
    this->CalculateRigidBodyProperties(RigidBody);

    //block m(1,1) of the mass matrix

    MatrixType m11 = ZeroMatrix(3,3);

    double TotalMass = 0;
    TotalMass = RigidBody.Mass;


    //block m(2,2) of the mass matrix
    MatrixType m22 = ZeroMatrix(3,3);

    Matrix InertiaDyadic = ZeroMatrix(3,3);
    InertiaDyadic = RigidBody.InertiaTensor;


    for( SizeType i=0; i < number_of_nodes; i++ )
      {

        double temp = TotalMass;

    	int RowIndex = i * (dofs_size);

    	for( SizeType k=0; k < dimension; k++ )
    	  {
    	    m11(k,k) = temp;
    	  }

    	m22 = InertiaDyadic * temp;

    	//Building the Local Tangent Inertia Matrix
    	BeamMathUtilsType::AddMatrix( rMassMatrix, m11, RowIndex, RowIndex );

	if(dimension == 2)
	  rMassMatrix(RowIndex+3,RowIndex+3) = m22(2,2);
	else
	  BeamMathUtilsType::AddMatrix( rMassMatrix, m22, RowIndex+3, RowIndex+3 );

      }


    KRATOS_CATCH("")

}



//************************************************************************************
//************************************************************************************

void RigidBodyElement::CalculateRotationLinearPartTensor(Vector& rRotationVector,
                                                         Matrix& rRotationTensor)

{
    KRATOS_TRY

    if( rRotationTensor.size1() != 3 )
      rRotationTensor.resize(3, 3, false);

    rRotationTensor = ZeroMatrix(3,3);

    // adding term 3
    BeamMathUtilsType::VectorToSkewSymmetricTensor( rRotationVector, rRotationTensor );

    rRotationTensor *= 0.5;

    double NormRotation =  norm_2(rRotationVector);
    Matrix RotationxRotation = outer_prod( rRotationVector, rRotationVector );

    if( NormRotation != 0 )
      RotationxRotation *= ( 1.0/( NormRotation * NormRotation ) );

    // adding term 1
    rRotationTensor += RotationxRotation;

    double Coefficient = 0;
    if( NormRotation != 0 )
      Coefficient = 0.5 * NormRotation / std::tan( 0.5 * NormRotation );

    for(SizeType i=0; i<3; i++)
      RotationxRotation(i,i) -= 1.0;

    RotationxRotation *= (-1)*Coefficient;

    // adding term 2
    rRotationTensor += RotationxRotation;


    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyElement::UpdateRigidBodyNodes(ProcessInfo& rCurrentProcessInfo)
{

     KRATOS_TRY

     Matrix SkewSymVariable(3,3);
     noalias(SkewSymVariable) = ZeroMatrix(3,3);
     Vector RadiusVector(3);
     noalias(RadiusVector) = ZeroVector(3);
     Vector Variable(3);
     noalias(Variable) = ZeroVector(3);
     Vector AngularVariable(3);
     noalias(AngularVariable) = ZeroVector(3);
     array_1d<double,3>     VariableArray;

     array_1d<double,3> Radius;

     Node<3>::Pointer rCenterOfGravity = this->GetGeometry()(0);

     if( rCenterOfGravity->Is(SLAVE) ){
       Element& MasterElement = this->GetGeometry()[0].GetValue(MASTER_ELEMENTS).back();
       rCenterOfGravity = MasterElement.GetGeometry()(0);
     }

     array_1d<double, 3 >&  Center              = rCenterOfGravity->GetInitialPosition();

     array_1d<double, 3 >&  Displacement        = rCenterOfGravity->FastGetSolutionStepValue(DISPLACEMENT);
     array_1d<double, 3 >&  Rotation            = rCenterOfGravity->FastGetSolutionStepValue(ROTATION);
     array_1d<double, 3 >&  StepRotation        = rCenterOfGravity->FastGetSolutionStepValue(STEP_ROTATION);
     array_1d<double, 3 >&  DeltaRotation       = rCenterOfGravity->FastGetSolutionStepValue(DELTA_ROTATION);

     array_1d<double, 3 >&  Velocity            = rCenterOfGravity->FastGetSolutionStepValue(VELOCITY);
     array_1d<double, 3 >&  Acceleration        = rCenterOfGravity->FastGetSolutionStepValue(ACCELERATION);
     array_1d<double, 3 >&  AngularVelocity     = rCenterOfGravity->FastGetSolutionStepValue(ANGULAR_VELOCITY);
     array_1d<double, 3 >&  AngularAcceleration = rCenterOfGravity->FastGetSolutionStepValue(ANGULAR_ACCELERATION);

     //std::cout<<" [  MasterElement "<<this->Id() ];
     //std::cout<<" [ Rotation:"<<Rotation<<",StepRotation:"<<StepRotation<<",DeltaRotation:"<<DeltaRotation<<"]"<<std::endl;
     //std::cout<<" [ Velocity:"<<Velocity<<",Acceleration:"<<Acceleration<<",Displacement:"<<Displacement<<",DeltaDisplacement"<<Displacement-rCenterOfGravity->FastGetSolutionStepValue(DISPLACEMENT,1)<<"]"<<std::endl;

     for (NodesContainerType::iterator i = mpNodes->begin(); i != mpNodes->end(); ++i)
     {
       //Get rotation matrix
       QuaternionType TotalQuaternion = QuaternionType::FromRotationVector<array_1d<double,3> >(Rotation);

       Radius = (i)->GetInitialPosition() - Center;

       Matrix RotationMatrix;
       TotalQuaternion.ToRotationMatrix(RotationMatrix);

       for(int j=0; j<3; j++)
         RadiusVector[j] = Radius[j];

       RadiusVector = prod( RotationMatrix, RadiusVector );

       for(int j=0; j<3; j++)
         Radius[j] = RadiusVector[j];

       //TotalQuaternion.RotateVector3<array_1d<double,3> >(Radius);

       array_1d<double, 3 >&  NodeDisplacement  = (i)->FastGetSolutionStepValue(DISPLACEMENT);
       array_1d<double, 3 >&  NodeRotation      = (i)->FastGetSolutionStepValue(ROTATION);
       array_1d<double, 3 >&  NodeStepRotation  = (i)->FastGetSolutionStepValue(STEP_ROTATION);
       array_1d<double, 3 >&  NodeDeltaRotation = (i)->FastGetSolutionStepValue(DELTA_ROTATION);

       noalias(NodeDisplacement)  = ( (Center + Displacement)  + Radius ) - (i)->GetInitialPosition();
       noalias(NodeRotation)      = Rotation;
       noalias(NodeStepRotation)  = StepRotation;
       noalias(NodeDeltaRotation) = DeltaRotation;


       for(int j=0; j<3; j++)
         RadiusVector[j] = Radius[j];

       //********************
       for(int j=0; j<3; j++)
         Variable[j] = AngularVelocity[j];

       //compute the skewsymmmetric tensor of the angular velocity
       BeamMathUtilsType::VectorToSkewSymmetricTensor(Variable, SkewSymVariable);

       //compute the contribution of the angular velocity to the velocity v = Wxr
       Variable = prod(SkewSymVariable,RadiusVector);

       for(int j=0; j<3; j++)
         VariableArray[j] = Variable[j];

       (i)->FastGetSolutionStepValue(VELOCITY)               = Velocity + VariableArray;


       //********************

       //centripetal acceleration:
       for(int j=0; j<3; j++)
         AngularVariable[j] = AngularVelocity[j];

       //compute the skewsymmmetric tensor of the angular velocity
       BeamMathUtilsType::VectorToSkewSymmetricTensor(AngularVariable, SkewSymVariable);

       AngularVariable = prod(SkewSymVariable,Variable); //ac = Wx(Wxr)


       for(int j=0; j<3; j++)
         Variable[j] = AngularAcceleration[j];

       //compute the skewsymmmetric tensor of the angular acceleration
       BeamMathUtilsType::VectorToSkewSymmetricTensor(Variable, SkewSymVariable);

       //compute the contribution of the angular velocity to the velocity a = Axr
       Variable = prod(SkewSymVariable,RadiusVector);

       for(int j=0; j<3; j++)
         VariableArray[j] = Variable[j] + AngularVariable[j];

       (i)->FastGetSolutionStepValue(ACCELERATION)           = Acceleration + VariableArray;


       //********************
       (i)->FastGetSolutionStepValue(ANGULAR_VELOCITY)       = AngularVelocity;
       (i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION)   = AngularAcceleration;

       // 	std::cout<<"  [ Finalize Rigid Body Link Point : [Id:"<<(i)->Id()<<"] "<<std::endl;
       // 	std::cout<<"  [ Displacement:"<<NodeDisplacement<<" / StepRotation"<<NodeStepRotation<<" ] "<<std::endl;
       // 	std::cout<<"  [ Rotation:"<<NodeRotation<<" / Angular Acceleration"<<AngularAcceleration<<" ] "<<std::endl;

     }

     KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

RigidBodyElement::SizeType RigidBodyElement::GetDofsSize()
{
  KRATOS_TRY

  const SizeType dimension = GetGeometry().WorkingSpaceDimension();
  const SizeType number_of_nodes  = GetGeometry().PointsNumber();

  return (number_of_nodes*(dimension*(dimension + 1)*0.5));

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int  RigidBodyElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (GetGeometry().WorkingSpaceDimension() != 3 || GetGeometry().size()!=1 )
    {
      KRATOS_THROW_ERROR( std::invalid_argument, "This element works only in 3D and with 1 noded geometry", "")
    }

    //verify that the variables are correctly initialized
    if(VELOCITY.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"VELOCITY has Key zero! (check if the application is correctly registered", "" )
    if(DISPLACEMENT.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered", "" )
    if(ACCELERATION.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"ACCELERATION has Key zero! (check if the application is correctly registered", "" )
    if(DENSITY.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"DENSITY has Key zero! (check if the application is correctly registered", "" )
     if(VOLUME_ACCELERATION.Key() == 0)
        // KRATOS_THROW_ERROR( std::invalid_argument,"VOLUME_ACCELERATION has Key zero! (check if the application is correctly registered", "" )
    if(NODAL_MASS.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"NODAL_MASS has Key zero! (check if the application is correctly registered", "" )
    if(LOCAL_INERTIA_TENSOR.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"LOCAL_INERTIA_TENSOR has Key zero! (check if the application is correctly registered", "" )
    if(ROTATION.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"ROTATION has Key zero! (check if the application is correctly registered", "" )

    //verify that the dofs exist
    for(SizeType i=0; i<this->GetGeometry().size(); i++)
    {
        if(this->GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false)
            KRATOS_THROW_ERROR( std::invalid_argument,"missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id() )
        if(this->GetGeometry()[i].HasDofFor(DISPLACEMENT_X) == false || this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false || this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Z) == false)
            KRATOS_THROW_ERROR( std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id() )
    }

    //verify that the area is given by properties
    if (this->GetProperties().Has(NODAL_MASS)==false)
    {
        if( GetValue(NODAL_MASS) == 0.0 )
            KRATOS_THROW_ERROR( std::logic_error,"NODAL_MASS not provided for this element", this->Id() )
    }

    //verify that the inertia is given by properties
    if (this->GetProperties().Has(LOCAL_INERTIA_TENSOR)==false)
    {
        if( GetValue(LOCAL_INERTIA_TENSOR)(0,0) == 0.0 )
	  KRATOS_THROW_ERROR( std::logic_error,"LOCAL_INERTIA_TENSOR not provided for this element ", this->Id() )
    }

    return 0;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    rSerializer.save("InitialLocalQuaternion",mInitialLocalQuaternion);
    rSerializer.save("RigidBodyNodes",mpNodes);
}

void RigidBodyElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    rSerializer.load("InitialLocalQuaternion",mInitialLocalQuaternion);
    rSerializer.load("RigidBodyNodes",mpNodes);
}


} // Namespace Kratos
