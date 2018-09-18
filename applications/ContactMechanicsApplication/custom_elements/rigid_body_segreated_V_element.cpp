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

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

RigidBodySegregatedVElement::RigidBodySegregatedVElement(IndexType NewId,GeometryType::Pointer pGeometry)
    :RigidBodyElement(NewId, pGeometry)
{
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

RigidBodySegregatedVElement::RigidBodySegregatedVElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    :RigidBodyElement(NewId, pGeometry, pProperties)
{
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

RigidBodySegregatedVElement::RigidBodySegregatedVElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, NodesContainerType::Pointer pNodes)
    :RigidBodyElement(NewId, pGeometry, pProperties)
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

RigidBodySegregatedVElement::RigidBodySegregatedVElement(RigidBodySegregatedVElement const& rOther)
    :RigidBodyElement(rOther)
{
}

//*********************************CREATE*********************************************
//************************************************************************************

Element::Pointer RigidBodySegregatedVElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
  return Kratos::make_shared<RigidBodySegregatedVElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}


//*********************************CLONE**********************************************
//************************************************************************************

Element::Pointer RigidBodySegregatedVElement::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
{
  RigidBodySegregatedVElement NewElement( NewId, GetGeometry().Create(ThisNodes), pGetProperties(), mpNodes );

  NewElement.mInitialLocalQuaternion = this->mInitialLocalQuaternion;
  NewElement.SetData(this->GetData());
  NewElement.SetFlags(this->GetFlags());

  return Kratos::make_shared<RigidBodySegregatedVElement>(NewElement);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

RigidBodySegregatedVElement::~RigidBodySegregatedVElement()
{
}

//************************************************************************************
//************************************************************************************

void RigidBodySegregatedVElement::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
{
  rElementalDofList.resize(0);

  switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
  {
    case VELOCITY_STEP:
      {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();
        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
          ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
          ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
          if( dimension == 2 ){
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));
          }
          else{
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Z));
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));
          }
        }
        break;
      }
    case PRESSURE_STEP:
      {
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }

}

//************************************************************************************
//************************************************************************************

void RigidBodySegregatedVElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
  SizeType dofs_size = this->GetDofsSize();

  if ( rResult.size() != dofs_size )
    rResult.resize(dofs_size, false);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();
        dofs_size = dimension * (dimension + 1) * 0.5;
        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
          SizeType index = i * (dofs_size);
          rResult[index]   = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
          rResult[index+1] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
          if( dimension == 2 ){
            rResult[index+2] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
          }
          else{
            rResult[index+2] = GetGeometry()[i].GetDof(VELOCITY_Z).EquationId();
            rResult[index+3] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
            rResult[index+4] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
            rResult[index+5] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
          }
        }
        break;
      }
    case PRESSURE_STEP:
      {
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }
}


//*********************************DISPLACEMENT***************************************
//************************************************************************************

void RigidBodySegregatedVElement::GetValuesVector(Vector& rValues, int Step)
{
  KRATOS_TRY

  SizeType dofs_size = this->GetDofsSize();

  if ( rValues.size() != dofs_size )
    rValues.resize( dofs_size, false );

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();
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
        break;
      }
    case PRESSURE_STEP:
      {
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }
  KRATOS_CATCH("")
}

//************************************VELOCITY****************************************
//************************************************************************************

void RigidBodySegregatedVElement::GetFirstDerivativesVector(Vector& rValues, int Step)
{
  KRATOS_TRY

  SizeType dofs_size = this->GetDofsSize();

  if ( rValues.size() != dofs_size )
    rValues.resize( dofs_size, false );

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        const SizeType number_of_nodes  = GetGeometry().size();
        const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
        dofs_size = dimension * (dimension + 1) * 0.5;
        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
          SizeType index = i * dimension * (dimension + 1) * 0.5;
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
        break;
      }
    case PRESSURE_STEP:
      {
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
  }
  KRATOS_CATCH("")
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void RigidBodySegregatedVElement::GetSecondDerivativesVector(Vector& rValues, int Step)
{
  KRATOS_TRY

  SizeType dofs_size = this->GetDofsSize();

  if ( rValues.size() != dofs_size )
    rValues.resize( dofs_size, false );

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();
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
        break;
      }
    case PRESSURE_STEP:
      {
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
  }
  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodySegregatedVElement::Initialize()
{
  KRATOS_TRY

  RigidBodyElement::Initialize();

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodySegregatedVElement::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyElement::InitializeSolutionStep(rCurrentProcessInfo);
        break;
      }
    case PRESSURE_STEP:
      {
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodySegregatedVElement::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  this->SetProcessInformation(rCurrentProcessInfo);

  RigidBodyElement::InitializeNonLinearIteration(rCurrentProcessInfo);

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodySegregatedVElement::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyElement::FinalizeNonLinearIteration(rCurrentProcessInfo);
        break;
      }
    case PRESSURE_STEP:
      {
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }


  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodySegregatedVElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyElement::FinalizeSolutionStep(rCurrentProcessInfo);
        break;
      }
    case PRESSURE_STEP:
      {
        //set as VELOCITY STEP for gauss point calculations:
        mStepVariable = VELOCITY_STEP;
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void RigidBodySegregatedVElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                         ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  SolidElement::CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodySegregatedVElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                        ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  RigidBodyElement::CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodySegregatedVElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                       VectorType& rRightHandSideVector,
                                                       ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  RigidBodyElement::CalculateLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodySegregatedVElement::CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
                                                                          VectorType& rRightHandSideVector,
                                                                          ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  RigidBodyElement::CalculateSecondDerivativesContributions(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyElement::CalculateElementalSystem(LocalSystemComponents& rLocalSystem,
                                                ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    //create and initialize element variables:
    ElementVariables Variables;
    Variables.Initialize(dimension,rCurrentProcessInfo);

    // std::cout<<" RigidBodyElement "<<this->Id()<<std::endl;
    // std::cout<<" [Displacement "<<GetGeometry()[0].FastGetSolutionStepValue( DISPLACEMENT )<<"]"<<std::endl;
    // std::cout<<" [Rotation     "<<GetGeometry()[0].FastGetSolutionStepValue( ROTATION )<<"]"<<std::endl;

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

void RigidBodySegregatedVElement::CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                                                ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  RigidBodyElement::CalculateSecondDerivativesLHS(rLeftHandSideMatrix, rCurrentProcessInfo);

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodySegregatedVElement::CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
                                                                ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  RigidBodyElement::CalculateSecondDerivativesRHS(rRightHandSideVector, rCurrentProcessInfo);

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

//Inertia in the SPATIAL configuration
void RigidBodySegregatedVElement::CalculateAndAddInertiaLHS(MatrixType& rLeftHandSideMatrix,
                                                            ElementVariables& rVariables)
{
  KRATOS_TRY

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
void RigidBodySegregatedVElement::CalculateAndAddInertiaRHS(VectorType& rRightHandSideVector,
                                                            ElementVariables& rVariables)
{
  KRATOS_TRY

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

RigidBodyElement::Sizetype RigidBodySegregatedVElement::GetDofsSize()
{
  KRATOS_TRY

  SizeType size = 0;
  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      const SizeType dimension = GetGeometry().WorkingSpaceDimension();
      const SizeType number_of_nodes  = GetGeometry().PointsNumber();
      size = number_of_nodes*dimension*(dimension + 1)*0.5; //size for velocity
      break;
    case PRESSURE_STEP:
      size = 0;
      break;
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;

  }
  return size;

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

int RigidBodySegregatedVElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  // Perform base element checks
  int ErrorCode = 0;
  ErrorCode = RigidBodyElement::Check(rCurrentProcessInfo);

  return ErrorCode;

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodySegregatedVElement::save( Serializer& rSerializer ) const
{
  KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, RigidBodyElement )
}

void RigidBodySegregatedVElement::load( Serializer& rSerializer )
{
  KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, RigidBodyElement )
}


} // Namespace Kratos
