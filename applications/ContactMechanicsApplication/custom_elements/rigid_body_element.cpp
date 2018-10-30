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
    const SizeType dimension       = GetGeometry().WorkingSpaceDimension();

    ElementalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_X));
    ElementalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_Y));
    if( dimension == 2 ){
      ElementalDofList.push_back(GetGeometry()[0].pGetDof(ROTATION_Z));
    }
    else{
      ElementalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_Z));
      ElementalDofList.push_back(GetGeometry()[0].pGetDof(ROTATION_X));
      ElementalDofList.push_back(GetGeometry()[0].pGetDof(ROTATION_Y));
      ElementalDofList.push_back(GetGeometry()[0].pGetDof(ROTATION_Z));
    }
}

//************************************************************************************
//************************************************************************************

void RigidBodyElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size = this->GetDofsSize();

    if ( rResult.size() != dofs_size )
        rResult.resize(dofs_size, false);

    rResult[0] = GetGeometry()[0].GetDof(DISPLACEMENT_X).EquationId();
    rResult[1] = GetGeometry()[0].GetDof(DISPLACEMENT_Y).EquationId();
    if( dimension == 2 ){
      rResult[2] = GetGeometry()[0].GetDof(ROTATION_Z).EquationId();
    }
    else{
      rResult[2] = GetGeometry()[0].GetDof(DISPLACEMENT_Z).EquationId();
      rResult[3] = GetGeometry()[0].GetDof(ROTATION_X).EquationId();
      rResult[4] = GetGeometry()[0].GetDof(ROTATION_Y).EquationId();
      rResult[5] = GetGeometry()[0].GetDof(ROTATION_Z).EquationId();
    }
}


//*********************************DISPLACEMENT***************************************
//************************************************************************************

void RigidBodyElement::GetValuesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size = this->GetDofsSize();

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    rValues[0] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_X, Step );
    rValues[1] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Y, Step );
    if( dimension == 2 ){
      rValues[2] = GetGeometry()[0].GetSolutionStepValue( ROTATION_Z, Step );
    }
    else{
      rValues[2] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Z, Step );
      rValues[3] = GetGeometry()[0].GetSolutionStepValue( ROTATION_X, Step );
      rValues[4] = GetGeometry()[0].GetSolutionStepValue( ROTATION_Y, Step );
      rValues[5] = GetGeometry()[0].GetSolutionStepValue( ROTATION_Z, Step );
    }

    KRATOS_CATCH("")
}

//************************************VELOCITY****************************************
//************************************************************************************

void RigidBodyElement::GetFirstDerivativesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size = this->GetDofsSize();

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    rValues[0] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_X, Step );
    rValues[1] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_Y, Step );
    if( dimension == 2 ){
      rValues[2] = GetGeometry()[0].GetSolutionStepValue( ANGULAR_VELOCITY_Z, Step );
    }
    else{
      rValues[2] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_Z, Step );
      rValues[3] = GetGeometry()[0].GetSolutionStepValue( ANGULAR_VELOCITY_X, Step );
      rValues[4] = GetGeometry()[0].GetSolutionStepValue( ANGULAR_VELOCITY_Y, Step );
      rValues[5] = GetGeometry()[0].GetSolutionStepValue( ANGULAR_VELOCITY_Z, Step );
    }

    KRATOS_CATCH("")
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void RigidBodyElement::GetSecondDerivativesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size = this->GetDofsSize();

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    rValues[0] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_X, Step );
    rValues[1] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_Y, Step );
    if( dimension == 2 ){
      rValues[2] = GetGeometry()[0].GetSolutionStepValue( ANGULAR_ACCELERATION_Z, Step );
    }
    else{
      rValues[2] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_Z, Step );
      rValues[3] = GetGeometry()[0].GetSolutionStepValue( ANGULAR_ACCELERATION_X, Step );
      rValues[4] = GetGeometry()[0].GetSolutionStepValue( ANGULAR_ACCELERATION_Y, Step );
      rValues[5] = GetGeometry()[0].GetSolutionStepValue( ANGULAR_ACCELERATION_Z, Step );
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
    const SizeType dofs_size = this->GetDofsSize();

    if ( rCalculationFlags.Is(RigidBodyElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != dofs_size )
            rLeftHandSideMatrix.resize( dofs_size, dofs_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( dofs_size, dofs_size ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( rCalculationFlags.Is(RigidBodyElement::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != dofs_size )
	    rRightHandSideVector.resize( dofs_size, false );

	noalias( rRightHandSideVector ) = ZeroVector( dofs_size ); //resetting RHS
    }
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


    this->MapLocalToGlobalSystem(rLocalSystem);

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyElement::MapLocalToGlobalSystem(LocalSystemComponents& rLocalSystem)
{
    KRATOS_TRY

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

RigidBodyElement::ArrayType& RigidBodyElement::MapToInitialLocalFrame(ArrayType& rVariable)
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

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    rVariables.VolumeForce = CalculateVolumeForce(rVariables.VolumeForce);

    double DomainSize = rVariables.RigidBody.Mass;

    //gravity load
    Vector GravityLoad(dimension);
    noalias(GravityLoad) = ZeroVector(dimension);

    for ( SizeType j = 0; j < dimension; j++ )
    {
      GravityLoad[j] = rVariables.VolumeForce[j] * DomainSize;
    }

    //substract because is added as a component of the InertiaRHS and is substracted again later in the scheme
    BeamMathUtilsType::SubstractVector( GravityLoad, rRightHandSideVector, 0 );

    // std::cout<<" Rigid Element Gravity x Mass"<<GravityLoad<<std::endl;
    // std::cout<<" Rigid Body External Force : "<<rRightHandSideVector<<std::endl;

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

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    if( (rRHSVariable == RESIDUAL_VECTOR) ){

      if ( rDestinationVariable == FORCE_RESIDUAL )
      {

        GetGeometry()[0].SetLock();

        array_1d<double, 3 > &ForceResidual = GetGeometry()[0].FastGetSolutionStepValue(FORCE_RESIDUAL);

        for(SizeType j=0; j<dimension; j++)
        {
          ForceResidual[j] += rRHSVector[j];
        }

        GetGeometry()[0].UnSetLock();
      }
      else if( rDestinationVariable == MOMENT_RESIDUAL )
      {
        int index = dimension;

        GetGeometry()[0].SetLock();

        array_1d<double, 3 > &MomentResidual = GetGeometry()[0].FastGetSolutionStepValue(MOMENT_RESIDUAL);

        if( dimension == 2 ){
          MomentResidual[2] += rRHSVector[index];
        }
        else{
          for(SizeType j=0; j<dimension; j++)
          {
            MomentResidual[j] += rRHSVector[index + j];
          }
        }
        GetGeometry()[0].UnSetLock();
      }
    }

    KRATOS_CATCH("")
}

//************************************CALCULATE VOLUME ACCELERATION*******************
//************************************************************************************

RigidBodyElement::ArrayType& RigidBodyElement::CalculateVolumeForce(ArrayType& rVolumeForce)
{
    KRATOS_TRY

    const SizeType dimension       = GetGeometry().WorkingSpaceDimension();

    noalias(rVolumeForce) = ZeroVector(3);

    if( GetGeometry()[0].SolutionStepsDataHas(VOLUME_ACCELERATION) ) //temporary, will be checked once at the beginning only
      for(SizeType i=0; i<dimension; ++i)
        rVolumeForce[i] += GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION)[i];

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

      double dofs_size = rLeftHandSideMatrix.size1();

      //2.-Calculate Inertial Forces:
      if ( rRightHandSideVector.size() != dofs_size )
	rRightHandSideVector.resize( dofs_size, false );

      noalias(rRightHandSideVector) = ZeroVector( dofs_size ); //resetting RHS
    }

    // std::cout<<" RIGID BODY RHS "<<rRightHandSideVector<<std::endl;
    // std::cout<<" RIGID BODY LHS "<<rLeftHandSideMatrix<<std::endl;

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

  // std::cout<<" Rigid body LHS "<<rLeftHandSideMatrix<<std::endl;
}

//************************************************************************************
//************************************************************************************

void RigidBodyElement::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                          ElementVariables& rVariables)
{
  // add inertia RHS
  this->CalculateAndAddInertiaRHS(rRightHandSideVector, rVariables);

  // calculate and add external forces
  this->CalculateAndAddExternalForces(rRightHandSideVector, rVariables);

  // std::cout<<" Rigid body RHS "<<rRightHandSideVector<<std::endl;

}

//************************************************************************************
//************************************************************************************

//Inertia in the SPATIAL configuration
void RigidBodyElement::CalculateAndAddInertiaLHS(MatrixType& rLeftHandSideMatrix,
                                                 ElementVariables& rVariables)
{
    KRATOS_TRY

    const ProcessInfo& rCurrentProcessInfo = rVariables.GetProcessInfo();
    const SizeType dimension       = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size       = this->GetDofsSize();

    if(rLeftHandSideMatrix.size1() != dofs_size)
      rLeftHandSideMatrix.resize (dofs_size, dofs_size, false);

    rLeftHandSideMatrix = ZeroMatrix( dofs_size, dofs_size );

    //rCurrentProcessInfo must give it:
    double Newmark0 = 0;
    double Newmark1 = 0;
    double Newmark2 = 0;

    this->GetTimeIntegrationParameters(Newmark0,Newmark1,Newmark2,rCurrentProcessInfo);

    //Current Compound Rotation Vector
    ArrayType CurrentCompoundRotationVector = GetGeometry()[0].FastGetSolutionStepValue(ROTATION);
    CurrentCompoundRotationVector = MapToInitialLocalFrame(CurrentCompoundRotationVector);

    //Previous Compound Rotation Vector
    // ArrayType PreviousCompoundRotationVector = GetGeometry()[0].FastGetSolutionStepValue(ROTATION,1);
    // PreviousCompoundRotationVector = MapToInitialLocalFrame(PreviousCompoundRotationVector);

    //Current Step Rotation Vector
    ArrayType CurrentStepRotationVector = GetGeometry()[0].FastGetSolutionStepValue(STEP_ROTATION);
    CurrentStepRotationVector = MapToInitialLocalFrame(CurrentStepRotationVector);

    //Angular Velocity Vector
    ArrayType AngularVelocityVector = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
    AngularVelocityVector = MapToInitialLocalFrame(AngularVelocityVector);

    //CurrentAngular Acceleration Vector
    ArrayType CurrentAngularAccelerationVector = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_ACCELERATION);
    CurrentAngularAccelerationVector = MapToInitialLocalFrame(CurrentAngularAccelerationVector);

    //PreviousAngular Acceleration Vector
    ArrayType PreviousAngularAccelerationVector = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_ACCELERATION,1);
    PreviousAngularAccelerationVector = MapToInitialLocalFrame(PreviousAngularAccelerationVector);


    double AlphaM = rCurrentProcessInfo[BOSSAK_ALPHA];
    ArrayType AngularAccelerationVector = (1.0-AlphaM)*CurrentAngularAccelerationVector + AlphaM*(PreviousAngularAccelerationVector);

    //Set step variables to local frame (current Frame is the local frame)

    //1.-Get rotation matrix
    QuaternionType TotalQuaternion;

    //1.1.-Current Rotation Matrix
    Matrix CurrentRotationMatrix(3,3);
    TotalQuaternion = QuaternionType::FromRotationVector(CurrentCompoundRotationVector);
    TotalQuaternion.ToRotationMatrix( CurrentRotationMatrix );

    //1.2.-Previous Rotation Matrix

    //option 1:
    // Matrix PreviousRotationMatrix(3,3);
    // TotalQuaternion = QuaternionType::FromRotationVector(PreviousCompoundRotationVector);
    // TotalQuaternion.ToRotationMatrix( PreviousRotationMatrix );

    //2.-Get inertia dyadic
    Matrix InertiaDyadic = ZeroMatrix(3,3);
    InertiaDyadic = rVariables.RigidBody.InertiaTensor;
    //Inertia dyadic expressed in the spatial frame
    InertiaDyadic = prod(CurrentRotationMatrix,InertiaDyadic);
    InertiaDyadic = prod(InertiaDyadic,trans(CurrentRotationMatrix));


    // SIMO -----------------------------
    //2.- Compute Term 1:

    Matrix MassTerm1(3,3);

    Vector InertiaxAngularVelocity     = prod( InertiaDyadic, AngularVelocityVector );
    Vector InertiaxAngularAcceleration = prod( InertiaDyadic, AngularAccelerationVector );

    // previous implementation
    // Vector VectorTerm1 = MathUtils<double>::CrossProduct( AngularVelocityVector, InertiaxAngularVelocity);
    // VectorTerm1 += InertiaxAngularAcceleration;

    Matrix TensorAngularVelocity(3,3);
    BeamMathUtilsType::VectorToSkewSymmetricTensor( AngularVelocityVector, TensorAngularVelocity );

    Vector VectorTerm1(3);
    noalias(VectorTerm1) = prod( TensorAngularVelocity, InertiaxAngularVelocity );
    VectorTerm1 += InertiaxAngularAcceleration;

    BeamMathUtilsType::VectorToSkewSymmetricTensor(VectorTerm1, MassTerm1);

    //3.- Compute Term 2:

     // Matrix TensorAngularVelocity(3,3);
    // BeamMathUtilsType::VectorToSkewSymmetricTensor( AngularVelocityVector, TensorAngularVelocity );

    Matrix InertiaxAngularVelocityTensor(3,3);
    BeamMathUtilsType::VectorToSkewSymmetricTensor( InertiaxAngularVelocity, InertiaxAngularVelocityTensor );

    Matrix TensorAngularVelocityxInertia(3,3);
    noalias(TensorAngularVelocityxInertia) = prod( TensorAngularVelocity, InertiaDyadic );

    Matrix InertiaxTensorAngularVelocity(3,3);
    noalias(InertiaxTensorAngularVelocity) = prod( InertiaDyadic, TensorAngularVelocity );

    //MassTerm2 = InertiaDyadic - Newmark2 * InertiaxAngularVelocityTensor + Newmark2 * TensorAngularVelocityxInertia + Newmark2 * InertiaxTensorAngularVelocity; //Cardona

    Matrix MassTerm2(3,3);
    noalias(MassTerm2) = (1.0-AlphaM) * InertiaDyadic - Newmark2 * InertiaxAngularVelocityTensor + Newmark2 * TensorAngularVelocityxInertia; //Simo
    MassTerm2 *= Newmark1;

    Matrix MassMatrixBlock2(3,3);
    noalias(MassMatrixBlock2) = (-1) * MassTerm1 + MassTerm2;

    // Compute Linear Part of the Step Rotation
    Matrix LinearPartRotationTensor(3,3);

    double NormCurrentStepRotation =  norm_2(CurrentStepRotationVector);
    if( NormCurrentStepRotation != 0 ){
      this->CalculateRotationLinearPartTensor( CurrentStepRotationVector , LinearPartRotationTensor );
    }
    else{
      //std::cout<<" Attention... problem sure "<<std::endl;
      noalias(LinearPartRotationTensor) = IdentityMatrix(3);
    }

    //this->CalculateRotationLinearPartTensor( CurrentStepRotationVector , LinearPartRotationTensor );
    // std::cout<<" MassMatrixBlock 2 SPA "<<MassMatrixBlock2<<std::endl;
    // std::cout<<" LinearPartRotationTensor SPA "<<LinearPartRotationTensor<<std::endl;

    //block 1 of the mass matrix
    MatrixType m11(dimension,dimension);
    noalias(m11) = IdentityMatrix(dimension);
    m11 *= (1.0-AlphaM) * Newmark1 * rVariables.RigidBody.Mass;

    //block 2 of the mass matrix
    MatrixType m22(3,3);
    noalias(m22) = prod(MassMatrixBlock2, LinearPartRotationTensor);

    //Building the Local Tangent Inertia Matrix
    BeamMathUtilsType::AddMatrix(rLeftHandSideMatrix, m11, 0, 0);
    if(dimension == 2)
      rLeftHandSideMatrix(2,2) = m22(2,2);
    else
      BeamMathUtilsType::AddMatrix(rLeftHandSideMatrix, m22, 3, 3);


    // std::cout<<" Rigid Body LeftHandSideMatrix "<<rLeftHandSideMatrix<<std::endl;

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

void RigidBodyElement::GetTimeIntegrationParameters(double& rP0,double& rP1,double& rP2, const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  double DeltaTime = rCurrentProcessInfo[DELTA_TIME];
  rP0 = 1.0;
  rP1 = (1.0/ ( DeltaTime * DeltaTime * rCurrentProcessInfo[NEWMARK_BETA] ));
  rP2 = ( DeltaTime * rCurrentProcessInfo[NEWMARK_GAMMA] );

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

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size = this->GetDofsSize();

    if(rRightHandSideVector.size() != dofs_size)
      rRightHandSideVector.resize(dofs_size, false);

    noalias(rRightHandSideVector) = ZeroVector( dofs_size );

    //Current Compound Rotation Vector
    ArrayType CurrentCompoundRotationVector = GetGeometry()[0].FastGetSolutionStepValue(ROTATION);
    CurrentCompoundRotationVector = MapToInitialLocalFrame(CurrentCompoundRotationVector);

    //CurrentLinear Acceleration Vector
    ArrayType CurrentLinearAccelerationVector = GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION);
    CurrentLinearAccelerationVector = MapToInitialLocalFrame(CurrentLinearAccelerationVector);

    //PreviousLinear Acceleration Vector
    ArrayType PreviousLinearAccelerationVector = GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION,1);
    PreviousLinearAccelerationVector = MapToInitialLocalFrame(PreviousLinearAccelerationVector);

    //Angular Velocity Vector
    ArrayType AngularVelocityVector = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
    AngularVelocityVector = MapToInitialLocalFrame(AngularVelocityVector);

    //CurrentAngular Acceleration Vector
    ArrayType CurrentAngularAccelerationVector = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_ACCELERATION);
    CurrentAngularAccelerationVector = MapToInitialLocalFrame(CurrentAngularAccelerationVector);

    //PreviousAngular Acceleration Vector
    ArrayType PreviousAngularAccelerationVector = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_ACCELERATION,1);
    PreviousAngularAccelerationVector = MapToInitialLocalFrame(PreviousAngularAccelerationVector);

    double AlphaM = rCurrentProcessInfo[BOSSAK_ALPHA];
    ArrayType LinearAccelerationVector  = (1.0-AlphaM) * CurrentLinearAccelerationVector + AlphaM * (PreviousLinearAccelerationVector);
    ArrayType AngularAccelerationVector = (1.0-AlphaM) * CurrentAngularAccelerationVector + AlphaM * (PreviousAngularAccelerationVector);

    QuaternionType TotalQuaternion = QuaternionType::FromRotationVector(CurrentCompoundRotationVector);
    Matrix CurrentRotationMatrix   = ZeroMatrix(3,3);
    TotalQuaternion.ToRotationMatrix( CurrentRotationMatrix );

    //-----------------
    //block 1 of the inertial force vector

    //Compute Linear Term:

    Vector LinearInertialForceVector(dimension);

    for(SizeType i=0; i<dimension; ++i)
      LinearInertialForceVector[i] = rVariables.RigidBody.Mass * LinearAccelerationVector[i];

    //-----------------
    //block 2 of the inertial force vector (rotations part::to be defined)

    //Get inertia dyadic
    Matrix InertiaDyadic(3,3);
    noalias(InertiaDyadic) = rVariables.RigidBody.InertiaTensor;
    //Inertia dyadic expressed in the initial frame
    InertiaDyadic = prod(CurrentRotationMatrix,InertiaDyadic);
    InertiaDyadic = prod(InertiaDyadic,trans(CurrentRotationMatrix));

    //Compute Angular Term:

    Vector InertiaxAngularVelocity(3);
    noalias(InertiaxAngularVelocity) = prod(InertiaDyadic, AngularVelocityVector);
    Vector InertiaxAngularAcceleration(3);
    noalias(InertiaxAngularAcceleration) = prod(InertiaDyadic, AngularAccelerationVector);

    Matrix TensorAngularVelocity(3,3);
    BeamMathUtilsType::VectorToSkewSymmetricTensor(AngularVelocityVector, TensorAngularVelocity);

    Vector AngularInertialForceVector(3);
    noalias(AngularInertialForceVector) = prod(TensorAngularVelocity, InertiaxAngularVelocity);

    // CROSS PRODUCT of AxB = prod( skewA, B )  where (skewA) = [Ax] = [A]^v =>  (skewA)^T = hat(A) (nomenclature)
    AngularInertialForceVector += InertiaxAngularAcceleration;

    //compose total acceleration integral function:

    Vector TotalInertialForceVector(dofs_size);
    noalias(TotalInertialForceVector) = ZeroVector(dofs_size);

    BeamMathUtilsType::AddVector(LinearInertialForceVector, TotalInertialForceVector, 0);

    if(dimension == 2)
      TotalInertialForceVector[dofs_size-1] += AngularInertialForceVector[2];
    else
      BeamMathUtilsType::AddVector(AngularInertialForceVector, TotalInertialForceVector, 3);


    BeamMathUtilsType::AddVector(TotalInertialForceVector, rRightHandSideVector, 0);

    //std::cout<<" Rigid Body: rRightHandSideVector "<<rRightHandSideVector<<std::endl;

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyElement::CalculateMassMatrix(MatrixType& rMassMatrix,
                                           ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size = this->GetDofsSize();

    if(rMassMatrix.size1() != dofs_size || rMassMatrix.size2() != dofs_size)
      rMassMatrix.resize(dofs_size,dofs_size,false);

    noalias(rMassMatrix) = ZeroMatrix(dofs_size,dofs_size);

    // Rigid Body Properties
    RigidBodyProperties RigidBody;
    this->CalculateRigidBodyProperties(RigidBody);

    //block 1 of the mass matrix
    MatrixType m11(dimension,dimension);
    noalias(m11) = ZeroMatrix(dimension,dimension);
    for(SizeType i=0; i<dimension; ++i)
    {
      m11(i,i) = RigidBody.Mass;
    }

    //block 2 of the mass matrix
    MatrixType m22(3,3);
    noalias(m22) = RigidBody.InertiaTensor;
    m22 *= RigidBody.Mass;

    //Building the Local Tangent Inertia Matrix
    BeamMathUtilsType::AddMatrix(rMassMatrix, m11, 0, 0);

    if(dimension == 2)
      rMassMatrix(2,2) = m22(2,2);
    else
      BeamMathUtilsType::AddMatrix(rMassMatrix, m22, 3, 3);

    KRATOS_CATCH("")
}



//************************************************************************************
//************************************************************************************

void RigidBodyElement::CalculateRotationLinearPartTensor(ArrayType& rRotationVector,
                                                         Matrix& rRotationTensor)

{
    KRATOS_TRY

    if( rRotationTensor.size1() != 3 )
      rRotationTensor.resize(3, 3, false);

    // adding term 3
    BeamMathUtilsType::VectorToSkewSymmetricTensor(rRotationVector, rRotationTensor);

    rRotationTensor *= 0.5;

    double NormRotation =  norm_2(rRotationVector);
    Matrix RotationxRotation(3,3);
    noalias(RotationxRotation) = outer_prod(rRotationVector, rRotationVector);

    if( NormRotation != 0 )
      RotationxRotation *= ( 1.0/(NormRotation * NormRotation) );

    // adding term 1
    rRotationTensor = RotationxRotation;

    double Coefficient = 0;
    if( NormRotation != 0 )
      Coefficient = 0.5 * NormRotation / std::tan(0.5 * NormRotation);

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

     Node<3>::Pointer rCenterOfGravity = this->GetGeometry()(0);

     if( rCenterOfGravity->Is(SLAVE) ){
       Element& MasterElement = this->GetGeometry()[0].GetValue(MASTER_ELEMENTS).back();
       rCenterOfGravity = MasterElement.GetGeometry()(0);
     }

     ArrayType&  Center              = rCenterOfGravity->GetInitialPosition();

     ArrayType&  Displacement        = rCenterOfGravity->FastGetSolutionStepValue(DISPLACEMENT);
     ArrayType&  Rotation            = rCenterOfGravity->FastGetSolutionStepValue(ROTATION);
     ArrayType&  StepRotation        = rCenterOfGravity->FastGetSolutionStepValue(STEP_ROTATION);

     ArrayType&  Velocity            = rCenterOfGravity->FastGetSolutionStepValue(VELOCITY);
     ArrayType&  Acceleration        = rCenterOfGravity->FastGetSolutionStepValue(ACCELERATION);
     ArrayType&  AngularVelocity     = rCenterOfGravity->FastGetSolutionStepValue(ANGULAR_VELOCITY);
     ArrayType&  AngularAcceleration = rCenterOfGravity->FastGetSolutionStepValue(ANGULAR_ACCELERATION);

     std::cout<<" [ MasterElement "<<this->Id()<<" ]"<<std::endl;
     std::cout<<" [ Nodes_size "<<mpNodes->size()<<" ]"<<std::endl;
     std::cout<<" [ Fixed DisplacementY DOF: "<<this->GetGeometry()[0].IsFixed(DISPLACEMENT_Y)<<"]"<<std::endl;
     std::cout<<" [ Rotation:"<<Rotation<<",StepRotation:"<<StepRotation<<"]"<<std::endl;
     std::cout<<" [ Velocity:"<<Velocity<<",Acceleration:"<<Acceleration<<",Displacement:"<<Displacement<<",DeltaDisplacement"<<Displacement-rCenterOfGravity->FastGetSolutionStepValue(DISPLACEMENT,1)<<"]"<<std::endl;
     std::cout<<" [ AngularVelocity:"<<AngularVelocity<<",AngularAcceleration:"<<AngularAcceleration<<"]"<<std::endl;


     ArrayType Radius;
     ArrayType Variable;
     Matrix SkewSymVariable(3,3);

     for (NodesContainerType::iterator i = mpNodes->begin(); i != mpNodes->end(); ++i)
     {
       //Get rotation matrix
       QuaternionType TotalQuaternion = QuaternionType::FromRotationVector<ArrayType>(Rotation);

       Radius = (i)->GetInitialPosition() - Center;

       Matrix RotationMatrix(3,3);
       TotalQuaternion.ToRotationMatrix(RotationMatrix);
       Radius = prod(RotationMatrix, Radius);

       noalias((i)->FastGetSolutionStepValue(DISPLACEMENT)) = ( (Center + Displacement) + Radius ) - (i)->GetInitialPosition();
       noalias((i)->FastGetSolutionStepValue(ROTATION)) = Rotation;
       noalias((i)->FastGetSolutionStepValue(STEP_ROTATION)) = StepRotation;
       noalias((i)->FastGetSolutionStepValue(ANGULAR_VELOCITY)) = AngularVelocity;
       noalias((i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION)) = AngularAcceleration;
       noalias((i)->FastGetSolutionStepValue(VELOCITY)) = Velocity;
       noalias((i)->FastGetSolutionStepValue(ACCELERATION)) = Acceleration;

       // Update velocity:

       //compute the skewsymmmetric tensor of the angular velocity
       BeamMathUtilsType::VectorToSkewSymmetricTensor(AngularVelocity, SkewSymVariable);

       //compute the contribution of the angular velocity to the velocity v = Wxr
       noalias(Variable) = prod(SkewSymVariable,Radius);

       noalias((i)->FastGetSolutionStepValue(VELOCITY)) += Variable;

       // Update Acceleration:

       //centripetal acceleration:
       Variable = prod(SkewSymVariable,Variable); //ac = Wx(Wxr)

       noalias((i)->FastGetSolutionStepValue(ACCELERATION)) += Variable;

       //compute the skewsymmmetric tensor of the angular acceleration
       BeamMathUtilsType::VectorToSkewSymmetricTensor(AngularAcceleration, SkewSymVariable);

       //compute the contribution of the angular velocity to the velocity a = Axr
       noalias(Variable) = prod(SkewSymVariable,Radius);

       noalias((i)->FastGetSolutionStepValue(ACCELERATION)) += Variable;

       // std::cout<<"  [ Finalize Rigid Body Link Point : [Id:"<<(i)->Id()<<"] "<<std::endl;
       // std::cout<<"  [ Displacement:"<<NodeDisplacement<<" / StepRotation"<<NodeStepRotation<<" ] "<<std::endl;
       // std::cout<<"  [ Rotation:"<<NodeRotation<<" / Angular Acceleration"<<AngularAcceleration<<" ] "<<std::endl;

     }

     KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

RigidBodyElement::SizeType RigidBodyElement::GetDofsSize()
{
  KRATOS_TRY

  const SizeType dimension = GetGeometry().WorkingSpaceDimension();

  return (dimension * (dimension + 1) * 0.5);

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

    if(GetGeometry().size()!=1)
    {
      KRATOS_THROW_ERROR( std::invalid_argument, "This element works only with 1 noded geometry", "")
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
