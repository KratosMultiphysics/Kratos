//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:             September 2018 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/rigid_body_segregated_V_element.hpp"
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

void RigidBodySegregatedVElement::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
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
          rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
          rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
          if( dimension == 2 ){
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));
          }
          else{
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Z));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));
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

void RigidBodySegregatedVElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
  SizeType dofs_size = this->GetDofsSize();

  if ( rResult.size() != dofs_size )
    rResult.resize(dofs_size, false);

  switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
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
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
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

  RigidBodyElement::CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);

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

void RigidBodySegregatedVElement::GetTimeIntegrationParameters(double& rP0,double& rP1,double& rP2, const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  double DeltaTime = rCurrentProcessInfo[DELTA_TIME];
  rP0 = rCurrentProcessInfo[NEWMARK_GAMMA] / DeltaTime *rCurrentProcessInfo[NEWMARK_BETA];
  rP1 = (1.0/ ( DeltaTime * DeltaTime * rCurrentProcessInfo[NEWMARK_BETA] ));
  rP2 = ( DeltaTime * rCurrentProcessInfo[NEWMARK_GAMMA] );

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

RigidBodySegregatedVElement::SizeType RigidBodySegregatedVElement::GetDofsSize()
{
  KRATOS_TRY

  SizeType size = 0;
  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();
        const SizeType number_of_nodes  = GetGeometry().PointsNumber();
        size = number_of_nodes*dimension*(dimension + 1)*0.5; //size for velocity
        break;
      }
    case PRESSURE_STEP:
      {
        size = 0;
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;

  }
  return size;

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void RigidBodySegregatedVElement::SetProcessInformation(const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  mStepVariable = StepType(rCurrentProcessInfo[SEGREGATED_STEP]);

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
