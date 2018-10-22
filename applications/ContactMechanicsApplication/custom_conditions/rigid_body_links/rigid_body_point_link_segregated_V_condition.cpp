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
#include "custom_conditions/rigid_body_links/rigid_body_point_link_segregated_V_condition.hpp"

#include "contact_mechanics_application_variables.h"

namespace Kratos
{


//***********************************************************************************
//***********************************************************************************
RigidBodyPointLinkSegregatedVCondition::RigidBodyPointLinkSegregatedVCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : RigidBodyPointLinkCondition(NewId, pGeometry)
{
}

//***********************************************************************************
//***********************************************************************************
RigidBodyPointLinkSegregatedVCondition::RigidBodyPointLinkSegregatedVCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : RigidBodyPointLinkCondition(NewId, pGeometry, pProperties)
{
  for ( SizeType i = 0; i < GetGeometry().size(); i++ )
  {
    GetGeometry()[i].Set(SLAVE); //Flag to set MASTER_ELEMENTS in that nodes (if is SLAVE, a MASTER is required)
  }
}

//************************************************************************************
//************************************************************************************
RigidBodyPointLinkSegregatedVCondition::RigidBodyPointLinkSegregatedVCondition( RigidBodyPointLinkSegregatedVCondition const& rOther )
    : RigidBodyPointLinkCondition(rOther)
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer RigidBodyPointLinkSegregatedVCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
  return Kratos::make_shared<RigidBodyPointLinkSegregatedVCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}


//************************************CLONE*******************************************
//************************************************************************************

Condition::Pointer RigidBodyPointLinkSegregatedVCondition::Clone(IndexType NewId, NodesArrayType const& rThisNodes) const
{
  return this->Create( NewId, rThisNodes, pGetProperties() );
}

//***********************************************************************************
//***********************************************************************************
RigidBodyPointLinkSegregatedVCondition::~RigidBodyPointLinkSegregatedVCondition()
{
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkSegregatedVCondition::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  rConditionDofList.resize(0);

  switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
  {
    case VELOCITY_STEP:
      {
        RigidBodyPointLinkCondition::GetDofList(rConditionDofList, rCurrentProcessInfo);
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

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkSegregatedVCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  rResult.resize( 0, false );

  switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
  {
    case VELOCITY_STEP:
      {
        RigidBodyPointLinkCondition::EquationIdVector(rResult, rCurrentProcessInfo);
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

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkSegregatedVCondition::GetValuesVector(Vector& rValues, int Step)
{
  KRATOS_TRY

  rValues.resize(0);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyPointLinkCondition::GetValuesVector(rValues, Step);
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

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkSegregatedVCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
{
  KRATOS_TRY

  rValues.resize(0);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyPointLinkCondition::GetFirstDerivativesVector(rValues, Step);
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


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkSegregatedVCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
{
  KRATOS_TRY

  rValues.resize(0);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyPointLinkCondition::GetSecondDerivativesVector(rValues, Step);
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

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkSegregatedVCondition::Initialize()
{
  KRATOS_TRY

  RigidBodyPointLinkCondition::Initialize();

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkSegregatedVCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyPointLinkCondition::InitializeSolutionStep(rCurrentProcessInfo);
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

void RigidBodyPointLinkSegregatedVCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                                                   VectorType& rRightHandSideVector,
                                                                   ProcessInfo& rCurrentProcessInfo )
{
  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyPointLinkCondition::CalculateLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
        break;
      }
    case PRESSURE_STEP:
      {
        rLeftHandSideMatrix.resize(0,0,false);
        rRightHandSideVector.resize(0);
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }

  //KRATOS_INFO("")<<mStepVariable<<" LHS:"<<rLeftHandSideMatrix<<" RHS:"<<rRightHandSideVector<<std::endl;

}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkSegregatedVCondition::CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
                                                                                     VectorType& rRightHandSideVector,
                                                                                     ProcessInfo& rCurrentProcessInfo)
{
  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyPointLinkCondition::CalculateSecondDerivativesContributions(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
        break;
      }
    case PRESSURE_STEP:
      {
        rLeftHandSideMatrix.resize(0,0,false);
        rRightHandSideVector.resize(0);
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }

  //KRATOS_INFO("")<<mStepVariable<<" 2LHS:"<<rLeftHandSideMatrix<<" 2RHS:"<<rRightHandSideVector<<std::endl;
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkSegregatedVCondition::CalculateRightHandSide( VectorType& rRightHandSideVector,
                                                                     ProcessInfo& rCurrentProcessInfo )
{

  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyPointLinkCondition::CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);
        break;
      }
    case PRESSURE_STEP:
      {
        rRightHandSideVector.resize(0);
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }

}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkSegregatedVCondition::CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                                                           ProcessInfo& rCurrentProcessInfo)
{

  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyPointLinkCondition::CalculateSecondDerivativesLHS(rLeftHandSideMatrix,rCurrentProcessInfo);
        break;
      }
    case PRESSURE_STEP:
      {
        rLeftHandSideMatrix.resize(0,0,false);
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }
}


//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkSegregatedVCondition::CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
                                                                           ProcessInfo& rCurrentProcessInfo)
{
  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyPointLinkCondition::CalculateSecondDerivativesRHS(rRightHandSideVector,rCurrentProcessInfo);
        break;
      }
    case PRESSURE_STEP:
      {
        rRightHandSideVector.resize(0);
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }
}


//************************************************************************************
//************************************************************************************

RigidBodyPointLinkSegregatedVCondition::SizeType RigidBodyPointLinkSegregatedVCondition::GetDofsSize()
{
  KRATOS_TRY

  SizeType size = 0;
  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        size = RigidBodyPointLinkCondition::GetDofsSize();
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

void RigidBodyPointLinkSegregatedVCondition::SetProcessInformation(const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  mStepVariable = StepType(rCurrentProcessInfo[SEGREGATED_STEP]);

  KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

int RigidBodyPointLinkSegregatedVCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
  return 0;
}

} // Namespace Kratos.
