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
  mStepVariable = VELOCITY_STEP;
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

RigidBodySegregatedVElement::RigidBodySegregatedVElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, NodesContainerType::Pointer pNodes)
    :RigidBodyElement(NewId, pGeometry, pProperties)
{
  mStepVariable = VELOCITY_STEP;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

RigidBodySegregatedVElement::RigidBodySegregatedVElement(RigidBodySegregatedVElement const& rOther)
    :RigidBodyElement(rOther)
    ,mStepVariable(rOther.mStepVariable)
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
  NewElement.mStepVariable = mStepVariable;
  
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
  this->SetProcessInformation(rCurrentProcessInfo);
  
  switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
  {
    case VELOCITY_STEP:
      {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();
        const SizeType dofs_size = this->GetDofsSize();

        if ( rResult.size() != dofs_size )
          rResult.resize(dofs_size, false);

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
        const SizeType dofs_size = this->GetDofsSize();
        if ( rResult.size() != dofs_size )
          rResult.resize(dofs_size, false);
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

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyElement::GetValuesVector(rValues, Step);
        break;
      }
    case PRESSURE_STEP:
      {
        SizeType dofs_size = this->GetDofsSize();
        if ( rValues.size() != dofs_size )
          rValues.resize( dofs_size, false );
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

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyElement::GetFirstDerivativesVector(rValues, Step);
        break;
      }
    case PRESSURE_STEP:
      {
        SizeType dofs_size = this->GetDofsSize();
        if ( rValues.size() != dofs_size )
          rValues.resize( dofs_size, false );
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

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyElement::GetSecondDerivativesVector(rValues, Step);
        break;
      }
    case PRESSURE_STEP:
      {
        SizeType dofs_size = this->GetDofsSize();
        if ( rValues.size() != dofs_size )
          rValues.resize( dofs_size, false );
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

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyElement::CalculateSecondDerivativesContributions(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
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

void RigidBodySegregatedVElement::CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                                                ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyElement::CalculateSecondDerivativesLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
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

void RigidBodySegregatedVElement::CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
                                                                ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyElement::CalculateSecondDerivativesRHS(rRightHandSideVector, rCurrentProcessInfo);
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
        size = number_of_nodes * dimension * (dimension + 1) * 0.5; //size for velocity
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
    if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
      KRATOS_THROW_ERROR( std::invalid_argument,"missing variable VELOCITY on node ", this->GetGeometry()[i].Id() )
    if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false || this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false || this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
      KRATOS_THROW_ERROR( std::invalid_argument,"missing one of the dofs for the variable VELOCITY on node ", GetGeometry()[i].Id() )
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

void RigidBodySegregatedVElement::save( Serializer& rSerializer ) const
{
  KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, RigidBodyElement )
}

void RigidBodySegregatedVElement::load( Serializer& rSerializer )
{
  KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, RigidBodyElement )
}


} // Namespace Kratos
