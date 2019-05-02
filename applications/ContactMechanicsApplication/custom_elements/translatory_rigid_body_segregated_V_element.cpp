//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:               October 2018 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/translatory_rigid_body_segregated_V_element.hpp"

#include "contact_mechanics_application_variables.h"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TranslatoryRigidBodySegregatedVElement::TranslatoryRigidBodySegregatedVElement(IndexType NewId,GeometryType::Pointer pGeometry)
    :TranslatoryRigidBodyElement(NewId, pGeometry)
{
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TranslatoryRigidBodySegregatedVElement::TranslatoryRigidBodySegregatedVElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    :TranslatoryRigidBodyElement(NewId, pGeometry, pProperties)
{
  mStepVariable = VELOCITY_STEP;
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TranslatoryRigidBodySegregatedVElement::TranslatoryRigidBodySegregatedVElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, NodesContainerType::Pointer pNodes)
    :TranslatoryRigidBodyElement(NewId, pGeometry, pProperties, pNodes)
{
  mStepVariable = VELOCITY_STEP;
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

TranslatoryRigidBodySegregatedVElement::TranslatoryRigidBodySegregatedVElement(TranslatoryRigidBodySegregatedVElement const& rOther)
    :TranslatoryRigidBodyElement(rOther)
    ,mStepVariable(rOther.mStepVariable)
{
}

//*********************************CREATE*********************************************
//************************************************************************************

Element::Pointer TranslatoryRigidBodySegregatedVElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
  return Kratos::make_shared<TranslatoryRigidBodySegregatedVElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}


//*********************************CLONE**********************************************
//************************************************************************************

Element::Pointer TranslatoryRigidBodySegregatedVElement::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
{

  TranslatoryRigidBodySegregatedVElement NewElement( NewId, GetGeometry().Create(ThisNodes), pGetProperties(), mpNodes );

  NewElement.mInitialLocalQuaternion = this->mInitialLocalQuaternion;
  NewElement.SetData(this->GetData());
  NewElement.SetFlags(this->GetFlags());
  NewElement.mStepVariable = mStepVariable;

  return Kratos::make_shared<TranslatoryRigidBodySegregatedVElement>(NewElement);

}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

TranslatoryRigidBodySegregatedVElement::~TranslatoryRigidBodySegregatedVElement()
{
}

//************************************************************************************
//************************************************************************************

void TranslatoryRigidBodySegregatedVElement::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
{
  rElementalDofList.resize(0);

  switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
  {
    case VELOCITY_STEP:
      {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
          rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
          rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
          if( dimension ==3 )
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Z));
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

void TranslatoryRigidBodySegregatedVElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
  this->SetProcessInformation(rCurrentProcessInfo);

  switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
  {
    case VELOCITY_STEP:
      {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
        const SizeType dofs_size = this->GetDofsSize();

        if ( rResult.size() != dofs_size )
          rResult.resize(dofs_size, false);

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
          int index = i * ( dimension );
          rResult[index]   = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
          rResult[index+1] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
          if( dimension ==3 )
            rResult[index+2] = GetGeometry()[i].GetDof(VELOCITY_Z).EquationId();
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

void TranslatoryRigidBodySegregatedVElement::GetValuesVector(Vector& rValues, int Step)
{
  KRATOS_TRY

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        TranslatoryRigidBodyElement::GetValuesVector(rValues, Step);
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

  KRATOS_CATCH( "" )
}

//************************************VELOCITY****************************************
//************************************************************************************

//************************************************************************************
//************************************************************************************
void TranslatoryRigidBodySegregatedVElement::GetFirstDerivativesVector(Vector& rValues, int Step)
{
  KRATOS_TRY

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        TranslatoryRigidBodyElement::GetFirstDerivativesVector(rValues, Step);
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

void TranslatoryRigidBodySegregatedVElement::GetSecondDerivativesVector(Vector& rValues, int Step)
{
  KRATOS_TRY

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        TranslatoryRigidBodyElement::GetSecondDerivativesVector(rValues, Step);
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

void TranslatoryRigidBodySegregatedVElement::Initialize()
{
    KRATOS_TRY

    TranslatoryRigidBodyElement::Initialize();

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void TranslatoryRigidBodySegregatedVElement::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        TranslatoryRigidBodyElement::InitializeSolutionStep(rCurrentProcessInfo);
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

void TranslatoryRigidBodySegregatedVElement::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  this->SetProcessInformation(rCurrentProcessInfo);

  TranslatoryRigidBodyElement::InitializeNonLinearIteration(rCurrentProcessInfo);

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void TranslatoryRigidBodySegregatedVElement::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        TranslatoryRigidBodyElement::FinalizeNonLinearIteration(rCurrentProcessInfo);
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

void TranslatoryRigidBodySegregatedVElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        TranslatoryRigidBodyElement::FinalizeSolutionStep(rCurrentProcessInfo);
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

void TranslatoryRigidBodySegregatedVElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                                    ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  TranslatoryRigidBodyElement::CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void TranslatoryRigidBodySegregatedVElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                                   ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  TranslatoryRigidBodyElement::CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void TranslatoryRigidBodySegregatedVElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                                  VectorType& rRightHandSideVector,
                                                                  ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  TranslatoryRigidBodyElement::CalculateLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);

  //KRATOS_INFO("")<<mStepVariable<<" LHS:"<<rLeftHandSideMatrix<<" RHS:"<<rRightHandSideVector<<std::endl;

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void TranslatoryRigidBodySegregatedVElement::CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
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
        TranslatoryRigidBodyElement::CalculateSecondDerivativesContributions(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
        break;
      }
    case PRESSURE_STEP:
      {
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }

  //KRATOS_INFO("")<<mStepVariable<<" 2LHS:"<<rLeftHandSideMatrix<<" 2RHS:"<<rRightHandSideVector<<std::endl;

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void TranslatoryRigidBodySegregatedVElement::CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                                                           ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        TranslatoryRigidBodyElement::CalculateSecondDerivativesLHS(rLeftHandSideMatrix, rCurrentProcessInfo);
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

void TranslatoryRigidBodySegregatedVElement::CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
                                                                           ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        TranslatoryRigidBodyElement::CalculateSecondDerivativesRHS(rRightHandSideVector, rCurrentProcessInfo);
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

void TranslatoryRigidBodySegregatedVElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //process information
  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        TranslatoryRigidBodyElement::CalculateMassMatrix(rMassMatrix, rCurrentProcessInfo);
        break;
      }
    case PRESSURE_STEP:
      {
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }

  //KRATOS_INFO("")<<mStepVariable<<" MassM:"<<rMassMatrix<<std::endl;

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void TranslatoryRigidBodySegregatedVElement::GetTimeIntegrationParameters(double& rP0,double& rP1,double& rP2, const ProcessInfo& rCurrentProcessInfo)
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

TranslatoryRigidBodySegregatedVElement::SizeType TranslatoryRigidBodySegregatedVElement::GetDofsSize()
{
  KRATOS_TRY

  SizeType size = 0;
  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        size = TranslatoryRigidBodyElement::GetDofsSize();
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

//***********************************************************************************
//***********************************************************************************

void TranslatoryRigidBodySegregatedVElement::UpdateRigidBodyNodes(ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        TranslatoryRigidBodyElement::UpdateRigidBodyNodes(rCurrentProcessInfo);
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

void TranslatoryRigidBodySegregatedVElement::SetProcessInformation(const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  mStepVariable = StepType(rCurrentProcessInfo[SEGREGATED_STEP]);

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

int TranslatoryRigidBodySegregatedVElement::Check(const ProcessInfo& rCurrentProcessInfo)
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

    return 0;

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void TranslatoryRigidBodySegregatedVElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, TranslatoryRigidBodyElement )
}

void TranslatoryRigidBodySegregatedVElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, TranslatoryRigidBodyElement )
}


} // Namespace Kratos
