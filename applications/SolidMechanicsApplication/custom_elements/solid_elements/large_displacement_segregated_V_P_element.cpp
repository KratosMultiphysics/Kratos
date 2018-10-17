//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:               April 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/solid_elements/large_displacement_segregated_V_P_element.hpp"
#include "solid_mechanics_application_variables.h"


namespace Kratos
{


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementSegregatedVPElement::LargeDisplacementSegregatedVPElement()
    : LargeDisplacementVElement()
{
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementSegregatedVPElement::LargeDisplacementSegregatedVPElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : LargeDisplacementVElement( NewId, pGeometry )
{
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementSegregatedVPElement::LargeDisplacementSegregatedVPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : LargeDisplacementVElement( NewId, pGeometry, pProperties )
{
  mStepVariable = VELOCITY_STEP;
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LargeDisplacementSegregatedVPElement::LargeDisplacementSegregatedVPElement( LargeDisplacementSegregatedVPElement const& rOther)
    :LargeDisplacementVElement(rOther)
    ,mStepVariable(rOther.mStepVariable)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

LargeDisplacementSegregatedVPElement&  LargeDisplacementSegregatedVPElement::operator=(LargeDisplacementSegregatedVPElement const& rOther)
{
    LargeDisplacementVElement::operator=(rOther);

    mStepVariable = rOther.mStepVariable;

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer LargeDisplacementSegregatedVPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
  return Kratos::make_shared< LargeDisplacementSegregatedVPElement >(NewId, GetGeometry().Create(rThisNodes), pProperties);
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer LargeDisplacementSegregatedVPElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    KRATOS_ERROR << "Calling Clone for a large displacement segregated VP 3D element ... illegal operation!!" << std::endl;

    LargeDisplacementSegregatedVPElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    //-----------//

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;


    if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
	NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

	if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
	  KRATOS_ERROR << "constitutive law not has the correct size " << NewElement.mConstitutiveLawVector.size() << std::endl;
      }

    for(SizeType i=0; i<mConstitutiveLawVector.size(); i++)
    {
      NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
    }

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Kratos::make_shared< LargeDisplacementSegregatedVPElement >(NewElement);
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementSegregatedVPElement::~LargeDisplacementSegregatedVPElement()
{
}


//************* GETTING METHODS
//************************************************************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
    rElementalDofList.resize(0);

    switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
    {
      case VELOCITY_STEP:
        {
          const SizeType dimension  = GetGeometry().WorkingSpaceDimension();
          for ( SizeType i = 0; i < GetGeometry().size(); i++ )
          {
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_X ) );
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_Y ) );

            if( dimension == 3 )
              rElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_Z ) );
          }
          break;
        }
      case PRESSURE_STEP:
        {
          for ( SizeType i = 0; i < GetGeometry().size(); i++ )
          {
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( PRESSURE ) );
          }
          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
    }
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{

    this->SetProcessInformation(rCurrentProcessInfo);

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size = this->GetDofsSize();

    if ( rResult.size() != dofs_size )
        rResult.resize( dofs_size, false );

    switch(mStepVariable)
    {
      case VELOCITY_STEP:
        {
          for ( SizeType i = 0; i < number_of_nodes; i++ )
          {
            SizeType index = i * dimension;
            rResult[index]     = GetGeometry()[i].GetDof( VELOCITY_X ).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof( VELOCITY_Y ).EquationId();

            if( dimension == 3)
              rResult[index + 2] = GetGeometry()[i].GetDof( VELOCITY_Z ).EquationId();
          }
          break;
        }
      case PRESSURE_STEP:
        {
          for ( SizeType i = 0; i < number_of_nodes; i++ )
          {
            rResult[i] = GetGeometry()[i].GetDof( PRESSURE ).EquationId();
          }
          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
    }

}


//************************************************************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    this->SetProcessInformation(rCurrentProcessInfo);

    SolidElement::InitializeExplicitContributions();

    switch(mStepVariable)
    {
      case VELOCITY_STEP:
        {

          for ( SizeType i = 0; i < mConstitutiveLawVector.size(); i++ )
            mConstitutiveLawVector[i]->InitializeSolutionStep( GetProperties(),
                                                               GetGeometry(),
                                                               row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ),
                                                               rCurrentProcessInfo );
          break;
        }
      case PRESSURE_STEP:
        {
          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
    }

    this->Set(SolidElement::FINALIZED_STEP,false);

    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  this->SetProcessInformation(rCurrentProcessInfo);

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY


  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    this->SetProcessInformation(rCurrentProcessInfo);

    switch(mStepVariable)
    {
      case VELOCITY_STEP:
        {
          SolidElement::FinalizeSolutionStep(rCurrentProcessInfo);
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

void LargeDisplacementSegregatedVPElement::SetProcessInformation(const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  mStepVariable = StepType(rCurrentProcessInfo[SEGREGATED_STEP]);

  KRATOS_CATCH( "" )
}

//*********************************DISPLACEMENT***************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::GetValuesVector( Vector& rValues, int Step )
{
    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size = this->GetDofsSize();

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    switch(mStepVariable)
    {
      case VELOCITY_STEP:
        {
          SizeType index = 0;
          for ( SizeType i = 0; i < number_of_nodes; i++ )
          {
            index = i * dimension;
            rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
            rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

            if ( dimension == 3 )
              rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );

          }
          break;
        }
      case PRESSURE_STEP:
        {
          for ( SizeType i = 0; i < number_of_nodes; i++ )
          {
            rValues[i]     = GetGeometry()[i].GetSolutionStepValue( PRESSURE, Step );
          }
          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
    }

}


//************************************VELOCITY****************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size = this->GetDofsSize();

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    switch(mStepVariable)
    {
      case VELOCITY_STEP:
        {
          SizeType index = 0;
          for ( SizeType i = 0; i < number_of_nodes; i++ )
          {
            index = i * dimension;
            rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
            rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

            if ( dimension == 3 )
              rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
          }
          break;
        }
      case PRESSURE_STEP:
        {
          for ( SizeType i = 0; i < number_of_nodes; i++ )
          {
            rValues[i]     = GetGeometry()[i].GetSolutionStepValue( PRESSURE_VELOCITY, Step );
          }
          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
    }

}

//*********************************ACCELERATION***************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size = this->GetDofsSize();

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    switch(mStepVariable)
    {
      case VELOCITY_STEP:
        {
          SizeType index = 0;
          for ( SizeType i = 0; i < number_of_nodes; i++ )
          {
            index = i * dimension;
            rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
            rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

            if ( dimension == 3 )
              rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
          }
          break;
        }
      case PRESSURE_STEP:
        {
          for ( SizeType i = 0; i < number_of_nodes; i++ )
          {
            rValues[i]     = GetGeometry()[i].GetSolutionStepValue( PRESSURE_ACCELERATION, Step );
          }
          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
    }

}


//************************************************************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //process information
    this->SetProcessInformation(rCurrentProcessInfo);

    SolidElement::CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //process information
    this->SetProcessInformation(rCurrentProcessInfo);

    SolidElement::CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //process information
    this->SetProcessInformation(rCurrentProcessInfo);

    SolidElement::CalculateLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
{
    KRATOS_TRY

    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

    switch(mStepVariable)
    {
      case VELOCITY_STEP:
        {
          // operation performed: add Km to the rLefsHandSideMatrix
          this->CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

          // operation performed: add Kg to the rLefsHandSideMatrix
          this->CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

          rLeftHandSideMatrix *= rVariables.GetProcessInfo()[DELTA_TIME]; // backward Euler Approach (BDF order 1)

          break;
        }
      case PRESSURE_STEP:
        {

          // operation performed: add Kpp to the rLefsHandSideMatrix
          this->CalculateAndAddKpp( rLeftHandSideMatrix, rVariables );

          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
    }

    //KRATOS_WATCH( rLeftHandSideMatrix )

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
{
  KRATOS_TRY

  VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        // operation performed: add InternalForces to the rRightHandSideVector
        this->CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight );

        // operation performed: add ExternalForces to the rRightHandSideVector
        this->CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

        break;
      }
    case PRESSURE_STEP:
      {
        // operation performed: add PressureForces to the rRightHandSideVector
        this->CalculateAndAddPressureForces( rRightHandSideVector, rVariables );

        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
  }

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::CalculateAndAddKpp(MatrixType& rLeftHandSideMatrix,
                                                              ElementDataType& rVariables)

{
  KRATOS_TRY

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
                                                                         ElementDataType & rVariables)
{
  KRATOS_TRY

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  // the internal step variable must be set because InitializeNonLinearIteration is not called before this method
  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        SolidElement::CalculateMassMatrix( rMassMatrix, rCurrentProcessInfo );
        break;
      }
    case PRESSURE_STEP:
      {
        const SizeType MatSize = this->GetDofsSize();
        if ( rMassMatrix.size1() != MatSize )
          rMassMatrix.resize( MatSize, MatSize, false );

        noalias(rMassMatrix) = ZeroMatrix( MatSize, MatSize );

        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  // the internal step variable must be set because InitializeNonLinearIteration is not called before this method
  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        SolidElement::CalculateDampingMatrix( rDampingMatrix, rCurrentProcessInfo );
        break;
      }
    case PRESSURE_STEP:
      {
        const SizeType MatSize = this->GetDofsSize();
        if ( rDampingMatrix.size1() != MatSize )
          rDampingMatrix.resize( MatSize, MatSize, false );

        noalias(rDampingMatrix) = ZeroMatrix( MatSize, MatSize );

        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

LargeDisplacementSegregatedVPElement::SizeType LargeDisplacementSegregatedVPElement::GetDofsSize()
{
  KRATOS_TRY

  const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
  const SizeType number_of_nodes  = GetGeometry().PointsNumber();

  SizeType size = 0;
  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      size = number_of_nodes * dimension; //size for velocity
      break;
    case PRESSURE_STEP:
      size = number_of_nodes; //size for pressure
      break;
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;

  }
  return size;

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

int  LargeDisplacementSegregatedVPElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Perform base element checks
    int ErrorCode = 0;
    ErrorCode = LargeDisplacementVElement::Check(rCurrentProcessInfo);

    // Check that the element nodes contain all required SolutionStepData and Degrees of freedom
    for(SizeType i=0; i<this->GetGeometry().size(); ++i)
      {
	// Nodal data
	Node<3> &rNode = this->GetGeometry()[i];
	KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE,rNode);
        //KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE_VELOCITY,rNode);
        //KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE_ACCELERATION,rNode);
	//KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUME_ACCELERATION,rNode);

	// Nodal dofs
	KRATOS_CHECK_DOF_IN_NODE(PRESSURE,rNode);
      }
    // Check compatibility with the constitutive law

    // Check that all required variables have been registered

    return ErrorCode;

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************


void LargeDisplacementSegregatedVPElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementVElement )
    int IntStepType = int(mStepVariable);
    rSerializer.save("StepVariable",IntStepType);
}

void LargeDisplacementSegregatedVPElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementVElement )
    int IntStepType;
    rSerializer.load("StepVariable",IntStepType);
    mStepVariable = StepType(IntStepType);
}


} // Namespace Kratos
