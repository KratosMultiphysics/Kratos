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
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LargeDisplacementSegregatedVPElement::LargeDisplacementSegregatedVPElement( LargeDisplacementSegregatedVPElement const& rOther)
    :LargeDisplacementVElement(rOther)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

LargeDisplacementSegregatedVPElement&  LargeDisplacementSegregatedVPElement::operator=(LargeDisplacementSegregatedVPElement const& rOther)
{
    LargeDisplacementVElement::operator=(rOther);

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer LargeDisplacementSegregatedVPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
  return Kratos::make_shared< LargeDisplacementSegregatedVPElement >(NewId, GetGeometry().Create(ThisNodes), pProperties);    
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
    rElementalDofList.resize( 0 );

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();
    
    switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
    {
      case VELOCITY_STEP:
        {
          for ( SizeType i = 0; i < GetGeometry().size(); i++ )
          {
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_X ) );
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_Y ) );

            if( dimension == 3 )
              rElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_Z ) );
          }
        }
      case PRESSURE_STEP:
        {
          for ( SizeType i = 0; i < GetGeometry().size(); i++ )
          {
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( PRESSURE ) );
          }
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
    }
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = GetDofsSize();

    if ( rResult.size() != dofs_size )
        rResult.resize( dofs_size, false );

    switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
    {
      case VELOCITY_STEP:
        {
          for ( SizeType i = 0; i < number_of_nodes; i++ )
          {
            int index = i * dimension;
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


//*********************************DISPLACEMENT***************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::GetValuesVector( Vector& rValues, int Step )
{
    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int   dofs_size        = GetDofsSize();
    
    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
    {
      case VELOCITY_STEP:
        {
          unsigned int index = 0;
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
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
    }

}


//************************************VELOCITY****************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int   dofs_size        = GetDofsSize();
    
    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
    {
      case VELOCITY_STEP:
        {
          unsigned int index = 0;
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
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
    }  
 
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = GetDofsSize();
    
    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
    {
      case VELOCITY_STEP:
        {
          unsigned int index = 0;
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
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
    }

}

//************************************************************************************
//************************************************************************************

void LargeDisplacementSegregatedVPElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
{
    KRATOS_TRY
       
    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix(); 

    switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
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
          this->CalculateAndAddKpp( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
    }
    
    //KRATOS_WATCH( rLeftHandSideMatrix )
  
    KRATOS_CATCH( "" )  
}




//************************************************************************************
//************************************************************************************

unsigned int LargeDisplacementSegregatedVPElement::GetDofsSize()
{
  KRATOS_TRY
     
  const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
  const SizeType number_of_nodes  = GetGeometry().PointsNumber();
  
  unsigned int size = 0;
  switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
  {
    case VELOCITY_STEP:
      size = number_of_nodes * dimension; //size for velocity
      break;
    case PRESSURE_STEP:
      size = number_of_nodes; //size for pressure
      break;
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
      
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
}

void LargeDisplacementSegregatedVPElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementVElement )
}


} // Namespace Kratos


