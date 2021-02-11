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
#include "custom_elements/solid_elements/large_displacement_V_element.hpp"
#include "solid_mechanics_application_variables.h"


namespace Kratos
{


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementVElement::LargeDisplacementVElement()
    : LargeDisplacementElement()
{
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementVElement::LargeDisplacementVElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : LargeDisplacementElement( NewId, pGeometry )
{
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementVElement::LargeDisplacementVElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : LargeDisplacementElement( NewId, pGeometry, pProperties )
{
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LargeDisplacementVElement::LargeDisplacementVElement( LargeDisplacementVElement const& rOther)
    :LargeDisplacementElement(rOther)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

LargeDisplacementVElement&  LargeDisplacementVElement::operator=(LargeDisplacementVElement const& rOther)
{
    LargeDisplacementElement::operator=(rOther);

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer LargeDisplacementVElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
  return Kratos::make_intrusive< LargeDisplacementVElement >(NewId, GetGeometry().Create(rThisNodes), pProperties);
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer LargeDisplacementVElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    KRATOS_THROW_ERROR( std::logic_error, "calling the default constructor for a large displacement 3D element ... illegal operation!!", "" )

    LargeDisplacementVElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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

    return Kratos::make_intrusive< LargeDisplacementVElement >(NewElement);
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementVElement::~LargeDisplacementVElement()
{
}


//************* GETTING METHODS
//************************************************************************************
//************************************************************************************



void LargeDisplacementVElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
    rElementalDofList.resize( 0 );

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    for ( SizeType i = 0; i < GetGeometry().size(); i++ )
    {
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_X ) );
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_Y ) );

        if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_Z ) );
    }
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementVElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size        = GetDofsSize();

    if ( rResult.size() != dofs_size )
        rResult.resize( dofs_size, false );

    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dimension;
        rResult[index]     = GetGeometry()[i].GetDof( VELOCITY_X ).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof( VELOCITY_Y ).EquationId();

        if( dimension == 3)
          rResult[index + 2] = GetGeometry()[i].GetDof( VELOCITY_Z ).EquationId();

    }

}


//************************************************************************************
//************************************************************************************

void LargeDisplacementVElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
{
    KRATOS_TRY

    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

    // operation performed: add Km to the rLefsHandSideMatrix
    this->CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // operation performed: add Kg to the rLefsHandSideMatrix
    this->CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    rLeftHandSideMatrix *= rVariables.GetProcessInfo()[DELTA_TIME]; // backward Euler Approach (BDF order 1)

    //KRATOS_WATCH( rLeftHandSideMatrix )

    KRATOS_CATCH( "" )
}




//************************************************************************************
//************************************************************************************

LargeDisplacementVElement::SizeType LargeDisplacementVElement::GetDofsSize() const
{
  KRATOS_TRY

  const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
  const SizeType number_of_nodes  = GetGeometry().PointsNumber();

  SizeType size = number_of_nodes * dimension; //usual size for velocity based elements

  return size;

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementVElement::SetElementData(ElementDataType& rVariables,
                                               ConstitutiveLaw::Parameters& rValues,
                                               const int & rPointNumber)
{

    //to take in account previous step for output print purposes
    unsigned int Alpha = 1; //current step
    if( this->Is(SolidElement::FINALIZED_STEP) ){
      Alpha = 0; //previous step
      this->GetHistoricalVariables(rVariables,rPointNumber);
    }

    if(rVariables.detF<0){

      KRATOS_WARNING(" [Element Ignored]") << "LargeDispVElement["<<this->Id()<<"] (|F|=" << rVariables.detF <<")  (Iter:"<<rVariables.GetProcessInfo()[NL_ITERATION_NUMBER]<<")"<<std::endl;
      rVariables.detJ = 0;

     SizeType number_of_nodes  = GetGeometry().PointsNumber();

      for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
        const array_1d<double, 3> & CurrentPosition      = GetGeometry()[i].Coordinates();
        const array_1d<double, 3> & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3> & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        array_1d<double, 3> PreviousPosition       = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
        KRATOS_WARNING("")<<" Node["<<GetGeometry()[i].Id()<<"]: (Position: (pre)"<<PreviousPosition<<",(cur)"<<CurrentPosition<<")"<<std::endl;
        //KRATOS_WARNING("")<<" (Displacement: (pre)"<<CurrentDisplacement<<",(cur)"<<PreviousDisplacement<<")"<<std::endl;
      }
      for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
        if( GetGeometry()[i].SolutionStepsDataHas(CONTACT_FORCE) ){
          const array_1d<double, 3 > & PreContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE,1);
          const array_1d<double, 3 > & ContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE);
          KRATOS_WARNING("")<<" (Contact: (pre)"<<PreContactForce<<",(cur)"<<ContactForce<<")["<<GetGeometry()[i].Id()<<"]"<<std::endl;
        }

      }

      this->Set(SELECTED,true);

      KRATOS_ERROR<<" [Element Failed] ["<<this->Id()<<"]"<<std::endl;

    }
    else{

      if(this->Is(SELECTED) && this->Is(ACTIVE)){
        this->Set(SELECTED,false);
        KRATOS_WARNING("")<<" Undo SELECTED LargeDispVElement "<<this->Id()<<std::endl;
      }

    }


    //Compute strain rate measures if they are required by the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    mConstitutiveLawVector[rPointNumber]->GetLawFeatures(LawFeatures);

    bool strain_rate_measure = false;
    for(SizeType i=0; i<LawFeatures.mStrainMeasures.size(); i++)
    {
      if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Velocity_Gradient)
	strain_rate_measure = true;
    }

    if( strain_rate_measure ){
      //Compute symmetric spatial velocity gradient [DN_DX = dN/dx_n*1] stored in a vector
      GeometryType& rGeometry = GetGeometry();
      ElementUtilities::CalculateVelocityGradientVector( rVariables.StrainVector, rGeometry, rVariables.DN_DX, Alpha );
      Flags &ConstitutiveLawOptions=rValues.GetOptions();
      ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    }

    //Compute F and detF (from 0 to n+1) : store it in H variable and detH
    rVariables.detH = rVariables.detF * rVariables.detF0;
    noalias(rVariables.H) = prod( rVariables.F, rVariables.F0 );

    rValues.SetDeterminantF(rVariables.detH);
    rValues.SetDeformationGradientF(rVariables.H);
    rValues.SetStrainVector(rVariables.StrainVector);
    rValues.SetStressVector(rVariables.StressVector);
    rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
    rValues.SetShapeFunctionsValues(rVariables.N);

}

//************************************************************************************
//************************************************************************************

int  LargeDisplacementVElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Perform base element checks
    int ErrorCode = 0;
    ErrorCode = LargeDisplacementElement::Check(rCurrentProcessInfo);

    // Check that the element nodes contain all required SolutionStepData and Degrees of freedom
    for(SizeType i=0; i<this->GetGeometry().size(); ++i)
      {
	// Nodal data
	Node<3> &rNode = this->GetGeometry()[i];
	KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rNode);
	//KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUME_ACCELERATION,rNode);

	// Nodal dofs
	KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X,rNode);
	KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y,rNode);
	if( rCurrentProcessInfo[SPACE_DIMENSION] == 3)
	  KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Z,rNode);
      }
    // Check compatibility with the constitutive law

    // Check that all required variables have been registered

    return ErrorCode;

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************


void LargeDisplacementVElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement )
}

void LargeDisplacementVElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement )
}


} // Namespace Kratos
