//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/solid_elements/large_displacement_element.hpp"
#include "solid_mechanics_application_variables.h"


namespace Kratos
{


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementElement::LargeDisplacementElement( )
    :SolidElement( )
{
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementElement::LargeDisplacementElement( IndexType NewId, GeometryType::Pointer pGeometry )
    :SolidElement( NewId, pGeometry )
{
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementElement::LargeDisplacementElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    :SolidElement( NewId, pGeometry, pProperties )
{
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LargeDisplacementElement::LargeDisplacementElement( LargeDisplacementElement const& rOther)
    :SolidElement(rOther)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

LargeDisplacementElement&  LargeDisplacementElement::operator=(LargeDisplacementElement const& rOther)
{
    SolidElement::operator=(rOther);

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer LargeDisplacementElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    KRATOS_ERROR << " calling the default method Create for a large displacement element " << std::endl;
    return Kratos::make_shared< LargeDisplacementElement >(NewId, GetGeometry().Create(rThisNodes), pProperties);
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer LargeDisplacementElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    KRATOS_ERROR << " calling the default method Clone for a large displacement element " << std::endl;

    LargeDisplacementElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;


    if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
	NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

	if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
	  KRATOS_ERROR << " constitutive law not has the correct size large displacement element " << std::endl;
      }

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Kratos::make_shared< LargeDisplacementElement >(NewElement);
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementElement::~LargeDisplacementElement()
{
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementElement::GetHistoricalVariables( ElementDataType& rVariables, const double& rPointNumber )
{
    //Deformation Gradient F ( set to identity )
    unsigned int size =  rVariables.F.size1();

    rVariables.detF = 1;
    noalias(rVariables.F) = IdentityMatrix(size);

}

//************************************************************************************
//************************************************************************************

void LargeDisplacementElement::SetElementData(ElementDataType& rVariables,
                                              ConstitutiveLaw::Parameters& rValues,
                                              const int & rPointNumber)
{

    //to take in account previous step for output print purposes
    if( this->Is(SolidElement::FINALIZED_STEP) ){
      this->GetHistoricalVariables(rVariables,rPointNumber);
    }

    //to be accurate calculus must stop
    if(rVariables.detF<0){

      KRATOS_WARNING(" [Element Ignored]") << "LargeDispElement["<<this->Id()<<"] (|F|=" << rVariables.detF <<")  (Iter:"<<rVariables.GetProcessInfo()[NL_ITERATION_NUMBER]<<")"<<std::endl;
      rVariables.detJ = 0;

      SizeType number_of_nodes  = GetGeometry().PointsNumber();

      for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
        array_1d<double, 3> & CurrentPosition      = GetGeometry()[i].Coordinates();
        array_1d<double, 3> & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double, 3> & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        array_1d<double, 3> PreviousPosition       = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
        KRATOS_WARNING("")<<" Node["<<GetGeometry()[i].Id()<<"]: (Position: (pre)"<<PreviousPosition<<",(cur)"<<CurrentPosition<<")"<<std::endl;
        //KRATOS_WARNING("")<<" (Displacement: (pre)"<<CurrentDisplacement<<",(cur)"<<PreviousDisplacement<<")"<<std::endl;
      }
      for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
        if( GetGeometry()[i].SolutionStepsDataHas(CONTACT_FORCE) ){
          array_1d<double, 3 > & PreContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE,1);
          array_1d<double, 3 > & ContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE);
          KRATOS_WARNING("")<<" (Contact: (pre)"<<PreContactForce<<",(cur)"<<ContactForce<<")["<<GetGeometry()[i].Id()<<"]"<<std::endl;
        }

      }

      this->Set(SELECTED,true);

      KRATOS_ERROR<<" [Element Failed] ["<<this->Id()<<"]"<<std::endl;

    }
    else{

      if(this->Is(SELECTED) && this->Is(ACTIVE)){
        this->Set(SELECTED,false);
        KRATOS_WARNING("")<<" Undo SELECTED LargeDispElement "<<this->Id()<<std::endl;
      }

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

void LargeDisplacementElement::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    SolidElement::InitializeSolutionStep(rCurrentProcessInfo);
    this->Set(SolidElement::FINALIZED_STEP, false);

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void LargeDisplacementElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    SolidElement::FinalizeSolutionStep(rCurrentProcessInfo);

    this->Set(SolidElement::FINALIZED_STEP, true);

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementElement::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
						   ElementDataType& rVariables,
						   double& rIntegrationWeight)

{
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();
    Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rVariables.StressVector );
    Matrix ReducedKg = prod( rVariables.DN_DX, rIntegrationWeight * Matrix( prod( StressTensor, trans( rVariables.DN_DX ) ) ) ); //to be optimized
    MathUtils<double>::ExpandAndAddReducedMatrix( rLeftHandSideMatrix, ReducedKg, dimension );

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementElement::CalculateGreenLagrangeStrain(const Matrix& rF, Vector& rStrainVector )
{
    KRATOS_TRY

    const SizeType dimension   = GetGeometry().WorkingSpaceDimension();

    //Right Cauchy-Green Calculation
    Matrix C ( dimension, dimension );
    noalias( C ) = prod( trans( rF ), rF );

    if( dimension == 2 )
    {

        //Green Lagrange Strain Calculation
        if ( rStrainVector.size() != 3 ) rStrainVector.resize( 3, false );

        rStrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

        rStrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

        rStrainVector[2] = C( 0, 1 ); // xy

    }
    else if( dimension == 3 )
    {

        //Green Lagrange Strain Calculation
        if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

        rStrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

        rStrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

        rStrainVector[2] = 0.5 * ( C( 2, 2 ) - 1.00 );

        rStrainVector[3] = C( 0, 1 ); // xy

        rStrainVector[4] = C( 1, 2 ); // yz

        rStrainVector[5] = C( 0, 2 ); // xz

    }
    else
    {
        KRATOS_ERROR << " something is wrong with the dimension green-lagrange strain " << std::endl;
    }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementElement::CalculateAlmansiStrain(const Matrix& rF, Vector& rStrainVector )
{
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    //Left Cauchy-Green Calculation
    Matrix LeftCauchyGreen(dimension, dimension);
    noalias(LeftCauchyGreen) = prod( rF, trans( rF ) );

    //Calculating the inverse of the jacobian
    Matrix InverseLeftCauchyGreen ( dimension, dimension );
    double det_b=0;
    MathUtils<double>::InvertMatrix( LeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    if( dimension == 2 )
    {

        //Almansi Strain Calculation
        rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );

        rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );

        rStrainVector[2] = - InverseLeftCauchyGreen( 0, 1 ); // xy

    }
    else if( dimension == 3 )
    {

        //Almansi Strain Calculation
        if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

        rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );

        rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );

        rStrainVector[2] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 2, 2 ) );

        rStrainVector[3] = - InverseLeftCauchyGreen( 0, 1 ); // xy

        rStrainVector[4] = - InverseLeftCauchyGreen( 1, 2 ); // yz

        rStrainVector[5] = - InverseLeftCauchyGreen( 0, 2 ); // xz

    }
    else
    {
        KRATOS_ERROR << " something is wrong with the dimension almansi strain " << std::endl;
    }


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementElement::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    SolidElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{

    KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number );

    if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR  || rVariable == ALMANSI_STRAIN_VECTOR )
    {
        //create and initialize element variables:
        ElementDataType Variables;
        this->InitializeElementData(Variables,rCurrentProcessInfo);

        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

	    //to take in account previous step writing
	    if( this->Is(SolidElement::FINALIZED_STEP) ){
	      this->GetHistoricalVariables(Variables,PointNumber);
	      noalias(Variables.H) = prod(Variables.F,Variables.F0);
	    }

            //Compute Green-Lagrange Strain
	    if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR )
	      this->CalculateGreenLagrangeStrain( Variables.H, Variables.StrainVector );
            else
	      this->CalculateAlmansiStrain( Variables.H, Variables.StrainVector );


            if ( rOutput[PointNumber].size() != Variables.StrainVector.size() )
	      rOutput[PointNumber].resize( Variables.StrainVector.size(), false );

            rOutput[PointNumber] = Variables.StrainVector;

        }

    }
    else
    {
      SolidElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementElement::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    SolidElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

int LargeDisplacementElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Perform base element checks
    int ErrorCode = 0;
    ErrorCode = SolidElement::Check(rCurrentProcessInfo);

    // Check compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

    // Check that the constitutive law has the correct dimension
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();
    if( dimension == 2 )
    {
      if( LawFeatures.mOptions.IsNot(ConstitutiveLaw::PLANE_STRAIN_LAW) && LawFeatures.mOptions.IsNot(ConstitutiveLaw::PLANE_STRESS_LAW) && LawFeatures.mOptions.IsNot(ConstitutiveLaw::AXISYMMETRIC_LAW) )
	KRATOS_ERROR <<  "wrong constitutive law used. This is a 2D element. Expected plane state or axisymmetric :: element id = " << this->Id() << std::endl;
    }

    return ErrorCode;

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SolidElement )
}

void LargeDisplacementElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SolidElement )
}


} // Namespace Kratos
