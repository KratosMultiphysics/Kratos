//
//   Project Name:        KratosDamApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/small_displacement_element.hpp"
#include "includes/constitutive_law.h"

#include "dam_application_variables.h"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallDisplacementElement::SmallDisplacementElement( )
    :SolidElement( )
{
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallDisplacementElement::SmallDisplacementElement( IndexType NewId, GeometryType::Pointer pGeometry )
    :SolidElement( NewId, pGeometry )
{
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallDisplacementElement::SmallDisplacementElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    :SolidElement( NewId, pGeometry, pProperties )
{
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

SmallDisplacementElement::SmallDisplacementElement( SmallDisplacementElement const& rOther)
    :SolidElement(rOther)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

SmallDisplacementElement&  SmallDisplacementElement::operator=(SmallDisplacementElement const& rOther)
{
    SolidElement::operator=(rOther);

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer SmallDisplacementElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
  return Kratos::make_intrusive< SmallDisplacementElement >(NewId, GetGeometry().Create(rThisNodes), pProperties);
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer SmallDisplacementElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    SmallDisplacementElement NewElement(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

    if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
	NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

	if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
	  KRATOS_ERROR << " constitutive law not has the correct size small displacement element " << std::endl;
      }


    for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
      {
	NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
      }

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Kratos::make_intrusive< SmallDisplacementElement >(NewElement);
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

SmallDisplacementElement::~SmallDisplacementElement()
{
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::SetElementData(ElementDataType& rVariables,
                                              ConstitutiveLaw::Parameters& rValues,
                                              const int & rPointNumber)
{
    KRATOS_TRY

    rValues.SetStrainVector(rVariables.StrainVector);
    rValues.SetStressVector(rVariables.StressVector);
    rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
    rValues.SetShapeFunctionsValues(rVariables.N);

    if(rVariables.detJ<0)
      {
	KRATOS_ERROR << " (small displacement) ELEMENT INVERTED |J|<0 : " << rVariables.detJ << std::endl;
      }

    rValues.SetDeterminantF(rVariables.detF);
    rValues.SetDeformationGradientF(rVariables.F);

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::InitializeElementData (ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    SolidElement::InitializeElementData(rVariables, rCurrentProcessInfo);

    //set variables including all integration points values

    //Calculate Delta Position
    rVariables.DeltaPosition = this->CalculateTotalDeltaPosition(rVariables.DeltaPosition);

    //calculating the reference jacobian from initial cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
						   ElementDataType& rVariables,
						   double& rIntegrationWeight)

{
  // small displacement element is linear no geometric stiffness present
}


//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void SmallDisplacementElement::CalculateKinematics(ElementDataType& rVariables, const double& rPointNumber)
{
    KRATOS_TRY

    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();
    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //Parent to reference configuration
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ;
    MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

    //Compute cartesian derivatives  [dN/dx_n]
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber] , InvJ );

    //Displacement Gradient H  [dU/dx_n]
    this->CalculateDisplacementGradient( rVariables.H, rVariables.DN_DX );

    //Set Shape Functions Values for this integration point
    noalias(rVariables.N) = matrix_row<const Matrix>( Ncontainer, rPointNumber);

    //Compute the deformation matrix B
    const GeometryType& rGeometry = GetGeometry();
    ElementUtilities::CalculateLinearDeformationMatrix(rVariables.B,rGeometry,rVariables.DN_DX);

    //Compute infinitessimal strain
    this->CalculateInfinitesimalStrain( rVariables.H, rVariables.StrainVector );


    KRATOS_CATCH( "" )
}


//*************************COMPUTE DISPLACEMENT GRADIENT******************************
//************************************************************************************

void SmallDisplacementElement::CalculateDisplacementGradient(Matrix& rH, const Matrix& rDN_DX)
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    noalias(rH) = ZeroMatrix(dimension, dimension);

    if( dimension == 2 )
    {

        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {

            const array_1d<double, 3 > & Displacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

            rH ( 0 , 0 ) += Displacement[0]*rDN_DX ( i , 0 );
            rH ( 0 , 1 ) += Displacement[0]*rDN_DX ( i , 1 );
            rH ( 1 , 0 ) += Displacement[1]*rDN_DX ( i , 0 );
            rH ( 1 , 1 ) += Displacement[1]*rDN_DX ( i , 1 );
        }
    }
    else if( dimension == 3 )
    {

        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {

            const array_1d<double, 3 > & Displacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

            rH ( 0 , 0 ) += Displacement[0]*rDN_DX ( i , 0 );
            rH ( 0 , 1 ) += Displacement[0]*rDN_DX ( i , 1 );
            rH ( 0 , 2 ) += Displacement[0]*rDN_DX ( i , 2 );
            rH ( 1 , 0 ) += Displacement[1]*rDN_DX ( i , 0 );
            rH ( 1 , 1 ) += Displacement[1]*rDN_DX ( i , 1 );
            rH ( 1 , 2 ) += Displacement[1]*rDN_DX ( i , 2 );
            rH ( 2 , 0 ) += Displacement[2]*rDN_DX ( i , 0 );
            rH ( 2 , 1 ) += Displacement[2]*rDN_DX ( i , 1 );
            rH ( 2 , 2 ) += Displacement[2]*rDN_DX ( i , 2 );
        }
    }
    else
    {
        KRATOS_ERROR << " something is wrong with the dimension displacement gradient " << std::endl;
    }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::CalculateInfinitesimalStrain(const Matrix& rH, Vector& rStrainVector )
{
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 )
    {

        //Infinitesimal Strain Calculation
        if ( rStrainVector.size() != 3 ) rStrainVector.resize( 3, false );

        rStrainVector[0] = rH( 0, 0 );

        rStrainVector[1] = rH( 1, 1 );

        rStrainVector[2] = (rH( 0, 1 ) + rH( 1, 0 )); // xy

    }
    else if( dimension == 3 )
    {

        //Infinitesimal Strain Calculation
        if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

        rStrainVector[0] = rH( 0, 0 );

        rStrainVector[1] = rH( 1, 1 );

        rStrainVector[2] = rH( 2, 2 );

        rStrainVector[3] = ( rH( 0, 1 ) + rH( 1, 0 ) ); // xy

        rStrainVector[4] = ( rH( 1, 2 ) + rH( 2, 1 ) ); // yz

        rStrainVector[5] = ( rH( 0, 2 ) + rH( 2, 0 ) ); // xz

    }
    else
    {

        KRATOS_ERROR << " something is wrong with the dimension infinitesimal strain " << std::endl;

    }

    KRATOS_CATCH( "" )

}


//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    SolidElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{

    KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if( rOutput.size() != integration_points_number )
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

void SmallDisplacementElement::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    SolidElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

int SmallDisplacementElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Perform base element checks
    int ErrorCode = 0;
    ErrorCode = SolidElement::Check(rCurrentProcessInfo);

    // Check that the element nodes contain all required SolutionStepData and Degrees of freedom
    for(SizeType i=0; i<this->GetGeometry().size(); ++i)
      {
	// Nodal data
	Node<3> &rNode = this->GetGeometry()[i];
	KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rNode);
	//KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUME_ACCELERATION,rNode);

	// Nodal dofs
	KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X,rNode);
	KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y,rNode);
	if( rCurrentProcessInfo[SPACE_DIMENSION] == 3)
	  KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z,rNode);
      }

    // Check compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

    bool correct_strain_measure = false;
    for(unsigned int i=0; i<LawFeatures.mStrainMeasures.size(); i++)
    {
      if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Infinitesimal)
	correct_strain_measure = true;
    }

    if( correct_strain_measure == false )
      KRATOS_ERROR <<  "constitutive law is not compatible with the small displacements element type" << std::endl;

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


void SmallDisplacementElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SolidElement )
}

void SmallDisplacementElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SolidElement )
}


} // Namespace Kratos
