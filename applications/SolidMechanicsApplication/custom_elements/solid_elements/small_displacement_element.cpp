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
#include "includes/define.h"
#include "custom_elements/solid_elements/small_displacement_element.hpp"
#include "includes/constitutive_law.h"

#include "solid_mechanics_application_variables.h"


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
    return Element::Pointer( new SmallDisplacementElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
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

    return Element::Pointer( new SmallDisplacementElement(NewElement) );
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

SmallDisplacementElement::~SmallDisplacementElement()
{
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::SetElementData(ElementDataPointerType& pVariables,
                                              ConstitutiveLaw::Parameters& rValues,
                                              const int & rPointNumber)
{
    KRATOS_TRY

    rValues.SetStrainVector(pVariables->StrainVector);
    rValues.SetStressVector(pVariables->StressVector);
    rValues.SetConstitutiveMatrix(pVariables->ConstitutiveMatrix);
    rValues.SetShapeFunctionsDerivatives(pVariables->DN_DX);
    rValues.SetShapeFunctionsValues(pVariables->N);

    if(pVariables->detJ<0)
      {
	KRATOS_ERROR << " (small displacement) ELEMENT INVERTED |J|<0 : " << pVariables->detJ << std::endl;
      }

    rValues.SetDeterminantF(pVariables->detF);
    rValues.SetDeformationGradientF(pVariables->F);

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::InitializeElementData (ElementDataPointerType & pVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    SolidElement::InitializeElementData(pVariables, rCurrentProcessInfo);

    //set variables including all integration points values

    //Calculate Delta Position
    pVariables->DeltaPosition = this->CalculateTotalDeltaPosition(pVariables->DeltaPosition);

    //calculating the reference jacobian from initial cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    pVariables->J = GetGeometry().Jacobian( pVariables->J, mThisIntegrationMethod, pVariables->DeltaPosition );

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
						   ElementDataPointerType& pVariables,
						   double& rIntegrationWeight)

{
  // small displacement element is linear no geometric stiffness present
}


//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void SmallDisplacementElement::CalculateKinematics(ElementDataPointerType& pVariables, const double& rPointNumber)
{
    KRATOS_TRY

    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = pVariables->GetShapeFunctionsGradients();
    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = pVariables->GetShapeFunctions();

    //Parent to reference configuration
    pVariables->StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ;
    MathUtils<double>::InvertMatrix( pVariables->J[rPointNumber], InvJ, pVariables->detJ);

    //Compute cartesian derivatives  [dN/dx_n]
    noalias( pVariables->DN_DX ) = prod( DN_De[rPointNumber] , InvJ );

    //Displacement Gradient H  [dU/dx_n]
    this->CalculateDisplacementGradient( pVariables->H, pVariables->DN_DX );

    //Set Shape Functions Values for this integration point
    noalias(pVariables->N) = matrix_row<const Matrix>( Ncontainer, rPointNumber);

    //Compute the deformation matrix B
    this->CalculateDeformationMatrix( pVariables->B, pVariables->DN_DX );

    //Compute infinitessimal strain
    this->CalculateInfinitesimalStrain( pVariables->H, pVariables->StrainVector );


    KRATOS_CATCH( "" )
}


//*************************COMPUTE DISPLACEMENT GRADIENT******************************
//************************************************************************************

void SmallDisplacementElement::CalculateDisplacementGradient(Matrix& rH, const Matrix& rDN_DX)
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType& dimension       = this->Dimension();

    noalias(rH) = ZeroMatrix(dimension, dimension);

    if( dimension == 2 )
    {

        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {

            array_1d<double, 3 > & Displacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

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

            array_1d<double, 3 > & Displacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

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

    const SizeType& dimension = this->Dimension();

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
void SmallDisplacementElement::CalculateDeformationMatrix(Matrix& rB, const Matrix& rDN_DX)
{
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().PointsNumber();
    const SizeType& dimension       = this->Dimension();
    unsigned int voigt_size            = dimension * (dimension +1) * 0.5;

    if ( rB.size1() != voigt_size || rB.size2() != dimension*number_of_nodes )
      rB.resize(voigt_size, dimension*number_of_nodes, false );

    if( dimension == 2 )
    {
        unsigned int index = 0;
        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
            index = 2 * i;

            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 0, index + 1 ) = 0.0;
            rB( 1, index + 0 ) = 0.0;
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rDN_DX( i, 1 );
            rB( 2, index + 1 ) = rDN_DX( i, 0 );

        }

    }
    else if( dimension == 3 )
    {
      unsigned int index = 0;
      for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
	  index = 3 * i;

	  rB( 0, index + 0 ) = rDN_DX( i, 0 );
	  rB( 1, index + 1 ) = rDN_DX( i, 1 );
	  rB( 2, index + 2 ) = rDN_DX( i, 2 );

	  rB( 3, index + 0 ) = rDN_DX( i, 1 );
	  rB( 3, index + 1 ) = rDN_DX( i, 0 );

	  rB( 4, index + 1 ) = rDN_DX( i, 2 );
	  rB( 4, index + 2 ) = rDN_DX( i, 1 );

	  rB( 5, index + 0 ) = rDN_DX( i, 2 );
	  rB( 5, index + 2 ) = rDN_DX( i, 0 );

        }

    }
    else
    {
        KRATOS_ERROR << " something is wrong with the dimension strain matrix " << std::endl;

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
        ElementDataPointerType Variables(make_unique<ElementDataType>());
        this->InitializeElementData(Variables,rCurrentProcessInfo);

        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            if ( rOutput[PointNumber].size() != Variables->StrainVector.size() )
                rOutput[PointNumber].resize( Variables->StrainVector.size(), false );

            rOutput[PointNumber] = Variables->StrainVector;

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
    const SizeType& dimension = this->Dimension();
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
