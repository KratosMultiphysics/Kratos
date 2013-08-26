//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/total_lagrangian_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "solid_mechanics_application.h"


namespace Kratos
{


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TotalLagrangianElement::TotalLagrangianElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : LargeDisplacementElement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TotalLagrangianElement::TotalLagrangianElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : LargeDisplacementElement( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

TotalLagrangianElement::TotalLagrangianElement( TotalLagrangianElement const& rOther)
    :LargeDisplacementElement(rOther)
    ,mTotalDomainInitialSize(rOther.mTotalDomainInitialSize)
    ,mInvJ0(rOther.mInvJ0)
    ,mDetJ0(rOther.mDetJ0)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

TotalLagrangianElement&  TotalLagrangianElement::operator=(TotalLagrangianElement const& rOther)
{
    LargeDisplacementElement::operator=(rOther);

    mInvJ0.clear();
    mInvJ0.resize( rOther.mInvJ0.size());

    for(unsigned int i=0; i<<mConstitutiveLawVector.size(); i++)
    {
        mInvJ0[i]=rOther.mInvJ0[i];
    }

    mTotalDomainInitialSize = rOther.mTotalDomainInitialSize;
    mDetJ0 = rOther.mDetJ0;

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer TotalLagrangianElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new TotalLagrangianElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer TotalLagrangianElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    TotalLagrangianElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    //-----------//

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

    if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
	NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());
	
	if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
	  KRATOS_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() );
      }
    

    for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
      {
	NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
      }


    //-----------//

    if ( NewElement.mInvJ0.size() != mInvJ0.size() )
      NewElement.mInvJ0.resize(mInvJ0.size());

    for(unsigned int i=0; i<<mInvJ0.size(); i++)
    {
        NewElement.mInvJ0[i] = mInvJ0[i];
    }

    NewElement.mTotalDomainInitialSize = mTotalDomainInitialSize;
    NewElement.mDetJ0 = mDetJ0;

        
    return Element::Pointer( new TotalLagrangianElement(NewElement) );
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

TotalLagrangianElement::~TotalLagrangianElement()
{
}



//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************


void TotalLagrangianElement::Initialize()
{
    KRATOS_TRY

    LargeDisplacementElement::Initialize();

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //Resizing jacobian inverses container
    mInvJ0.resize( integration_points.size() );
    mDetJ0.resize( integration_points.size(), false );


    //Compute jacobian inverses and set the domain initial size:
    GeometryType::JacobiansType J0;
    J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod );
    mTotalDomainInitialSize = 0.00;

    //calculating the inverse J0

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //getting informations for integration
        double IntegrationWeight = integration_points[PointNumber].Weight();

        //calculating and storing inverse of the jacobian and the parameters needed
        MathUtils<double>::InvertMatrix( J0[PointNumber], mInvJ0[PointNumber], mDetJ0[PointNumber] );

        //calculating the total area
        mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;
    }


    KRATOS_CATCH( "" )
}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//************************************************************************************
//************************************************************************************

void TotalLagrangianElement::SetGeneralVariables(GeneralVariables& rVariables,
        ConstitutiveLaw::Parameters& rValues,
        const int & rPointNumber)
{
    LargeDisplacementElement::SetGeneralVariables(rVariables,rValues,rPointNumber);

    //Set extra options for the contitutive law
    Flags &ConstitutiveLawOptions=rValues.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::INITIAL_CONFIGURATION);

}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void TotalLagrangianElement::CalculateKinematics(GeneralVariables& rVariables,
        const double& rPointNumber)

{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //Parent to reference configuration
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;

    //Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], mInvJ0[rPointNumber] );

    //Deformation Gradient F
    noalias( rVariables.F ) = prod( rVariables.J[rPointNumber], mInvJ0[rPointNumber] );

    //Jacobian Determinant for the isoparametric and numerical integration
    rVariables.detJ = mDetJ0[rPointNumber];

    //Determinant of the Deformation Gradient F0
    // (in this element F = F0, then the F0 is set to the identity for coherence in the constitutive law)
    rVariables.detF0 = 1;
    rVariables.F0    = identity_matrix<double> ( dimension );

    //Set Shape Functions Values for this integration point
    rVariables.N=row( Ncontainer, rPointNumber);

    //Compute the deformation matrix B
    this->CalculateDeformationMatrix(rVariables.B,rVariables.F,rVariables.DN_DX);


    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void TotalLagrangianElement::CalculateDeformationMatrix(Matrix& rB,
        Matrix& rF,
        Matrix& rDN_DX)
{
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 2 * i;

            rB( 0, index + 0 ) = rF( 0, 0 ) * rDN_DX( i, 0 );
            rB( 0, index + 1 ) = rF( 1, 0 ) * rDN_DX( i, 0 );
            rB( 1, index + 0 ) = rF( 0, 1 ) * rDN_DX( i, 1 );
            rB( 1, index + 1 ) = rF( 1, 1 ) * rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rF( 0, 0 ) * rDN_DX( i, 1 ) + rF( 0, 1 ) * rDN_DX( i, 0 );
            rB( 2, index + 1 ) = rF( 1, 0 ) * rDN_DX( i, 1 ) + rF( 1, 1 ) * rDN_DX( i, 0 );

        }

    }
    else if( dimension == 3 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 3 * i;

            rB( 0, index + 0 ) = rF( 0, 0 ) * rDN_DX( i, 0 );
            rB( 0, index + 1 ) = rF( 1, 0 ) * rDN_DX( i, 0 );
            rB( 0, index + 2 ) = rF( 2, 0 ) * rDN_DX( i, 0 );
            rB( 1, index + 0 ) = rF( 0, 1 ) * rDN_DX( i, 1 );
            rB( 1, index + 1 ) = rF( 1, 1 ) * rDN_DX( i, 1 );
            rB( 1, index + 2 ) = rF( 2, 1 ) * rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rF( 0, 2 ) * rDN_DX( i, 2 );
            rB( 2, index + 1 ) = rF( 1, 2 ) * rDN_DX( i, 2 );
            rB( 2, index + 2 ) = rF( 2, 2 ) * rDN_DX( i, 2 );
            rB( 3, index + 0 ) = rF( 0, 0 ) * rDN_DX( i, 1 ) + rF( 0, 1 ) * rDN_DX( i, 0 );
            rB( 3, index + 1 ) = rF( 1, 0 ) * rDN_DX( i, 1 ) + rF( 1, 1 ) * rDN_DX( i, 0 );
            rB( 3, index + 2 ) = rF( 2, 0 ) * rDN_DX( i, 1 ) + rF( 2, 1 ) * rDN_DX( i, 0 );
            rB( 4, index + 0 ) = rF( 0, 1 ) * rDN_DX( i, 2 ) + rF( 0, 2 ) * rDN_DX( i, 1 );
            rB( 4, index + 1 ) = rF( 1, 1 ) * rDN_DX( i, 2 ) + rF( 1, 2 ) * rDN_DX( i, 1 );
            rB( 4, index + 2 ) = rF( 2, 1 ) * rDN_DX( i, 2 ) + rF( 2, 2 ) * rDN_DX( i, 1 );
            rB( 5, index + 0 ) = rF( 0, 2 ) * rDN_DX( i, 0 ) + rF( 0, 0 ) * rDN_DX( i, 2 );
            rB( 5, index + 1 ) = rF( 1, 2 ) * rDN_DX( i, 0 ) + rF( 1, 0 ) * rDN_DX( i, 2 );
            rB( 5, index + 2 ) = rF( 2, 2 ) * rDN_DX( i, 0 ) + rF( 2, 0 ) * rDN_DX( i, 2 );

        }

    }
    else
    {

        KRATOS_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

    }

    KRATOS_CATCH( "" )
}


//************************************CALCULATE TOTAL MASS****************************
//************************************************************************************

double& TotalLagrangianElement::CalculateTotalMass( double& rTotalMass )
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rTotalMass = mTotalDomainInitialSize * GetProperties()[DENSITY];

    if( dimension == 2 )
        rTotalMass *= GetProperties()[THICKNESS];

    return rTotalMass;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************



void TotalLagrangianElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement );
    rSerializer.save("mTotalDomainInitialSize",mTotalDomainInitialSize);
    rSerializer.save("InvJ0",mInvJ0);
    rSerializer.save("DetJ0",mDetJ0);
}

void TotalLagrangianElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement );
    rSerializer.load("mTotalDomainInitialSize",mTotalDomainInitialSize);
    rSerializer.load("InvJ0",mInvJ0);
    rSerializer.load("DetJ0",mDetJ0);

}


} // Namespace Kratos


