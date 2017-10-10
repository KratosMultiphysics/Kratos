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
#include "custom_elements/solid_elements/total_lagrangian_element.hpp"
#include "solid_mechanics_application_variables.h"


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

    for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
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
	  KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() )
      }
    

    for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
      {
	NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
      }


    //-----------//

    if ( NewElement.mInvJ0.size() != mInvJ0.size() )
      NewElement.mInvJ0.resize(mInvJ0.size());

    for(unsigned int i=0; i<mInvJ0.size(); i++)
    {
        NewElement.mInvJ0[i] = mInvJ0[i];
    }

    NewElement.mTotalDomainInitialSize = mTotalDomainInitialSize;
    NewElement.mDetJ0 = mDetJ0;

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

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


//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void TotalLagrangianElement::CalculateKinematics(ElementVariables& rVariables,
        const double& rPointNumber)

{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();
    
    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //Parent to reference configuration
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;

    //Jacobian Determinant for the isoparametric and numerical integration
    //
    rVariables.detJ = mDetJ0[rPointNumber];

    //Calculating the cartesian derivatives [dN/dx_n] = [dN/d£][d£/dx_0]
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], mInvJ0[rPointNumber] );

    //Deformation Gradient F [dx_n+1/dx_0] = [dx_n+1/d£] [d£/dx_0]
    noalias( rVariables.F ) = prod( rVariables.j[rPointNumber], mInvJ0[rPointNumber] );

    //Determinant of the deformation gradient F
    rVariables.detF  = MathUtils<double>::Det(rVariables.F);

    //
    //
    //

    //
    //

    //Determinant of the Deformation Gradient F0
    // (in this element F = F0, then F0 is set to the identity for coherence in the constitutive law)
    rVariables.detF0 = 1;
    rVariables.F0    = identity_matrix<double> ( dimension );

    //Set Shape Functions Values for this integration point
    rVariables.N=row( Ncontainer, rPointNumber);

    //Compute the deformation matrix B
    CalculateDeformationMatrix(rVariables.B,rVariables.F,rVariables.DN_DX);

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

    rB.clear(); //set all components to zero

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

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

    }

    KRATOS_CATCH( "" )
}


//************************************CALCULATE TOTAL MASS****************************
//************************************************************************************

double& TotalLagrangianElement::CalculateTotalMass( double& rTotalMass, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rTotalMass = mTotalDomainInitialSize * GetProperties()[DENSITY];

    if( dimension == 2 ){
      if ( this->GetProperties().Has( THICKNESS ) )
	rTotalMass *= GetProperties()[THICKNESS];
    }

    return rTotalMass;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

    
void TotalLagrangianElement::GetHistoricalVariables( ElementVariables& rVariables, const double& rPointNumber )
{

}


//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

double& TotalLagrangianElement::CalculateVolumeChange( double& rVolumeChange, ElementVariables& rVariables )
{
    KRATOS_TRY
      
    rVolumeChange = 1.0;

    return rVolumeChange;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

int TotalLagrangianElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Perform base element checks
    int ErrorCode = 0;
    ErrorCode = LargeDisplacementElement::Check(rCurrentProcessInfo);
    
    return ErrorCode;

    KRATOS_CATCH( "" );
}
  
//************************************************************************************
//************************************************************************************


void TotalLagrangianElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement )
    rSerializer.save("mTotalDomainInitialSize",mTotalDomainInitialSize);
    rSerializer.save("InvJ0",mInvJ0);
    rSerializer.save("DetJ0",mDetJ0);
}

void TotalLagrangianElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement )
    rSerializer.load("mTotalDomainInitialSize",mTotalDomainInitialSize);
    rSerializer.load("InvJ0",mInvJ0);
    rSerializer.load("DetJ0",mDetJ0);

}


} // Namespace Kratos


