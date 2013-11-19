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
#include "custom_elements/spatial_lagrangian_U_P_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "solid_mechanics_application.h"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SpatialLagrangianUPElement::SpatialLagrangianUPElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : LargeDisplacementUPElement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SpatialLagrangianUPElement::SpatialLagrangianUPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : LargeDisplacementUPElement( NewId, pGeometry, pProperties )
{
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

SpatialLagrangianUPElement::SpatialLagrangianUPElement( SpatialLagrangianUPElement const& rOther)
    :LargeDisplacementUPElement(rOther)
    ,mDeformationGradientF0(rOther.mDeformationGradientF0)
    ,mDeterminantF0(rOther.mDeterminantF0)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

SpatialLagrangianUPElement&  SpatialLagrangianUPElement::operator=(SpatialLagrangianUPElement const& rOther)
{
    LargeDisplacementUPElement::operator=(rOther);

    mDeformationGradientF0.clear();
    mDeformationGradientF0.resize(rOther.mDeformationGradientF0.size());

    for(unsigned int i=0; i<<mConstitutiveLawVector.size(); i++)
    {
        mDeformationGradientF0[i] = rOther.mDeformationGradientF0[i];
    }

    mDeterminantF0 = rOther.mDeterminantF0;

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer SpatialLagrangianUPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new SpatialLagrangianUPElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer SpatialLagrangianUPElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    SpatialLagrangianUPElement NewElement(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    //-----------//

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

    if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
	NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());
	
	if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
	  KRATOS_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() )
      }
    

    for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
      {
	NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
      }


    //-----------//

    if ( NewElement.mDeformationGradientF0.size() != mDeformationGradientF0.size() )
      NewElement.mDeformationGradientF0.resize(mDeformationGradientF0.size());

    for(unsigned int i=0; i<<mDeformationGradientF0.size(); i++)
    {
        NewElement.mDeformationGradientF0[i] = mDeformationGradientF0[i];
    }

    NewElement.mDeterminantF0 = mDeterminantF0;

        
    return Element::Pointer( new SpatialLagrangianUPElement(NewElement) );
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

SpatialLagrangianUPElement::~SpatialLagrangianUPElement()
{
}

//************************************************************************************
//************************************************************************************


//*********************************SET DOUBLE VALUE***********************************
//************************************************************************************

void SpatialLagrangianUPElement::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{

  if (rVariable == DETERMINANT_F){

    const unsigned int& integration_points_number = mConstitutiveLawVector.size();

    
    for ( unsigned int PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
      {
	mDeterminantF0[PointNumber] = rValues[PointNumber];
      }

  }
  else{

    LargeDisplacementUPElement::SetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

  }


}


//************************************************************************************
//************************************************************************************

//**********************************GET DOUBLE VALUE**********************************
//************************************************************************************


void SpatialLagrangianUPElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{

  if (rVariable == DETERMINANT_F){

    const unsigned int& integration_points_number = mConstitutiveLawVector.size();

    if ( rValues.size() != integration_points_number )
      rValues.resize( integration_points_number );
    
    for ( unsigned int PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
      {
	rValues[PointNumber] = mDeterminantF0[PointNumber];
      }

  }
  else{

    LargeDisplacementUPElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

  }

}

//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void SpatialLagrangianUPElement::Initialize()
{
    KRATOS_TRY

    LargeDisplacementUPElement::Initialize();

    SizeType integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //Resize historic deformation gradient
    if ( mDeformationGradientF0.size() != integration_points_number )
        mDeformationGradientF0.resize( integration_points_number );

    if ( mDeterminantF0.size() != integration_points_number )
        mDeterminantF0.resize( integration_points_number, false );

    for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
    {
        mDeterminantF0[PointNumber] = 1;
        mDeformationGradientF0[PointNumber] = identity_matrix<double> (dimension);
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void SpatialLagrangianUPElement::InitializeGeneralVariables (GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    LargeDisplacementElement::InitializeGeneralVariables(rVariables,rCurrentProcessInfo);

    //Calculate Delta Position
    rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);

    //set variables including all integration points values

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );


}


////************************************************************************************
////************************************************************************************

void SpatialLagrangianUPElement::FinalizeStepVariables( GeneralVariables & rVariables, const double& rPointNumber )
{ 
    //update internal (historical) variables
    mDeterminantF0[rPointNumber]         = rVariables.detF0 ;
    mDeformationGradientF0[rPointNumber] = rVariables.F0;
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

void SpatialLagrangianUPElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
{

    //contributions to stiffness matrix calculated on the reference config
    rVariables.detF0   *= rVariables.detF;
    double DeterminantF = rVariables.detF;
    rVariables.detF = 1; //in order to simplify updated and spatial lagrangian

    LargeDisplacementUPElement::CalculateAndAddLHS( rLocalSystem, rVariables, rIntegrationWeight );

    rVariables.detF     = DeterminantF;
    rVariables.detF0   /= rVariables.detF;
    //KRATOS_WATCH( rLeftHandSideMatrix )
}


//************************************************************************************
//************************************************************************************

void SpatialLagrangianUPElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
{
    //contribution to external forces
    rVariables.detF0   *= rVariables.detF;
    double DeterminantF = rVariables.detF;
    rVariables.detF = 1; //in order to simplify updated and spatial lagrangian

    LargeDisplacementUPElement::CalculateAndAddRHS( rLocalSystem, rVariables, rVolumeForce, rIntegrationWeight );

    rVariables.detF     = DeterminantF;
    rVariables.detF0   /= rVariables.detF;
    //KRATOS_WATCH( rRightHandSideVector )
}


//************************************************************************************
//************************************************************************************

void SpatialLagrangianUPElement::SetGeneralVariables(GeneralVariables& rVariables,
        ConstitutiveLaw::Parameters& rValues,
        const int & rPointNumber)
{
    LargeDisplacementElement::SetGeneralVariables(rVariables,rValues,rPointNumber);

    //Set extra options for the contitutive law
    Flags &ConstitutiveLawOptions=rValues.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::FINAL_CONFIGURATION);

}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void SpatialLagrangianUPElement::CalculateKinematics(GeneralVariables& rVariables,
        const double& rPointNumber)

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

    //Step domain size
    rVariables.DomainSize = rVariables.detJ;

    //Compute cartesian derivatives [dN/dx_n]
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], InvJ );

    //Current Deformation Gradient [dx_n+1/dx_n]
    //this->CalculateDeformationGradient (rVariables.DN_DX, rVariables.F, rVariables.DeltaPosition);

    //Deformation Gradient F [dx_n+1/dx_n] to be updated
    noalias( rVariables.F ) = prod( rVariables.j[rPointNumber], InvJ );

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n+1]
    Matrix Invj;
    MathUtils<double>::InvertMatrix( rVariables.j[rPointNumber], Invj, rVariables.detJ); //overwrites detJ

    //Compute cartesian derivatives [dN/dx_n+1]
    rVariables.DN_DX = prod( DN_De[rPointNumber], Invj ); //overwrites DX now is the current position dx

    //Determinant of the Deformation Gradient F0
    rVariables.detF0 = mDeterminantF0[rPointNumber];
    rVariables.F0    = mDeformationGradientF0[rPointNumber];

    //Set Shape Functions Values for this integration point
    rVariables.N=row( Ncontainer, rPointNumber);

    //Compute the deformation matrix B
    this->CalculateDeformationMatrix(rVariables.B, rVariables.F, rVariables.DN_DX);


    KRATOS_CATCH( "" )
}


//*************************COMPUTE DEFORMATION GRADIENT*******************************
//************************************************************************************

void SpatialLagrangianUPElement::CalculateDeformationGradient(const Matrix& rDN_DX,
        Matrix& rF,
        Matrix& rDeltaPosition)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rF = identity_matrix<double> ( dimension );

    if( dimension == 2 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            rF ( 0 , 0 ) += rDeltaPosition(i,0)*rDN_DX ( i , 0 );
            rF ( 0 , 1 ) += rDeltaPosition(i,0)*rDN_DX ( i , 1 );
            rF ( 1 , 0 ) += rDeltaPosition(i,1)*rDN_DX ( i , 0 );
            rF ( 1 , 1 ) += rDeltaPosition(i,1)*rDN_DX ( i , 1 );
        }

    }
    else if( dimension == 3)
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {

            rF ( 0 , 0 ) += rDeltaPosition(i,0)*rDN_DX ( i , 0 );
            rF ( 0 , 1 ) += rDeltaPosition(i,0)*rDN_DX ( i , 1 );
            rF ( 0 , 2 ) += rDeltaPosition(i,0)*rDN_DX ( i , 2 );
            rF ( 1 , 0 ) += rDeltaPosition(i,1)*rDN_DX ( i , 0 );
            rF ( 1 , 1 ) += rDeltaPosition(i,1)*rDN_DX ( i , 1 );
            rF ( 1 , 2 ) += rDeltaPosition(i,1)*rDN_DX ( i , 2 );
            rF ( 2 , 0 ) += rDeltaPosition(i,2)*rDN_DX ( i , 0 );
            rF ( 2 , 1 ) += rDeltaPosition(i,2)*rDN_DX ( i , 1 );
            rF ( 2 , 2 ) += rDeltaPosition(i,2)*rDN_DX ( i , 2 );
        }

    }
    else
    {

        KRATOS_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

    }

    KRATOS_CATCH( "" )
}




//************************************************************************************
//************************************************************************************


void SpatialLagrangianUPElement::CalculateDeformationMatrix(Matrix& rB,
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

            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rDN_DX( i, 1 );
            rB( 2, index + 1 ) = rDN_DX( i, 0 );

        }

    }
    else if( dimension == 3 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 3 * i;

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

        KRATOS_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

    }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************


void SpatialLagrangianUPElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementUPElement )
    rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.save("DeterminantF0",mDeterminantF0);
}

void SpatialLagrangianUPElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementUPElement )
    rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.load("DeterminantF0",mDeterminantF0);
}


} // Namespace Kratos


