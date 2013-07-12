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
#include "custom_elements/total_lagrangian_3D_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "solid_mechanics_application.h"


namespace Kratos
{
 

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  TotalLagrangian3DElement::TotalLagrangian3DElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : LargeDisplacement3DElement( NewId, pGeometry )
  {
    //DO NOT ADD DOFS HERE!!!
  }


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  TotalLagrangian3DElement::TotalLagrangian3DElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : LargeDisplacement3DElement( NewId, pGeometry, pProperties )
  {
    //DO NOT ADD DOFS HERE!!!
  }


  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  TotalLagrangian3DElement::TotalLagrangian3DElement( TotalLagrangian3DElement const& rOther)
    :LargeDisplacement3DElement(rOther)
    ,mTotalDomainInitialSize(rOther.mTotalDomainInitialSize)
    ,mInvJ0(rOther.mInvJ0)
    ,mDetJ0(rOther.mDetJ0)
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  TotalLagrangian3DElement&  TotalLagrangian3DElement::operator=(TotalLagrangian3DElement const& rOther)
  {
    LargeDisplacement3DElement::operator=(rOther);

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

  Element::Pointer TotalLagrangian3DElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
  {
    return Element::Pointer( new TotalLagrangian3DElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
  }


  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  TotalLagrangian3DElement::~TotalLagrangian3DElement()
  {
  }



  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************


  void TotalLagrangian3DElement::Initialize()
  {
    KRATOS_TRY

    LargeDisplacement3DElement::Initialize();

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



  //************************************************************************************
  //************************************************************************************

  double& TotalLagrangian3DElement::CalculateIntegrationWeight(double& rIntegrationWeight)
  {
     return rIntegrationWeight;
  }



  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************


  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************


  void TotalLagrangian3DElement::CalculateKinematics(Standard& rVariables,
						     const double& rPointNumber)

  {
    KRATOS_TRY
      
    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //Parent to reference configuration
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;

    // Parent to reference configuration
    Matrix J ( dimension , dimension);
    J = GetGeometry().Jacobian( J, rPointNumber , mThisIntegrationMethod );
    
    //Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], mInvJ0[rPointNumber] );

    //Deformation Gradient F
    noalias( rVariables.F ) = prod( J, mInvJ0[rPointNumber] );

    //Jacobian Determinant for the isoparametric and numerical integration
    rVariables.detJ = mDetJ0[rPointNumber];

    //Determinant of the Deformation Gradient F0
    // (in this element F = F0, then the F0 is set to the identity for coherence in the constitutive law)
    rVariables.detF0 = 1;
    rVariables.F0    = identity_matrix<double>(rVariables.F0.size1());

    //Compute the deformation matrix B
    this->CalculateDeformationMatrix(rVariables.B,rVariables.F,rVariables.DN_DX);


    KRATOS_CATCH( "" )
      }



  //************************************************************************************
  //************************************************************************************

  void TotalLagrangian3DElement::CalculateDeformationMatrix(Matrix& rB,
							    Matrix& rF,
							    Matrix& rDN_DX)
  {
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
 
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	unsigned int index = dimension * i;

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

    KRATOS_CATCH( "" )
      }


  //************************************CALCULATE TOTAL MASS****************************
  //************************************************************************************

  double& TotalLagrangian3DElement::CalculateTotalMass( double& rTotalMass )
  {
    KRATOS_TRY

    rTotalMass = mTotalDomainInitialSize * GetProperties()[DENSITY];

    return rTotalMass;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************



  void TotalLagrangian3DElement::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacement3DElement );
    rSerializer.save("mTotalDomainInitialSize",mTotalDomainInitialSize);
    rSerializer.save("InvJ0",mInvJ0);
    rSerializer.save("DetJ0",mDetJ0);
   }

  void TotalLagrangian3DElement::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacement3DElement );
    rSerializer.load("mTotalDomainInitialSize",mTotalDomainInitialSize);
    rSerializer.load("InvJ0",mInvJ0);
    rSerializer.load("DetJ0",mDetJ0);

  }


} // Namespace Kratos


