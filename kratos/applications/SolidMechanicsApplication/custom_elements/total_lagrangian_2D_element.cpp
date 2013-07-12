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
#include "custom_elements/total_lagrangian_2D_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "solid_mechanics_application.h"


namespace Kratos
{


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  TotalLagrangian2DElement::TotalLagrangian2DElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : TotalLagrangian3DElement( NewId, pGeometry )
  {
    //DO NOT ADD DOFS HERE!!!
  }


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  TotalLagrangian2DElement::TotalLagrangian2DElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : TotalLagrangian3DElement( NewId, pGeometry, pProperties )
  {
  }


  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  TotalLagrangian2DElement::TotalLagrangian2DElement( TotalLagrangian2DElement const& rOther)
    :TotalLagrangian3DElement(rOther)
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  TotalLagrangian2DElement&  TotalLagrangian2DElement::operator=(TotalLagrangian2DElement const& rOther)
  {
    TotalLagrangian3DElement::operator=(rOther);
    return *this;
  }


  //*********************************OPERATIONS*****************************************
  //************************************************************************************

  Element::Pointer TotalLagrangian2DElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
  {
    return Element::Pointer( new TotalLagrangian2DElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
  }


  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  TotalLagrangian2DElement::~TotalLagrangian2DElement()
  {
  }


  //************* GETTING METHODS
  //************************************************************************************
  //************************************************************************************


  void TotalLagrangian2DElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
  {
    rElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
      {
	rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
	rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
      }
  }


  //************************************************************************************
  //************************************************************************************

  void TotalLagrangian2DElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
  {
    int number_of_nodes = GetGeometry().size();
    unsigned int element_size = number_of_nodes * 2;

    if ( rResult.size() != element_size )
      rResult.resize( element_size, false );

    for ( int i = 0; i < number_of_nodes; i++ )
      {
	int index = i * 2;
	rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
	rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
      }

  }


  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************


  void TotalLagrangian2DElement::InitializeStandardVariables (Standard & rVariables, const ProcessInfo& rCurrentProcessInfo)
  {
  
    const unsigned int number_of_nodes = GetGeometry().size();
 
    rVariables.B.resize( 3 , number_of_nodes * 2 );
  
    rVariables.F.resize( 2, 2 );

    rVariables.F0.resize( 2, 2 );
  
    rVariables.ConstitutiveMatrix.resize( 3, 3 );
  
    rVariables.StrainVector.resize( 3 );
  
    rVariables.StressVector.resize( 3 );
  
    rVariables.DN_DX.resize( number_of_nodes, 2 );
  
  }

  //************************************************************************************
  //************************************************************************************

  double& TotalLagrangian2DElement::CalculateIntegrationWeight(double& rIntegrationWeight)
  {
    rIntegrationWeight *= GetProperties()[THICKNESS];
    
    return rIntegrationWeight;
  }



  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************


  void TotalLagrangian2DElement::CalculateGreenLagrangeStrain(const Matrix& rF,
							      Vector& rStrainVector )
  {
    KRATOS_TRY

    //Right Cauchy-Green Calculation
    Matrix C ( 2, 2 );
    noalias( C ) = prod( trans( rF ), rF );

    //Green Lagrange Strain Calculation
    if ( rStrainVector.size() != 3 ) rStrainVector.resize( 3, false );

    rStrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

    rStrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

    rStrainVector[2] = C( 0, 1 ); // xy


    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************

  void TotalLagrangian2DElement::CalculateAlmansiStrain(const Matrix& rF,
						      Vector& rStrainVector )
  {
    KRATOS_TRY

      //Left Cauchy-Green Calculation
      Matrix LeftCauchyGreen = prod( rF, trans( rF ) );

    //Calculating the inverse of the jacobian 
    Matrix InverseLeftCauchyGreen ( 2, 2 );
    double det_b=0;
    MathUtils<double>::InvertMatrix( LeftCauchyGreen, InverseLeftCauchyGreen, det_b);


    //Almansi Strain Calculation
    if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

    rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );

    rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );

    rStrainVector[2] = - InverseLeftCauchyGreen( 0, 1 ); // xy



    KRATOS_CATCH( "" )
      }

  //************************************************************************************
  //************************************************************************************

  void TotalLagrangian2DElement::CalculateDeformationMatrix(Matrix& rB,
							    Matrix& rF,
							    Matrix& rDN_DX)
  {
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
 
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

    KRATOS_CATCH( "" )
      }


  //************************************CALCULATE TOTAL MASS****************************
  //************************************************************************************

  double& TotalLagrangian2DElement::CalculateTotalMass( double& rTotalMass )
  {
    KRATOS_TRY

    rTotalMass = mTotalDomainInitialSize * GetProperties()[DENSITY] * GetProperties()[THICKNESS];

    return rTotalMass;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************
  /**
   * This function provides the place to perform checks on the completeness of the input.
   * It is designed to be called only once (or anyway, not often) typically at the beginning
   * of the calculations, so to verify that nothing is missing from the input
   * or that no common error is found.
   * @param rCurrentProcessInfo
   */
  int  TotalLagrangian2DElement::Check( const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    LargeDisplacement3DElement::Check(rCurrentProcessInfo);

    if ( THICKNESS.Key() == 0 )
      KRATOS_ERROR( std::invalid_argument, "THICKNESS has Key zero! (check if the application is correctly registered", "" );
    
    if ( this->GetProperties().Has( THICKNESS ) == false )
      KRATOS_ERROR( std::logic_error, "THICKNESS not provided for element ", this->Id() );   

    return 0;

    KRATOS_CATCH( "" );
  }


  void TotalLagrangian2DElement::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, TotalLagrangian3DElement );
 }

  void TotalLagrangian2DElement::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, TotalLagrangian3DElement );

  }


} // Namespace Kratos


