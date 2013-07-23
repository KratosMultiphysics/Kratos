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
#include "custom_elements/updated_lagrangian_U_P_2D_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "solid_mechanics_application.h"


namespace Kratos
{
 

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  UpdatedLagrangianUP2DElement::UpdatedLagrangianUP2DElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : UpdatedLagrangianUP3DElement( NewId, pGeometry )
  {
    //DO NOT ADD DOFS HERE!!!
  }


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  UpdatedLagrangianUP2DElement::UpdatedLagrangianUP2DElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : UpdatedLagrangianUP3DElement( NewId, pGeometry, pProperties )
  {
  }


  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  UpdatedLagrangianUP2DElement::UpdatedLagrangianUP2DElement( UpdatedLagrangianUP2DElement const& rOther)
    : UpdatedLagrangianUP3DElement(rOther)
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  UpdatedLagrangianUP2DElement&  UpdatedLagrangianUP2DElement::operator=(UpdatedLagrangianUP2DElement const& rOther)
  {
    UpdatedLagrangianUP3DElement::operator=(rOther);

    return *this;
  }


  //*********************************OPERATIONS*****************************************
  //************************************************************************************

  Element::Pointer UpdatedLagrangianUP2DElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
  {
    return Element::Pointer( new UpdatedLagrangianUP2DElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
  }


  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  UpdatedLagrangianUP2DElement::~UpdatedLagrangianUP2DElement()
  {
  }

  //************* GETTING METHODS
  //************************************************************************************
  //************************************************************************************


  void UpdatedLagrangianUP2DElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
  {
    rElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
      {
	rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
	rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
	rElementalDofList.push_back( GetGeometry()[i].pGetDof( PRESSURE ));
      }

    
  }


  //************************************************************************************
  //************************************************************************************

  void UpdatedLagrangianUP2DElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
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
	rResult[index + 2] = GetGeometry()[i].GetDof( PRESSURE ).EquationId();
      }

  }

  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************


  void UpdatedLagrangianUP2DElement::InitializeGeneralVariables (GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
  {
  
    const unsigned int number_of_nodes = GetGeometry().size();
 
    rVariables.detF  = 1;

    rVariables.detF0 = 1;

    rVariables.B.resize( 3 , number_of_nodes * 2 );
  
    rVariables.F.resize( 2, 2 );
  
    rVariables.ConstitutiveMatrix.resize( 3, 3 );
  
    rVariables.StrainVector.resize( 3 );
  
    rVariables.StressVector.resize( 3 );
  
    rVariables.ElasticLeftCGVector.resize( 4 );
 
    rVariables.DN_DX.resize( number_of_nodes, 2 );

    //set variables including all integration points values

    //reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));
 
    //reading shape functions local gradients
    rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));
    
    //calculating the jacobian from cartesian coordinates to parent coordinates for all integration points
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod );

    //Calculate Delta Position
    rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);
   
  }

  //************************************************************************************
  //************************************************************************************

  double& UpdatedLagrangianUP2DElement::CalculateIntegrationWeight(double& rIntegrationWeight)
  {
     rIntegrationWeight *= GetProperties()[THICKNESS];

     return rIntegrationWeight;
  }

  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************


  void UpdatedLagrangianUP2DElement::CalculateGreenLagrangeStrain(const Matrix& rF,
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

  void UpdatedLagrangianUP2DElement::CalculateAlmansiStrain(const Matrix& rF,
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
    if ( rStrainVector.size() != 3 ) rStrainVector.resize( 3, false );

    rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );

    rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );

    rStrainVector[2] = - InverseLeftCauchyGreen( 0, 1 ); // xy



    KRATOS_CATCH( "" )
      }

  //*************************COMPUTE DEFORMATION GRADIENT*******************************
  //************************************************************************************

  void UpdatedLagrangianUP2DElement::CalculateDeformationGradient(const Matrix& rDN_DX,
								Matrix& rF,
								Matrix& DeltaPosition)
  {
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    rF = identity_matrix<double> ( 2 );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	rF ( 0 , 0 ) += DeltaPosition(i,0)*rDN_DX ( i , 0 );
	rF ( 0 , 1 ) += DeltaPosition(i,0)*rDN_DX ( i , 1 );
	rF ( 1 , 0 ) += DeltaPosition(i,1)*rDN_DX ( i , 0 );
	rF ( 1 , 1 ) += DeltaPosition(i,1)*rDN_DX ( i , 1 );
      }

      

    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************

  void UpdatedLagrangianUP2DElement::CalculateDeformationMatrix(Matrix& rB,
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

  double& UpdatedLagrangianUP2DElement::CalculateTotalMass( double& rTotalMass )
  {
    KRATOS_TRY

    rTotalMass = GetGeometry().DomainSize() * GetProperties()[DENSITY] * GetProperties()[THICKNESS];

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
  int  UpdatedLagrangianUP2DElement::Check( const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    LargeDisplacementUP3DElement::Check(rCurrentProcessInfo);

    if ( THICKNESS.Key() == 0 )
      KRATOS_ERROR( std::invalid_argument, "THICKNESS has Key zero! (check if the application is correctly registered", "" );
    
    if ( this->GetProperties().Has( THICKNESS ) == false )
      KRATOS_ERROR( std::logic_error, "THICKNESS not provided for element ", this->Id() );   

    return 0;

    KRATOS_CATCH( "" );
  }


  //************************************************************************************
  //************************************************************************************



  void UpdatedLagrangianUP2DElement::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, UpdatedLagrangianUP3DElement );

   }

  void UpdatedLagrangianUP2DElement::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, UpdatedLagrangianUP3DElement );

  }


} // Namespace Kratos


