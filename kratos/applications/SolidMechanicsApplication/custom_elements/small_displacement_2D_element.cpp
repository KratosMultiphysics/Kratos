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
#include "custom_elements/small_displacement_2D_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "solid_mechanics_application.h"


namespace Kratos
{
 

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  SmallDisplacement2DElement::SmallDisplacement2DElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : SmallDisplacement3DElement( NewId, pGeometry )
  {
    //DO NOT ADD DOFS HERE!!!
  }


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  SmallDisplacement2DElement::SmallDisplacement2DElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : SmallDisplacement3DElement( NewId, pGeometry, pProperties )
  {
  }


  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  SmallDisplacement2DElement::SmallDisplacement2DElement( SmallDisplacement2DElement const& rOther)
    :SmallDisplacement3DElement(rOther)
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  SmallDisplacement2DElement&  SmallDisplacement2DElement::operator=(SmallDisplacement2DElement const& rOther)
  {
    SmallDisplacement3DElement::operator=(rOther);

    return *this;
  }


  //*********************************OPERATIONS*****************************************
  //************************************************************************************

  Element::Pointer SmallDisplacement2DElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
  {
    return Element::Pointer( new SmallDisplacement2DElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
  }


  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  SmallDisplacement2DElement::~SmallDisplacement2DElement()
  {
  }


  //************* GETTING METHODS
  //************************************************************************************
  //************************************************************************************


  void SmallDisplacement2DElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
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

  void SmallDisplacement2DElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
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

  void SmallDisplacement2DElement::InitializeStandardVariables (Standard & rVariables, const ProcessInfo& rCurrentProcessInfo)
  {

    const unsigned int number_of_nodes = GetGeometry().size();
 
    rVariables.B.resize( 3 , number_of_nodes * 2 );

    rVariables.H.resize( 2, 2 );

    rVariables.ConstitutiveMatrix.resize( 3, 3 );
  
    rVariables.StrainVector.resize( 3 );
  
    rVariables.StressVector.resize( 3 );

    rVariables.DN_DX.resize( number_of_nodes, 2 );
 
    //needed parameters for consistency with the general constitutive law: small displacements 
    rVariables.detF0 = 1;
    rVariables.detF  = 1;
    rVariables.F0    = identity_matrix<double>(2);
    rVariables.F     = identity_matrix<double>(2);

    //set variables including all integration points values

    //reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));
 
    //reading shape functions local gradients
    rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));
    
    //calculating the jacobian from cartesian coordinates to parent coordinates for all integration points
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod );


  }

  //************************************************************************************
  //************************************************************************************

  double& SmallDisplacement2DElement::CalculateIntegrationWeight(double& rIntegrationWeight)
  {

    rIntegrationWeight *= GetProperties()[THICKNESS];

    return rIntegrationWeight;
  }



  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************


  //*************************COMPUTE DISPLACEMENT GRADIENT******************************
  //************************************************************************************

  void SmallDisplacement2DElement::CalculateDisplacementGradient(const Matrix& rDN_DX,
								 Matrix& rH)
  {
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    rH = zero_matrix<double> ( 2 );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {

	array_1d<double, 3 > & Displacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

	rH ( 0 , 0 ) += Displacement[0]*rDN_DX ( i , 0 );
	rH ( 0 , 1 ) += Displacement[0]*rDN_DX ( i , 1 );
	rH ( 1 , 0 ) += Displacement[1]*rDN_DX ( i , 0 );
	rH ( 1 , 1 ) += Displacement[1]*rDN_DX ( i , 1 );
     }


    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************

  void SmallDisplacement2DElement::CalculateInfinitesimalStrain(const Matrix& rH,
								 Vector& rStrainVector )
  {
    KRATOS_TRY


      //Green Lagrange Strain Calculation
      if ( rStrainVector.size() != 3 ) rStrainVector.resize( 3, false );

    rStrainVector[0] = rH( 0, 0 );

    rStrainVector[1] = rH( 1, 1 );

    rStrainVector[2] = (rH( 0, 1 ) + rH( 1, 0 )); // xy

    KRATOS_CATCH( "" )

      }


  //************************************************************************************
  //************************************************************************************

  void SmallDisplacement2DElement::CalculateDeformationMatrix(Matrix& rB,
							      Matrix& rDN_DX)
  {
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
 
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	unsigned int index = 2 * i;

	rB( 0, index + 0 ) = rDN_DX( i, 0 );
	rB( 1, index + 1 ) = rDN_DX( i, 1 );
	rB( 2, index + 0 ) = rDN_DX( i, 1 );
	rB( 2, index + 1 ) = rDN_DX( i, 0 );

      }

    KRATOS_CATCH( "" )
      }


  //************************************CALCULATE TOTAL MASS****************************
  //************************************************************************************

  double& SmallDisplacement2DElement::CalculateTotalMass( double& TotalMass )
  {
    KRATOS_TRY

    TotalMass = GetGeometry().DomainSize() * GetProperties()[DENSITY] * GetProperties()[THICKNESS];

    return TotalMass;

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
  int  SmallDisplacement2DElement::Check( const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    SmallDisplacement3DElement::Check(rCurrentProcessInfo);

    if ( THICKNESS.Key() == 0 )
      KRATOS_ERROR( std::invalid_argument, "THICKNESS has Key zero! (check if the application is correctly registered", "" );
    
    if ( this->GetProperties().Has( THICKNESS ) == false )
      KRATOS_ERROR( std::logic_error, "THICKNESS not provided for element ", this->Id() );   

    return 0;

    KRATOS_CATCH( "" );
  }


  void SmallDisplacement2DElement::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer,  SmallDisplacement3DElement );
   }

  void SmallDisplacement2DElement::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SmallDisplacement3DElement );
  }


} // Namespace Kratos


