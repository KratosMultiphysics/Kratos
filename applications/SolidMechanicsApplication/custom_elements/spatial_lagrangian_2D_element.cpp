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
#include "custom_elements/spatial_lagrangian_2D_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "solid_mechanics_application.h"


namespace Kratos
{
 

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  SpatialLagrangian2DElement::SpatialLagrangian2DElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : SpatialLagrangian3DElement( NewId, pGeometry )
  {
    //DO NOT ADD DOFS HERE!!!
  }


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  SpatialLagrangian2DElement::SpatialLagrangian2DElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : SpatialLagrangian3DElement( NewId, pGeometry, pProperties )
  {
    //DO NOT ADD DOFS HERE!!!
  }


  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  SpatialLagrangian2DElement::SpatialLagrangian2DElement( SpatialLagrangian2DElement const& rOther)
    :SpatialLagrangian3DElement(rOther)
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  SpatialLagrangian2DElement&  SpatialLagrangian2DElement::operator=(SpatialLagrangian2DElement const& rOther)
  {
    SpatialLagrangian3DElement::operator=(rOther);

    return *this;
  }


  //*********************************OPERATIONS*****************************************
  //************************************************************************************

  Element::Pointer SpatialLagrangian2DElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
  {
    return Element::Pointer( new SpatialLagrangian2DElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
  }


  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  SpatialLagrangian2DElement::~SpatialLagrangian2DElement()
  {
  }



  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************


  void SpatialLagrangian2DElement::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
  {
  
    const unsigned int number_of_nodes = GetGeometry().size();
 
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
 
    //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points
    rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod, rVariables.DeltaPosition );

  }

  //************************************************************************************
  //************************************************************************************

  double& SpatialLagrangian2DElement::CalculateIntegrationWeight(double& rIntegrationWeight)
  {
     rIntegrationWeight *= GetProperties()[THICKNESS];

     return rIntegrationWeight;
  }



  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************


  //*************************COMPUTE DEFORMATION GRADIENT*******************************
  //************************************************************************************

  void SpatialLagrangian2DElement::CalculateDeformationGradient(const Matrix& rDN_DX,
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

  void SpatialLagrangian2DElement::CalculateDeformationMatrix(Matrix& rB,
							      Matrix& rF,
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

  double& SpatialLagrangian2DElement::CalculateTotalMass( double& rTotalMass )
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
  int  SpatialLagrangian2DElement::Check( const ProcessInfo& rCurrentProcessInfo )
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

  //************************************************************************************
  //************************************************************************************


  void SpatialLagrangian2DElement::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SpatialLagrangian3DElement );
   }

  void SpatialLagrangian2DElement::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SpatialLagrangian3DElement );
  }


} // Namespace Kratos


