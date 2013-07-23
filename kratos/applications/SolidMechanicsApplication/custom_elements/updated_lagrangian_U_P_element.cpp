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
#include "custom_elements/updated_lagrangian_U_P_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "solid_mechanics_application.h"


namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  UpdatedLagrangianUPElement::UpdatedLagrangianUPElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : SpatialLagrangianUPElement( NewId, pGeometry )
  {
    //DO NOT ADD DOFS HERE!!!
  }


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  UpdatedLagrangianUPElement::UpdatedLagrangianUPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : SpatialLagrangianUPElement( NewId, pGeometry, pProperties )
  {
  }


  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  UpdatedLagrangianUPElement::UpdatedLagrangianUPElement( UpdatedLagrangianUPElement const& rOther)
    :SpatialLagrangianUPElement(rOther)
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  UpdatedLagrangianUPElement&  UpdatedLagrangianUPElement::operator=(UpdatedLagrangianUPElement const& rOther)
  {
    SpatialLagrangianUPElement::operator=(rOther);

    return *this;
  }


  //*********************************OPERATIONS*****************************************
  //************************************************************************************

  Element::Pointer UpdatedLagrangianUPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
  {
    return Element::Pointer( new UpdatedLagrangianUPElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
  }


  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  UpdatedLagrangianUPElement::~UpdatedLagrangianUPElement()
  {
  }


  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************
  

  void UpdatedLagrangianUPElement::InitializeGeneralVariables (GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
  {
    LargeDisplacementUPElement::InitializeGeneralVariables(rVariables,rCurrentProcessInfo);

    //Calculate Delta Position
    rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);
 
  }


  //************************************************************************************
  //************************************************************************************

  void UpdatedLagrangianUPElement::SetGeneralVariables (GeneralVariables& rVariables,
							 ConstitutiveLaw::Parameters& rValues,
							 const int & rPointNumber)
  {
    LargeDisplacementUPElement::SetGeneralVariables(rVariables,rValues,rPointNumber);

    //Set extra options for the contitutive law
    Flags &ConstitutiveLawOptions=rValues.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::LAST_KNOWN_CONFIGURATION);

  }

  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************


  void UpdatedLagrangianUPElement::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, GeneralVariables& rVariables, double& rIntegrationWeight)
  {

    //contributions to stiffness matrix calculated on the reference config

    CalculatePushForwardDN_DX( rVariables ); //to be compatible with the updated lagrangian configuration

    LargeDisplacementUPElement::CalculateAndAddLHS( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    //KRATOS_WATCH(rLeftHandSideMatrix)
  }


  //************************************************************************************
  //************************************************************************************

  void UpdatedLagrangianUPElement::CalculateAndAddRHS(VectorType& rRightHandSideVector, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
  {

    //contribution to external forces

    LargeDisplacementUPElement::CalculateAndAddRHS( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );
    
    //KRATOS_WATCH(rRightHandSideVector)
  }


  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************


  void UpdatedLagrangianUPElement::CalculateKinematics(GeneralVariables& rVariables,
						     const double& rPointNumber)

  {
    KRATOS_TRY
      
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //Parent to reference configuration
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;
    
    //Calculating the inverse of the jacobian and the parameters needed
    Matrix InvJ;
    MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

    //Compute cartesian derivatives
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], InvJ );

    //Current Deformation Gradient
    this->CalculateDeformationGradient (rVariables.DN_DX, rVariables.F, rVariables.DeltaPosition);

    //Determinant of the Deformation Gradient F0
    rVariables.detF0 = mDeterminantF0[rPointNumber];
    rVariables.F0    = mDeformationGradientF0[rPointNumber];

    //Set Shape Functions Values for this integration point
    rVariables.N=row( Ncontainer, rPointNumber);

    //Compute the deformation matrix B
    this->CalculateDeformationMatrix(rVariables.B, rVariables.F, rVariables.DN_DX);


    KRATOS_CATCH( "" )
      }



  //************************************************************************************
  //************************************************************************************


  void UpdatedLagrangianUPElement::CalculateDeformationMatrix(Matrix& rB,
								   Matrix& rF,
								   Matrix& rDN_DX)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 ){

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
    else if( dimension == 3 ){
     
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
    else{

      KRATOS_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

    }

    KRATOS_CATCH( "" )
      }
  


  //************************************************************************************
  //************************************************************************************

  void UpdatedLagrangianUPElement::CalculatePushForwardDN_DX(GeneralVariables& rVariables)
  {
    //SET DN_DX WITH THE NON LINEAR PART  = F^-1 * DN_DX

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //Inverse of Deformation Gradient:
    Matrix InverseF ( dimension , dimension );
    double detF=0;
    MathUtils<double>::InvertMatrix( trans(rVariables.F), InverseF, detF);
    
    Matrix DN_DX_new( number_of_nodes, dimension  );
    DN_DX_new.clear();

    for ( unsigned int l= 0; l< number_of_nodes; l++ )
      {  
	for ( unsigned int m = 0; m< dimension; m++ )
	  {		  
	    for ( unsigned int n = 0; n< dimension; n++ )
	      {
		DN_DX_new ( l , m ) += rVariables.DN_DX ( l , n ) * InverseF ( m , n );
	      }
	  }
      }

    rVariables.DN_DX = DN_DX_new;

  }

  //************************************************************************************
  //************************************************************************************


  void UpdatedLagrangianUPElement::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SpatialLagrangianUPElement );
   }

  void UpdatedLagrangianUPElement::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SpatialLagrangianUPElement );
  }


} // Namespace Kratos


