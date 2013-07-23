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
#include "custom_elements/updated_lagrangian_U_P_3D_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "solid_mechanics_application.h"


namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  UpdatedLagrangianUP3DElement::UpdatedLagrangianUP3DElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : SpatialLagrangianUP3DElement( NewId, pGeometry )
  {
    //DO NOT ADD DOFS HERE!!!
  }


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  UpdatedLagrangianUP3DElement::UpdatedLagrangianUP3DElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : SpatialLagrangianUP3DElement( NewId, pGeometry, pProperties )
  {
  }


  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  UpdatedLagrangianUP3DElement::UpdatedLagrangianUP3DElement( UpdatedLagrangianUP3DElement const& rOther)
    :SpatialLagrangianUP3DElement(rOther)
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  UpdatedLagrangianUP3DElement&  UpdatedLagrangianUP3DElement::operator=(UpdatedLagrangianUP3DElement const& rOther)
  {
    SpatialLagrangianUP3DElement::operator=(rOther);

    return *this;
  }


  //*********************************OPERATIONS*****************************************
  //************************************************************************************

  Element::Pointer UpdatedLagrangianUP3DElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
  {
    return Element::Pointer( new UpdatedLagrangianUP3DElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
  }


  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  UpdatedLagrangianUP3DElement::~UpdatedLagrangianUP3DElement()
  {
  }


  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************
  

  void UpdatedLagrangianUP3DElement::InitializeGeneralVariables (GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
  {
    LargeDisplacement3DElement::InitializeGeneralVariables(rVariables,rCurrentProcessInfo);

    //Calculate Delta Position
    rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);
 
  }


  //************************************************************************************
  //************************************************************************************

  void UpdatedLagrangianUP3DElement::SetGeneralVariables (GeneralVariables& rVariables,
							 ConstitutiveLaw::Parameters& rValues,
							 const int & rPointNumber)
  {
    LargeDisplacement3DElement::SetGeneralVariables(rVariables,rValues,rPointNumber);

    //Set extra options for the contitutive law
    Flags &ConstitutiveLawOptions=rValues.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::LAST_KNOWN_CONFIGURATION);

  }

  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************


  void UpdatedLagrangianUP3DElement::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, GeneralVariables& rVariables, double& rIntegrationWeight)
  {

    //contributions to stiffness matrix calculated on the reference config

    CalculatePushForwardDN_DX( rVariables ); //to be compatible with the updated lagrangian configuration

    LargeDisplacementUP3DElement::CalculateAndAddLHS( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    //KRATOS_WATCH(rLeftHandSideMatrix)
  }


  //************************************************************************************
  //************************************************************************************

  void UpdatedLagrangianUP3DElement::CalculateAndAddRHS(VectorType& rRightHandSideVector, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
  {

    //contribution to external forces

    LargeDisplacementUP3DElement::CalculateAndAddRHS( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );
    
    //KRATOS_WATCH(rRightHandSideVector)
  }


  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************


  void UpdatedLagrangianUP3DElement::CalculateKinematics(GeneralVariables& rVariables,
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
    rVariables.ElasticLeftCGVector = mElasticLeftCauchyGreenVector[rPointNumber];

    //Set Shape Functions Values for this integration point
    rVariables.N=row( Ncontainer, rPointNumber);

    //Compute the deformation matrix B
    this->CalculateDeformationMatrix(rVariables.B, rVariables.F, rVariables.DN_DX);


    KRATOS_CATCH( "" )
      }



  //************************************************************************************
  //************************************************************************************


  void UpdatedLagrangianUP3DElement::CalculateDeformationMatrix(Matrix& rB,
								   Matrix& rF,
								   Matrix& rDN_DX)
  {
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
 
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

    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************

  void UpdatedLagrangianUP3DElement::CalculatePushForwardDN_DX(GeneralVariables& rVariables)
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


  void UpdatedLagrangianUP3DElement::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SpatialLagrangianUP3DElement );
   }

  void UpdatedLagrangianUP3DElement::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SpatialLagrangianUP3DElement );
  }


} // Namespace Kratos


