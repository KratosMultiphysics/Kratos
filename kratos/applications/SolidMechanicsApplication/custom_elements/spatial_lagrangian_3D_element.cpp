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
#include "custom_elements/spatial_lagrangian_3D_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "solid_mechanics_application.h"


namespace Kratos
{
 

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  SpatialLagrangian3DElement::SpatialLagrangian3DElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : LargeDisplacement3DElement( NewId, pGeometry )
  {
    //DO NOT ADD DOFS HERE!!!
  }


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  SpatialLagrangian3DElement::SpatialLagrangian3DElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : LargeDisplacement3DElement( NewId, pGeometry, pProperties )
  {
    //DO NOT ADD DOFS HERE!!!
  }


  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  SpatialLagrangian3DElement::SpatialLagrangian3DElement( SpatialLagrangian3DElement const& rOther)
    :LargeDisplacement3DElement(rOther)
    ,mElasticLeftCauchyGreenVector(rOther.mElasticLeftCauchyGreenVector)
    ,mDeterminantF0(rOther.mDeterminantF0)
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  SpatialLagrangian3DElement&  SpatialLagrangian3DElement::operator=(SpatialLagrangian3DElement const& rOther)
  {
    LargeDisplacement3DElement::operator=(rOther);

    mElasticLeftCauchyGreenVector.clear();
    mElasticLeftCauchyGreenVector.resize(rOther.mElasticLeftCauchyGreenVector.size());

    for(unsigned int i=0; i<<mConstitutiveLawVector.size(); i++)
      {
	mElasticLeftCauchyGreenVector [i] = rOther.mElasticLeftCauchyGreenVector[i];
      }

    mDeterminantF0 = rOther.mDeterminantF0;


    return *this;
  }


  //*********************************OPERATIONS*****************************************
  //************************************************************************************

  Element::Pointer SpatialLagrangian3DElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
  {
    return Element::Pointer( new SpatialLagrangian3DElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
  }


  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  SpatialLagrangian3DElement::~SpatialLagrangian3DElement()
  {
  }



  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************


  void SpatialLagrangian3DElement::Initialize()
  {
    KRATOS_TRY

    LargeDisplacement3DElement::Initialize();

    SizeType integration_points_number=GetGeometry().IntegrationPointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //Resize historic deformation gradient
    mElasticLeftCauchyGreenVector.resize( integration_points_number );
    mDeterminantF0.resize( integration_points_number, false );
    
    for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
      {
	mDeterminantF0[PointNumber] = 1;
	mElasticLeftCauchyGreenVector[PointNumber] = ZeroVector(dimension * 2);
	for( unsigned int i=0; i<3; i++)
	  mElasticLeftCauchyGreenVector[PointNumber][i] = 1;
	
      }

    KRATOS_CATCH( "" )
  }


 //************************************************************************************
  //************************************************************************************

  void SpatialLagrangian3DElement::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
  {
    LargeDisplacement3DElement::InitializeGeneralVariables(rVariables,rCurrentProcessInfo);

    //Calculate Delta Position
    rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);
 
    //set variables including all integration points values
   
    //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points
    rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod, rVariables.DeltaPosition );


  }


  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************
  //************************************************************************************
  //************************************************************************************

  //************************************************************************************

  void SpatialLagrangian3DElement::SetGeneralVariables(GeneralVariables& rVariables,
						       ConstitutiveLaw::Parameters& rValues,
						       const int & rPointNumber)
  {
    LargeDisplacement3DElement::SetGeneralVariables(rVariables,rValues,rPointNumber);

    //Set extra options for the contitutive law
    Flags &ConstitutiveLawOptions=rValues.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::FINAL_CONFIGURATION);

  }

  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************


  void SpatialLagrangian3DElement::CalculateKinematics(GeneralVariables& rVariables,
						     const double& rPointNumber)

  {
    KRATOS_TRY
      
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //Parent to reference configuration
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;
    
    //Calculating the inverse of the jacobian and the parameters needed
    Matrix InvJ;
    MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

    //Compute cartesian derivatives
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], InvJ );

    //Current Deformation Gradient
    this->CalculateDeformationGradient (rVariables.DN_DX, rVariables.F, rVariables.DeltaPosition);

    //Calculating the inverse of the jacobian and the parameters needed
    Matrix Invj;
    MathUtils<double>::InvertMatrix( rVariables.j[rPointNumber], Invj, rVariables.detJ); //overwrites detJ 

    //Compute cartesian derivatives
    rVariables.DN_DX = prod( DN_De[rPointNumber], Invj ); //overwrites DX now is the current position dx

    //Determinant of the Deformation Gradient F0
    rVariables.detF0                 = mDeterminantF0[rPointNumber];
    rVariables.ElasticLeftCGVector   = mElasticLeftCauchyGreenVector[rPointNumber];

    //Set Shape Functions Values for this integration point
    rVariables.N=row( Ncontainer, rPointNumber);

    //Compute the deformation matrix B
    this->CalculateDeformationMatrix(rVariables.B, rVariables.F, rVariables.DN_DX);


    KRATOS_CATCH( "" )
      }


  //*************************COMPUTE DELTA POSITION*************************************
  //************************************************************************************

  Matrix& SpatialLagrangian3DElement::CalculateDeltaPosition(Matrix & DeltaPosition)
  {
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    
    DeltaPosition = zero_matrix<double>( number_of_nodes , dimension);
   
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {	    
	array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
	array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
	    	    
	//std::cout<<" CurrentDisplacement "<<CurrentDisplacement<<" PreviousDisplacement "<<PreviousDisplacement<<std::endl;
	for ( unsigned int j = 0; j < dimension; j++ )
	  {	    
	    DeltaPosition(i,j) = CurrentDisplacement[j]-PreviousDisplacement[j];
	  }
      }


    return DeltaPosition;

    KRATOS_CATCH( "" )
      }


  //*************************COMPUTE DEFORMATION GRADIENT*******************************
  //************************************************************************************

  void SpatialLagrangian3DElement::CalculateDeformationGradient(const Matrix& rDN_DX,
								Matrix& rF,
								Matrix& DeltaPosition)
  {
    KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    rF = identity_matrix<double> ( 3 );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {

	rF ( 0 , 0 ) += DeltaPosition(i,0)*rDN_DX ( i , 0 );
	rF ( 0 , 1 ) += DeltaPosition(i,0)*rDN_DX ( i , 1 );
	rF ( 0 , 2 ) += DeltaPosition(i,0)*rDN_DX ( i , 2 );
	rF ( 1 , 0 ) += DeltaPosition(i,1)*rDN_DX ( i , 0 );
	rF ( 1 , 1 ) += DeltaPosition(i,1)*rDN_DX ( i , 1 );
	rF ( 1 , 2 ) += DeltaPosition(i,1)*rDN_DX ( i , 2 );
	rF ( 2 , 0 ) += DeltaPosition(i,2)*rDN_DX ( i , 0 );
	rF ( 2 , 1 ) += DeltaPosition(i,2)*rDN_DX ( i , 1 );
	rF ( 2 , 2 ) += DeltaPosition(i,2)*rDN_DX ( i , 2 );
      }

      

    KRATOS_CATCH( "" )
      }




  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangian3DElement::CalculateDeformationMatrix(Matrix& rB,
							      Matrix& rF,
							      Matrix& rDN_DX)
  {
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
 
    rB.clear();
    
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

    KRATOS_CATCH( "" )
      }




  //************************************************************************************
  //************************************************************************************


  void SpatialLagrangian3DElement::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacement3DElement );
    rSerializer.save("ElasticLeftCauchyGreenVector",mElasticLeftCauchyGreenVector);
    rSerializer.save("DeterminantF0",mDeterminantF0);

   }

  void SpatialLagrangian3DElement::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacement3DElement );
    rSerializer.load("ElasticLeftCauchyGreenVector",mElasticLeftCauchyGreenVector);
    rSerializer.load("DeterminantF0",mDeterminantF0);

  }


} // Namespace Kratos


