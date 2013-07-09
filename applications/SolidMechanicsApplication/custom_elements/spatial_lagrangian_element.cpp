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
#include "custom_elements/spatial_lagrangian_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "solid_mechanics_application.h"


namespace Kratos
{


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  SpatialLagrangianElement::SpatialLagrangianElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
  {
    //DO NOT ADD DOFS HERE!!!
  }


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  SpatialLagrangianElement::SpatialLagrangianElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
  {
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
  }


  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  SpatialLagrangianElement::SpatialLagrangianElement( SpatialLagrangianElement const& rOther)
    :Element(rOther)
    ,mThisIntegrationMethod(rOther.mThisIntegrationMethod)
    ,mConstitutiveLawVector(rOther.mConstitutiveLawVector)
    ,mDeformationGradientF0(rOther.mDeformationGradientF0)
    ,mDeterminantF0(rOther.mDeterminantF0)
  {
  }



  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  SpatialLagrangianElement&  SpatialLagrangianElement::operator=(SpatialLagrangianElement const& rOther)
  {
    Element::operator=(rOther);

    mThisIntegrationMethod = rOther.mThisIntegrationMethod;

    mConstitutiveLawVector.clear();
    mConstitutiveLawVector.resize(rOther.mConstitutiveLawVector.size());

    mDeformationGradientF0.clear();
    mDeformationGradientF0.resize(rOther.mDeformationGradientF0.size());


    for(unsigned int i=0; i<<mConstitutiveLawVector.size(); i++)
      {
        mConstitutiveLawVector[i] = rOther.mConstitutiveLawVector[i];
	mDeformationGradientF0 [i] = rOther.mDeformationGradientF0[i];
      }

    mDeterminantF0 = rOther.mDeterminantF0;

    return *this;
  }



  //*********************************OPERATIONS*****************************************
  //************************************************************************************

  Element::Pointer SpatialLagrangianElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
  {
    return Element::Pointer(new SpatialLagrangianElement( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
  }


  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************


  SpatialLagrangianElement::~SpatialLagrangianElement()
  {
  }



  //************* GETTING METHODS
  //************************************************************************************
  //************************************************************************************

  SpatialLagrangianElement::IntegrationMethod SpatialLagrangianElement::GetIntegrationMethod() const
  {
    return mThisIntegrationMethod;
  }


  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangianElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
  {
    rElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
      {
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

        if ( GetGeometry().WorkingSpaceDimension() == 3 )
	  {
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
	  }
      }
  }

  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangianElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
  {
    int number_of_nodes = GetGeometry().size();
    int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int dim2 = number_of_nodes * dimension;

    if ( rResult.size() != dim2 )
      rResult.resize( dim2, false );

    for ( int i = 0; i < number_of_nodes; i++ )
      {
        int index = i * dimension;
        rResult[index] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

        if ( dimension == 3 )
	  rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
      }
  }

  //*********************************DISPLACEMENT***************************************
  //************************************************************************************

  void SpatialLagrangianElement::GetValuesVector( Vector& rValues, int Step )
  {
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * dimension;

    if ( rValues.size() != MatSize ) rValues.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        unsigned int index = i * dimension;
        rValues[index] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 )
	  rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
      }
  }


  //************************************VELOCITY****************************************
  //************************************************************************************

  void SpatialLagrangianElement::GetFirstDerivativesVector( Vector& rValues, int Step )
  {
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * dimension;

    if ( rValues.size() != MatSize ) rValues.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        unsigned int index = i * dimension;
        rValues[index] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

        if ( dimension == 3 )
	  rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
      }
  }

  //*********************************ACCELERATION***************************************
  //************************************************************************************

  void SpatialLagrangianElement::GetSecondDerivativesVector( Vector& rValues, int Step )
  {
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * dimension;

    if ( rValues.size() != MatSize ) rValues.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        unsigned int index = i * dimension;
        rValues[index] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dimension == 3 )
	  rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
      }
  }


  //*********************************SET DOUBLE VALUE***********************************
  //************************************************************************************

  void SpatialLagrangianElement::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
							      std::vector<double>& rValues,
							      const ProcessInfo& rCurrentProcessInfo )
  {
    for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
      {
	mConstitutiveLawVector[PointNumber]->SetValue( rVariable,
						       rValues[PointNumber], rCurrentProcessInfo );
      }
  }

  //*********************************SET VECTOR VALUE***********************************
  //************************************************************************************

  void SpatialLagrangianElement::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
  {

    for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
      {
	mConstitutiveLawVector[PointNumber]->SetValue( rVariable,
						       rValues[PointNumber], rCurrentProcessInfo );
      }
    

  }


  //*********************************SET MATRIX VALUE***********************************
  //************************************************************************************

  void SpatialLagrangianElement::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
							      std::vector<Matrix>& rValues,
							      const ProcessInfo& rCurrentProcessInfo )
  {
    for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
      {
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable,
						       rValues[PointNumber], rCurrentProcessInfo );
      }

  }

  //********************************SET CONSTITUTIVE VALUE******************************
  //************************************************************************************

  void SpatialLagrangianElement::SetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
							      std::vector<ConstitutiveLaw::Pointer>& rValues,
							      const ProcessInfo& rCurrentProcessInfo )
  {
    if(rVariable == CONSTITUTIVE_LAW)
      {
	if ( mConstitutiveLawVector.size() != rValues.size() )
	  {
	    mConstitutiveLawVector.resize(rValues.size());

	    if( mConstitutiveLawVector.size() != GetGeometry().IntegrationPointsNumber() )
	      KRATOS_ERROR( std::logic_error, "constitutive law not has the correct size ", mConstitutiveLawVector.size() );
	  }
     
	for(unsigned int i=0; i<rValues.size(); i++)
	  {
	    mConstitutiveLawVector[i] = rValues[i]->Clone();

	  }
      }

    // if(rVariable == CONSTITUTIVE_LAW_POINTER)
    //   {
    //     if ( mConstitutiveLawVector.size() != rValues.size() )
    // 	{
    // 	  mConstitutiveLawVector.resize(rValues.size());

    // 	  if( mConstitutiveLawVector.size() != GetGeometry().IntegrationPointsNumber() )
    // 	    KRATOS_ERROR( std::logic_error, "constitutive law not has the correct size ", mConstitutiveLawVector.size() );
    // 	}
     
    //     for(unsigned int i=0; i<rValues.size(); i++)
    // 	{
    // 	  mConstitutiveLawVector[i] = rValues[i];

    // 	}
    //   }
    
  }

  //*********************************GET DOUBLE VALUE***********************************
  //************************************************************************************

  void SpatialLagrangianElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
							      std::vector<double>& rValues,
							      const ProcessInfo& rCurrentProcessInfo )
  {
    if ( rValues.size() != GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
      rValues.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(), false );

    if ( rVariable == VON_MISES_STRESS ){
      Vector Values;
      CalculateOnIntegrationPoints( rVariable, Values, rCurrentProcessInfo );
      for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
	rValues[ii] = Values[ii];
    }
    else{
      for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
        rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
    }
  }


  //**********************************GET VECTOR VALUE**********************************
  //************************************************************************************

  void SpatialLagrangianElement::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
							      std::vector<Vector>& rValues,
							      const ProcessInfo& rCurrentProcessInfo )
  {
    const unsigned int& size = GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size();

    if ( rValues.size() != size )
      rValues.resize( size );



    for ( unsigned int PointNumber = 0;  PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
      {
	rValues[PointNumber] =
	  mConstitutiveLawVector[PointNumber]->GetValue( rVariable, rValues[PointNumber] );
      }

    if ( rVariable == PK2_STRESS_VECTOR ||  rVariable == CAUCHY_STRESS_VECTOR )
      CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

  }

  //***********************************GET MATRIX VALUE*********************************
  //************************************************************************************

  void SpatialLagrangianElement::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
							      std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
  {

    CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
 
  }


  //********************************GET CONSTITUTIVE VALUE******************************
  //************************************************************************************

  void SpatialLagrangianElement::GetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
							      std::vector<ConstitutiveLaw::Pointer>& rValues,
							      const ProcessInfo& rCurrentProcessInfo )
  {
    // if(rVariable == CONSTITUTIVE_LAW || rVariable == CONSTITUTIVE_LAW_POINTER)
    //   {
    //     if ( rValues.size() != mConstitutiveLawVector.size() )
    // 	{
    // 	  rValues.resize(mConstitutiveLawVector.size());
    // 	}
     
    //     for(unsigned int i=0; i<rValues.size(); i++)
    // 	{
    // 	  rValues[i] = mConstitutiveLawVector[i];
    // 	}
    //   }

  }

  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangianElement::Initialize()
  {
    KRATOS_TRY

      //Number of integration points
      SizeType integration_points_number=GetGeometry().IntegrationPointsNumber();

    //Constitutive Law initialisation
    if ( mConstitutiveLawVector.size() != integration_points_number )
      {
        mConstitutiveLawVector.resize( integration_points_number );
      }


    //Material initialisation
    InitializeMaterial();

    //Resize historic deformation gradient
    mDeformationGradientF0.resize( integration_points_number );
    mDeterminantF0.resize( integration_points_number, false );



    KRATOS_CATCH( "" )
      }



  ////************************************************************************************
  ////************************************************************************************

  void SpatialLagrangianElement::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
  {

    ClearNodalForces();
    
    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
      mConstitutiveLawVector[i]->InitializeSolutionStep( GetProperties(),
							 GetGeometry(),
							 row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ),
							 CurrentProcessInfo );
  }

  ////************************************************************************************
  ////************************************************************************************

  void SpatialLagrangianElement::InitializeNonLinearIteration( ProcessInfo& CurrentProcessInfo )
  {
    ClearNodalForces();
  }


  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangianElement::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
  {
    KRATOS_TRY

    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),CurrentProcessInfo);

    Flags &Options=Values.GetOptions();

    Options.Set(ConstitutiveLaw::COMPUTE_STRESS);
    Options.Reset(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    Standard Variables;
    InitializeStandardVariables(Variables,CurrentProcessInfo);

    //reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {

        //COMPUTE kinematics B,F,DN_DX ...
        CalculateKinematics(Variables,PointNumber);

        //set standart parameters
        SetStandardParameters(Variables,Values,PointNumber);

        //call the constitutive law to update material variables
        //returns the variables increment (stresses, strains and internal)
        mConstitutiveLawVector[PointNumber]->FinalizeMaterialResponseCauchy (Values);

	mDeformationGradientF0[PointNumber] = Values.GetDeformationGradientF0();
	mDeterminantF0[PointNumber]         = Values.GetDeterminantF0();

        mConstitutiveLawVector[PointNumber]->FinalizeSolutionStep( GetProperties(),
								   GetGeometry(),
								   Variables.N,
								   CurrentProcessInfo );
      }

    KRATOS_CATCH( "" )
      }

  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangianElement::InitializeMaterial()
  {
    KRATOS_TRY

      if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
	{
	  for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
	    {
	      mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
	      mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(),
							     row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
	    }
	}
      else
	{
	  KRATOS_ERROR( std::logic_error, "a constitutive law needs to be specified for the element with ID ", this->Id() )
	    }

    KRATOS_CATCH( "" )
      }

  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangianElement::ResetConstitutiveLaw()
  {
    KRATOS_TRY

      if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
        {
	  for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
	    mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }

    KRATOS_CATCH( "" )
      }

  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangianElement::ClearNodalForces()
  {
    KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	
	array_1d<double, 3 > & ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(FORCE_EXTERNAL);
	array_1d<double, 3 > & InternalForce = GetGeometry()[i].FastGetSolutionStepValue(FORCE_INTERNAL);
	array_1d<double, 3 > & DynamicForce  = GetGeometry()[i].FastGetSolutionStepValue(FORCE_DYNAMIC);
	
	ExternalForce.clear();
	InternalForce.clear();
	DynamicForce.clear();

      }

    KRATOS_CATCH( "" )
      }




  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************


  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************


  void SpatialLagrangianElement::CalculateKinematics(Standard& rVariables,
						     const double& rPointNumber)

  {
    KRATOS_TRY

      const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Parent to current configuration
    Matrix j ( dimension , dimension);
    j = CalculateJacobian ( j ); //Must be generalized


    // Parent to reference configuration
    Matrix J ( dimension , dimension);
    J = GetGeometry().Jacobian( J, rPointNumber , mThisIntegrationMethod );


    Matrix InvJ;
    //Calculating the inverse of the jacobian and the parameters needed
    MathUtils<double>::InvertMatrix( J, InvJ, rVariables.detJ);

    //Compute cartesian derivatives
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber] , InvJ );

    //Current Deformation Gradient
    CalculateDeformationGradient (rVariables.DN_DX, rVariables.F );


    Matrix Invj;
    //Calculating the inverse of the jacobian and the parameters needed
    MathUtils<double>::InvertMatrix( j, Invj, rVariables.detJ); //overwrites detJ 

    //Compute cartesian derivatives
    rVariables.DN_DX = prod( DN_De[rPointNumber] , Invj ); //overwrites DX now is the current position dx

    //Compute the deformation matrix
    CalculateDeformationMatrix(rVariables.B,rVariables.F,rVariables.DN_DX);

    //Compute Almansi Strain Increment == Simetric Gradient of Displacements
    CalculateAlmansiStrainIncrement(rVariables.F,rVariables.DN_DX,rVariables.StrainVector);

    //Left Cauchy-Green Calculation
    //noalias( rMatrix ) = prod( rVariables.F, trans( rVariables.F ) );


    KRATOS_CATCH( "" )
      }

  //*************************COMPUTE JACOBIAN GRADIENT*******************************
  //************************************************************************************

  Matrix & SpatialLagrangianElement::CalculateJacobian(Matrix& rJ)
  {
    KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    std::vector<array_1d<double, 3 > > CurrentPosition (number_of_nodes);

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if ( dimension == 2 )
      {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
	    
	    array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
	    array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
	    array_1d<double, 3 > & ReferencePosition    = GetGeometry()[i].Coordinates();
	    
	    CurrentPosition[i] = ReferencePosition + (CurrentDisplacement-PreviousDisplacement);

	  }
    
	rJ.resize( dimension , dimension );

	rJ( 0, 0 ) = -( CurrentPosition[0][0] ) + ( CurrentPosition[1][0] );
	rJ( 1, 0 ) = -( CurrentPosition[0][1] ) + ( CurrentPosition[1][1] );
	rJ( 0, 1 ) = -( CurrentPosition[0][0] ) + ( CurrentPosition[2][0] );
	rJ( 1, 1 ) = -( CurrentPosition[0][1] ) + ( CurrentPosition[2][1] );
      }    

    return rJ;

    KRATOS_CATCH( "" )
      }


  //*************************COMPUTE DEFORMATION GRADIENT*******************************
  //************************************************************************************

  void SpatialLagrangianElement::CalculateDeformationGradient(const Matrix& rDN_DX,
							      Matrix& rF )
  {
    KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();


    if ( dimension == 2 )
      {
        rF=identity_matrix<double> ( 2 );

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
            //Displacement from the reference to the current configuration
            array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            array_1d<double, 3 > DeltaDisplacement=CurrentDisplacement-PreviousDisplacement;
	    
	    //std::cout<<" Node "<<GetGeometry()[i].Id()<<" DeltaDisplacement "<<DeltaDisplacement<<std::endl;

            rF ( 0 , 0 ) += DeltaDisplacement[0]*rDN_DX ( i , 0 );
            rF ( 0 , 1 ) += DeltaDisplacement[0]*rDN_DX ( i , 1 );
            rF ( 1 , 0 ) += DeltaDisplacement[1]*rDN_DX ( i , 0 );
            rF ( 1 , 1 ) += DeltaDisplacement[1]*rDN_DX ( i , 1 );

	  }

	//std::cout<<" rF "<<rF<<std::endl;
      }


    if ( dimension == 3 )
      {

        rF=identity_matrix<double> ( 3 );

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
            //Displacement from the reference to the current configuration
            array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            array_1d<double, 3 > DeltaDisplacement=CurrentDisplacement-PreviousDisplacement;

            rF ( 0 , 0 ) += DeltaDisplacement[0]*rDN_DX ( i , 0 );
            rF ( 0 , 1 ) += DeltaDisplacement[0]*rDN_DX ( i , 1 );
            rF ( 0 , 2 ) += DeltaDisplacement[0]*rDN_DX ( i , 2 );
            rF ( 1 , 0 ) += DeltaDisplacement[1]*rDN_DX ( i , 0 );
            rF ( 1 , 1 ) += DeltaDisplacement[1]*rDN_DX ( i , 1 );
            rF ( 1 , 2 ) += DeltaDisplacement[1]*rDN_DX ( i , 2 );
            rF ( 2 , 0 ) += DeltaDisplacement[0]*rDN_DX ( i , 0 );
            rF ( 2 , 1 ) += DeltaDisplacement[0]*rDN_DX ( i , 1 );
            rF ( 2 , 2 ) += DeltaDisplacement[0]*rDN_DX ( i , 2 );
	  }

      }

    KRATOS_CATCH( "" )
      }



  //****************************COMPUTE VELOCITY GRADIENT*******************************
  //************************************************************************************

  void SpatialLagrangianElement::CalculateVelocityGradient(const Matrix& rDN_DX,
							   Matrix& rDF )
  {
    KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();


    if ( dimension == 2 )
      {
        rDF=zero_matrix<double> ( 2 );

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
            //Displacement from the reference to the current configuration
            array_1d<double, 3 > & CurrentVelocity  = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);

            rDF ( 0 , 0 ) += CurrentVelocity[0]*rDN_DX ( i , 0 );
            rDF ( 0 , 1 ) += CurrentVelocity[0]*rDN_DX ( i , 1 );
            rDF ( 1 , 0 ) += CurrentVelocity[1]*rDN_DX ( i , 0 );
            rDF ( 1 , 1 ) += CurrentVelocity[1]*rDN_DX ( i , 1 );
	  }

      }


    if ( dimension == 3 )
      {

        rDF=zero_matrix<double> ( 3 );

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
            //Displacement from the reference to the current configuration
            array_1d<double, 3 > & CurrentVelocity  = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);

            rDF ( 0 , 0 ) += CurrentVelocity[0]*rDN_DX ( i , 0 );
            rDF ( 0 , 1 ) += CurrentVelocity[0]*rDN_DX ( i , 1 );
            rDF ( 0 , 2 ) += CurrentVelocity[0]*rDN_DX ( i , 2 );
            rDF ( 1 , 0 ) += CurrentVelocity[1]*rDN_DX ( i , 0 );
            rDF ( 1 , 1 ) += CurrentVelocity[1]*rDN_DX ( i , 1 );
            rDF ( 1 , 2 ) += CurrentVelocity[1]*rDN_DX ( i , 2 );
            rDF ( 2 , 0 ) += CurrentVelocity[0]*rDN_DX ( i , 0 );
            rDF ( 2 , 1 ) += CurrentVelocity[0]*rDN_DX ( i , 1 );
            rDF ( 2 , 2 ) += CurrentVelocity[0]*rDN_DX ( i , 2 );
	  }

      }

    KRATOS_CATCH( "" )
      }

  //*********************************COMPUTE STRAIN*************************************
  //************************************************************************************

  void SpatialLagrangianElement::CalculateAlmansiStrain(const Matrix& rF,
							Vector& rStrainVector ) //with rF=F0 -> total strain
  {
    KRATOS_TRY

      unsigned int dimension = GetGeometry().WorkingSpaceDimension();


    //Left Cauchy-Green Calculation
    Matrix LeftCauchyGreen = prod( rF, trans( rF ) );
    Matrix InverseLeftCauchyGreen;
    double det_b=0;
    //Calculating the inverse of the jacobian and the parameters needed
    MathUtils<double>::InvertMatrix( LeftCauchyGreen, InverseLeftCauchyGreen, det_b);


    if ( dimension == 2 )
      {
        if ( rStrainVector.size() != 3 ) rStrainVector.resize( 3, false );

        rStrainVector[0] = 0.5 * ( 1.00 - InverseLeftCauchyGreen( 0, 0 ) );

        rStrainVector[1] = 0.5 * ( 1.00 - InverseLeftCauchyGreen( 1, 1 ) );

        rStrainVector[2] = - InverseLeftCauchyGreen( 0, 1 );

      }

    if ( dimension == 3 )
      {
        if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

        rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );

        rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );

        rStrainVector[2] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 2, 2 ) );

        rStrainVector[3] = - InverseLeftCauchyGreen( 0, 1 ); // xy

        rStrainVector[4] = - InverseLeftCauchyGreen( 1, 2 ); // yz

        rStrainVector[5] = - InverseLeftCauchyGreen( 0, 2 ); // xz
      }

    KRATOS_CATCH( "" )
      }


  //*********************************COMPUTE STRAIN INCREMENT***************************
  //************************************************************************************

  void SpatialLagrangianElement::CalculateAlmansiStrainIncrement(const Matrix& rF,
								 const Matrix& rDN_DX,
								 Vector& rStrainVector ) //current increment
  {
    KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
   
    if ( dimension == 2 )
      {
	//Displacement Gradient
	Matrix H=zero_matrix<double> ( 2 );
	
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
            //Displacement from the reference to the current configuration
            array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            array_1d<double, 3 > DeltaDisplacement=CurrentDisplacement-PreviousDisplacement;

            H ( 0 , 0 ) += DeltaDisplacement[0]*rDN_DX ( i , 0 );
            H ( 0 , 1 ) += DeltaDisplacement[0]*rDN_DX ( i , 1 );
            H ( 1 , 0 ) += DeltaDisplacement[1]*rDN_DX ( i , 0 );
            H ( 1 , 1 ) += DeltaDisplacement[1]*rDN_DX ( i , 1 );

            //std::cout<<" node "<<i<<" -> DeltaDisplacement : "<<DeltaDisplacement<<std::endl;

	  }

        if ( rStrainVector.size() != 3 ) rStrainVector.resize( 3, false );

        rStrainVector[0] = 0.5 * ( H( 0 , 0 ) + H( 0 , 0 ) );

        rStrainVector[1] = 0.5 * ( H( 1 , 1 ) + H( 1 , 1 ) );

        rStrainVector[2] = ( H( 0 , 1 ) + H( 1 , 0 ) );

      }


    if ( dimension == 3 )
      {

        Matrix H=zero_matrix<double> ( 3 );

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
            //Displacement from the reference to the current configuration
            array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            array_1d<double, 3 > DeltaDisplacement=CurrentDisplacement-PreviousDisplacement;

            H ( 0 , 0 ) += DeltaDisplacement[0]*rDN_DX ( i , 0 );
            H ( 0 , 1 ) += DeltaDisplacement[0]*rDN_DX ( i , 1 );
            H ( 0 , 2 ) += DeltaDisplacement[0]*rDN_DX ( i , 2 );
            H ( 1 , 0 ) += DeltaDisplacement[1]*rDN_DX ( i , 0 );
            H ( 1 , 1 ) += DeltaDisplacement[1]*rDN_DX ( i , 1 );
            H ( 1 , 2 ) += DeltaDisplacement[1]*rDN_DX ( i , 2 );
            H ( 2 , 0 ) += DeltaDisplacement[0]*rDN_DX ( i , 0 );
            H ( 2 , 1 ) += DeltaDisplacement[0]*rDN_DX ( i , 1 );
            H ( 2 , 2 ) += DeltaDisplacement[0]*rDN_DX ( i , 2 );
	  }

        if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

        rStrainVector[0] = 0.5 * ( H( 0 , 0 ) + H( 0 , 0 ) );

        rStrainVector[1] = 0.5 * ( H( 1 , 1 ) + H( 1 , 1 ) );

        rStrainVector[2] = 0.5 * ( H( 2 , 2 ) + H( 2 , 2 ) );

        rStrainVector[3] = ( H( 0 , 1 ) + H( 1 , 0 ) ); //xy

        rStrainVector[4] = ( H( 1 , 2 ) + H( 2 , 1 ) ); //yz

        rStrainVector[5] = ( H( 0 , 2 ) + H( 2 , 0 ) ); //xz

      }


    KRATOS_CATCH( "" )
      }


  //*************************DECIMAL CORRECTION OF STRAINS******************************
  //************************************************************************************

  void SpatialLagrangianElement::DecimalCorrection(Vector& rStrainVector)
  { 
    KRATOS_TRY
    
      for ( unsigned int i = 0; i < rStrainVector.size(); i++ )
	{
	  if( rStrainVector[i]*rStrainVector[i]<1e-24 )
	    {
	      rStrainVector[i]=0;
	    }
	
	}

    KRATOS_CATCH( "" )
      }

  //*************************COMPUTE DEFORMATION MATRIX B*******************************
  //************************************************************************************

  void SpatialLagrangianElement::CalculateDeformationMatrix(Matrix& rB,
							    Matrix& rF,
							    Matrix& rDN_DX)
  {
    KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rB = zero_matrix<double> ( rB.size1() , rB.size2() );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        unsigned int index = dimension * i;

        if ( dimension == 2 )
	  {
            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rDN_DX( i, 1 );
            rB( 2, index + 1 ) = rDN_DX( i, 0 );
	  }
        else
	  {
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

    KRATOS_CATCH( "" )
      }


  //***********************COMPUTE LOCAL SYSTEM CONTRIBUTIONS***************************
  //************************************************************************************


  void SpatialLagrangianElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
  {
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag  = true;


    CalculateElementalSystem( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );


  }

  //************************************************************************************
  //************************************************************************************


  void SpatialLagrangianElement::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
  {
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateElementalSystem( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
  }


  //************************************************************************************
  //************************************************************************************


  void SpatialLagrangianElement::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
  {
    if (rLeftHandSideMatrix.size1() != 0)
      rLeftHandSideMatrix.resize(0, 0);
  }


  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangianElement::SetStandardParameters(Standard& rVariables,
						       ConstitutiveLaw::Parameters& rValues,
						       const int & rPointNumber)
  {
    rVariables.detF =MathUtils<double>::Det(rVariables.F);

    if(rVariables.detF<0)
      KRATOS_ERROR(std::invalid_argument,"DeterminantF < 0","");

    rValues.SetDeformationGradientF0(mDeformationGradientF0[rPointNumber]);
    rValues.SetDeterminantF0(mDeterminantF0[rPointNumber]);

    rValues.SetDeterminantF(rVariables.detF);
    rValues.SetDeformationGradientF(rVariables.F);
    rValues.SetStrainVector(rVariables.StrainVector);
    rValues.SetStressVector(rVariables.StressVector);
    rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rValues.SetShapeFunctionsDevivatives(rVariables.DN_DX);

    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );
    rVariables.N=row(Ncontainer , rPointNumber);
    rValues.SetShapeFunctionsValues(rVariables.N);
  }

  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangianElement::InitializeStandardVariables (Standard & rVariables, const ProcessInfo& rCurrentProcessInfo)
  {
  
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
 
    unsigned int StrainSize;
  
    if ( dimension == 2 )
      {
	StrainSize = 3;
      }
    else
      {
	StrainSize = 6;
      } 

    rVariables.B.resize( StrainSize, number_of_nodes * dimension );
  
    rVariables.F.resize( dimension, dimension );
  
    rVariables.F0.resize( dimension, dimension );

    rVariables.ConstitutiveMatrix.resize( StrainSize, StrainSize );
  
    rVariables.StrainVector.resize( StrainSize );
  
    rVariables.StressVector.resize( StrainSize );
  
    rVariables.DN_DX.resize( number_of_nodes, dimension );
  
    
  }


  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangianElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
							  VectorType& rRightHandSideVector,
							  bool CalculateStiffnessMatrixFlag,
							  bool CalculateResidualVectorFlag )
  {

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dimension;
      
    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
      {
	if ( rLeftHandSideMatrix.size1() != MatSize )
	  rLeftHandSideMatrix.resize( MatSize, MatSize, false );
	  
	noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
      }
      
      
    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
      {
	if ( rRightHandSideVector.size() != MatSize )
	  rRightHandSideVector.resize( MatSize, false );
	  
	rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
      }
  }



  //************************************************************************************
  //************************************************************************************


  void SpatialLagrangianElement::CalculateElementalSystem( MatrixType& rLeftHandSideMatrix,
							   VectorType& rRightHandSideVector,
							   ProcessInfo& rCurrentProcessInfo,
							   bool CalculateStiffnessMatrixFlag,
							   bool CalculateResidualVectorFlag )
  {
    KRATOS_TRY

    //Create Constitutive Law Parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    //Set Constitutive Law Flags:
    Flags &Options=Values.GetOptions();
    Options.Set(ConstitutiveLaw::COMPUTE_STRESS);
    Options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    //Create and Initialize Element Variables:
    Standard Variables;
    InitializeStandardVariables(Variables,rCurrentProcessInfo);

    //Initialize sizes for the system components:
    InitializeSystemMatrices(rLeftHandSideMatrix,rRightHandSideVector,CalculateStiffnessMatrixFlag,CalculateResidualVectorFlag);


    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    
    //auxiliary terms
    Vector BodyForce;

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {

        //COMPUTE kinematics B,F,DN_DX ...
        CalculateKinematics(Variables,PointNumber);

        //set standart parameters
        SetStandardParameters(Variables,Values,PointNumber);

        //CALL the constitutive law
        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy (Values);

        double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;


        if ( dimension == 2 ) IntegrationWeight *= GetProperties()[THICKNESS];

        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
	  {
            //contributions to stiffness matrix calculated on the reference config

            // operation performed: add Km to the rLefsHandSideMatrix
            CalculateAndAddKm( rLeftHandSideMatrix, Variables.B, Variables.ConstitutiveMatrix, IntegrationWeight );

            // operation performed: add Kg to the rLefsHandSideMatrix
            CalculateAndAddKg( rLeftHandSideMatrix, Variables.DN_DX, Variables.StressVector, IntegrationWeight );
	  }

        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
	  {
            //contribution to external forces
            BodyForce = GetProperties()[BODY_FORCE]*GetProperties()[DENSITY];

            // operation performed: rRightHandSideVector += ExtForce*IntegrationWeight
            CalculateAndAddExternalForces(Variables.N,rCurrentProcessInfo, BodyForce, rRightHandSideVector, IntegrationWeight );

            // operation performed: rRightHandSideVector -= IntForce*IntegrationWeight
            CalculateAndAddInternalForces(Variables.B,Variables.StressVector,rRightHandSideVector,IntegrationWeight);

	  }
      }


    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangianElement::MassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

      //lumped
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int NumberOfNodes = GetGeometry().size();
    unsigned int MatSize = dimension * NumberOfNodes;

    if ( rMassMatrix.size1() != MatSize )
      rMassMatrix.resize( MatSize, MatSize, false );

    rMassMatrix = ZeroMatrix( MatSize, MatSize );

    double TotalMass = GetGeometry().DomainSize() * GetProperties()[DENSITY];

    if ( dimension == 2 ) TotalMass *= GetProperties()[THICKNESS];

    Vector LumpFact;

    LumpFact = GetGeometry().LumpingFactors( LumpFact );

    for ( unsigned int i = 0; i < NumberOfNodes; i++ )
      {
        double temp = LumpFact[i] * TotalMass;

        for ( unsigned int j = 0; j < dimension; j++ )
	  {
            unsigned int index = i * dimension + j;
            rMassMatrix( index, index ) = temp;
	  }
      }

    // std::cout<<std::endl;
    // std::cout<<" Mass Matrix "<<rMassMatrix<<std::endl;


    KRATOS_CATCH( "" )
      }

  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangianElement::DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY
      unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dimension;

    if ( rDampMatrix.size1() != MatSize )
      rDampMatrix.resize( MatSize, MatSize, false );

    noalias( rDampMatrix ) = ZeroMatrix( MatSize, MatSize );

    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************


  double SpatialLagrangianElement::CalculateIntegrationWeight( double& GaussPointWeight, double& DetJ0 )
  {
    //to permorm the integration over the reference domain we need to include
    // the thickness in 2D
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    double weight = GaussPointWeight;

    weight *= DetJ0;

    if ( dimension == 2 ) weight *= GetProperties()[THICKNESS];

    return weight;
  }

  //************************************************************************************
  //************************************************************************************

  inline void SpatialLagrangianElement::CalculateAndAddExternalForces(const Vector& N,
								      const ProcessInfo& CurrentProcessInfo,
								      Vector& BodyForce,
								      VectorType& rRightHandSideVector,
								      double& rIntegrationWeight)
								    
  {
    KRATOS_TRY
      unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    
    double Fext=0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        int index = dimension * i;
	
	array_1d<double, 3 > & ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(FORCE_EXTERNAL);

	Fext=0;
        for ( unsigned int j = 0; j < dimension; j++ )
	  {
	    Fext = rIntegrationWeight * N[i] * BodyForce[j];
            rRightHandSideVector[index + j] +=Fext;
	    ExternalForce[j] += Fext;
	  }
      }
    
    // std::cout<<std::endl;
    // std::cout<<" Fext "<<rRightHandSideVector<<std::endl;


    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************

  inline void SpatialLagrangianElement::CalculateAndAddInternalForces(Matrix & rB,
								      Vector& rStressVector,
								      VectorType& rRightHandSideVector,
								      double& rIntegrationWeight
								      )
  {
    KRATOS_TRY

      VectorType InternalForces = rIntegrationWeight * prod( trans( rB ), rStressVector );
    noalias( rRightHandSideVector ) -= InternalForces;
      
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	unsigned int indexu  = dimension * i;
	array_1d<double, 3 > & InternalForce = GetGeometry()[i].FastGetSolutionStepValue(FORCE_INTERNAL);

	for ( unsigned int j = 0; j < dimension; j++ )
	  {
  	    InternalForce[j] -= InternalForces [indexu+j];
	  }
      }

    
    // std::cout<<std::endl;
    // std::cout<<" Fint "<<InternalForces<<std::endl;

    KRATOS_CATCH( "" )
      }



  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangianElement::CalculateAndAddKm(MatrixType& rK,
						   Matrix& rB,
						   Matrix& rD,
						   double& rIntegrationWeight
						   )
  {
    KRATOS_TRY

      //contributions to stiffness matrix calculated on the reference config
      noalias( rK ) += prod( trans( rB ),  rIntegrationWeight * Matrix( prod( rD, rB ) ) ); //to be optimized to remove the temporary

    // std::cout<<std::endl;
    // std::cout<<" Kmat "<<rK<<std::endl;

    KRATOS_CATCH( "" )
      }



  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangianElement::CalculateAndAddKg(MatrixType& rK,
						   Matrix& rDN_DX,
						   Vector& rStressVector,
						   double& rIntegrationWeight
						   )
  {
    KRATOS_TRY
      // unsigned int dimension = mpReferenceGeometry->WorkingSpaceDimension();
      //Matrix<double> StressTensor = MathUtils<double>::StressVectorToTensor(StressVector);
      //Matrix<double> ReducedKg(DN_Dx.RowsNumber(),DN_Dx.RowsNumber());
      //Matrix<double>::MatMulAndAdd_B_D_Btrans(ReducedKg,weight,DN_Dx,StressTensor);
      //MathUtils<double>::ExpandAndAddReducedMatrix(K,ReducedKg,dimension);

      MatrixType Kh=rK;

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rStressVector );
    Matrix ReducedKg = prod( rDN_DX,  rIntegrationWeight * Matrix( prod( StressTensor, trans( rDN_DX ) ) ) ); //to be optimized

    // std::cout<<std::endl;
    // std::cout<<" ReducedKg "<<ReducedKg<<std::endl;

    MathUtils<double>::ExpandAndAddReducedMatrix( rK, ReducedKg, dimension );

    // std::cout<<" DN_DX "<<rDN_DX<<std::endl;

    // std::cout<<std::endl;
    // std::cout<<" Kgeo "<<rK-Kh<<std::endl;


    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangianElement::CalculateOnIntegrationPoints( const Variable<double>& rVariable, Vector& rOutput, const ProcessInfo& rCurrentProcessInfo )
  {

    KRATOS_TRY

    if ( rOutput.size() != GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
      rOutput.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(), false );

    if ( rVariable == VON_MISES_STRESS )
      {

	ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
	Flags &Options=Values.GetOptions();
	Options.Set(ConstitutiveLaw::COMPUTE_STRESS);

	Standard Variables;
	InitializeStandardVariables (Variables,rCurrentProcessInfo);
 
	//reading integration points
	const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

	for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
	  {
	    //COMPUTE kinematics B,F,DN_DX ...
	    CalculateKinematics(Variables,PointNumber);
	
	    //set standart parameters
	    SetStandardParameters(Variables,Values,PointNumber);
	    
	    //CALL the constitutive law
	    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

	    ComparisonUtils EquivalentStress;
	    rOutput[PointNumber] =  EquivalentStress.CalculateVonMises(Variables.StressVector);
	  }
      }
    else{
      
      for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
        rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rOutput[ii] );
    }

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangianElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo )
  {

    KRATOS_TRY
    unsigned int StrainSize;

    if ( GetGeometry().WorkingSpaceDimension() == 2 )
      {
        StrainSize = 3;
      }
    else
      {
        StrainSize = 6;
      }


    if ( rVariable == PK2_STRESS_VECTOR )
      {

	ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
	Flags &Options=Values.GetOptions();
	Options.Set(ConstitutiveLaw::COMPUTE_STRESS);

	Standard Variables;
	InitializeStandardVariables (Variables,rCurrentProcessInfo);
      
	//reading integration points
	const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

	for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
	  {
	  
	    //COMPUTE kinematics B,F,DN_DX ...
	    CalculateKinematics(Variables,PointNumber);
	
	    //set standart parameters
	    SetStandardParameters(Variables,Values,PointNumber);
	    
	    //CALL the constitutive law
	    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

	    Variables.StressVector=Values.GetStressVector(Variables.StressVector);

	    mConstitutiveLawVector[PointNumber]->TransformStresses(Variables.StressVector,Variables.F,Variables.detF,ConstitutiveLaw::StressMeasure_Cauchy,ConstitutiveLaw::StressMeasure_PK2);

          if ( rOutput[PointNumber].size() != StrainSize )
	      rOutput[PointNumber].resize( StrainSize, false );

	    rOutput[PointNumber] = Variables.StressVector;
	  }
      }
    else
      {
	for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
	  {
            if ( rOutput[ii].size() != StrainSize )
	      rOutput[ii].resize( StrainSize, false );

            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable , rOutput[ii] );
	  }
      }
  
    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangianElement::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo )
  {

    KRATOS_TRY

    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
    
    unsigned int StrainSize = 0;
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    if ( dimension == 2 )
      StrainSize = 3;
    else
      StrainSize = 6;


    Standard Variables;
    InitializeStandardVariables(Variables,rCurrentProcessInfo);
    
    //reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );


    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {

        //COMPUTE kinematics B,F,DN_DX ...
        CalculateKinematics(Variables,PointNumber);

        if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
	  {

            Vector TotalStrainVector ( StrainSize );
            Matrix F0;
            //I need to solve it asking the F0 explicitly.
            mConstitutiveLawVector[PointNumber]->GetValue(DEFORMATION_GRADIENT,F0);
            CalculateAlmansiStrain( F0, TotalStrainVector );

	    TotalStrainVector = mConstitutiveLawVector[PointNumber]->TransformStrains(TotalStrainVector,Variables.F,ConstitutiveLaw::StrainMeasure_Almansi,ConstitutiveLaw::StrainMeasure_GreenLagrange);

            if ( rOutput[PointNumber].size2() != TotalStrainVector.size() )
	      rOutput[PointNumber].resize( 1, TotalStrainVector.size(), false );

            for ( unsigned int ii = 0; ii < TotalStrainVector.size(); ii++ )
	      rOutput[PointNumber]( 0, ii ) = TotalStrainVector[ii];

	  }
        else if ( rVariable == PK2_STRESS_TENSOR )   // in fact is a vector due to a mismatch on kratos-gid printing
	  { 

	    Flags &Options=Values.GetOptions();
	    Options.Set(ConstitutiveLaw::COMPUTE_STRESS);
     
	    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
	      {
	  
		//COMPUTE kinematics B,F,DN_DX ...
		CalculateKinematics(Variables,PointNumber);
	
		//set standart parameters
		SetStandardParameters(Variables,Values,PointNumber);
	    
		//CALL the constitutive law
		mConstitutiveLawVector[PointNumber]->CalculateMaterialResponsePK2(Values);
	
		if ( rOutput[PointNumber].size2() != Variables.StrainVector.size() )
		  rOutput[PointNumber].resize( 1, Variables.StrainVector.size(), false );
	    
		for ( unsigned int ii = 0; ii < Variables.StrainVector.size(); ii++ )
		  {
		    rOutput[PointNumber]( 0, ii ) = Variables.StressVector[ii];
		  }

	      }


	  }
        else if ( rVariable == CAUCHY_STRESS_TENSOR )  // in fact is a vector due to a mismatch on kratos-gid printing
	  {

	    Flags &Options=Values.GetOptions();
	    Options.Set(ConstitutiveLaw::COMPUTE_STRESS);
     
	    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
	      {
	  
		//COMPUTE kinematics B,F,DN_DX ...
		CalculateKinematics(Variables,PointNumber);
	
		//set standart parameters
		SetStandardParameters(Variables,Values,PointNumber);
	    
		//CALL the constitutive law
		mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);


		if ( rOutput[PointNumber].size2() != Variables.StrainVector.size() )
		  rOutput[PointNumber].resize( 1, Variables.StrainVector.size(), false );
	    
		for ( unsigned int ii = 0; ii < Variables.StrainVector.size(); ii++ )
		  {
		    rOutput[PointNumber]( 0, ii ) = Variables.StressVector[ii];
		  }


	      }
	  }

      }

    KRATOS_CATCH( "" )
      }




  //************************************************************************************
  //************************************************************************************

  void SpatialLagrangianElement::Calculate( const Variable<double>& rVariable, double& rOutput, const ProcessInfo& rCurrentProcessInfo )
  {

    double lamda = 1.00; // parametro que depende del tipo de problema y del elemento pag 308 libro dinamica de Barbat
    double c1 = 0.00; //sqrt(GetProperties()[YOUNG_MODULUS]/GetProperties()[DENSITY]); velocidad del sonido en el medio
    double c2 = 0.00; // norma de la velocidad actual dentro del elemento
    double c = 0.00;
    double wmax = 0.00;
    Vector Values( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() );
    Vector Velocities;

    GetFirstDerivativesVector( Velocities, 0 );

    if ( rVariable == DELTA_TIME )
      {
        for ( unsigned int PointNumber = 0;
	      PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size();
	      PointNumber++ )
	  {
            mConstitutiveLawVector[PointNumber]-> GetValue( DELTA_TIME, c1 );
            Values[PointNumber] = c1;
	  }
      }

    c1 = ( *std::max_element( Values.begin(), Values.end() ) );

    c2 = norm_2( Velocities );

    c = ( c1 > c2 ) ? c1 : c2;


    double le = GetGeometry().Length();
    //KRATOS_WATCH(le)

    /// maxima frecuencia de un elemento
    wmax = ( lamda * c ) / le;
    rOutput = 2.0 / wmax;
    //KRATOS_WATCH(rOutput)

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
  int  SpatialLagrangianElement::Check( const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

      unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();



    //verify that the variables are correctly initialized

    if ( VELOCITY.Key() == 0 )
      KRATOS_ERROR( std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "" );

    if ( DISPLACEMENT.Key() == 0 )
      KRATOS_ERROR( std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "" );

    if ( ACCELERATION.Key() == 0 )
      KRATOS_ERROR( std::invalid_argument, "ACCELERATION has Key zero! (check if the application is correctly registered", "" );

    if ( DENSITY.Key() == 0 )
      KRATOS_ERROR( std::invalid_argument, "DENSITY has Key zero! (check if the application is correctly registered", "" );

    if ( BODY_FORCE.Key() == 0 )
      KRATOS_ERROR( std::invalid_argument, "BODY_FORCE has Key zero! (check if the application is correctly registered", "" );

    if ( THICKNESS.Key() == 0 )
      KRATOS_ERROR( std::invalid_argument, "THICKNESS has Key zero! (check if the application is correctly registered", "" );

    //verify that the dofs exist
    for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
      {
        if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
	  KRATOS_ERROR( std::invalid_argument, "missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id() );

        if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )
	  KRATOS_ERROR( std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id() );
      }

    //verify that the constitutive law exists
    if ( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false )
      {
        KRATOS_ERROR( std::logic_error, "constitutive law not provided for property ", this->GetProperties().Id() );
      }

    //Verify that the body force is defined
    if ( this->GetProperties().Has( BODY_FORCE ) == false )
      {
        KRATOS_ERROR( std::logic_error, "BODY_FORCE not provided for property ", this->GetProperties().Id() )
	  }

    //verify that the constitutive law has the correct dimension
    if ( dimension == 2 )
      {
        if ( this->GetProperties().Has( THICKNESS ) == false )
	  KRATOS_ERROR( std::logic_error, "THICKNESS not provided for element ", this->Id() );

        if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 3 )
	  KRATOS_ERROR( std::logic_error, "wrong constitutive law used. This is a 2D element! expected strain size is 3 (el id = ) ", this->Id() );
      }
    else
      {
        if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 6 )
	  KRATOS_ERROR( std::logic_error, "wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) ", this->Id() );
      }

    //check constitutive law
    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
      {
        return mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo );
      }

    //check if it is in the XY plane for 2D case


    return 0;

    KRATOS_CATCH( "" );
  }


  void SpatialLagrangianElement::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    int IntMethod = int(mThisIntegrationMethod);
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.save("DeterminantF0",mDeterminantF0);
  }

  void SpatialLagrangianElement::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.load("DeterminantF0",mDeterminantF0);

  }



} // Namespace Kratos


