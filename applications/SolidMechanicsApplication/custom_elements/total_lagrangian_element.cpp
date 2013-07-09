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
#include "custom_elements/total_lagrangian_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "solid_mechanics_application.h"


namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  TotalLagrangianElement::TotalLagrangianElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
  {
    //DO NOT ADD DOFS HERE!!!
  }


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  TotalLagrangianElement::TotalLagrangianElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
  {
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
  }


  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  TotalLagrangianElement::TotalLagrangianElement( TotalLagrangianElement const& rOther)
    :Element(rOther)
    ,mThisIntegrationMethod(rOther.mThisIntegrationMethod)
    ,mConstitutiveLawVector(rOther.mConstitutiveLawVector)
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  TotalLagrangianElement&  TotalLagrangianElement::operator=(TotalLagrangianElement const& rOther)
  {
    Element::operator=(rOther);

    mThisIntegrationMethod = rOther.mThisIntegrationMethod;

    mConstitutiveLawVector.clear();
    mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

    for(unsigned int i=0; i<<mConstitutiveLawVector.size(); i++)
      {
        mConstitutiveLawVector[i] = rOther.mConstitutiveLawVector[i];
      }

    return *this;
  }


  //*********************************OPERATIONS*****************************************
  //************************************************************************************

  Element::Pointer TotalLagrangianElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
  {
    return Element::Pointer( new TotalLagrangianElement( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
  }


  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  TotalLagrangianElement::~TotalLagrangianElement()
  {
  }


  //************* GETTING METHODS
  //************************************************************************************
  //************************************************************************************

  TotalLagrangianElement::IntegrationMethod TotalLagrangianElement::GetIntegrationMethod() const
  {
    return mThisIntegrationMethod;
  }

  //************************************************************************************
  //************************************************************************************

  void TotalLagrangianElement::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo )
  {
    ElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
      {
	ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
	ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

	if ( GetGeometry().WorkingSpaceDimension() == 3 )
	  {
	    ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
	  }
      }
  }


  //************************************************************************************
  //************************************************************************************

  void TotalLagrangianElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
  {
    int number_of_nodes = GetGeometry().size();
    int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int dim2 = number_of_nodes * dim;

    if ( rResult.size() != dim2 )
      rResult.resize( dim2, false );

    for ( int i = 0; i < number_of_nodes; i++ )
      {
	int index = i * dim;
	rResult[index] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
	rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

	if ( dim == 3 )
	  rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
      }

  }

  //*********************************DISPLACEMENT***************************************
  //************************************************************************************

  void TotalLagrangianElement::GetValuesVector( Vector& values, int Step )
  {
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * dim;

    if ( values.size() != MatSize ) values.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	unsigned int index = i * dim;
	values[index] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
	values[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

	if ( dim == 3 )
	  values[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
      }
  }


  //************************************VELOCITY****************************************
  //************************************************************************************

  void TotalLagrangianElement::GetFirstDerivativesVector( Vector& values, int Step )
  {
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * dim;

    if ( values.size() != MatSize ) values.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	unsigned int index = i * dim;
	values[index] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
	values[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

	if ( dim == 3 )
	  values[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
      }
  }

  //*********************************ACCELERATION***************************************
  //************************************************************************************

  void TotalLagrangianElement::GetSecondDerivativesVector( Vector& values, int Step )
  {
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * dim;

    if ( values.size() != MatSize ) values.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	unsigned int index = i * dim;
	values[index] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
	values[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

	if ( dim == 3 )
	  values[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
      }
  }


  //*********************************SET DOUBLE VALUE***********************************
  //************************************************************************************

  void TotalLagrangianElement::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
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

  void TotalLagrangianElement::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, 
							    std::vector<Vector>& rValues, 
							    const ProcessInfo& rCurrentProcessInfo )
  {

    for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
      {
	mConstitutiveLawVector[PointNumber]->SetValue( rVariable,
						       rValues[PointNumber], rCurrentProcessInfo );
      }

  }

  //*********************************SET MATRIX VALUE***********************************
  //************************************************************************************


  void TotalLagrangianElement::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, 
							    std::vector<Matrix>& rValues, 
							    const ProcessInfo& rCurrentProcessInfo )
  {
    for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
      {
	mConstitutiveLawVector[PointNumber]->SetValue( rVariable,
						       rValues[PointNumber], rCurrentProcessInfo );
      }

  }

  //*********************************GET DOUBLE VALUE***********************************
  //************************************************************************************


  void TotalLagrangianElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
							    std::vector<double>& rValues,
							    const ProcessInfo& rCurrentProcessInfo )
  {
    if ( rValues.size() != GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
      rValues.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(), false );

    for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
      rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
  }


  //**********************************GET VECTOR VALUE**********************************
  //************************************************************************************


  void TotalLagrangianElement::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, 
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

  void TotalLagrangianElement::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
							    std::vector<Matrix>& rValues, 
							    const ProcessInfo& rCurrentProcessInfo )
  {

	CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
  }


  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************


  void TotalLagrangianElement::Initialize()
  {
    KRATOS_TRY

      const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //Constitutive Law initialisation
    if ( mConstitutiveLawVector.size() != integration_points.size() )
      {
	mConstitutiveLawVector.resize( integration_points.size() );
      }

    //Material initialisation
    InitializeMaterial();


    //resizing jacobian inverses containers
    mInvJ0.resize( integration_points.size() );
    mDetJ0.resize( integration_points.size(), false );


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

  void TotalLagrangianElement::SetStandardParameters(Standard& rVariables,
						     ConstitutiveLaw::Parameters& rValues,
						     const int & rPointNumber)
  {
    rVariables.detF  = MathUtils<double>::Det(rVariables.F);
    
    if(rVariables.detF<0)
      KRATOS_ERROR(std::invalid_argument,"DeterminantF < 0","");

    //in this element F = F0, then the F0 is set to the identity for coherence in the constitutive law.
    double detF0 = 1;
    Matrix F0    = identity_matrix<double> (rVariables.F.size1());
    rValues.SetDeterminantF0(detF0);
    rValues.SetDeformationGradientF0(F0);

    rValues.SetDeterminantF(rVariables.detF);
    rValues.SetDeformationGradientF(rVariables.F);
    rValues.SetStrainVector(rVariables.StrainVector);
    rValues.SetStressVector(rVariables.StressVector);
    rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rValues.SetShapeFunctionsDevivatives(rVariables.DN_DX);

    Flags &Options=rValues.GetOptions();
    Options.Reset(ConstitutiveLaw::COMPUTE_STRAIN); //to use the already computed strain in the linear_elastic_law
   
    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );
    rVariables.N=row(Ncontainer , rPointNumber);
    rValues.SetShapeFunctionsValues(rVariables.N);
  }

  //************************************************************************************
  //************************************************************************************

  void TotalLagrangianElement::InitializeStandardVariables (Standard & rVariables, const ProcessInfo& rCurrentProcessInfo)
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
  
    rVariables.ConstitutiveMatrix.resize( StrainSize, StrainSize );
  
    rVariables.StrainVector.resize( StrainSize );
  
    rVariables.StressVector.resize( StrainSize );
  
    rVariables.DN_DX.resize( number_of_nodes, dimension );
  
  }


  //************************************************************************************
  //************************************************************************************

  void TotalLagrangianElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
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

  void TotalLagrangianElement::CalculateElementalSystem( MatrixType& rLeftHandSideMatrix,
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

    //KRATOS_WATCH(J)

    //auxiliary terms
    Vector BodyForce;

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {

	//COMPUTE kinematics B,F,DN_DX ...
        CalculateKinematics(Variables,PointNumber);

        //set standart parameters
        SetStandardParameters(Variables,Values,PointNumber);

	//compute stresses and constitutive parameters
	mConstitutiveLawVector[PointNumber]->CalculateMaterialResponsePK2(Values);


	//calculating weights for integration on the "reference configuration"
	double IntegrationWeight = integration_points[PointNumber].Weight() * mDetJ0[PointNumber];

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

	    // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
	    CalculateAndAddExternalForces(Variables.N, rCurrentProcessInfo, BodyForce, rRightHandSideVector, IntegrationWeight );

	    // operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
	    CalculateAndAddInternalForces(Variables.B,Variables.StressVector,rRightHandSideVector,IntegrationWeight);
	    
	  }
      }


    KRATOS_CATCH( "" )
      }

  //************************************************************************************
  //************************************************************************************

  void TotalLagrangianElement::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
  {
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateElementalSystem( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
  }

  //************************************************************************************
  //************************************************************************************

  void TotalLagrangianElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
  {
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateElementalSystem( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );


  }


  ////************************************************************************************
  ////************************************************************************************

  void TotalLagrangianElement::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
  {
    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
      mConstitutiveLawVector[i]->InitializeSolutionStep( GetProperties(),
							 GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ),
							 CurrentProcessInfo );
  }

  ////************************************************************************************
  ////************************************************************************************

  void TotalLagrangianElement::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
  {
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
        mConstitutiveLawVector[PointNumber]->FinalizeMaterialResponsePK2 (Values);

	
	mConstitutiveLawVector[PointNumber]->FinalizeSolutionStep( GetProperties(),
							 GetGeometry(),
							 Variables.N,
							 CurrentProcessInfo );
      }
  }

  //************************************************************************************
  //************************************************************************************

  void TotalLagrangianElement::InitializeMaterial()
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
	KRATOS_ERROR( std::logic_error, "a constitutive law needs to be specified for the element with ID ", this->Id() )
	  KRATOS_CATCH( "" )
	  }


  //************************************************************************************
  //************************************************************************************

  void TotalLagrangianElement::ResetConstitutiveLaw()
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

  inline void TotalLagrangianElement::CalculateAndAddExternalForces(const Vector& rN,
								   const ProcessInfo& rCurrentProcessInfo,
								   Vector& rBodyForce,
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
	
	Fext = 0;
	for ( unsigned int j = 0; j < dimension; j++ )
	  {
	    Fext = rIntegrationWeight * rN[i] * rBodyForce[j];
	    rRightHandSideVector[index + j] += Fext;
	    ExternalForce[j] +=Fext;
	  }
      }

    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************

  inline void TotalLagrangianElement::CalculateAndAddInternalForces(Matrix & rB,
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
 //************************************************************************************
  //************************************************************************************

  void TotalLagrangianElement::CalculateAndAddKm(MatrixType& rK,
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

  void TotalLagrangianElement::CalculateAndAddKg(MatrixType& rK,
						   Matrix& rDN_DX,
						   Vector& rStressVector,
						   double& rIntegrationWeight
						   )

  {
    KRATOS_TRY
 
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rStressVector );
    Matrix ReducedKg = prod( rDN_DX, rIntegrationWeight * Matrix( prod( StressTensor, trans( rDN_DX ) ) ) ); //to be optimized
    MathUtils<double>::ExpandAndAddReducedMatrix( rK, ReducedKg, dimension );

    KRATOS_CATCH( "" )
      }



  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************


  void TotalLagrangianElement::CalculateKinematics(Standard& rVariables,
						     const double& rPointNumber)

  {
    KRATOS_TRY
      
    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Parent to reference configuration
    Matrix J ( dimension , dimension);
    J = GetGeometry().Jacobian( J, rPointNumber , mThisIntegrationMethod );
    
    //Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], mInvJ0[rPointNumber] );

    //Deformation Gradient F0
    noalias( rVariables.F ) = prod( J, mInvJ0[rPointNumber] );

    //Compute the deformation matrix B
    CalculateDeformationMatrix(rVariables.B,rVariables.F,rVariables.DN_DX);

    //Right Cauchy-Green Calculation
    Matrix C ( dimension, dimension );
    noalias( C ) = prod( trans( rVariables.F ), rVariables.F );
    
    //Compute Green-Lagrange Strain 
    CalculateGreenLagrangeStrain( C, rVariables.StrainVector );
    
    //Correction:
    //DecimalCorrection( rVariables.StrainVector );


    KRATOS_CATCH( "" )
      }

  //************************************************************************************
  //************************************************************************************

  void TotalLagrangianElement::CalculateGreenLagrangeStrain(const Matrix& C,
							    Vector& StrainVector )
  {
    KRATOS_TRY
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();
 
    if ( dimension == 2 )
      {
	if ( StrainVector.size() != 3 ) StrainVector.resize( 3, false );

	StrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

	StrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

	StrainVector[2] = C( 0, 1 );
      }

    if ( dimension == 3 )
      {
	if ( StrainVector.size() != 6 ) StrainVector.resize( 6, false );

	StrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

	StrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

	StrainVector[2] = 0.5 * ( C( 2, 2 ) - 1.00 );

	StrainVector[3] = C( 0, 1 ); // xy

	StrainVector[4] = C( 1, 2 ); // yz

	StrainVector[5] = C( 0, 2 ); // xz
      }

    KRATOS_CATCH( "" )
      }

  //************************************************************************************
  //************************************************************************************

  void TotalLagrangianElement::CalculateDeformationMatrix(Matrix& B,
							  Matrix& F,
							  Matrix& DN_DX)
  {
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
 
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	unsigned int index = dimension * i;

	if ( dimension == 2 )
	  {
	    B( 0, index + 0 ) = F( 0, 0 ) * DN_DX( i, 0 );
	    B( 0, index + 1 ) = F( 1, 0 ) * DN_DX( i, 0 );
	    B( 1, index + 0 ) = F( 0, 1 ) * DN_DX( i, 1 );
	    B( 1, index + 1 ) = F( 1, 1 ) * DN_DX( i, 1 );
	    B( 2, index + 0 ) = F( 0, 0 ) * DN_DX( i, 1 ) + F( 0, 1 ) * DN_DX( i, 0 );
	    B( 2, index + 1 ) = F( 1, 0 ) * DN_DX( i, 1 ) + F( 1, 1 ) * DN_DX( i, 0 );
	  }
	else
	  {
	    B( 0, index + 0 ) = F( 0, 0 ) * DN_DX( i, 0 );
	    B( 0, index + 1 ) = F( 1, 0 ) * DN_DX( i, 0 );
	    B( 0, index + 2 ) = F( 2, 0 ) * DN_DX( i, 0 );
	    B( 1, index + 0 ) = F( 0, 1 ) * DN_DX( i, 1 );
	    B( 1, index + 1 ) = F( 1, 1 ) * DN_DX( i, 1 );
	    B( 1, index + 2 ) = F( 2, 1 ) * DN_DX( i, 1 );
	    B( 2, index + 0 ) = F( 0, 2 ) * DN_DX( i, 2 );
	    B( 2, index + 1 ) = F( 1, 2 ) * DN_DX( i, 2 );
	    B( 2, index + 2 ) = F( 2, 2 ) * DN_DX( i, 2 );
	    B( 3, index + 0 ) = F( 0, 0 ) * DN_DX( i, 1 ) + F( 0, 1 ) * DN_DX( i, 0 );
	    B( 3, index + 1 ) = F( 1, 0 ) * DN_DX( i, 1 ) + F( 1, 1 ) * DN_DX( i, 0 );
	    B( 3, index + 2 ) = F( 2, 0 ) * DN_DX( i, 1 ) + F( 2, 1 ) * DN_DX( i, 0 );
	    B( 4, index + 0 ) = F( 0, 1 ) * DN_DX( i, 2 ) + F( 0, 2 ) * DN_DX( i, 1 );
	    B( 4, index + 1 ) = F( 1, 1 ) * DN_DX( i, 2 ) + F( 1, 2 ) * DN_DX( i, 1 );
	    B( 4, index + 2 ) = F( 2, 1 ) * DN_DX( i, 2 ) + F( 2, 2 ) * DN_DX( i, 1 );
	    B( 5, index + 0 ) = F( 0, 2 ) * DN_DX( i, 0 ) + F( 0, 0 ) * DN_DX( i, 2 );
	    B( 5, index + 1 ) = F( 1, 2 ) * DN_DX( i, 0 ) + F( 1, 0 ) * DN_DX( i, 2 );
	    B( 5, index + 2 ) = F( 2, 2 ) * DN_DX( i, 0 ) + F( 2, 0 ) * DN_DX( i, 2 );
	  }

      }

    KRATOS_CATCH( "" )
      }



  //************************************************************************************
  //************************************************************************************

  void TotalLagrangianElement::MassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

      //lumped
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int NumberOfNodes = GetGeometry().size();
    unsigned int MatSize = dimension * NumberOfNodes;

    if ( rMassMatrix.size1() != MatSize )
      rMassMatrix.resize( MatSize, MatSize, false );

    rMassMatrix = ZeroMatrix( MatSize, MatSize );

    double TotalMass = mTotalDomainInitialSize * GetProperties()[DENSITY];

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

    KRATOS_CATCH( "" )
      }

  //************************************************************************************
  //************************************************************************************

  void TotalLagrangianElement::DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY
      unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dim;

    if ( rDampMatrix.size1() != MatSize )
      rDampMatrix.resize( MatSize, MatSize, false );

    noalias( rDampMatrix ) = ZeroMatrix( MatSize, MatSize );

    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************

  void TotalLagrangianElement::CalculateOnIntegrationPoints( const Variable<double>& rVariable, Vector& rOutput, const ProcessInfo& rCurrentProcessInfo )
  {
    if ( rOutput.size() != GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
      rOutput.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(), false );


    if ( rVariable == VON_MISES_STRESS )
      {
	double StressSize = 3;
	if ( GetGeometry().WorkingSpaceDimension() == 2 )
	  StressSize = 3;
	else
	  StressSize = 6;
      
        Vector StressVector( StressSize );

        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
	  {
            StressVector = mConstitutiveLawVector[ii]->GetValue(CAUCHY_STRESS_VECTOR,StressVector);

	    ComparisonUtils EquivalentStress;
	    rOutput[ii] =  EquivalentStress.CalculateVonMises(StressVector);
	  }
      }
    else{
      
      for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
        rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rOutput[ii] );
    }

  }

  //************************************************************************************
  //************************************************************************************

  void TotalLagrangianElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo )
  {
    unsigned int StrainSize;

    if ( GetGeometry().WorkingSpaceDimension() == 2 ) 
      StrainSize = 3;
    else 
      StrainSize = 6;

    
    if ( rVariable == CAUCHY_STRESS_VECTOR )
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
	    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponsePK2(Values);

	    Variables.StressVector=Values.GetStressVector(Variables.StressVector);

	    mConstitutiveLawVector[PointNumber]->TransformStresses(Variables.StressVector,Variables.F,Variables.detF,ConstitutiveLaw::StressMeasure_PK2,ConstitutiveLaw::StressMeasure_Cauchy);

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

  }

  //************************************************************************************
  //************************************************************************************

  void TotalLagrangianElement::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo )
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

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

 
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {

	//COMPUTE kinematics B,F,DN_DX ...
        CalculateKinematics(Variables,PointNumber);

	if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
	  {

	    if ( rOutput[PointNumber].size2() != Variables.StrainVector.size() )
	      rOutput[PointNumber].resize( 1, Variables.StrainVector.size(), false );

	    for ( unsigned int ii = 0; ii < Variables.StrainVector.size(); ii++ )
	      rOutput[PointNumber]( 0, ii ) = Variables.StrainVector[ii];
	  }
	else if ( rVariable == PK2_STRESS_TENSOR )
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

		Variables.StressVector=Values.GetStressVector(Variables.StressVector);

		if ( rOutput[PointNumber].size2() != Variables.StrainVector.size() )
		  rOutput[PointNumber].resize( 1, Variables.StrainVector.size(), false );
	    
		for ( unsigned int ii = 0; ii < Variables.StrainVector.size(); ii++ )
		  {
		    rOutput[PointNumber]( 0, ii ) = Variables.StressVector[ii];
		  }

	      }

	  }
	else if ( rVariable == CAUCHY_STRESS_TENSOR )
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

		Variables.StressVector=Values.GetStressVector(Variables.StressVector);

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

  void TotalLagrangianElement::Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo )
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
    Output = 2.0 / wmax;
    //KRATOS_WATCH(Output)

  }

 //*************************DECIMAL CORRECTION OF STRAINS******************************
  //************************************************************************************

  void TotalLagrangianElement::DecimalCorrection(Vector& rStrainVector)
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


  //************************************************************************************
  //************************************************************************************
  /**
   * This function provides the place to perform checks on the completeness of the input.
   * It is designed to be called only once (or anyway, not often) typically at the beginning
   * of the calculations, so to verify that nothing is missing from the input
   * or that no common error is found.
   * @param rCurrentProcessInfo
   */
  int  TotalLagrangianElement::Check( const ProcessInfo& rCurrentProcessInfo )
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


  void TotalLagrangianElement::save( Serializer& rSerializer ) const
  {
    //  std::cout << "Saving the TotalLagrangianElement #" << Id() << std::endl;
    rSerializer.save( "Name", "TotalLagrangianElement" );
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
  }

  void TotalLagrangianElement::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    //  std::cout << "Loading the TotalLagrangianElement #" << Id() << std::endl;
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
  }


} // Namespace Kratos


