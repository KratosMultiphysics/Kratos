//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/beam_elements/large_displacement_beam_element.hpp"

#include "solid_mechanics_application_variables.h"


//NOTE:
//ALL TRANSFORMATION MATRICES Q ARE GIVEN IN COLUMNS
//the transformation of coordinates is done via (e'= QT * e) and (e = Q * e')
//vector value :  v' = QT * v
//matrix value :  A' = QT * A * Q
//these transformations change the reference axes
//they do not rotate quantities


namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  LargeDisplacementBeamElement::LargeDisplacementBeamElement(IndexType NewId,GeometryType::Pointer pGeometry)
    : BeamElement(NewId, pGeometry)
  {
  }

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************


  LargeDisplacementBeamElement::LargeDisplacementBeamElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : BeamElement(NewId, pGeometry, pProperties)
  {
    mReducedIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
    mFullIntegrationMethod = mReducedIntegrationMethod;
    this->IncreaseIntegrationMethod(mFullIntegrationMethod,1);
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  LargeDisplacementBeamElement::LargeDisplacementBeamElement(LargeDisplacementBeamElement const& rOther)
    :BeamElement(rOther)
    ,mReducedIntegrationMethod(rOther.mReducedIntegrationMethod)
    ,mFullIntegrationMethod(rOther.mFullIntegrationMethod)
    ,mIterationCounter(rOther.mIterationCounter)
    ,mInvJ0(rOther.mInvJ0)
    ,mCurrentCurvatureVectors(rOther.mCurrentCurvatureVectors)
    ,mPreviousCurvatureVectors(rOther.mPreviousCurvatureVectors)
    ,mFrameQuaternionsReduced(rOther.mFrameQuaternionsReduced)
    ,mFrameQuaternionsFull(rOther.mFrameQuaternionsFull)
  {
  }

  //*********************************OPERATIONS*****************************************
  //************************************************************************************

  Element::Pointer LargeDisplacementBeamElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
  {
    return Kratos::make_shared<LargeDisplacementBeamElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  LargeDisplacementBeamElement::~LargeDisplacementBeamElement()
  {
  }


  //************************************************************************************
  //************************************************************************************

  LargeDisplacementBeamElement::IntegrationMethod LargeDisplacementBeamElement::GetReducedIntegrationMethod() const
  {
    return mReducedIntegrationMethod;
  }

  //************************************************************************************
  //************************************************************************************

  LargeDisplacementBeamElement::IntegrationMethod LargeDisplacementBeamElement::GetFullIntegrationMethod() const
  {
    return mFullIntegrationMethod;
  }

  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::Initialize()
  {
    KRATOS_TRY

    BeamElement::Initialize();

    //------------- REDUCED QUADRATURE INTEGRATION

    IntegrationMethod ReducedIntegrationMethod = this->GetReducedIntegrationMethod();

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( ReducedIntegrationMethod );

    //Curvature Initialization
    if ( mCurrentCurvatureVectors.size() != integration_points.size() )
      {
        mCurrentCurvatureVectors.resize( integration_points.size() );
      }

    for ( unsigned int i = 0; i < mCurrentCurvatureVectors.size(); i++ )
      {
	mCurrentCurvatureVectors[i].resize(3,false);
	noalias(mCurrentCurvatureVectors[i]) = ZeroVector(3);
      }


    if ( mPreviousCurvatureVectors.size() != integration_points.size() )
      {
        mPreviousCurvatureVectors.resize( integration_points.size() );
      }

    for ( unsigned int i = 0; i < mPreviousCurvatureVectors.size(); i++ )
      {
	mPreviousCurvatureVectors[i].resize(3,false);
	noalias(mPreviousCurvatureVectors[i]) = ZeroVector(3);
      }

    //-------------

    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    Matrix Identity(dimension,dimension);
    noalias(Identity) = IdentityMatrix(dimension);

    //Local Quaternions Initialization
    if ( mFrameQuaternionsReduced.size() != integration_points.size() )
      {
        mFrameQuaternionsReduced.resize( integration_points.size() );
      }


    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {
	mFrameQuaternionsReduced[PointNumber]  = QuaternionType::FromRotationMatrix( Identity );
      }


    //-------------


    //Compute jacobian inverses and set the domain initial size:
    GeometryType::JacobiansType J0;
    J0 = GetGeometry().Jacobian( J0, ReducedIntegrationMethod );

    //calculating the inverse J0
    Vector Jacobian(3);
    noalias(Jacobian) = ZeroVector(3);

    //Calculating the jacobian [dx_0/d£]
    for( SizeType i=0; i<dimension; i++)
      {
	Jacobian[i] = J0[0](i,0);
      }

    mInvJ0 = 1.0/norm_2(Jacobian);


    //------------- FULL QUADRATURE INTEGRATION

    IntegrationMethod FullIntegrationMethod = this->GetFullIntegrationMethod();

    const GeometryType::IntegrationPointsArrayType& integration_points_full = GetGeometry().IntegrationPoints( FullIntegrationMethod );


    //Local Quaternions Initialization
    if ( mFrameQuaternionsFull.size() != integration_points_full.size() )
      {
        mFrameQuaternionsFull.resize( integration_points_full.size() );
      }

    for ( unsigned int PointNumber = 0; PointNumber < integration_points_full.size(); PointNumber++ )
      {
       	mFrameQuaternionsFull[PointNumber]  = QuaternionType::FromRotationMatrix( Identity );
      }

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    BeamElement::InitializeSolutionStep(rCurrentProcessInfo);

    //predict is done after initialize solution step -- perform this operations in initialize to write results correctly --


    // Update Frame Quaternions:

    //create and initialize element variables:
    ElementDataType Variables;
    this->InitializeElementData(Variables,rCurrentProcessInfo);

    IntegrationMethod ThisIntegrationMethod = mThisIntegrationMethod;

    //reduced quadrature integration:
    mThisIntegrationMethod = mReducedIntegrationMethod;

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //reading shape functions
    Variables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    //get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = Variables.GetShapeFunctions();

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {
	//set shape functions values for this integration point
	noalias(Variables.N) = matrix_row<const Matrix>( Ncontainer, PointNumber);

	//compute local to global frame
	this->CalculateFrameMapping( Variables, PointNumber );

	//update frame quaternion
	mFrameQuaternionsReduced[PointNumber] = QuaternionType::FromRotationMatrix(Variables.CurrentRotationMatrix);
      }


    //full quadrature integration:
    mThisIntegrationMethod = mFullIntegrationMethod;

    const GeometryType::IntegrationPointsArrayType& full_integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //reading shape functions
    Variables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    //get the shape functions for the order of the integration method [N]
    const Matrix& NFcontainer = Variables.GetShapeFunctions();

    for ( unsigned int PointNumber = 0; PointNumber < full_integration_points.size(); PointNumber++ )
      {
	//set shape functions values for this integration point
	noalias(Variables.N) = matrix_row<const Matrix>( NFcontainer, PointNumber);

	//compute local to global frame
	this->CalculateFrameMapping( Variables, PointNumber );

	//update frame quaternion
	mFrameQuaternionsFull[PointNumber] = QuaternionType::FromRotationMatrix(Variables.CurrentRotationMatrix);
      }

    mThisIntegrationMethod = ThisIntegrationMethod;

    mIterationCounter = 0;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY


    KRATOS_CATCH( "" )
  }

  ////************************************************************************************
  ////************************************************************************************

  void LargeDisplacementBeamElement::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    mIterationCounter++;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    // Update Curvature Vectors:

    for ( unsigned int i = 0; i < mCurrentCurvatureVectors.size(); i++ )
      {
	mPreviousCurvatureVectors[i] = mCurrentCurvatureVectors[i] ;
      }

    BeamElement::FinalizeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateCurrentCurvature(ElementDataType& rVariables, const Variable<array_1d<double, 3 > >& rVariable)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();

    Vector CurrentStepRotationDerivativesVector(3);
    noalias(CurrentStepRotationDerivativesVector) = ZeroVector(3);
    Vector CurrentStepRotationVector(3);
    noalias(CurrentStepRotationVector) = ZeroVector(3);
    Vector CurrentValueVector(3);
    noalias(CurrentValueVector) = ZeroVector(3);

    //jelenic:crisfield [j-c]:start
    QuaternionType QuaternionValue;
    std::vector<QuaternionType> CurrentStepRotationQuaternions;
    //------------------[j-c]:end

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
    	//Current Rotation Derivatives
    	CurrentValueVector = GetNodalCurrentValue( rVariable, CurrentValueVector, i );
    	CurrentValueVector = this->MapToInitialLocalFrame( CurrentValueVector, rVariables.PointNumber );

	CurrentStepRotationVector += rVariables.N[i] * (CurrentValueVector);
    	CurrentStepRotationDerivativesVector += rVariables.DN_DX(i,0) * ( CurrentValueVector );

	//------------------[j-c]:start
	QuaternionValue = QuaternionType::FromRotationVector(CurrentValueVector);
	CurrentStepRotationQuaternions.push_back(QuaternionValue);
	//------------------[j-c]:end
      }

    //------------------[j-c]:start
    QuaternionValue = (CurrentStepRotationQuaternions.front()).conjugate() * (CurrentStepRotationQuaternions.back());
    QuaternionValue.ToRotationVector(CurrentValueVector);
    CurrentValueVector *= rVariables.N[0];
    QuaternionValue = QuaternionType::FromRotationVector(CurrentValueVector);
    QuaternionValue = (CurrentStepRotationQuaternions.front()) * QuaternionValue;
    QuaternionValue.ToRotationVector( CurrentStepRotationVector );
    //------------------[j-c]:end


    QuaternionType CurrentStepRotationVectorQuaternion = QuaternionType::FromRotationVector( CurrentStepRotationVector );
    rVariables.PreviousCurvatureVector = BeamMathUtilsType::MapToReferenceLocalFrame( CurrentStepRotationVectorQuaternion, rVariables.PreviousCurvatureVector);

    double CurrentStepRotationVectorModulus = norm_2(CurrentStepRotationVector);

    double CurrentStepRotationVectorSinusI   = 0;
    double CurrentStepRotationVectorCosinusI = 0;

    //------------------[option 0]:start
    // double CurrentStepRotationVectorSinusII  = 0;
    // Vector EquationVectorPartA(3); noalias(EquationVectorPartA) = ZeroVector(3);
    // Vector EquationVectorPartB(3); noalias(EquationVectorPartB) = ZeroVector(3);
    // Vector EquationVectorPartC(3); noalias(EquationVectorPartC) = ZeroVector(3);
    //------------------[option 0]:end

    Matrix EquationMatrixPartA(3,3);
    noalias(EquationMatrixPartA) = ZeroMatrix(3,3);
    Matrix EquationMatrixPartB(3,3);
    noalias(EquationMatrixPartB) = ZeroMatrix(3,3);
    Matrix EquationMatrixPartC(3,3);
    noalias(EquationMatrixPartC) = ZeroMatrix(3,3);

    Matrix DiagonalMatrix(3,3);
    noalias(DiagonalMatrix) = IdentityMatrix(3);

    Matrix SkewSymCurrentStepRotation(3,3);
    noalias(SkewSymCurrentStepRotation) = ZeroMatrix(3,3);
    BeamMathUtilsType::VectorToSkewSymmetricTensor(CurrentStepRotationVector, SkewSymCurrentStepRotation);

    if( CurrentStepRotationVectorModulus != 0.00 ){

      CurrentStepRotationVectorSinusI   = std::sin(CurrentStepRotationVectorModulus) / CurrentStepRotationVectorModulus;
      CurrentStepRotationVectorCosinusI = std::cos(CurrentStepRotationVectorModulus);

      //------------------[option 0]:start
      // CurrentStepRotationVectorSinusII  = std::sin(0.5 * CurrentStepRotationVectorModulus) / (0.5 * CurrentStepRotationVectorModulus);

      // EquationVectorPartA = (CurrentStepRotationVectorSinusI) * CurrentStepRotationDerivativesVector;
      // EquationVectorPartB = ( 1.0 - CurrentStepRotationVectorSinusI ) * inner_prod( CurrentStepRotationVector, CurrentStepRotationDerivativesVector ) * ( 1.0 / ( CurrentStepRotationVectorModulus * CurrentStepRotationVectorModulus ) ) * CurrentStepRotationVector;
      // EquationVectorPartC =  0.5 * ( CurrentStepRotationVectorSinusII * CurrentStepRotationVectorSinusII ) * MathUtils<double>::CrossProduct(CurrentStepRotationVector, CurrentStepRotationDerivativesVector);

      //------------------[j-c]:start
      // EquationVectorPartC =  (1 - CurrentStepRotationVectorCosinusI) * ( 1.0 / ( CurrentStepRotationVectorModulus * CurrentStepRotationVectorModulus ) ) * MathUtils<double>::CrossProduct(CurrentStepRotationVector, CurrentStepRotationDerivativesVector);
      //------------------[j-c]:end

      //------------------[option 0]:end

      EquationMatrixPartA = (CurrentStepRotationVectorSinusI) * DiagonalMatrix;
      EquationMatrixPartB = (( 1.0 - CurrentStepRotationVectorSinusI ) / ( CurrentStepRotationVectorModulus * CurrentStepRotationVectorModulus )) * outer_prod(CurrentStepRotationVector, CurrentStepRotationVector);
      EquationMatrixPartC = ((1 - CurrentStepRotationVectorCosinusI) / ( CurrentStepRotationVectorModulus * CurrentStepRotationVectorModulus )) * SkewSymCurrentStepRotation;

    }

    //------------------[option 0]:start
    //Vector EquationVector = EquationVectorPartA + EquationVectorPartB + EquationVectorPartC;
    // rVariables.CurrentCurvatureVector  = rVariables.PreviousCurvatureVector;
    // rVariables.CurrentCurvatureVector += EquationVector;
    //------------------[option 0]:end

    Matrix EquationMatrix = EquationMatrixPartA + EquationMatrixPartB + EquationMatrixPartC;

    rVariables.CurrentCurvatureVector  = rVariables.PreviousCurvatureVector;
    rVariables.CurrentCurvatureVector += prod(EquationMatrix, CurrentStepRotationDerivativesVector);

    //------------------------

    KRATOS_CATCH( "" )
  }


  //*****************************************************************************
  //*****************************************************************************

  void LargeDisplacementBeamElement::MapToSpatialFrame(const ElementDataType& rVariables, Matrix& rVariable)
  {
    KRATOS_TRY

    BeamMathUtilsType::MapLocalToGlobal3D(rVariables.CurrentRotationMatrix, rVariable);

    KRATOS_CATCH( "" )
  }


  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateKinematics(ElementDataType& rVariables, const unsigned int& rPointNumber)
  {
    KRATOS_TRY

    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    const SizeType number_of_nodes  = GetGeometry().size();

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //point number
    rVariables.PointNumber = rPointNumber;

    //Set Shape Functions Values for this integration point
    noalias(rVariables.N) = matrix_row<const Matrix>( Ncontainer, rPointNumber);

    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();

    //TOTAL LAGRANGIAN
    //Compute cartesian derivatives [dN/dx_0]
    rVariables.DN_DX = mInvJ0 * DN_De[rPointNumber];
    rVariables.detJ  = 1.0/mInvJ0;

    //compute local to global frame
    this->CalculateFrameMapping(rVariables,rPointNumber);

    //Compute current centroid DISPLACEMENT DERIVATIVES:
    noalias(rVariables.CurrentAxisPositionDerivatives) = ZeroVector(dimension);

    Vector CurrentValueVector(3);
    noalias(CurrentValueVector) = ZeroVector(3);


    if( this->Is(BeamElement::FINALIZED_STEP) ){

      rVariables.DeltaPosition = this->CalculateDeltaPosition(rVariables.DeltaPosition);

      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  //A: Current Nodes Position
	  CurrentValueVector = GetGeometry()[i].Coordinates();

	  for ( SizeType j = 0; j < dimension; j++ )
	    {
	      CurrentValueVector[j] -= rVariables.DeltaPosition(i,j);
	    }

	  CurrentValueVector = this->MapToInitialLocalFrame( CurrentValueVector, rPointNumber );

	  //Current Frame Axis Position derivative
	  for( SizeType j = 0; j < dimension; j++ )
	    {
	      rVariables.CurrentAxisPositionDerivatives[j] +=  rVariables.DN_DX(i,0) * ( CurrentValueVector[j] );
	    }
	}
    }
    else{

      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  //A: Current Nodes Position
	  CurrentValueVector = GetGeometry()[i].Coordinates();
	  CurrentValueVector = this->MapToInitialLocalFrame( CurrentValueVector, rPointNumber );

	  //Current Frame Axis Position derivative
	  for( SizeType j = 0; j < dimension; j++ )
	    {
	      rVariables.CurrentAxisPositionDerivatives[j] +=  rVariables.DN_DX(i,0) * ( CurrentValueVector[j] );
	    }
	}
    }


    //Compute current CURVATURES
    if( this->Is(BeamElement::FINALIZED_STEP) ){

      rVariables.CurrentCurvatureVector = mCurrentCurvatureVectors[rPointNumber];

    }
    else{

      // if( mIterationCounter == 0 ){

	rVariables.PreviousCurvatureVector = mPreviousCurvatureVectors[rPointNumber];
      	rVariables.CurrentCurvatureVector = mPreviousCurvatureVectors[rPointNumber];

	this->CalculateCurrentCurvature(rVariables, STEP_ROTATION);

      // }
      // else{

      // 	rVariables.PreviousCurvatureVector = mCurrentCurvatureVectors[rPointNumber];
      // 	rVariables.CurrentCurvatureVector = mCurrentCurvatureVectors[rPointNumber];

      // 	this->CalculateCurrentCurvature(rVariables, DELTA_ROTATION);
      // }

    }

    KRATOS_CATCH( "" )
  }


  //*************************COMPUTE FRAME MAPPING POSITION*****************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateFrameMapping(ElementDataType& rVariables,const unsigned int& rPointNumber)
  {
    KRATOS_TRY

    if( mThisIntegrationMethod == mReducedIntegrationMethod ){
      mFrameQuaternionsReduced[rPointNumber].ToRotationMatrix(rVariables.PreviousRotationMatrix);
    }

    if( mThisIntegrationMethod == mFullIntegrationMethod ){
      mFrameQuaternionsFull[rPointNumber].ToRotationMatrix(rVariables.PreviousRotationMatrix);
    }

    //option frame 1:
    Vector CurrentStepRotationVector(3);
    noalias(CurrentStepRotationVector) = ZeroVector(3);
    this->GetLocalCurrentValue(STEP_ROTATION, CurrentStepRotationVector, rVariables.N);

    Matrix ExponentialRotationMatrix(3,3);
    noalias(ExponentialRotationMatrix) = ZeroMatrix(3,3);

    BeamMathUtilsType::ExponentialTransform( CurrentStepRotationVector, ExponentialRotationMatrix );

    rVariables.CurrentRotationMatrix = prod(ExponentialRotationMatrix, rVariables.PreviousRotationMatrix);


    //option frame 2:

    // const SizeType number_of_nodes  = GetGeometry().size();
    // Vector CurrentValueVector(3);
    // noalias(CurrentValueVector) = ZeroVector(3);

    // std::vector<QuaternionType> NodeQuaternions;
    // QuaternionType QuaternionValue;

    // //strains due to displacements and rotations
    // for ( SizeType i = 0; i < number_of_nodes; i++ )
    //   {

    //     //Rotations
    //     CurrentValueVector = GetNodalCurrentValue( ROTATION, CurrentValueVector, i );

    //     //Current Frame is the Local Frame
    //     CurrentValueVector = this->MapToInitialLocalFrame( CurrentValueVector, rPointNumber );

    //     QuaternionValue = QuaternionType::FromRotationVector(CurrentValueVector);

    //     NodeQuaternions.push_back(QuaternionValue);

    //   }

    // QuaternionValue = (NodeQuaternions.front()).conjugate() * (NodeQuaternions.back());

    // QuaternionValue.ToRotationVector(CurrentValueVector);

    // CurrentValueVector *= rVariables.N[0];

    // QuaternionValue = QuaternionType::FromRotationVector(CurrentValueVector);

    // QuaternionValue = (NodeQuaternions.front()) * QuaternionValue;

    // // Current Rotation Matrix
    // QuaternionValue.ToRotationMatrix(rVariables.CurrentRotationMatrix);


    KRATOS_CATCH( "" )
  }


  //*********************************SET STRAIN VARIABLES*******************************
  //************************************************************************************

  void LargeDisplacementBeamElement::UpdateStrainVariables(ElementDataType& rVariables, const unsigned int& rPointNumber)
  {
    KRATOS_TRY

    mCurrentCurvatureVectors[rPointNumber] = rVariables.CurrentCurvatureVector;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************


  void LargeDisplacementBeamElement::CalculateDynamicSystem( LocalSystemComponents& rLocalSystem,
							     ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    //create and initialize element variables:
    ElementDataType Variables;

    IntegrationMethod ThisIntegrationMethod = mThisIntegrationMethod;
    //full quadrature integration:
    mThisIntegrationMethod = mFullIntegrationMethod;

    this->InitializeElementData(Variables,rCurrentProcessInfo);

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    // initialize variables short version:

    //Compute Section Properties:
    this->CalculateSectionProperties(Variables.Section);
    Variables.Length = GetGeometry().Length();

    //Set equilibrium point initial:0/mid:0.5/final:1
    if( rCurrentProcessInfo.Has(EQUILIBRIUM_POINT) )
      Variables.Alpha = rCurrentProcessInfo[EQUILIBRIUM_POINT];
    else
      Variables.Alpha = 1;

    //reading shape functions
    Variables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    //get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = Variables.GetShapeFunctions();

    // initialize variables short version;

    MatrixType  LocalLeftHandSideMatrix;
    VectorType  LocalRightHandSideVector;
    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( LocalLeftHandSideMatrix, LocalRightHandSideVector, rLocalSystem.CalculationFlags );

    //(in fact is the two nodes beam element, full quadrature 2 integration points)
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {
	//integration point number
	Variables.PointNumber = PointNumber;

	//set shape functions values for this integration point
	noalias(Variables.N) = matrix_row<const Matrix>( Ncontainer, PointNumber);

	//compute local to global frame
	this->CalculateFrameMapping( Variables, PointNumber );

	Variables.detJ = 1.0/mInvJ0;

	double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
	IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );

	if ( rLocalSystem.CalculationFlags.Is(BeamElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
	  {
	    LocalLeftHandSideMatrix.clear();

    	    this->CalculateAndAddInertiaLHS( LocalLeftHandSideMatrix, Variables, rCurrentProcessInfo, IntegrationWeight ); // (R_N+1, R_N)

	    BeamMathUtilsType::MapLocalToGlobal3D(mInitialLocalQuaternion,LocalLeftHandSideMatrix);

	    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();
	    rLeftHandSideMatrix += LocalLeftHandSideMatrix;

	    //std::cout<<"["<<this->Id()<<"] Beam RotatedDynamic rLeftHandSideMatrix "<<rLeftHandSideMatrix<<std::endl;

	  }

	if ( rLocalSystem.CalculationFlags.Is(BeamElement::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
	  {
	    LocalRightHandSideVector.clear();

	    this->CalculateAndAddInertiaRHS( LocalRightHandSideVector, Variables, rCurrentProcessInfo, IntegrationWeight );

	    BeamMathUtilsType::MapLocalToGlobal3D(mInitialLocalQuaternion,LocalRightHandSideVector);

	    VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
	    rRightHandSideVector += LocalRightHandSideVector;

	    //std::cout<<"["<<this->Id()<<"] Beam RotatedDynamic rRightHandSideVector "<<rRightHandSideVector<<std::endl;
	  }

      }

    // VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
    // std::cout<<"["<<this->Id()<<"] Beam RotatedDynamic rRightHandSideVector "<<rRightHandSideVector<<std::endl;
    // MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();
    // std::cout<<"["<<this->Id()<<"] Beam RotatedDynamic rLeftHandSideMatrix "<<rLeftHandSideMatrix<<std::endl;

    mThisIntegrationMethod = ThisIntegrationMethod;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateConstitutiveMatrix(ElementDataType& rVariables)
  {
    KRATOS_TRY

    //Material Elastic constitutive matrix
    this->CalculateMaterialConstitutiveMatrix(rVariables.ConstitutiveMatrix, rVariables);

    //Spatial Elastic constitutive matrix
    this->MapToSpatialFrame( rVariables, rVariables.ConstitutiveMatrix);

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateStressResultants(ElementDataType& rVariables, const unsigned int& rPointNumber)
  {
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    Vector StrainResultants = rVariables.CurrentAxisPositionDerivatives;
    Vector StrainCouples    = rVariables.CurrentCurvatureVector;

    //Reference frame given by the Frame Rotation
    StrainResultants = prod( trans(rVariables.CurrentRotationMatrix), StrainResultants );
    StrainCouples    = prod( trans(rVariables.CurrentRotationMatrix), StrainCouples );

    Vector E1(3);
    noalias(E1) = ZeroVector(3);
    E1[0] = 1.0;
    StrainResultants -= E1;


    for ( SizeType i = 0; i < dimension; i++ )
      {
	rVariables.StrainVector[i]   = StrainResultants[i];
	rVariables.StrainVector[i+3] = StrainCouples[i];
      }


    Matrix ConstitutiveMatrix(6,6);
    noalias(ConstitutiveMatrix) = ZeroMatrix(6,6);
    this->CalculateMaterialConstitutiveMatrix(ConstitutiveMatrix, rVariables);

    //Reference Stress Vector
    rVariables.StressVector = prod( ConstitutiveMatrix, rVariables.StrainVector );

    //std::cout<<" Stress "<<rVariables.StressVector<<" Strain "<<rVariables.StrainVector<<std::endl;

    Vector StressResultants(3);
    noalias(StressResultants) = ZeroVector(dimension);
    Vector StressCouples(3);
    noalias(StressCouples) = ZeroVector(dimension);

    for ( SizeType i = 0; i < dimension; i++ )
      {
	StressResultants[i] = rVariables.StressVector[i];
	StressCouples[i]    = rVariables.StressVector[i+3];
      }

    //Current frame given by the Frame Rotation
    StressResultants = prod( rVariables.CurrentRotationMatrix, StressResultants );
    StressCouples    = prod( rVariables.CurrentRotationMatrix, StressCouples );

    for ( SizeType i = 0; i < dimension; i++ )
      {
	rVariables.StressVector[i]   = StressResultants[i];
	rVariables.StressVector[i+3] = StressCouples[i];
      }

    //set variables for the initialization update
    this->UpdateStrainVariables(rVariables,rPointNumber);

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  //Strain Energy Calculation
  void LargeDisplacementBeamElement::CalculateStrainEnergy(double& rEnergy, ElementDataType& rVariables, const ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight)
  {
    KRATOS_TRY

    //Internal Energy Calculation: alpha = 1
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    //compute Strain Resultants and Couples
    Vector StrainResultants = rVariables.CurrentAxisPositionDerivatives;
    Vector StrainCouples    = rVariables.CurrentCurvatureVector;

    //Reference frame given by the Frame Rotation
    StrainResultants = prod( trans(rVariables.CurrentRotationMatrix), StrainResultants );
    StrainCouples    = prod( trans(rVariables.CurrentRotationMatrix), StrainCouples );

    Vector E1(3);
    noalias(E1) = ZeroVector(3);
    E1[0] = 1.0;
    StrainResultants -= E1;

    //----------------

    Matrix ConstitutiveMatrix(6,6);
    noalias(ConstitutiveMatrix) = ZeroMatrix(6,6);
    this->CalculateMaterialConstitutiveMatrix(ConstitutiveMatrix, rVariables);

    //----------------

    Vector StrainVector(6);
    noalias(StrainVector) = ZeroVector(6);
    for ( SizeType i = 0; i < dimension; i++ )
      {
    	StrainVector[i]   = StrainResultants[i];
    	StrainVector[i+3] = StrainCouples[i];
      }

    //Reference Stress Vector
    Vector StressVector = prod( ConstitutiveMatrix, StrainVector );

    rEnergy += 0.5 * (inner_prod(StressVector, StrainVector)) * rIntegrationWeight ;


    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
								   ElementDataType& rVariables,
								   Vector& rVolumeForce,
								   double& rIntegrationWeight)
  {
    KRATOS_TRY

    BeamElement::CalculateAndAddExternalForces(rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight);

    //follower load forces (to implement)
    this->CalculateAndAddFollowerForces( rRightHandSideVector, rVariables, rIntegrationWeight );

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateExternalForcesEnergy(double& rEnergy,
								   ElementDataType& rVariables,
								   Vector& rVolumeForce,
								   double& rIntegrationWeight)
  {
    KRATOS_TRY

    SizeType number_of_nodes  = GetGeometry().PointsNumber();

    double DomainSize = rVariables.Section.Area;

    //---------------------

    // External Energy Calculation:
    Vector CurrentValueVector(6);
    noalias(CurrentValueVector) = ZeroVector(3);
    Vector CurrentDisplacementVector(3);
    noalias(CurrentDisplacementVector) = ZeroVector(3);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	//Current Displacement Vector
	CurrentValueVector = GetNodalCurrentValue( DISPLACEMENT, CurrentValueVector, i );
	CurrentDisplacementVector += rVariables.N[i] * CurrentValueVector;
      }

    CurrentDisplacementVector = this->MapToInitialLocalFrame( CurrentDisplacementVector, rVariables.PointNumber  );

    //for more than one integration point
    Vector PotencialForceVector(3);
    noalias(PotencialForceVector) = ZeroVector(3);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
    	PotencialForceVector  += rIntegrationWeight * rVariables.N[i] * rVolumeForce * DomainSize;
      }


    rEnergy += inner_prod(PotencialForceVector,CurrentDisplacementVector);

    //---------------------

    // External Energy Calculation:
    noalias(CurrentValueVector) = ZeroVector(3);
    noalias(CurrentDisplacementVector) = ZeroVector(3);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	//Current Displacement Vector
	CurrentValueVector = GetNodalCurrentValue( ROTATION, CurrentValueVector, i );
	CurrentDisplacementVector += rVariables.N[i] * CurrentValueVector;
      }

    CurrentDisplacementVector     = this->MapToInitialLocalFrame( CurrentDisplacementVector, rVariables.PointNumber  );

    //for more than one integration point
    Vector PotencialMomentVector(3);
    noalias(PotencialMomentVector) = ZeroVector(3);

    Matrix SkewSymMatrix(3,3);
    noalias(SkewSymMatrix) = ZeroMatrix(3,3);
    Vector E1(3);
    noalias(E1) = ZeroVector(3);
    E1[0] = 1.0;

    BeamMathUtilsType::VectorToSkewSymmetricTensor(PotencialForceVector,SkewSymMatrix); // m = f x r = skewF · r
    PotencialMomentVector = prod(SkewSymMatrix,E1);

    rEnergy += inner_prod(PotencialMomentVector,CurrentDisplacementVector);


    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateAndAddFollowerForces(VectorType& rRightHandSideVector,
								   ElementDataType& rVariables,
								   double& rIntegrationWeight)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    //Initialize Local Matrices
    VectorType Fi(6);
    noalias(Fi) = ZeroVector(6);
    Vector ResultantsVector(6);
    noalias(ResultantsVector) = ZeroVector(6);

    //locate follower load resultants in skew-symmetric form
    Vector FollowerLoad(dimension);
    noalias(FollowerLoad) = ZeroVector(dimension);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	if( GetGeometry()[i].SolutionStepsDataHas( FOLLOWER_FORCE_LOAD ) )
	  FollowerLoad += rVariables.N[i] * GetGeometry()[i].FastGetSolutionStepValue( FOLLOWER_FORCE_LOAD );
      }

    //Current Frame given by the frame rotation
    FollowerLoad = prod( rVariables.CurrentRotationMatrix, FollowerLoad );

    BeamMathUtilsType::AddVector(FollowerLoad, ResultantsVector, 3);

    //Initialize Local Matrices
    MatrixType OperatorI(6,6);
    noalias(OperatorI) = ZeroMatrix(6,6);

    unsigned int RowIndex = 0;
    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	RowIndex = i * (dimension * 2);

	this->CalculateOperator( OperatorI, rVariables.N, i );

    	//nodal force vector
    	Fi  = prod( OperatorI, ResultantsVector );

	Fi *= rIntegrationWeight;

	BeamMathUtilsType::SubstractVector( Fi, rRightHandSideVector, RowIndex );

      }


    KRATOS_CATCH( "" )

  }

  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
								   ElementDataType& rVariables,
								   double& rIntegrationWeight)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = rRightHandSideVector.size();

    VectorType Fi(6);
    noalias(Fi) = ZeroVector(6);

    //Initialize Local Matrices
    MatrixType DifferentialOperatorI(MatSize,MatSize);
    noalias(DifferentialOperatorI) = ZeroMatrix(MatSize,MatSize);

    unsigned int RowIndex = 0;

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {

	RowIndex = i * (dimension * 2);

	noalias(Fi) = ZeroVector(6);

 	this->CalculateDifferentialOperator(DifferentialOperatorI, rVariables, i, rVariables.Alpha );

	//nodal force vector
	Fi = prod( DifferentialOperatorI, rVariables.StressVector );

	Fi *= rIntegrationWeight;

	BeamMathUtilsType::SubstractVector( Fi, rRightHandSideVector, RowIndex );

      }

    //std::cout<<" Fint "<<rRightHandSideVector<<std::endl;

    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateInternalForcesEnergy(double& rEnergy,
								   ElementDataType& rVariables,
								   double& rIntegrationWeight)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const unsigned int MatSize = this->GetDofsSize();

    VectorType Fi(6);
    noalias(Fi) = ZeroVector(6);

    //Initialize Local Matrices
    MatrixType DifferentialOperatorI(MatSize,MatSize);
    noalias(DifferentialOperatorI) = ZeroMatrix(MatSize,MatSize);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	noalias(Fi) = ZeroVector(6);

 	this->CalculateDifferentialOperator(DifferentialOperatorI, rVariables, i, rVariables.Alpha );

	//nodal force vector
	Fi  = prod( DifferentialOperatorI, rVariables.StressVector );
	Fi *= rIntegrationWeight;

	// Internal Energy Calculation:
	Vector Movements;
	this->GetCurrentNodalMovements( Movements, i, rVariables.PointNumber);

	rEnergy += inner_prod( Fi, Movements );
      }


    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateOperator(MatrixType& rOperator,
						       Vector& rN,
						       const int& rNode)
  {
    KRATOS_TRY

    //Initialize Local Matrices
    if( rOperator.size1() != 6 )
      rOperator.resize(6, 6, false);

    noalias(rOperator) = ZeroMatrix(6,6);

    rOperator( 0, 0 ) =  rN[rNode];
    rOperator( 1, 1 ) =  rN[rNode];
    rOperator( 2, 2 ) =  rN[rNode];
    rOperator( 3, 3 ) =  rN[rNode];
    rOperator( 4, 4 ) =  rN[rNode];
    rOperator( 5, 5 ) =  rN[rNode];

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateDifferentialOperator(MatrixType& rDifferentialOperator,
								   ElementDataType& rVariables,
								   const int& rNode,
								   double alpha)
  {
    KRATOS_TRY

    //Differencial operator transposed

    //Initialize Local Matrices
    if( rDifferentialOperator.size1() != 6 )
      rDifferentialOperator.resize(6, 6, false);

    noalias(rDifferentialOperator) = ZeroMatrix(6,6);

    rDifferentialOperator( 0, 0 ) =  rVariables.DN_DX( rNode, 0 );
    rDifferentialOperator( 1, 1 ) =  rVariables.DN_DX( rNode, 0 );
    rDifferentialOperator( 2, 2 ) =  rVariables.DN_DX( rNode, 0 );
    rDifferentialOperator( 3, 3 ) =  rVariables.DN_DX( rNode, 0 );
    rDifferentialOperator( 4, 4 ) =  rVariables.DN_DX( rNode, 0 );
    rDifferentialOperator( 5, 5 ) =  rVariables.DN_DX( rNode, 0 );

    //locate stress resultants in skew-symmetric "transposed" form
    Matrix SkewSymResultants(3,3);
    noalias(SkewSymResultants) = ZeroMatrix(3,3);

    BeamMathUtilsType::VectorToSkewSymmetricTensor(rVariables.CurrentAxisPositionDerivatives, SkewSymResultants);

    SkewSymResultants *= (-1) * rVariables.N[rNode];

    for ( unsigned int i = 0; i < 3; i++ )
      {
	for ( unsigned int j = 0; j < 3; j++ )
	  {
	    rDifferentialOperator( i+3, j ) = SkewSymResultants(i,j);
   	  }
      }

    KRATOS_CATCH( "" )
  }



  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateDiscreteOperator(MatrixType& rDiscreteOperator,
							       ElementDataType& rVariables,
							       const int& rNode)
  {
    KRATOS_TRY

    //Initialize Local Matrices
    if( rDiscreteOperator.size1() != 6 || rDiscreteOperator.size2() != 9)
      rDiscreteOperator.resize(6, 9, false);

    noalias(rDiscreteOperator) = ZeroMatrix(6,9);

    rDiscreteOperator( 0, 0 ) =  rVariables.DN_DX( rNode, 0 );
    rDiscreteOperator( 1, 1 ) =  rVariables.DN_DX( rNode, 0 );
    rDiscreteOperator( 2, 2 ) =  rVariables.DN_DX( rNode, 0 );
    rDiscreteOperator( 3, 3 ) =  rVariables.DN_DX( rNode, 0 );
    rDiscreteOperator( 4, 4 ) =  rVariables.DN_DX( rNode, 0 );
    rDiscreteOperator( 5, 5 ) =  rVariables.DN_DX( rNode, 0 );

    rDiscreteOperator( 3, 6 ) =  rVariables.N[ rNode ];
    rDiscreteOperator( 4, 7 ) =  rVariables.N[ rNode ];
    rDiscreteOperator( 5, 8 ) =  rVariables.N[ rNode ];

    //std::cout<<" DiscreteOperator "<<rDiscreteOperator<<std::endl;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  //Geometric operator
  void LargeDisplacementBeamElement::CalculateBmatrix(MatrixType& rBmatrix,
						      ElementDataType& rVariables)
  {
    KRATOS_TRY

    //Initialize Local Matrices
    if( rBmatrix.size1() != 9 )
      rBmatrix.resize(9, 9, false);

    noalias(rBmatrix) = ZeroMatrix(9,9);

    Vector StressResultants(3);
    Vector StressCouples(3);
    for ( unsigned int i = 0; i < 3; i++ )
      {
	StressResultants[i] = rVariables.StressVector[i];
        StressCouples[i] = rVariables.StressVector[i+3];
      }

    //NOTE: avoid Kuug noise in plane ploblems
    if( fabs(inner_prod(StressResultants,StressCouples)) < 1e-15 )
      noalias(StressResultants) = ZeroVector(3);

    //locate stress resultants in skew-symmetric form
    Matrix SkewSymStressResultants(3,3);
    noalias(SkewSymStressResultants) = ZeroMatrix(3,3);

    BeamMathUtilsType::VectorToSkewSymmetricTensor(StressResultants, SkewSymStressResultants);

    for ( unsigned int i = 0; i < 3; i++ )
      {
	for ( unsigned int j = 0; j < 3; j++ )
	  {
	    rBmatrix( i+6, j ) =  SkewSymStressResultants(i,j);
	    rBmatrix( i, j+6 ) = -SkewSymStressResultants(i,j);
   	  }
      }

    //locate stress couples in skew-symmetric form
    Matrix SkewSymStressCouples(3,3);
    noalias(SkewSymStressCouples) = ZeroMatrix(3,3);

    BeamMathUtilsType::VectorToSkewSymmetricTensor(StressCouples, SkewSymStressCouples);

    //locate stress couples in skew-symmetric form
    for ( unsigned int i = 0; i < 3; i++ )
      {
	for ( unsigned int j = 0; j < 3; j++ )
	  {
	    rBmatrix( i+3, j+6 ) = -SkewSymStressCouples(i,j);
	    //EXTRA :: Symmetrize the Bmatrix: (convergence performance if linear rotations)
	    //rBmatrix( i+6, j+3 ) =  SkewSymStressCouples(i,j);
	    //EXTRA ::

   	  }
      }

    Matrix AxialStressMatrix(3,3);
    noalias(AxialStressMatrix) = ZeroMatrix(3,3);

    //axial skew-symmetric matrix (Simo-Vu-Quoc)
    // AxialStressMatrix = outer_prod( StressResultants, rVariables.CurrentAxisPositionDerivatives );

    // std::cout<<" StressResultants "<<StressResultants<<std::endl;
    // std::cout<<" AxisPositionDerivatives "<<rVariables.CurrentAxisPositionDerivatives<<std::endl;
    // std::cout<<" AxialStressMatrix "<<AxialStressMatrix<<std::endl;

    // double AxialStressScalar = inner_prod( StressResultants, rVariables.CurrentAxisPositionDerivatives );

    // //EXTRA :: Symmetrize the Bmatrix: (convergence performance if linear rotations)
    // // Matrix ExtraAxialStressMatrix = outer_prod( rVariables.CurrentAxisPositionDerivatives, StressResultants );
    // // AxialStressMatrix += ExtraAxialStressMatrix;
    // //EXTRA ::

    // AxialStressMatrix( 0, 0 ) -= AxialStressScalar;
    // AxialStressMatrix( 1, 1 ) -= AxialStressScalar;
    // AxialStressMatrix( 2, 2 ) -= AxialStressScalar;


    //std::cout<<" AxialStressMatrix "<<AxialStressMatrix<<std::endl;

    //axial skew-symmetric matrix (Crisfield)
    Matrix AxialSkewSymMatrix(3,3);
    noalias(AxialSkewSymMatrix) = ZeroMatrix(3,3);
    BeamMathUtilsType::VectorToSkewSymmetricTensor(rVariables.CurrentAxisPositionDerivatives, AxialSkewSymMatrix);
    AxialStressMatrix = prod(AxialSkewSymMatrix,SkewSymStressResultants);

    //std::cout<<" AxialStressMatrix "<<AxialStressMatrix<<std::endl;

    // EXTRA MOD ::Improve symmetry of the Bmatrix: (convergence performance if linear rotations)
    // for ( unsigned int i = 0; i < rBmatrix.size1(); i++ )
    //   {
    // 	 for ( unsigned int j = 0; j < rBmatrix.size2(); j++ )
    // 	   {
    // 	     if(i!=j)
    // 	       rBmatrix(i,j) = 0 ;
    // 	   }
    //   }
    // EXTRA MOD ::Improve symmetry of the Bmatrix: (convergence performance if linear rotations)

    for ( unsigned int i = 0; i < 3; i++ )
      {
    	for ( unsigned int j = 0; j < 3; j++ )
    	  {
    	    rBmatrix( i+6, j+6 ) = AxialStressMatrix(i,j);
    	  }
      }

    // std::cout<<" Bmatrix : "<<std::endl;
    // for ( unsigned int i = 0; i < rBmatrix.size1(); i++ )
    //   {
    // 	std::cout<<"["<<i<<"] [";
    // 	for ( unsigned int j = 0; j < rBmatrix.size2()-1; j++ )
    // 	  {
    // 	    std::cout<<rBmatrix(i,j)<<", ";
    // 	  }
    // 	std::cout<<rBmatrix(i,rBmatrix.size2()-1)<<" ]"<<std::endl;

    //   }

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
							 ElementDataType& rVariables,
							 double& rIntegrationWeight)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = dimension * 2;

    //Initialize Local Matrices
    MatrixType DifferentialOperatorI(MatSize,MatSize);
    noalias(DifferentialOperatorI) = ZeroMatrix(MatSize,MatSize);
    MatrixType DifferentialOperatorJ(MatSize,MatSize);
    noalias(DifferentialOperatorJ) = ZeroMatrix(MatSize,MatSize);

    MatrixType Kij(MatSize,MatSize);
    noalias(Kij) = ZeroMatrix(MatSize,MatSize);

    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	RowIndex = i * (dimension * 2);

	this->CalculateDifferentialOperator( DifferentialOperatorI, rVariables, i, rVariables.Alpha );

	for ( SizeType j = 0; j < number_of_nodes; j++ )
	  {

	    noalias(Kij) = ZeroMatrix(6,6);

	    ColIndex = j * (dimension * 2);

	    this->CalculateDifferentialOperator( DifferentialOperatorJ, rVariables, j, 1 );

	    noalias(Kij) = prod( DifferentialOperatorI, Matrix(prod( rVariables.ConstitutiveMatrix, trans(DifferentialOperatorJ) )) );

	    Kij *= rIntegrationWeight;

	    //std::cout<<" Kij "<<Kij<<std::endl;

	    //Building the Local Stiffness Matrix
	    BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, Kij, RowIndex, ColIndex );
	  }
      }

    //std::cout<<" Kuum "<<rLeftHandSideMatrix<<std::endl;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
							 ElementDataType& rVariables,
							 double& rIntegrationWeight)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();


    //MatrixType Kuum = rLeftHandSideMatrix;

    //Initialize Local Matrices
    MatrixType Bmatrix(9, 9);
    noalias(Bmatrix) = ZeroMatrix(9, 9);
    this->CalculateBmatrix( Bmatrix, rVariables );

    MatrixType DiscreteOperatorI(6,9);
    noalias(DiscreteOperatorI) = ZeroMatrix(6,9);
    MatrixType DiscreteOperatorJ(6,9);
    noalias(DiscreteOperatorJ) = ZeroMatrix(6,9);

    MatrixType Kij(6,6);
    noalias(Kij) = ZeroMatrix(6,6);

    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {

	RowIndex = i * (dimension * 2);

 	this->CalculateDiscreteOperator( DiscreteOperatorI, rVariables, i );

	for ( SizeType j = 0; j < number_of_nodes; j++ )
	  {

	    noalias(Kij) = ZeroMatrix(6,6);

	    ColIndex = j * (dimension * 2);

	    this->CalculateDiscreteOperator( DiscreteOperatorJ, rVariables, j );

	    noalias(Kij) = prod( DiscreteOperatorI, Matrix( prod( Bmatrix, trans(DiscreteOperatorJ) ) ) );

	    Kij *= rIntegrationWeight;

	    //Building the Local Stiffness Matrix
	    BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, Kij, RowIndex, ColIndex );

	  }
      }

    // bool symmetric = true;
    // MatrixType Kuug = rLeftHandSideMatrix-Kuum;
    // std::cout<<" Kuug : "<<std::endl;
    // for ( unsigned int i = 0; i < Kuug.size1(); i++ )
    //   {
    // 	std::cout<<"["<<i<<"] [";
    // 	for ( unsigned int j = 0; j < Kuug.size2()-1; j++ )
    // 	  {
    // 	    std::cout<<std::scientific<<Kuug(i,j)<<", ";
    // 	    if( Kuug(i,j) != Kuug(j,i) )
    // 	      symmetric = false;

    // 	  }
    // 	std::cout<<Kuug(i,Kuug.size2()-1)<<" ]"<<std::endl;

    //   }

    // if( symmetric == true )
    //   std::cout<<" Kuug SYMMETRIC "<<std::endl;
    // else
    //   std::cout<<" Kuug NON SYMMETRIC "<<std::endl;

    //std::cout<<" Kuug "<<rLeftHandSideMatrix-Kuum<<std::endl;

    // Local geometrical follower load stiffness
    this->CalculateAndAddKuuf( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    KRATOS_CATCH( "" )

  }

  //************************************************************************************
  //************************************************************************************


  void LargeDisplacementBeamElement::CalculateAndAddKuug2(MatrixType& rLeftHandSideMatrix,
							  ElementDataType& rVariables,
							  double& rIntegrationWeight)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();


    // MatrixType Kuum = rLeftHandSideMatrix;

    //Initialize Local Matrices
    MatrixType Kij(6,6);
    noalias(Kij) = ZeroMatrix(6,6);
    Matrix GabK(3,3);
    noalias(GabK) = ZeroMatrix(3,3);
    Matrix DiagonalMatrix(3,3);
    noalias(DiagonalMatrix) = IdentityMatrix(3);

    Vector StressResultants(3);
    noalias(StressResultants) = ZeroVector(3);
    for ( unsigned int i = 0; i < 3; i++ )
      {
	StressResultants[i] = rVariables.StressVector[i];
      }

    Vector StressCouples(3);
    noalias(StressCouples) = ZeroVector(3);
    for ( unsigned int i = 0; i < 3; i++ )
      {
	StressCouples[i] = rVariables.StressVector[i+3];
      }

    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {

	RowIndex = i * (dimension * 2);

	for ( SizeType j = 0; j < number_of_nodes; j++ )
	  {

	    noalias(Kij) = ZeroMatrix(6,6);

	    ColIndex = j * (dimension * 2);


	    //term 11 -> 0
	    //term 12
	    noalias(GabK) = ZeroMatrix(3,3);
	    Matrix SkewSymStressResultants(3,3);
	    noalias(SkewSymStressResultants) = ZeroMatrix(3,3);
	    BeamMathUtilsType::VectorToSkewSymmetricTensor(StressResultants, SkewSymStressResultants);
	    GabK = (-1) * (rVariables.DN_DX(i, 0) * rVariables.N[j]) * SkewSymStressResultants;
	    //Building the Local Stiffness Matrix
	    BeamMathUtilsType::AddMatrix( Kij, GabK, 0, 3 );

	    //term 21
	    noalias(GabK) = ZeroMatrix(3,3);
	    GabK = (rVariables.N[i] * rVariables.DN_DX(j, 0) ) * SkewSymStressResultants;
	    //Building the Local Stiffness Matrix
	    BeamMathUtilsType::AddMatrix( Kij, GabK, 3, 0 );


	    //term 22
	    noalias(GabK) = ZeroMatrix(3,3);

	    Matrix SkewSymStressCouples(3,3);
	    noalias(SkewSymStressCouples) = ZeroMatrix(3,3);
	    BeamMathUtilsType::VectorToSkewSymmetricTensor(StressCouples, SkewSymStressCouples);

	    GabK  = (-1) * (rVariables.DN_DX(i, 0) * rVariables.N[j]) * SkewSymStressCouples;
	    GabK += ( rVariables.N[i] * rVariables.N[j]) * outer_prod( StressResultants, rVariables.CurrentAxisPositionDerivatives );

	    GabK -= ( rVariables.N[i] * rVariables.N[j]) * inner_prod( StressResultants, rVariables.CurrentAxisPositionDerivatives ) * DiagonalMatrix;

	    //Building the Local Stiffness Matrix
	    BeamMathUtilsType::AddMatrix( Kij, GabK, 3, 3 );

	    Kij *= rIntegrationWeight;

	    //Building the Local Stiffness Matrix
	    BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, Kij, RowIndex, ColIndex );

	  }
      }

    // bool symmetric = true;
    // MatrixType Kuug = rLeftHandSideMatrix-Kuum;
    // std::cout<<" Kuug : "<<std::endl;
    // for ( unsigned int i = 0; i < Kuug.size1(); i++ )
    //   {
    // 	std::cout<<"["<<i<<"] [";
    // 	for ( unsigned int j = 0; j < Kuug.size2()-1; j++ )
    // 	  {
    // 	    std::cout<<std::scientific<<Kuug(i,j)<<", ";
    // 	    if( Kuug(i,j) != Kuug(j,i) )
    // 	      symmetric = false;

    // 	  }
    // 	std::cout<<Kuug(i,Kuug.size2()-1)<<" ]"<<std::endl;

    //   }

    // if( symmetric == true )
    //   std::cout<<" Kuug SYMMETRIC "<<std::endl;
    // else
    //   std::cout<<" Kuug NON SYMMETRIC "<<std::endl;

    //std::cout<<" Kuug2 "<<rLeftHandSideMatrix-Kuum<<std::endl;

    // Local geometrical follower load stiffness
    this->CalculateAndAddKuuf( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    KRATOS_CATCH( "" )

  }

  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateAndAddKuuf(MatrixType& rLeftHandSideMatrix,
							 ElementDataType& rVariables,
							 double& rIntegrationWeight)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    //Initialize Local Matrices
    MatrixType Kij(6,6);
    noalias(Kij) = ZeroMatrix(6,6);
    MatrixType Lij(3,3);
    noalias(Lij) = ZeroMatrix(3,3);

    //locate follower load resultants in skew-symmetric form
    Vector FollowerLoad(dimension);
    noalias(FollowerLoad) = ZeroVector(dimension);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	if( GetGeometry()[i].SolutionStepsDataHas( FOLLOWER_FORCE_LOAD ) )
	  FollowerLoad += rVariables.N[i] *  GetGeometry()[i].FastGetSolutionStepValue( FOLLOWER_FORCE_LOAD );
      }

    //Current Frame given by the frame rotation
    FollowerLoad = prod( rVariables.CurrentRotationMatrix, FollowerLoad );

    Matrix SkewSymResultants(dimension,dimension);
    noalias(SkewSymResultants) = ZeroMatrix(dimension,dimension);

    BeamMathUtilsType::VectorToSkewSymmetricTensor(FollowerLoad, SkewSymResultants);

    SkewSymResultants *= (-1);

    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {

	noalias(Kij) = ZeroMatrix(6,6);

	RowIndex = i * (dimension * 2);


	for ( SizeType j = 0; j < number_of_nodes; j++ )
	  {

	    ColIndex = j * (dimension * 2);

	    Lij =  rIntegrationWeight * rVariables.N[i] * rVariables.N[j] * SkewSymResultants;

	    //Building the Local Stiffness Matrix
	    BeamMathUtilsType::AddMatrix( Kij, Lij, 0, 3 );

	    //Building the Local Stiffness Matrix
	    BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, Kij, RowIndex, ColIndex );

	  }
      }


    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************

  //Inertia in the SPATIAL configuration
  void LargeDisplacementBeamElement::CalculateAndAddInertiaLHS(MatrixType& rLeftHandSideMatrix, ElementDataType& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight )
  {

    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize               = number_of_nodes * ( dimension * 2 );

    if(rLeftHandSideMatrix.size1() != MatSize)
      rLeftHandSideMatrix.resize (MatSize, MatSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix( MatSize, MatSize );

    SectionProperties Section;
    this->CalculateSectionProperties(Section);


    //rCurrentProcessInfo must give it:
    double DeltaTime = rCurrentProcessInfo[DELTA_TIME];

    double Beta = 1;
    double Newmark1 = 1;
    if( rCurrentProcessInfo.Has(NEWMARK_BETA) ){
      Beta = rCurrentProcessInfo[NEWMARK_BETA];
      Newmark1 = (1.0/ ( DeltaTime * DeltaTime * Beta ));
    }

    double Gamma = 1;
    double Newmark2 = 1;
    if( rCurrentProcessInfo.Has(NEWMARK_GAMMA) ){
      Gamma = rCurrentProcessInfo[NEWMARK_GAMMA];
      Newmark2 = ( DeltaTime * Gamma );
    }


    //block m(1,1) of the mass matrix

    MatrixType m11(3,3);
    noalias(m11) = ZeroMatrix(3,3);

    double TotalMass = 0;
    TotalMass  = this->CalculateTotalMass( Section, TotalMass );

    //block m(2,2) of the mass matrix

    MatrixType m22(3,3);
    noalias(m22) = ZeroMatrix(3,3);

    Vector CurrentCompoundRotationVector(3);
    noalias(CurrentCompoundRotationVector) = ZeroVector(3);
    Vector ReferenceCompoundRotationVector(3);
    noalias(ReferenceCompoundRotationVector) = ZeroVector(3);

    Vector CurrentStepRotationVector(3);
    noalias(CurrentStepRotationVector) = ZeroVector(3);
    Vector PreviousStepRotationVector(3);
    noalias(PreviousStepRotationVector) = ZeroVector(3);

    Vector AngularVelocityVector(3);
    noalias(AngularVelocityVector) = ZeroVector(3);
    Vector AngularAccelerationVector(3);
    noalias(AngularAccelerationVector) = ZeroVector(3);
    Vector CurrentAngularAccelerationVector(3);
    noalias(CurrentAngularAccelerationVector) = ZeroVector(3);
    Vector PreviousAngularAccelerationVector(3);
    noalias(PreviousAngularAccelerationVector) = ZeroVector(3);

    Vector CurrentValueVector(3);
    noalias(CurrentValueVector) = ZeroVector(3);
    Vector PreviousValueVector(3);
    noalias(PreviousValueVector) = ZeroVector(3);

    std::vector<QuaternionType> ReferenceNodeQuaternions;
    QuaternionType QuaternionValue;

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {

	//Current Compound Rotation Vector
	CurrentValueVector = GetNodalCurrentValue( ROTATION, CurrentValueVector, i );
	CurrentValueVector = this->MapToInitialLocalFrame( CurrentValueVector );

	CurrentCompoundRotationVector   += rVariables.N[i] * CurrentValueVector;

	//Reference Compound Rotation Vector
	PreviousValueVector = GetNodalPreviousValue( ROTATION, PreviousValueVector, i );
	PreviousValueVector = this->MapToInitialLocalFrame( PreviousValueVector );

	ReferenceCompoundRotationVector += rVariables.N[i] * PreviousValueVector;

	//Reference Frame is the Local Frame
	QuaternionValue = QuaternionType::FromRotationVector(PreviousValueVector);

	ReferenceNodeQuaternions.push_back(QuaternionValue);


	//Current Step Rotation Vector
	CurrentValueVector = GetNodalCurrentValue( STEP_ROTATION, CurrentValueVector, i );
	noalias(CurrentStepRotationVector) += rVariables.N[i] * CurrentValueVector;

	//Previous Step Rotation Vector
	PreviousValueVector = GetNodalPreviousValue( STEP_ROTATION, PreviousValueVector, i );
	noalias(PreviousStepRotationVector) += rVariables.N[i] * PreviousValueVector;

	//Angular Velocity Vector
	CurrentValueVector = GetNodalCurrentValue( ANGULAR_VELOCITY, CurrentValueVector, i );
	noalias(AngularVelocityVector) += rVariables.N[i] * CurrentValueVector;

	//Current Angular Acceleration Vector
	CurrentValueVector = GetNodalCurrentValue( ANGULAR_ACCELERATION, CurrentValueVector, i );
	noalias(CurrentAngularAccelerationVector) += rVariables.N[i] * CurrentValueVector;

        //Previous Angular Acceleration Vector
	CurrentValueVector = GetNodalPreviousValue( ANGULAR_ACCELERATION, CurrentValueVector, i );
	noalias(PreviousAngularAccelerationVector) += rVariables.N[i] * CurrentValueVector;

      }


    //Set step variables to local frame (current Frame is the local frame)
    PreviousStepRotationVector        = this->MapToInitialLocalFrame( PreviousStepRotationVector );
    CurrentStepRotationVector         = this->MapToInitialLocalFrame( CurrentStepRotationVector );
    AngularVelocityVector             = this->MapToInitialLocalFrame( AngularVelocityVector );
    CurrentAngularAccelerationVector  = this->MapToInitialLocalFrame( CurrentAngularAccelerationVector );
    PreviousAngularAccelerationVector = this->MapToInitialLocalFrame( PreviousAngularAccelerationVector );

    double AlphaM = 0;
    if( rCurrentProcessInfo.Has(BOSSAK_ALPHA) ){
      AlphaM = rCurrentProcessInfo[BOSSAK_ALPHA];
    }

    AngularAccelerationVector = (1.0-AlphaM)*CurrentAngularAccelerationVector + AlphaM*(PreviousAngularAccelerationVector);

    //1.-Get inertia dyadic
    Matrix InertiaDyadic(3,3);
    noalias(InertiaDyadic) = ZeroMatrix(3,3);
    this->CalculateInertiaDyadic( Section, InertiaDyadic );
    InertiaDyadic = prod(rVariables.CurrentRotationMatrix,InertiaDyadic);
    InertiaDyadic = prod(InertiaDyadic,trans(rVariables.CurrentRotationMatrix));

    //2.- Compute Term 1:

    Matrix MassTerm1(3,3);
    noalias(MassTerm1) = ZeroMatrix(3,3);

    Vector InertiaxAngularVelocity     = prod( InertiaDyadic, AngularVelocityVector );
    Vector InertiaxAngularAcceleration = prod( InertiaDyadic, AngularAccelerationVector );

    Matrix TensorAngularVelocity(3,3);
    noalias(TensorAngularVelocity) = ZeroMatrix(3,3);

    BeamMathUtilsType::VectorToSkewSymmetricTensor( AngularVelocityVector, TensorAngularVelocity );

    Vector VectorTerm1(3);
    noalias(VectorTerm1) = ZeroVector(3);

    VectorTerm1  = prod( TensorAngularVelocity, InertiaxAngularVelocity );

    VectorTerm1 += InertiaxAngularAcceleration;

    BeamMathUtilsType::VectorToSkewSymmetricTensor(VectorTerm1, MassTerm1);


    //3.- Compute Term 2:

    Matrix MassTerm2(3,3);
    noalias(MassTerm2) = ZeroMatrix(3,3);

    Matrix InertiaxAngularVelocityTensor(3,3);
    noalias(InertiaxAngularVelocityTensor) = ZeroMatrix(3,3);
    BeamMathUtilsType::VectorToSkewSymmetricTensor( InertiaxAngularVelocity, InertiaxAngularVelocityTensor );

    Matrix TensorAngularVelocityxInertia(3,3);
    noalias(TensorAngularVelocityxInertia) = ZeroMatrix(3,3);
    noalias(TensorAngularVelocityxInertia) = prod( TensorAngularVelocity, InertiaDyadic );

    Matrix InertiaxTensorAngularVelocity(3,3);
    noalias(InertiaxTensorAngularVelocity) = ZeroMatrix(3,3);
    noalias(InertiaxTensorAngularVelocity) = prod( InertiaDyadic, TensorAngularVelocity );

    //MassTerm2 = (1.0-AlphaM) * InertiaDyadic - Newmark2 * InertiaxAngularVelocityTensor + Newmark2 * TensorAngularVelocityxInertia + Newmark2 * InertiaxTensorAngularVelocity; //Cardona

    MassTerm2 = (1.0-AlphaM) * InertiaDyadic - Newmark2 * InertiaxAngularVelocityTensor + Newmark2 * TensorAngularVelocityxInertia; //Simo

    MassTerm2 *= Newmark1;

    MatrixType MassMatrixBlock2(3,3);
    noalias(MassMatrixBlock2) = ZeroMatrix(3,3);

    MassMatrixBlock2 = (-1) * MassTerm1 + MassTerm2;

    // Compute Linear Part of the Step Rotation
    Matrix LinearPartRotationTensor(3,3);
    noalias(LinearPartRotationTensor) = ZeroMatrix(3,3);

    Vector StepRotationVector = PreviousStepRotationVector; //CurrentStepRotationVector;

    double NormStepRotation =  norm_2(StepRotationVector);
    if( NormStepRotation != 0 ){
      this->CalculateRotationLinearPartTensor( StepRotationVector , LinearPartRotationTensor );
    }
    else{
      noalias(LinearPartRotationTensor) = IdentityMatrix(3);
    }


    MassMatrixBlock2 = prod( MassMatrixBlock2, LinearPartRotationTensor );

    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    Matrix DiagonalMatrix(3,3);
    noalias(DiagonalMatrix) = IdentityMatrix(3);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	noalias(m11) = ZeroMatrix(3,3);
	noalias(m22) = ZeroMatrix(3,3);

	RowIndex = i * (dimension * 2);


	for ( SizeType j = 0; j < number_of_nodes; j++ )
	  {

	    ColIndex = j * (dimension * 2);

	    m11 = (1.0-AlphaM) * Newmark1 * TotalMass * rVariables.N[i] * rVariables.N[j] * rIntegrationWeight * DiagonalMatrix;

	    m22 = rVariables.N[i] * rVariables.N[j] * rIntegrationWeight * MassMatrixBlock2;


	    //Building the Local Tangent Inertia Matrix
	    BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, m11, RowIndex, ColIndex );
	    BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, m22, RowIndex+3, ColIndex+3 );

	  }
      }

    KRATOS_CATCH( "" )

  }

  //************************************************************************************
  //************************************************************************************

  //Inertia in the SPATIAL configuration
  void LargeDisplacementBeamElement::CalculateAndAddInertiaRHS(VectorType& rRightHandSideVector, ElementDataType& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize               = number_of_nodes * ( dimension * 2 );

    if(rRightHandSideVector.size() != MatSize)
      rRightHandSideVector.resize(MatSize, false);

    noalias(rRightHandSideVector) = ZeroVector( MatSize );

    SectionProperties Section;
    this->CalculateSectionProperties(Section);


    double TotalMass = 0;
    TotalMass = this->CalculateTotalMass( Section, TotalMass );

    //displacements and rotations vector

    Vector LinearAccelerationVector(3);
    noalias(LinearAccelerationVector) = ZeroVector(3);
    Vector CurrentLinearAccelerationVector(3);
    noalias(CurrentLinearAccelerationVector) = ZeroVector(3);
    Vector PreviousLinearAccelerationVector(3);
    noalias(PreviousLinearAccelerationVector) = ZeroVector(3);

    Vector AngularVelocityVector(3);
    noalias(AngularVelocityVector) = ZeroVector(3);
    Vector AngularAccelerationVector(3);
    noalias(AngularAccelerationVector) = ZeroVector(3);
    Vector CurrentAngularAccelerationVector(3);
    noalias(CurrentAngularAccelerationVector) = ZeroVector(3);
    Vector PreviousAngularAccelerationVector(3);
    noalias(PreviousAngularAccelerationVector) = ZeroVector(3);


    Vector CurrentValueVector(3);
    noalias(CurrentValueVector) = ZeroVector(3);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	//Current Linear Acceleration Vector
	CurrentValueVector = GetNodalCurrentValue( ACCELERATION, CurrentValueVector, i );
	CurrentLinearAccelerationVector        += rVariables.N[i] * CurrentValueVector;

	//Previous Linear Acceleration Vector
	CurrentValueVector = GetNodalPreviousValue( ACCELERATION, CurrentValueVector, i );
	PreviousLinearAccelerationVector       += rVariables.N[i] * CurrentValueVector;

	//Angular Velocity Vector
	CurrentValueVector = GetNodalCurrentValue( ANGULAR_VELOCITY, CurrentValueVector, i );
	AngularVelocityVector                  += rVariables.N[i] * CurrentValueVector;

	//CurrentAngular Acceleration Vector
	CurrentValueVector = GetNodalCurrentValue( ANGULAR_ACCELERATION, CurrentValueVector, i );
	CurrentAngularAccelerationVector       += rVariables.N[i] * CurrentValueVector;

        //Previous Angular Acceleration Vector
	CurrentValueVector = GetNodalPreviousValue( ANGULAR_ACCELERATION, CurrentValueVector, i );
	PreviousAngularAccelerationVector      += rVariables.N[i] * CurrentValueVector;

      }

    //Set step variables to local frame (current Frame is the local frame)
    CurrentLinearAccelerationVector   = this->MapToInitialLocalFrame( CurrentLinearAccelerationVector );
    PreviousLinearAccelerationVector  = this->MapToInitialLocalFrame( PreviousLinearAccelerationVector );

    AngularVelocityVector             = this->MapToInitialLocalFrame( AngularVelocityVector );

    CurrentAngularAccelerationVector  = this->MapToInitialLocalFrame( CurrentAngularAccelerationVector );
    PreviousAngularAccelerationVector = this->MapToInitialLocalFrame( PreviousAngularAccelerationVector );

    double AlphaM = 0;
    if( rCurrentProcessInfo.Has(BOSSAK_ALPHA) ){
      AlphaM = rCurrentProcessInfo[BOSSAK_ALPHA];
    }

    LinearAccelerationVector  = (1.0-AlphaM) * CurrentLinearAccelerationVector  + AlphaM * (PreviousLinearAccelerationVector);
    AngularAccelerationVector = (1.0-AlphaM) * CurrentAngularAccelerationVector + AlphaM * (PreviousAngularAccelerationVector);

    //-----------------
    //block m(1) of the inertial force vector

    //Compute Linear Term
    Vector LinearInertialForceVector(3);
    noalias(LinearInertialForceVector) = ZeroVector(3);
    LinearInertialForceVector = TotalMass * LinearAccelerationVector;


    //-----------------
    //block m(2,2) of the inertial force vector (rotations part::to be defined)

    //Get inertia dyadic
    Matrix InertiaDyadic(3,3);
    noalias(InertiaDyadic) = ZeroMatrix(3,3);
    this->CalculateInertiaDyadic( Section, InertiaDyadic );
    InertiaDyadic = prod(rVariables.CurrentRotationMatrix,InertiaDyadic);
    InertiaDyadic = prod(InertiaDyadic,trans(rVariables.CurrentRotationMatrix));


    //Compute Angular Term:
    Vector InertiaxAngularVelocity     = prod( InertiaDyadic, AngularVelocityVector );
    Vector InertiaxAngularAcceleration = prod( InertiaDyadic, AngularAccelerationVector );

    Matrix TensorAngularVelocity(3,3);
    noalias(TensorAngularVelocity) = ZeroMatrix(3,3);
    BeamMathUtilsType::VectorToSkewSymmetricTensor( AngularVelocityVector, TensorAngularVelocity );

    Vector AngularInertialForceVector(3);
    noalias(AngularInertialForceVector) = ZeroVector(3);
    AngularInertialForceVector  = prod( TensorAngularVelocity, InertiaxAngularVelocity );

    // CROSS PRODUCT of AxB = prod( skewA, B )  where (skewA) = [Ax] = [A]^v =>  (skewA)^T = hat(A) (nomenclature)
    AngularInertialForceVector += InertiaxAngularAcceleration;


    //compose total acceleration integral function:
    Vector TotalInertialForceVector(6);
    noalias(TotalInertialForceVector) = ZeroVector(6);

    BeamMathUtilsType::AddVector(LinearInertialForceVector, TotalInertialForceVector, 0);
    BeamMathUtilsType::AddVector(AngularInertialForceVector, TotalInertialForceVector, 3);

    //Initialize Local Matrices
    VectorType Fi(6);
    noalias(Fi) = ZeroVector(6);
    MatrixType OperatorI(6,6);
    noalias(OperatorI) = ZeroMatrix(6,6);
    unsigned int RowIndex = 0;

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	RowIndex = i * (dimension * 2);

	noalias(Fi) = ZeroVector(6);

 	this->CalculateOperator( OperatorI, rVariables.N, i );

    	//nodal force vector
    	Fi  = prod( OperatorI, TotalInertialForceVector );

	Fi *= rIntegrationWeight;

	BeamMathUtilsType::AddVector(Fi, rRightHandSideVector, RowIndex);
      }


    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  //Kinetic Energy Calculation
  void LargeDisplacementBeamElement::CalculateKineticEnergy(double& rEnergy,
							    ElementDataType& rVariables,
							    const ProcessInfo& rCurrentProcessInfo,
							    double& rIntegrationWeight)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();

    SectionProperties Section;
    this->CalculateSectionProperties(Section);

    //Get total mass
    double TotalMass = 0;
    TotalMass = this->CalculateTotalMass( Section, TotalMass );

    //Get inertia dyadic
    Matrix InertiaDyadic(3,3);
    noalias(InertiaDyadic) = ZeroMatrix(3,3);
    this->CalculateInertiaDyadic( Section, InertiaDyadic );

    Matrix CurrentInertiaDyadic = prod(rVariables.CurrentRotationMatrix,InertiaDyadic);
    CurrentInertiaDyadic = prod(CurrentInertiaDyadic,trans(rVariables.CurrentRotationMatrix));

    // Kinetic Energy Calculation:
    Vector CurrentNodalVelocities(6);
    noalias(CurrentNodalVelocities) = ZeroVector(6);
    Vector CurrentVelocitiesVector(6);
    noalias(CurrentVelocitiesVector) = ZeroVector(6);

    Matrix KineticMatrix;
    this->GetKineticMatrix( KineticMatrix, TotalMass, CurrentInertiaDyadic );


    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {

	noalias(CurrentNodalVelocities) = ZeroVector(6);
	this->GetCurrentNodalVelocities( CurrentNodalVelocities, i, rVariables.PointNumber );
	CurrentVelocitiesVector += rVariables.N[i] * CurrentNodalVelocities;

      }


    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	rEnergy += 0.5 * rVariables.N[i] * inner_prod( CurrentVelocitiesVector, prod( KineticMatrix, CurrentVelocitiesVector) ) * rIntegrationWeight;
      }


    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  //Linear Momentum and Angular Momentum Calculation
  void LargeDisplacementBeamElement::CalculateMomentumRelations(array_1d<double,3>& rLinearMomentum,
								array_1d<double,3>& rAngularMomentum,
								ElementDataType& rVariables,
								const ProcessInfo& rCurrentProcessInfo,
								double& rIntegrationWeight)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();

    SectionProperties Section;
    this->CalculateSectionProperties(Section);

    //Get total mass
    double TotalMass = 0;
    TotalMass = this->CalculateTotalMass( Section, TotalMass );

    //Get inertia dyadic
    Matrix InertiaDyadic(3,3);
    noalias(InertiaDyadic) = ZeroMatrix(3,3);
    this->CalculateInertiaDyadic( Section, InertiaDyadic );

    //spatial
    Matrix CurrentInertiaDyadic = prod(rVariables.CurrentRotationMatrix,InertiaDyadic);
    CurrentInertiaDyadic = prod(CurrentInertiaDyadic,trans(rVariables.CurrentRotationMatrix));

    Matrix DiagonalMatrix = identity_matrix<double> (3);

    Vector CurrentValueVector(3);
    noalias(CurrentValueVector) = ZeroVector(3);
    Vector CurrentPositionVector(3);
    noalias(CurrentPositionVector) = ZeroVector(3);
    Vector CurrentDisplacementVector(3);
    noalias(CurrentDisplacementVector) = ZeroVector(3);
    Vector CurrentLinearVelocityVector(3);
    noalias(CurrentLinearVelocityVector) = ZeroVector(3);
    Vector CurrentAngularVelocityVector(3);
    noalias(CurrentAngularVelocityVector) = ZeroVector(3);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	//Current Position Vector
	CurrentValueVector     = GetGeometry()[i].Coordinates();
	CurrentPositionVector += rVariables.N[i] * CurrentValueVector;

	//Current Displacement Vector
	CurrentValueVector = GetNodalCurrentValue( DISPLACEMENT, CurrentValueVector, i );
	CurrentDisplacementVector += rVariables.N[i] * CurrentValueVector;

	//Current Linear Velocity Vector
	CurrentValueVector = GetNodalCurrentValue( VELOCITY, CurrentValueVector, i );
	CurrentLinearVelocityVector += rVariables.N[i] * CurrentValueVector;

	//Current Angular Velocity Vector
	CurrentValueVector = GetNodalCurrentValue( ANGULAR_VELOCITY, CurrentValueVector, i );
	CurrentAngularVelocityVector += rVariables.N[i] * CurrentValueVector;
      }

    CurrentPositionVector         = this->MapToInitialLocalFrame( CurrentPositionVector, rVariables.PointNumber );
    CurrentDisplacementVector     = this->MapToInitialLocalFrame( CurrentDisplacementVector, rVariables.PointNumber  );
    CurrentLinearVelocityVector   = this->MapToInitialLocalFrame( CurrentLinearVelocityVector, rVariables.PointNumber  );
    CurrentAngularVelocityVector  = this->MapToInitialLocalFrame( CurrentAngularVelocityVector, rVariables.PointNumber );



    // Calculate Linear Momentum and Angular Momentum
    Vector LinearMomentumVector(3);
    noalias(LinearMomentumVector) = ZeroVector(3);
    Vector AngularMomentumVector(3);
    noalias(AngularMomentumVector) = ZeroVector(3);


    //for more than one integration point
    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
    	LinearMomentumVector  += TotalMass * rVariables.N[i] * prod( DiagonalMatrix, CurrentLinearVelocityVector ) * rIntegrationWeight;
      }

    MathUtils<double>::CrossProduct(AngularMomentumVector, CurrentPositionVector, LinearMomentumVector);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
    	AngularMomentumVector += prod( CurrentInertiaDyadic, CurrentAngularVelocityVector ) * rVariables.N[i] * rIntegrationWeight;
      }

    //for only one integration point
    //LinearMomentumVector   = TotalMass * prod( DiagonalMatrix, CurrentLinearVelocityVector ) * rIntegrationWeight;
    //MathUtils<double>::CrossProduct(AngularMomentumVector, CurrentPositionVector, LinearMomentumVector);
    //AngularMomentumVector += prod( CurrentInertiaDyadic, CurrentAngularVelocityVector ) * rIntegrationWeight;

    // Note:
    // LocalTransformationMatrix is a rotation matrix with new base in columns
    // That means that the standard rotation K = Q·K'·QT and F = Q·F' is the correct transformation

    Matrix LocalTransformationMatrix(3,3);
    noalias(LocalTransformationMatrix) = ZeroMatrix(3,3);
    mInitialLocalQuaternion.ToRotationMatrix(LocalTransformationMatrix);

    LinearMomentumVector  = prod(LinearMomentumVector,  trans(LocalTransformationMatrix) );
    AngularMomentumVector = prod(AngularMomentumVector, trans(LocalTransformationMatrix) );

    for( unsigned int i=0; i<3; i++ )
      {
	rLinearMomentum[i]  += LinearMomentumVector[i];
	rAngularMomentum[i] += AngularMomentumVector[i];
      }


    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::GetKineticMatrix(Matrix& rKineticMatrix, const double& rMass, const Matrix& rInertia)
  {
    KRATOS_TRY

    if( rKineticMatrix.size1() != 6 || rKineticMatrix.size2() != 6 )
      rKineticMatrix.resize(6,6);

    noalias(rKineticMatrix) = ZeroMatrix(6,6);

    //Building the rotation matrix for the local element matrix
    for (unsigned int i=0; i<3; i++)
      {
	rKineticMatrix(i,i) = rMass;
      }

    for (unsigned int i=0; i<3; i++)
      {
	for(unsigned int j=0; j<3; j++)
	  {
	    rKineticMatrix(i+3,j+3) = rInertia(i,j);
	  }
      }

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************
  void LargeDisplacementBeamElement::GetCurrentNodalVelocities(Vector& rValues, const int& rNode, unsigned int PointNumber)
  {
    if( rValues.size() != 6 )
      rValues.resize(6);

    noalias(rValues) = ZeroVector(6);

    Vector CurrentValueVector(3);
    noalias(CurrentValueVector) = ZeroVector(3);
    CurrentValueVector = GetNodalCurrentValue( VELOCITY, CurrentValueVector, rNode );
    CurrentValueVector = this->MapToInitialLocalFrame( CurrentValueVector, PointNumber );

    rValues[0] = CurrentValueVector[0];
    rValues[1] = CurrentValueVector[1];
    rValues[2] = CurrentValueVector[2];


    noalias(CurrentValueVector) = ZeroVector(3);
    CurrentValueVector = GetNodalCurrentValue( ANGULAR_VELOCITY, CurrentValueVector, rNode );
    CurrentValueVector = this->MapToInitialLocalFrame( CurrentValueVector, PointNumber );

    rValues[3] = CurrentValueVector[0];
    rValues[4] = CurrentValueVector[1];
    rValues[5] = CurrentValueVector[2];

  }

  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::GetCurrentNodalMovements(Vector& rValues, const int& rNode, unsigned int PointNumber)
  {
    if( rValues.size() != 6 )
      rValues.resize(6);

    noalias(rValues) = ZeroVector(6);

    Vector CurrentValueVector(3);
    noalias(CurrentValueVector) = ZeroVector(3);
    CurrentValueVector = GetNodalCurrentValue( DISPLACEMENT, CurrentValueVector, rNode );
    CurrentValueVector = this->MapToInitialLocalFrame( CurrentValueVector, PointNumber );

    rValues[0] = CurrentValueVector[0];
    rValues[1] = CurrentValueVector[1];
    rValues[2] = CurrentValueVector[2];

    noalias(CurrentValueVector) = ZeroVector(3);
    CurrentValueVector = GetNodalCurrentValue( ROTATION, CurrentValueVector, rNode );
    CurrentValueVector = this->MapToInitialLocalFrame( CurrentValueVector, PointNumber );

    rValues[3] = CurrentValueVector[0];
    rValues[4] = CurrentValueVector[1];
    rValues[5] = CurrentValueVector[2];

  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
  {

    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize               = number_of_nodes * ( dimension * 2 );

    if(rMassMatrix.size1() != MatSize)
      rMassMatrix.resize (MatSize, MatSize, false);

    noalias(rMassMatrix) = ZeroMatrix(MatSize, MatSize);

    SectionProperties Section;
    this->CalculateSectionProperties(Section);

    //block m(1,1) of the mass matrix

    MatrixType m11(3,3);
    noalias(m11) = ZeroMatrix(3,3);

    double TotalMass = 0;
    TotalMass = this->CalculateTotalMass( Section, TotalMass );
    TotalMass *= GetGeometry().Length();

    Vector LumpFact(number_of_nodes);
    noalias(LumpFact) = ZeroVector(number_of_nodes);

    LumpFact = GetGeometry().LumpingFactors(LumpFact);


    //block m(2,2) of the mass matrix

    MatrixType m22(3,3);
    noalias(m22) = ZeroMatrix(3,3);

    Matrix InertiaDyadic(3,3);
    noalias(InertiaDyadic) = ZeroMatrix(3,3);
    this->CalculateInertiaDyadic( Section, InertiaDyadic );


    for( SizeType i=0; i < number_of_nodes; i++ )
      {

        double temp = LumpFact[i] * TotalMass;

	int RowIndex = i * (dimension * 2);

	for( SizeType k=0; k < dimension; k++ )
	  {
	    m11(k,k) = temp;
	  }

	m22 = InertiaDyadic * temp;


	//Building the Local Tangent Inertia Matrix
	BeamMathUtilsType::AddMatrix( rMassMatrix, m11, RowIndex, RowIndex );
	BeamMathUtilsType::AddMatrix( rMassMatrix, m22, RowIndex+3, RowIndex+3 );

      }


    // Transform Local to Global LHSMatrix:
    BeamMathUtilsType::MapLocalToGlobal3D(mInitialLocalQuaternion,rMassMatrix);


    KRATOS_CATCH( "" )

  }


  //************************************CALCULATE TOTAL MASS****************************
  //************************************************************************************

  double& LargeDisplacementBeamElement::CalculateTotalMass( SectionProperties& Section, double& rTotalMass )
  {
    KRATOS_TRY

    //const SizeType dimension  = GetGeometry().WorkingSpaceDimension();
    rTotalMass = ( Section.Area ) * GetProperties()[DENSITY];

    return rTotalMass;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateInertiaDyadic(SectionProperties& rSection, Matrix& rInertiaDyadic)
  {
    KRATOS_TRY

    if( rInertiaDyadic.size1() != 3 )
      rInertiaDyadic.resize(3, 3, false);

    noalias(rInertiaDyadic) = ZeroMatrix(3,3);

    //if the local axes are the principal axes of the cross section

    //Axis Local E1
    rInertiaDyadic(0,0) = rSection.Rotational_Inertia; //beam axis
    rInertiaDyadic(1,1) = rSection.Inertia_y; //vertial axis
    rInertiaDyadic(2,2) = rSection.Inertia_z; //horizontal axis

    rInertiaDyadic *= GetProperties()[DENSITY];

    //std::cout<<" INERTIA DYADIC "<<rInertiaDyadic<<std::endl;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateRotationLinearPartTensor(Vector& rRotationVector, Matrix& rRotationTensor)

  {
    KRATOS_TRY

    if( rRotationTensor.size1() != 3 )
      rRotationTensor.resize(3, 3, false);

    noalias(rRotationTensor) = ZeroMatrix(3,3);

    // adding term 3
    BeamMathUtilsType::VectorToSkewSymmetricTensor( rRotationVector, rRotationTensor );

    rRotationTensor *= 0.5;

    double NormRotation =  norm_2(rRotationVector);
    Matrix RotationxRotation = outer_prod( rRotationVector, rRotationVector );

    if( NormRotation != 0 )
      RotationxRotation *= ( 1.0/( NormRotation * NormRotation ) );

    // adding term 1
    rRotationTensor += RotationxRotation;

    double Coefficient = 0;
    if( NormRotation != 0 )
      Coefficient = 0.5 * NormRotation / std::tan( 0.5 * NormRotation );

    for(unsigned int i=0; i<3; i++)
      RotationxRotation(i,i) -= 1.0;

    RotationxRotation *= (-1)*Coefficient;

    // adding term 2
    rRotationTensor += RotationxRotation;


    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo )
  {

    KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_points_number )
      rOutput.resize( integration_points_number, false );

    if ( rVariable == KINETIC_ENERGY ){

      //create and initialize element variables:
      ElementDataType Variables;

      IntegrationMethod ThisIntegrationMethod = mThisIntegrationMethod;
      //full quadrature integration:
      mThisIntegrationMethod = mFullIntegrationMethod;

      this->InitializeElementData(Variables,rCurrentProcessInfo);

      const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

      // initialize variables short version:

      //Compute Section Properties:
      this->CalculateSectionProperties(Variables.Section);
      Variables.Length = GetGeometry().Length();

      //Set equilibrium point initial:0/mid:0.5/final:1
      if( rCurrentProcessInfo.Has(EQUILIBRIUM_POINT) )
	Variables.Alpha = rCurrentProcessInfo[EQUILIBRIUM_POINT];
      else
	Variables.Alpha = 1;

      //reading shape functions
      Variables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

      //get the shape functions for the order of the integration method [N]
      const Matrix& Ncontainer = Variables.GetShapeFunctions();

      double IntegrationWeight = 0;
      double Energy = 0;

      for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
	{
	  //integration point number
	  Variables.PointNumber = PointNumber;

	  //set shape functions values for this integration point
	  noalias(Variables.N) = matrix_row<const Matrix>( Ncontainer, PointNumber);

	  //compute local to global frame
	  this->CalculateFrameMapping( Variables, PointNumber );

	  Variables.detJ = 1.0/mInvJ0;

	  IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
	  IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );


	  Energy = 0;

	  //Kinetic energy calculation
	  this->CalculateKineticEnergy(Energy, Variables, rCurrentProcessInfo, IntegrationWeight);

	  rOutput[PointNumber] = Energy;

	  //Momentum calculation
	  //this->CalculateMomentumRelations(Variables, rCurrentProcessInfo, IntegrationWeight);

	}

      mThisIntegrationMethod = ThisIntegrationMethod;

    }
    else if ( rVariable == INTERNAL_ENERGY ){
      //create and initialize element variables:
      ElementDataType Variables;
      this->InitializeElementData(Variables,rCurrentProcessInfo);

      //reading integration points
      const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

      double IntegrationWeight = 1.0;
      double Energy = 0;
      for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
	{
	  //compute element kinematics  ...
	  this->CalculateKinematics(Variables,PointNumber);

	  //compute element ConstitutiveTensor
	  this->CalculateConstitutiveMatrix(Variables);

	  //compute element Strain and Stress Resultants and Couples
	  this->CalculateStressResultants(Variables, PointNumber);

	  IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
	  IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );

	  Energy = 0;
	  this->CalculateInternalForcesEnergy(Energy, Variables, IntegrationWeight);

	  rOutput[PointNumber] = Energy;
	}
    }
    else if ( rVariable == EXTERNAL_ENERGY ){

      //create and initialize element variables:
      ElementDataType Variables;
      this->InitializeElementData(Variables,rCurrentProcessInfo);

      //reading integration points
      const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

      //auxiliary terms
      Vector VolumeForce(3);
      noalias(VolumeForce) = ZeroVector(3);

      double IntegrationWeight = 1.0;
      double Energy = 0;

      for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
	{
	  //compute element kinematics  ...
	  this->CalculateKinematics(Variables,PointNumber);

	  //compute element ConstitutiveTensor
	  this->CalculateConstitutiveMatrix(Variables);

	  //compute element Strain and Stress Resultants and Couples
	  this->CalculateStressResultants(Variables, PointNumber);

	  IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
	  IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );

	  //contribution to external forces
	  VolumeForce  = this->CalculateVolumeForce( VolumeForce, Variables.N );

	  Energy = 0;
	  this->CalculateExternalForcesEnergy(Energy, Variables, VolumeForce, IntegrationWeight);

	  rOutput[PointNumber] = Energy;
	}
    }
    else{
      BeamElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }


    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::CalculateOnIntegrationPoints(  const Variable<array_1d<double, 3 > >& rVariable,
								    std::vector< array_1d<double, 3 > >& rOutput,
								    const ProcessInfo& rCurrentProcessInfo )
  {

    KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_points_number )
      rOutput.resize( integration_points_number );


    if ( rVariable == LINEAR_MOMENTUM || rVariable == ANGULAR_MOMENTUM ){

      //create and initialize element variables:
      ElementDataType Variables;
      this->InitializeElementData(Variables,rCurrentProcessInfo);

      IntegrationMethod ThisIntegrationMethod = mThisIntegrationMethod;
      //full quadrature integration:
      mThisIntegrationMethod = mFullIntegrationMethod;

      const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

      // initialize variables short version:

      //Compute Section Properties:
      this->CalculateSectionProperties(Variables.Section);
      Variables.Length = GetGeometry().Length();

      //Set equilibrium point initial:0/mid:0.5/final:1
      if( rCurrentProcessInfo.Has(EQUILIBRIUM_POINT) )
	Variables.Alpha = rCurrentProcessInfo[EQUILIBRIUM_POINT];
      else
	Variables.Alpha = 1;

      //reading shape functions
      Variables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

      //get the shape functions for the order of the integration method [N]
      const Matrix& Ncontainer = Variables.GetShapeFunctions();

      double IntegrationWeight = 0;
      array_1d<double,3> LinearMomentum;
      array_1d<double,3> AngularMomentum;
      LinearMomentum.clear();
      AngularMomentum.clear();

      for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
	{
	  //integration point number
	  Variables.PointNumber = PointNumber;

	  //set shape functions values for this integration point
	  noalias(Variables.N) = matrix_row<const Matrix>( Ncontainer, PointNumber);

	  //compute local to global frame
	  this->CalculateFrameMapping( Variables, PointNumber );

	  Variables.detJ = 1.0/mInvJ0;

	  IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
	  IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );

	  //Momentum calculation
	  this->CalculateMomentumRelations(LinearMomentum, AngularMomentum, Variables, rCurrentProcessInfo, IntegrationWeight);

	  if( rVariable == LINEAR_MOMENTUM )
	    rOutput[PointNumber] = LinearMomentum;

	  if( rVariable == ANGULAR_MOMENTUM )
	    rOutput[PointNumber] = AngularMomentum;
	}

      mThisIntegrationMethod = ThisIntegrationMethod;

    }
    else{
      BeamElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
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
  int  LargeDisplacementBeamElement::Check(const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    // Perform base element checks
    int ErrorCode = 0;
    ErrorCode = BeamElement::Check(rCurrentProcessInfo);

    return ErrorCode;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamElement::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BeamElement )
    int IntMethod = int(mReducedIntegrationMethod);
    rSerializer.save("ReducedIntegrationMethod",IntMethod);
    IntMethod = int(mFullIntegrationMethod);
    rSerializer.save("FullIntegrationMethod",IntMethod);
    rSerializer.save("IterationCounter",mIterationCounter);
    rSerializer.save("InvJ0",mInvJ0);
    rSerializer.save("CurrentCurvatureVectors",mCurrentCurvatureVectors);
    rSerializer.save("FrameQuaternionsReduced",mFrameQuaternionsReduced);
    rSerializer.save("FrameQuaternionsFull",mFrameQuaternionsFull);
  }

  void LargeDisplacementBeamElement::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BeamElement )
    int IntMethod;
    rSerializer.load("ReducedIntegrationMethod",IntMethod);
    mReducedIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("FullIntegrationMethod",IntMethod);
    mFullIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("IterationCounter",mIterationCounter);
    rSerializer.load("InvJ0",mInvJ0);
    rSerializer.load("CurrentCurvatureVectors",mCurrentCurvatureVectors);
    rSerializer.load("FrameQuaternionsReduced",mFrameQuaternionsReduced);
    rSerializer.load("FrameQuaternionsFull",mFrameQuaternionsFull);
  }


} // Namespace Kratos
