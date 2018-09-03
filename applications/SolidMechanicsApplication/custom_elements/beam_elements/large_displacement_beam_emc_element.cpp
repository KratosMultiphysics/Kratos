//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/beam_elements/large_displacement_beam_emc_element.hpp"

#include "solid_mechanics_application_variables.h"


namespace Kratos
{


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  LargeDisplacementBeamEMCElement::LargeDisplacementBeamEMCElement(IndexType NewId,GeometryType::Pointer pGeometry)
    : LargeDisplacementBeamElement(NewId, pGeometry)
  {
  }

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************


  LargeDisplacementBeamEMCElement::LargeDisplacementBeamEMCElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : LargeDisplacementBeamElement(NewId, pGeometry, pProperties)
  {
    KRATOS_TRY

    KRATOS_CATCH( "" )
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  LargeDisplacementBeamEMCElement::LargeDisplacementBeamEMCElement(LargeDisplacementBeamEMCElement const& rOther)
    :LargeDisplacementBeamElement(rOther)
    ,mCurrentStrainResultantsVector(rOther.mCurrentStrainResultantsVector)
    ,mPreviousStrainResultantsVector(rOther.mPreviousStrainResultantsVector)
  {

  }

  //*********************************OPERATIONS*****************************************
  //************************************************************************************

  Element::Pointer LargeDisplacementBeamEMCElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
  {
    return Kratos::make_shared<LargeDisplacementBeamEMCElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  LargeDisplacementBeamEMCElement::~LargeDisplacementBeamEMCElement()
  {
  }

  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamEMCElement::Initialize()
  {
    KRATOS_TRY

    LargeDisplacementBeamElement::Initialize();

    //------------- REDUCED QUADRATURE INTEGRATION

    IntegrationMethod ReducedIntegrationMethod = this->GetReducedIntegrationMethod();

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( ReducedIntegrationMethod );


    //Resultants Initialization
    if ( mCurrentStrainResultantsVector.size() != integration_points.size() )
      {
        mCurrentStrainResultantsVector.resize( integration_points.size() );
      }

    for ( unsigned int i = 0; i < mCurrentStrainResultantsVector.size(); i++ )
      {
	mCurrentStrainResultantsVector[i].resize(3,false);
	noalias(mCurrentStrainResultantsVector[i]) = ZeroVector(3);
      }


    //Resultants Initialization
    if ( mPreviousStrainResultantsVector.size() != integration_points.size() )
      {
        mPreviousStrainResultantsVector.resize( integration_points.size() );
      }

    for ( unsigned int i = 0; i < mPreviousStrainResultantsVector.size(); i++ )
      {
	mPreviousStrainResultantsVector[i].resize(3,false);
	noalias(mPreviousStrainResultantsVector[i]) = ZeroVector(3);
      }


    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamEMCElement::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    LargeDisplacementBeamElement::InitializeSolutionStep(rCurrentProcessInfo);

    for ( unsigned int i = 0; i < mPreviousStrainResultantsVector.size(); i++ )
      {
	mPreviousStrainResultantsVector[i] = mCurrentStrainResultantsVector[i] ;
      }

    KRATOS_CATCH( "" )
  }


  //*****************************************************************************
  //*****************************************************************************

  void LargeDisplacementBeamEMCElement::MapToSpatialFrame(const ElementDataType& rVariables, Matrix& rVariable)
  {

    KRATOS_TRY

    //matrix value :  A = Q * A' * QT

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    unsigned int MatSize = rVariable.size1();

    Matrix AuxiliarRotationMatrix(MatSize,MatSize);
    noalias(AuxiliarRotationMatrix) = ZeroMatrix(MatSize,MatSize);

    //Building the rotation matrix for the local element matrix N+Alpha
    for (unsigned int i=0; i<dimension; i++)
      {
	for(unsigned int j=0; j<dimension; j++)
	  {
	    AuxiliarRotationMatrix(i,j) = rVariables.AlphaRotationMatrix(i,j);
	  }
      }

    for (unsigned int i=0; i<dimension; i++)
      {
	for(unsigned int j=0; j<dimension; j++)
	  {
	    AuxiliarRotationMatrix(i+dimension,j+dimension) = rVariables.AlphaRotationMatrixAsterisk(i,j);
	  }
      }

    //Rotate Local Stiffness Matrix
    Matrix aux_matrix(MatSize,MatSize);
    noalias(aux_matrix) = ZeroMatrix(MatSize,MatSize);
    noalias(aux_matrix) = prod(AuxiliarRotationMatrix, rVariable);


    noalias(AuxiliarRotationMatrix) = ZeroMatrix(MatSize,MatSize);
    //Building the rotation matrix for the local element matrix N+1
    for (unsigned int kk=0; kk < MatSize; kk += dimension)
      {
        for (unsigned int i=0; i<dimension; i++)
	  {
            for(unsigned int j=0; j<dimension; j++)
	      {
		AuxiliarRotationMatrix(i+kk,j+kk) = rVariables.CurrentRotationMatrix(i,j);
	      }
	  }
      }

    //Transformed Matrix
    noalias(rVariable) = ZeroMatrix(MatSize,MatSize);
    noalias(rVariable) = prod(aux_matrix,trans(AuxiliarRotationMatrix));


    KRATOS_CATCH( "" )

  }

  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamEMCElement::InitializeElementData(ElementDataType& rVariables, const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    BeamElement::InitializeElementData(rVariables,rCurrentProcessInfo);

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    //Set equilibrium point initial:0/mid:0.5/final:1
    if( rCurrentProcessInfo.Has(EQUILIBRIUM_POINT) )
      rVariables.Alpha = rCurrentProcessInfo[EQUILIBRIUM_POINT];
    else
      rVariables.Alpha = 1;


    rVariables.DeltaTime = rCurrentProcessInfo[DELTA_TIME];

    rVariables.PreviousAxisPositionDerivatives.resize( dimension );
    rVariables.PreviousRotationMatrix.resize( dimension, dimension );

    KRATOS_CATCH( "" )
  }

  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************

  void LargeDisplacementBeamEMCElement::CalculateKinematics(ElementDataType& rVariables, const unsigned int& rPointNumber)
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

    if(rVariables.CurrentAxisPositionDerivatives.size() != dimension)
      rVariables.CurrentAxisPositionDerivatives.resize(dimension,false);

    if(rVariables.PreviousAxisPositionDerivatives.size() != dimension)
      rVariables.PreviousAxisPositionDerivatives.resize(dimension,false);

    noalias(rVariables.CurrentAxisPositionDerivatives) = ZeroVector(dimension);
    noalias(rVariables.PreviousAxisPositionDerivatives) = ZeroVector(dimension);

    //compute local to global frame
    this->CalculateFrameMapping(rVariables, rPointNumber );

    Vector CurrentValueVector(3);
    noalias(CurrentValueVector) = ZeroVector(3);


    //strains due to displacements and rotations

    if( this->Is(BeamElement::FINALIZED_STEP) ){

        //rVariables.DeltaPosition = this->CalculateDeltaPosition(rVariables.DeltaPosition);

        Matrix PreviousDeltaPosition;
        PreviousDeltaPosition = CalculatePreviousDeltaPosition(PreviousDeltaPosition);


	for ( unsigned int i = 0; i < number_of_nodes; i++ )
	{
	    //A: Current Nodes Position
	    CurrentValueVector = GetGeometry()[i].Coordinates();

	    for ( unsigned int j = 0; j < dimension; j++ )
	    {
		CurrentValueVector[j] -= rVariables.DeltaPosition(i,j);
	    }

	    CurrentValueVector = this->MapToInitialLocalFrame( CurrentValueVector, rPointNumber );

	    //Current Frame Axis Position derivative
	    for( unsigned int j = 0; j < dimension; j++ )
	    {
		rVariables.CurrentAxisPositionDerivatives[j] +=  rVariables.DN_DX(i,0) * ( CurrentValueVector[j] );
	    }


	    //B: Previous Nodes Position
	    CurrentValueVector = GetGeometry()[i].Coordinates();

	    for ( unsigned int j = 0; j < 3; j++ )
	    {
              CurrentValueVector[j] -= (rVariables.DeltaPosition(i,j) + PreviousDeltaPosition(i,j));
	    }

	    CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rPointNumber );

	    //Previous Frame Axis Position derivative
	    rVariables.PreviousAxisPositionDerivatives +=  rVariables.DN_DX(i,0) * ( CurrentValueVector );

	}
    }
    else{

	for ( unsigned int i = 0; i < number_of_nodes; i++ )
	{

	    //A: Current Nodes Position
	    CurrentValueVector = GetGeometry()[i].Coordinates();
	    CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rPointNumber );

	    //Current Frame Axis Position derivative
	    rVariables.CurrentAxisPositionDerivatives +=  rVariables.DN_DX(i,0) * ( CurrentValueVector );


	    //B: Previous Nodes Position
	    CurrentValueVector = GetGeometry()[i].Coordinates();

	    for ( unsigned int j = 0; j < 3; j++ )
	    {
		CurrentValueVector[j] -= rVariables.DeltaPosition(i,j);
	    }

	    CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rPointNumber );

	    //Previous Frame Axis Position derivative
	    rVariables.PreviousAxisPositionDerivatives +=  rVariables.DN_DX(i,0) * ( CurrentValueVector );

	}
    }

    //*************************************//

    //Compute current CURVATURES
    if( this->Is(BeamElement::FINALIZED_STEP) ){

	//set current STRAIN RESULTANTS
	rVariables.CurrentStrainResultantsVector  = mCurrentStrainResultantsVector[rPointNumber];
	rVariables.PreviousStrainResultantsVector = mPreviousStrainResultantsVector[rPointNumber];

	//set current CURVATURES
	rVariables.CurrentCurvatureVector  = mCurrentCurvatureVectors[rPointNumber];
	rVariables.PreviousCurvatureVector = mCurrentCurvatureVectors[rPointNumber];

    }
    else{

	//set current STRAIN RESULTANTS
	rVariables.CurrentStrainResultantsVector  = mPreviousStrainResultantsVector[rPointNumber];
	rVariables.PreviousStrainResultantsVector = mPreviousStrainResultantsVector[rPointNumber];

	//set current CURVATURES
	rVariables.CurrentCurvatureVector  = mPreviousCurvatureVectors[rPointNumber];
	rVariables.PreviousCurvatureVector = mPreviousCurvatureVectors[rPointNumber];
    }

    KRATOS_CATCH( "" )
  }

  //*************************COMPUTE PREVIOUS DELTA POSITION****************************
  //************************************************************************************


  Matrix& LargeDisplacementBeamEMCElement::CalculatePreviousDeltaPosition(Matrix & rDeltaPosition)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().PointsNumber();
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    if(rDeltaPosition.size1() != number_of_nodes || rDeltaPosition.size2() !=  dimension)
      rDeltaPosition.resize(number_of_nodes,dimension,false);

    //noalias(rDeltaPosition) = ZeroMatrix(number_of_nodes,dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        array_1d<double, 3 > & CurrentStepDisplacement = GetGeometry()[i].FastGetSolutionStepValue(STEP_DISPLACEMENT,1);

       for ( unsigned int j = 0; j < dimension; j++ )
	  {
	    rDeltaPosition(i,j) = CurrentStepDisplacement[j];
	  }

      }

    return rDeltaPosition;

    KRATOS_CATCH( "" )

  }

  //*************************COMPUTE FRAME MAPPING*************************************
  //************************************************************************************

  void LargeDisplacementBeamEMCElement::CalculateFrameMapping(ElementDataType& rVariables,const unsigned int& rPointNumber)
  {

    if( mThisIntegrationMethod == this->mReducedIntegrationMethod ){
      mFrameQuaternionsReduced[rPointNumber].ToRotationMatrix(rVariables.PreviousRotationMatrix);
    }

    if( mThisIntegrationMethod == this->mFullIntegrationMethod ){
      mFrameQuaternionsFull[rPointNumber].ToRotationMatrix(rVariables.PreviousRotationMatrix);
    }

    Vector CurrentStepRotationVector(3);
    noalias(CurrentStepRotationVector) = ZeroVector(3);
    this->GetLocalCurrentValue(STEP_ROTATION, CurrentStepRotationVector, rVariables.N);

    Matrix CayleyRotationMatrix(3,3);
    noalias(CayleyRotationMatrix) = ZeroMatrix(3,3);

    if(rVariables.Alpha == 1){ //quasi-static case exponential update
      BeamMathUtilsType::ExponentialTransform( CurrentStepRotationVector, CayleyRotationMatrix );
    }
    else{
      BeamMathUtilsType::CayleyTransform( CurrentStepRotationVector, CayleyRotationMatrix );
    }

    rVariables.CurrentRotationMatrix = prod(CayleyRotationMatrix, rVariables.PreviousRotationMatrix);

    //*------------------------------*//

    CalculateAlphaRotationMatrix( rVariables.PreviousRotationMatrix, rVariables.CurrentRotationMatrix, rVariables.AlphaRotationMatrix, rVariables.AlphaRotationMatrixAsterisk, rVariables.Alpha);

    //*------------------------------*//

  }


  //*********************************SET STRAIN VARIABLES*******************************
  //************************************************************************************

  void LargeDisplacementBeamEMCElement::UpdateStrainVariables(ElementDataType& rVariables, const unsigned int& rPointNumber)
  {
    KRATOS_TRY

    LargeDisplacementBeamElement::UpdateStrainVariables(rVariables, rPointNumber);

    mCurrentStrainResultantsVector[rPointNumber] = rVariables.CurrentStrainResultantsVector;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamEMCElement::CalculateAlphaRotationMatrix( const Matrix& rPreviousRotationMatrix,
								      const Matrix& rCurrentRotationMatrix,
								      Matrix& rAlphaRotationMatrix,
								      Matrix& rAlphaRotationMatrixAsterisk,
								      double alpha)
  {
    KRATOS_TRY

    rAlphaRotationMatrix = (1-alpha) * (rPreviousRotationMatrix) + alpha * rCurrentRotationMatrix;

    double DetRotationMatrix = 0;
    MathUtils<double>::InvertMatrix3( rAlphaRotationMatrix, rAlphaRotationMatrixAsterisk, DetRotationMatrix);
    rAlphaRotationMatrixAsterisk = DetRotationMatrix *  trans(rAlphaRotationMatrixAsterisk);

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamEMCElement::CalculateConstitutiveMatrix(ElementDataType& rVariables)
  {
    KRATOS_TRY

    // Material Elastic constitutive matrix
    this->CalculateMaterialConstitutiveMatrix(rVariables.ConstitutiveMatrix, rVariables);

    //Spatial Elastic constitutive matrix
    this->MapToSpatialFrame( rVariables, rVariables.ConstitutiveMatrix);

    rVariables.ConstitutiveMatrix *= rVariables.Alpha;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamEMCElement::CalculateStrainResultants(Vector& rStrainResultants, ElementDataType& rVariables, double alpha)
  {
    KRATOS_TRY

    // current Strain Resultants  N+alpha
    if( rStrainResultants.size() != 3 )
      rStrainResultants.resize(3, false);

    noalias(rStrainResultants) = ZeroVector(3);

    //OPTION ENERGY CONSERVATION START
    Vector CurrentStrainResultantsVector(3);
    noalias(CurrentStrainResultantsVector) = ZeroVector(3);
    CalculateCurrentStrainResultantsVector( rVariables, CurrentStrainResultantsVector, alpha );

    rStrainResultants  = (1-alpha) * rVariables.PreviousStrainResultantsVector + alpha * CurrentStrainResultantsVector;
    //OPTION ENERGY CONSERVATION END


    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamEMCElement::CalculateStrainCouples(Vector& rStrainCouples, ElementDataType& rVariables, double alpha)
  {
    KRATOS_TRY

    // current Strain Couples  N+alpha
    if( rStrainCouples.size() != 3 )
      rStrainCouples.resize(3, false);

    noalias(rStrainCouples) = ZeroVector(3);

    //------------------------

    //OPTION ENERGY CONSERVATION START
    Vector CurrentCurvatureVector(3);
    noalias(CurrentCurvatureVector) = ZeroVector(3);
    CalculateCurrentCurvatureVector( rVariables, CurrentCurvatureVector, alpha );

    rStrainCouples = (1-alpha) * rVariables.PreviousCurvatureVector + alpha * CurrentCurvatureVector;
    //OPTION ENERGY CONSERVATION END


    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamEMCElement::CalculateCurrentStrainResultantsVector(ElementDataType& rVariables,
									       Vector& rCurrentStrainResultantsVector,
									       double Alpha)
  {
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    Vector CurrentStepDisplacementDerivativesVector(3);
    noalias(CurrentStepDisplacementDerivativesVector) = ZeroVector(3);
    Vector CurrentStepRotationDerivativesVector(3);
    noalias(CurrentStepRotationDerivativesVector) = ZeroVector(3);
    Vector CurrentStepRotationVector(3);
    noalias(CurrentStepRotationVector) = ZeroVector(3);
    Vector CurrentStepDisplacementVector(3);
    noalias(CurrentStepDisplacementVector) = ZeroVector(3);
    Vector CurrentValueVector(3);
    noalias(CurrentValueVector) = ZeroVector(3);

    const SizeType number_of_nodes  = GetGeometry().size();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
    	//Current Step Rotation Derivatives
    	CurrentValueVector = GetNodalCurrentValue( STEP_ROTATION, CurrentValueVector, i );
    	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rVariables.PointNumber );

	CurrentStepRotationVector +=  rVariables.N[i] * CurrentValueVector;

    	CurrentStepRotationDerivativesVector += rVariables.DN_DX(i,0) * ( CurrentValueVector );

    	//Current Step Displacement Derivatives
	CurrentValueVector = GetNodalCurrentValue( STEP_DISPLACEMENT, CurrentValueVector, i );

	for ( unsigned int j = 0; j < dimension; j++ )
	  {
	    CurrentValueVector[j] = rVariables.DeltaPosition(i,j);
	  }

    	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rVariables.PointNumber );

	CurrentStepDisplacementVector +=  rVariables.N[i] * CurrentValueVector;

    	CurrentStepDisplacementDerivativesVector += rVariables.DN_DX(i,0) * ( CurrentValueVector );
      }

    if( Alpha != 1 ){ //dynamic case approach by simo.

      Matrix CayleyRotationMatrix(3,3);
      noalias(CayleyRotationMatrix) = ZeroMatrix(3,3);
      noalias(CurrentStepDisplacementDerivativesVector) = ZeroVector(3);
      noalias(CurrentStepRotationVector) = ZeroVector(3);
      Vector PreviousValueVector(3);
      noalias(PreviousValueVector) = ZeroVector(3);

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        //Current Linear Velocity Vector
        CurrentValueVector = GetNodalCurrentValue( VELOCITY, CurrentValueVector, i );
        CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rVariables.PointNumber );

        //Previous Linear Velocity Vector
        PreviousValueVector = GetNodalPreviousValue( VELOCITY, PreviousValueVector, i );
        PreviousValueVector = MapToInitialLocalFrame( PreviousValueVector, rVariables.PointNumber );

        CurrentStepDisplacementDerivativesVector += rVariables.DN_DX(i,0) * (CurrentValueVector + PreviousValueVector);


        //Current Step Rotation Derivatives
    	CurrentValueVector = GetNodalCurrentValue( STEP_ROTATION, CurrentValueVector, i );
    	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rVariables.PointNumber );

        BeamMathUtilsType::CayleyTransform( CurrentValueVector, CayleyRotationMatrix );

        //Current Angular Velocity Vector
        CurrentValueVector = GetNodalCurrentValue( ANGULAR_VELOCITY, CurrentValueVector, i );
        CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rVariables.PointNumber );

        //Previous Angular Velocity Vector
        PreviousValueVector = GetNodalPreviousValue( ANGULAR_VELOCITY, PreviousValueVector, i );
        PreviousValueVector = MapToInitialLocalFrame( PreviousValueVector, rVariables.PointNumber );

        CurrentStepRotationVector += rVariables.N[i] * (CurrentValueVector + prod( CayleyRotationMatrix, PreviousValueVector ));
      }

      CurrentStepDisplacementDerivativesVector *= 0.5 * rVariables.DeltaTime;
      CurrentStepRotationVector *= 0.5 * rVariables.DeltaTime;
    }


    double alpha = 0.5; //Alpha;   //quasi-static

    Matrix AlphaRotationMatrix(3,3);
    noalias(AlphaRotationMatrix) = ZeroMatrix(3,3);
    Matrix AlphaRotationMatrixAsterisk(3,3);
    noalias(AlphaRotationMatrixAsterisk) = ZeroMatrix(3,3);

    CalculateAlphaRotationMatrix( rVariables.PreviousRotationMatrix, rVariables.CurrentRotationMatrix, AlphaRotationMatrix, AlphaRotationMatrixAsterisk, alpha);

    Vector AxisPositionDerivativesAlpha = (1-alpha) * rVariables.PreviousAxisPositionDerivatives + alpha * rVariables.CurrentAxisPositionDerivatives;

    Vector CurrentStepxAxisPosition;
    MathUtils<double>::CrossProduct(CurrentStepxAxisPosition,CurrentStepRotationVector, AxisPositionDerivativesAlpha);

    Vector AxisDisplacementDerivatives = CurrentStepDisplacementDerivativesVector - CurrentStepxAxisPosition;


    //std::cout<<" ID "<<this->Id()<<" Previous "<<rVariables.PreviousStrainResultantsVector<<std::endl;

    //std::cout<<" Displacement "<<CurrentStepDisplacementVector<<std::endl;
    //-----------
    // CASE A:
    Vector CurrentStrainResultantsVectorA(3);
    noalias(CurrentStrainResultantsVectorA) = ZeroVector(3);
    CurrentStrainResultantsVectorA  = rVariables.PreviousStrainResultantsVector;
    CurrentStrainResultantsVectorA += prod( trans(AlphaRotationMatrix), AxisDisplacementDerivatives );

    //std::cout<<" StrainResultant A: "<< CurrentStrainResultantsVectorA <<std::endl;

    //-----------
    // CASE B:
    Vector CurrentStrainResultantsVectorB(3);
    noalias(CurrentStrainResultantsVectorB) = ZeroVector(3);
    CurrentStrainResultantsVectorB = rVariables.PreviousStrainResultantsVector;

    Vector DeltaAxisPositionDerivatives = rVariables.CurrentAxisPositionDerivatives - rVariables.PreviousAxisPositionDerivatives;

    CurrentStrainResultantsVectorB += prod( trans(AlphaRotationMatrix), CurrentStepDisplacementDerivativesVector );

    Matrix DeltaRotationMatrix = rVariables.CurrentRotationMatrix - rVariables.PreviousRotationMatrix;
    CurrentStrainResultantsVectorB += prod( trans(DeltaRotationMatrix), AxisPositionDerivativesAlpha );

    //std::cout<<" StrainResultant B: "<< CurrentStrainResultantsVectorB <<std::endl;

    //-----------
    // CASE C:
    Vector CurrentStrainResultantsVectorC(3);
    noalias(CurrentStrainResultantsVectorC) = ZeroVector(3);

    // Axis vector in the REFERENCE initial frame
    Vector E1(3);
    noalias(E1) = ZeroVector(3);
    E1[0] = 1.0;

    CurrentStrainResultantsVectorC = prod( trans(rVariables.CurrentRotationMatrix), rVariables.CurrentAxisPositionDerivatives ) - E1;

    //std::cout<<" StrainResultant C: "<< CurrentStrainResultantsVectorC <<std::endl;

    //dynamic and energy cases compatible
    rCurrentStrainResultantsVector = CurrentStrainResultantsVectorA;

    if( rVariables.Alpha == Alpha && rVariables.Alpha == 1){ //quasi-static cases compatible
      rCurrentStrainResultantsVector = CurrentStrainResultantsVectorC;
      //std::cout<<" QUASI-STATIC CASE STRESS RESULTANTS "<<std::endl;
    }

    //std::cout<<" CurrentStrainResultants "<<rCurrentStrainResultantsVector<<std::endl;


    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamEMCElement::CalculateCurrentCurvatureVector(ElementDataType& rVariables,
									Vector& rCurrentCurvatureVector,
									double Alpha)
  {
    KRATOS_TRY

    Vector CurrentStepRotationDerivativesVector(3);
    noalias(CurrentStepRotationDerivativesVector) = ZeroVector(3);
    Vector CurrentValueVector(3);
    noalias(CurrentValueVector) = ZeroVector(3);

    const SizeType number_of_nodes  = GetGeometry().size();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
    	//Current Rotation Derivatives
    	CurrentValueVector = GetNodalCurrentValue( STEP_ROTATION, CurrentValueVector, i );
    	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rVariables.PointNumber );

    	CurrentStepRotationDerivativesVector += rVariables.DN_DX(i,0) * ( CurrentValueVector );
      }

    // energy calculation and static cases compatible
    double alpha = Alpha;
    if( rVariables.Alpha != Alpha )
      alpha = rVariables.Alpha;


    Matrix AlphaRotationMatrix(3,3);
    noalias(AlphaRotationMatrix) = ZeroMatrix(3,3);
    Matrix AlphaRotationMatrixAsterisk(3,3);
    noalias(AlphaRotationMatrixAsterisk) = ZeroMatrix(3,3);

    CalculateAlphaRotationMatrix( rVariables.PreviousRotationMatrix, rVariables.CurrentRotationMatrix, AlphaRotationMatrix, AlphaRotationMatrixAsterisk, alpha);

    rCurrentCurvatureVector  = rVariables.PreviousCurvatureVector;
    rCurrentCurvatureVector += prod( trans(AlphaRotationMatrixAsterisk), CurrentStepRotationDerivativesVector );

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamEMCElement::CalculateStressResultants(ElementDataType& rVariables, const unsigned int& rPointNumber)
  {
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    //compute Strain Resultants and Couples
    Vector StrainResultants(3);
    noalias(StrainResultants) = ZeroVector(3);
    Vector StrainCouples(3);
    noalias(StrainCouples) = ZeroVector(3);

    CalculateStrainResultants(StrainResultants, rVariables, rVariables.Alpha);
    CalculateStrainCouples(StrainCouples, rVariables, rVariables.Alpha);

    //std::cout<<" CurrentAxisPositionDerivatives "<<rVariables.CurrentAxisPositionDerivatives<<std::endl;

    //std::cout<<" CurrentAxisDer "<<rVariables.CurrentAxisPositionDerivatives<<" R "<<rVariables.CurrentRotationMatrix<<std::endl;

    for ( unsigned int i = 0; i < dimension; i++ )
      {
	rVariables.StrainVector[i]   = StrainResultants[i];
	rVariables.StrainVector[i+3] = StrainCouples[i];
      }

    //----------------

    Matrix ConstitutiveMatrix(6,6);
    noalias(ConstitutiveMatrix) = ZeroMatrix(6,6);
    this->CalculateMaterialConstitutiveMatrix(ConstitutiveMatrix, rVariables);

    //Reference Stress Vector
    rVariables.StressVector = prod( ConstitutiveMatrix, rVariables.StrainVector );

    //std::cout<<" Stress "<<rVariables.StressVector<<" Strain "<<rVariables.StrainVector<<std::endl;

    Vector StressResultants(3);
    noalias(StressResultants) = ZeroVector(3);
    Vector StressCouples(3);
    noalias(StressCouples) = ZeroVector(3);

    for ( unsigned int i = 0; i < dimension; i++ )
      {
    	StressResultants[i] = rVariables.StressVector[i];
    	StressCouples[i]    = rVariables.StressVector[i+3];
      }

    //----------------

    //Current frame given by the Frame Rotation
    StressResultants = prod(rVariables.AlphaRotationMatrix, StressResultants);
    StressCouples    = prod(rVariables.AlphaRotationMatrixAsterisk, StressCouples);

    for ( unsigned int i = 0; i < dimension; i++ )
      {
    	rVariables.StressVector[i]   = StressResultants[i];
    	rVariables.StressVector[i+3] = StressCouples[i];
      }


    //Set current curvatures and strain resultans:
    CalculateCurrentStrainResultantsVector(rVariables, rVariables.CurrentStrainResultantsVector, rVariables.Alpha);
    CalculateCurrentCurvatureVector(rVariables, rVariables.CurrentCurvatureVector, rVariables.Alpha);

    //set variables for the initialization update
    this->UpdateStrainVariables(rVariables,rPointNumber);


    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  //Strain Energy Calculation
  void LargeDisplacementBeamEMCElement::CalculateStrainEnergy(double& rEnergy, ElementDataType& rVariables, const ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight)
  {
    KRATOS_TRY

    //Internal Energy Calculation: alpha = 1
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    //compute Strain Resultants and Couples
    Vector StrainResultants(3);
    noalias(StrainResultants) = ZeroVector(3);
    Vector StrainCouples(3);
    noalias(StrainCouples) = ZeroVector(3);

    double Alpha = 1.0;
    CalculateStrainResultants(StrainResultants, rVariables, Alpha);
    CalculateStrainCouples(StrainCouples, rVariables, Alpha);

    //----------------

    Matrix ConstitutiveMatrix(6,6);
    noalias(ConstitutiveMatrix) = ZeroMatrix(6,6);
    this->CalculateMaterialConstitutiveMatrix(ConstitutiveMatrix, rVariables);

    //----------------

    Vector StrainVector(6);
    noalias(StrainVector) = ZeroVector(6);
    for ( unsigned int i = 0; i < dimension; i++ )
      {
    	StrainVector[i]   = StrainResultants[i];
    	StrainVector[i+3] = StrainCouples[i];
      }

    //Reference Stress Vector
    Vector StressVector = prod( ConstitutiveMatrix, StrainVector );

    //std::cout<<" Strain "<<StrainVector<<" Stress "<<StressVector<<std::endl;


    rEnergy += 0.5 * (inner_prod(StressVector, StrainVector)) * rIntegrationWeight ;

    //std::cout<<" StrainEnergy "<<rEnergy<<" rIntegrationWeight "<<rIntegrationWeight<<std::endl;

    KRATOS_CATCH( "" )
  }



  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamEMCElement::CalculateAndAddFollowerForces(VectorType& rRightHandSideVector,
								      ElementDataType & rVariables,
								      double& rIntegrationWeight)
  {
    KRATOS_TRY


    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamEMCElement::CalculateDifferentialOperator(MatrixType& rDifferentialOperator,
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

    Vector AxisPositionDerivativesAlpha = (1-alpha) * rVariables.PreviousAxisPositionDerivatives + alpha * rVariables.CurrentAxisPositionDerivatives;

    BeamMathUtilsType::VectorToSkewSymmetricTensor(AxisPositionDerivativesAlpha, SkewSymResultants);

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


  void LargeDisplacementBeamEMCElement::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
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
    Vector StressCouples(3);
    for ( unsigned int i = 0; i < 3; i++ )
      {
	StressResultants[i] = rVariables.StressVector[i];
        StressCouples[i] = rVariables.StressVector[i+3];
      }

    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    //NOTE: avoid Kuug noise in plane ploblems
    if( fabs(inner_prod(StressResultants,StressCouples)) < 1e-15 )
      noalias(StressResultants) = ZeroVector(3);

    //Get frame step rotation
    Vector CurrentStepRotation(3);
    noalias(CurrentStepRotation) = ZeroVector(3);

    this->GetLocalCurrentValue(STEP_ROTATION, CurrentStepRotation, rVariables.N);

    Matrix SkewSymStepRotation(3,3);
    noalias(SkewSymStepRotation) = ZeroMatrix(3,3);
    BeamMathUtilsType::VectorToSkewSymmetricTensor(CurrentStepRotation, SkewSymStepRotation);

    Vector AxisPositionDerivativesAlpha = (1-rVariables.Alpha) * rVariables.PreviousAxisPositionDerivatives + rVariables.Alpha * rVariables.CurrentAxisPositionDerivatives;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {

	RowIndex = i * (dimension * 2);

	for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {

	    noalias(Kij) = ZeroMatrix(6,6);

	    ColIndex = j * (dimension * 2);


	    //term 11 -> 0
	    //term 12
	    noalias(GabK) = ZeroMatrix(3,3);
	    Matrix SkewSymStressResultants(3,3);
	    noalias(SkewSymStressResultants) = ZeroMatrix(3,3);
	    BeamMathUtilsType::VectorToSkewSymmetricTensor(StressResultants, SkewSymStressResultants);
	    GabK = SkewSymStressResultants;
	    Vector CurrentValueVector = prod( SkewSymStepRotation, StressResultants );
	    Matrix SkewSymValueVector(3,3);
	    noalias(SkewSymValueVector)= ZeroMatrix(3,3);
	    BeamMathUtilsType::VectorToSkewSymmetricTensor(CurrentValueVector, SkewSymValueVector);
	    GabK += (1-rVariables.Alpha) * SkewSymValueVector;
	    GabK *= (-1) * (rVariables.DN_DX(i, 0) * rVariables.N[j]);

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

	    GabK  =  SkewSymStressCouples;

	    GabK +=  (1-rVariables.Alpha) * outer_prod( StressCouples, CurrentStepRotation );

	    GabK +=  (1-rVariables.Alpha) * inner_prod( StressCouples, CurrentStepRotation ) * DiagonalMatrix;

	    GabK *= (-1) * (rVariables.DN_DX(i, 0) * rVariables.N[j]);

	    CurrentValueVector = StressResultants + (1-rVariables.Alpha) * prod( SkewSymStepRotation, StressResultants );

	    GabK += ( rVariables.N[i] * rVariables.N[j]) * outer_prod( CurrentValueVector, AxisPositionDerivativesAlpha );

	    GabK -= ( rVariables.N[i] * rVariables.N[j]) * inner_prod( CurrentValueVector, AxisPositionDerivativesAlpha ) * DiagonalMatrix;

	    //Building the Local Stiffness Matrix
	    BeamMathUtilsType::AddMatrix( Kij, GabK, 3, 3 );

	    Kij *= rIntegrationWeight * rVariables.Alpha;

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

  void LargeDisplacementBeamEMCElement::CalculateAndAddKuuf(MatrixType& rLeftHandSideMatrix,
							    ElementDataType& rVariables,
							    double& rIntegrationWeight)
  {
    KRATOS_TRY


    KRATOS_CATCH( "" )

  }

  //************************************************************************************
  //************************************************************************************

  //Inertia in the SPATIAL configuration
  void LargeDisplacementBeamEMCElement::CalculateAndAddInertiaLHS(MatrixType& rLeftHandSideMatrix, ElementDataType& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight )
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

    Vector CurrentValueVector(3);
    noalias(CurrentValueVector) = ZeroVector(3);
    Vector CurrentAngularVelocityVector(3);
    noalias(CurrentAngularVelocityVector) = ZeroVector(3);
    Vector CurrentStepRotation(3);
    noalias(CurrentStepRotation) = ZeroVector(3);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	//Current Angular Velocity Vector
	CurrentValueVector = GetNodalCurrentValue( ANGULAR_VELOCITY, CurrentValueVector, i );
	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rVariables.PointNumber );
	CurrentAngularVelocityVector += rVariables.N[i] * CurrentValueVector;

	//Get frame step rotation
	CurrentValueVector = GetNodalCurrentValue( STEP_ROTATION, CurrentValueVector, i );
	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rVariables.PointNumber );
	CurrentStepRotation += rVariables.N[i] * CurrentValueVector;
      }

    //std::cout<<" Current Rotation "<< CurrentStepRotation<<std::endl;

    //block m(1,1) of the mass matrix

    MatrixType m11(3,3);
    noalias(m11) = ZeroMatrix(3,3);

    double TotalMass = 0;
    TotalMass  = this->CalculateTotalMass( Section, TotalMass );

    //block m(2,2) of the mass matrix

    MatrixType m22(3,3);
    noalias(m22) = ZeroMatrix(3,3);

    //2.-Get inertia dyadic
    Matrix InertiaDyadic(3,3);
    noalias(InertiaDyadic) = ZeroMatrix(3,3);

    this->CalculateInertiaDyadic( Section, InertiaDyadic );

    Matrix CurrentInertiaDyadic = prod(rVariables.CurrentRotationMatrix,InertiaDyadic);
    CurrentInertiaDyadic = prod(CurrentInertiaDyadic,trans(rVariables.CurrentRotationMatrix));

    //std::cout<<" InertiaDyadic "<<InertiaDyadic<<" TotalMass "<<TotalMass<<std::endl;


    // Compute Linear Part of the Step Rotation
    Matrix LinearPartRotationTensor(3,3);
    noalias(LinearPartRotationTensor) = ZeroMatrix(3,3);

    double NormCurrentStepRotation =  norm_2(CurrentStepRotation);
    if( NormCurrentStepRotation != 0 ){

      this->CalculateRotationLinearPartTensor( CurrentStepRotation , LinearPartRotationTensor );

    }
    else{
      noalias(LinearPartRotationTensor) = IdentityMatrix(3);
    }

    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    Matrix DiagonalMatrix(3,3);
    noalias(DiagonalMatrix) = IdentityMatrix(3);

    Matrix SkewSymAngularMomentum(3,3);
    noalias(SkewSymAngularMomentum) = ZeroMatrix(3,3);
    Vector AngularMomentumVector(3);
    noalias(AngularMomentumVector) = ZeroVector(3);

    AngularMomentumVector = prod( CurrentInertiaDyadic, CurrentAngularVelocityVector );
    BeamMathUtilsType::VectorToSkewSymmetricTensor( AngularMomentumVector, SkewSymAngularMomentum );


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	noalias(m11) = ZeroMatrix(3,3);
	noalias(m22) = ZeroMatrix(3,3);

	RowIndex = i * (dimension * 2);

	for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {

	    ColIndex = j * (dimension * 2);

	    //complete mass matrix integration
	    m11 = TotalMass * rVariables.N[i] * rVariables.N[j] * rIntegrationWeight * DiagonalMatrix;

	    m22  = prod( CurrentInertiaDyadic, trans(LinearPartRotationTensor) );

	    m22 -= 0.5 * DeltaTime * SkewSymAngularMomentum;
	    m22 *= rVariables.N[i] * rVariables.N[j] * rIntegrationWeight;

	    m11 *= 2.0 / (DeltaTime * DeltaTime);
	    m22 *= 2.0 / (DeltaTime * DeltaTime);


	    //Building the Local Tangent Inertia Matrix
	    BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, m11, RowIndex, ColIndex );
	    BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, m22, RowIndex+3, ColIndex+3 );

	  }


      }

    //std::cout<<" rLeftHandSideDynamic "<<rLeftHandSideMatrix<<std::endl;

    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************

  //Inertia in the SPATIAL configuration
  void LargeDisplacementBeamEMCElement::CalculateAndAddInertiaRHS(VectorType& rRightHandSideVector, ElementDataType& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize            = number_of_nodes * ( dimension * 2 );

    if(rRightHandSideVector.size() != MatSize)
      rRightHandSideVector.resize(MatSize, false);

    noalias(rRightHandSideVector) = ZeroVector( MatSize );

    SectionProperties Section;
    this->CalculateSectionProperties(Section);

    //rCurrentProcessInfo must give it:
    double DeltaTime = rCurrentProcessInfo[DELTA_TIME];

    double TotalMass = 0;
    TotalMass = this->CalculateTotalMass( Section, TotalMass );

    Vector CurrentValueVector(3);
    noalias(CurrentValueVector) = ZeroVector(3);
    Vector CurrentLinearVelocityVector(3);
    noalias(CurrentLinearVelocityVector) = ZeroVector(3);
    Vector PreviousLinearVelocityVector(3);
    noalias(PreviousLinearVelocityVector) = ZeroVector(3);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	//Current Linear Velocity Vector
	CurrentValueVector = GetNodalCurrentValue( VELOCITY, CurrentValueVector, i );
	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rVariables.PointNumber );
	CurrentLinearVelocityVector  += rVariables.N[i] * CurrentValueVector;

	//Previous Linear Velocity Vector
	CurrentValueVector = GetNodalPreviousValue( VELOCITY, CurrentValueVector, i );
	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rVariables.PointNumber );
	PreviousLinearVelocityVector += rVariables.N[i] * CurrentValueVector;
      }

    //Compute Linear Term:
    Vector LinearInertialForceVector(3);
    noalias(LinearInertialForceVector) = ZeroVector(3);
    LinearInertialForceVector  = TotalMass * (CurrentLinearVelocityVector - PreviousLinearVelocityVector);


    //-----------------

    Vector CurrentAngularVelocityVector(3);
    noalias(CurrentAngularVelocityVector) = ZeroVector(3);
    Vector PreviousAngularVelocityVector(3);
    noalias(PreviousAngularVelocityVector) = ZeroVector(3);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	//Current Linear Velocity Vector
	CurrentValueVector = GetNodalCurrentValue( ANGULAR_VELOCITY, CurrentValueVector, i );
	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rVariables.PointNumber );
	CurrentAngularVelocityVector  += rVariables.N[i] * CurrentValueVector;

	//Previous Linear Velocity Vector
	CurrentValueVector = GetNodalPreviousValue( ANGULAR_VELOCITY, CurrentValueVector, i );
	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rVariables.PointNumber );
	PreviousAngularVelocityVector += rVariables.N[i] * CurrentValueVector;
      }



    //Compute Angular Term:

    //Get inertia dyadic
    Matrix InertiaDyadic(3,3);
    noalias(InertiaDyadic) = ZeroMatrix(3,3);
    this->CalculateInertiaDyadic( Section, InertiaDyadic );

    //std::cout<<" I "<<InertiaDyadic<<std::endl;

    //Current Inertia Dyadic
    Matrix CurrentInertiaDyadic = prod(rVariables.CurrentRotationMatrix,InertiaDyadic);
    CurrentInertiaDyadic = prod(CurrentInertiaDyadic,trans(rVariables.CurrentRotationMatrix));

    //Previous Inertia Dyadic
    Matrix PreviousInertiaDyadic = prod(rVariables.PreviousRotationMatrix,InertiaDyadic);
    PreviousInertiaDyadic = prod(PreviousInertiaDyadic,trans(rVariables.PreviousRotationMatrix));


    Vector AngularInertialForceVector(3);
    noalias(AngularInertialForceVector) = ZeroVector(3);
    AngularInertialForceVector = prod( CurrentInertiaDyadic, CurrentAngularVelocityVector ) - prod( PreviousInertiaDyadic, PreviousAngularVelocityVector );


    // Build incremental momentum vector
    Vector TotalInertialForceVector(6);
    noalias(TotalInertialForceVector) = ZeroVector(6);

    BeamMathUtilsType::AddVector(LinearInertialForceVector, TotalInertialForceVector, 0);
    BeamMathUtilsType::AddVector(AngularInertialForceVector, TotalInertialForceVector, 3);


    //-----------------
    VectorType Fi(6);
    noalias(Fi) = ZeroVector(6);

    //Initialize Local Matrices
    unsigned int RowIndex = 0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {

	RowIndex = i * (dimension * 2);

	noalias(Fi) = ZeroVector(6);

	//complete matrix
	Fi = TotalInertialForceVector * rVariables.N[i] * rIntegrationWeight / DeltaTime;

	BeamMathUtilsType::AddVector(Fi, rRightHandSideVector, RowIndex);
      }

    //std::cout<<" rRightHandSideVector "<<rRightHandSideVector<<std::endl;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementBeamEMCElement::CalculateRotationLinearPartTensor(Vector& rRotationVector, Matrix& rRotationTensor)

  {
    KRATOS_TRY

    if( rRotationTensor.size1() != 3 )
      rRotationTensor.resize(3, 3, false);

    noalias(rRotationTensor) = ZeroMatrix(3,3);

    BeamMathUtilsType::VectorToSkewSymmetricTensor( rRotationVector, rRotationTensor );

    rRotationTensor *= (-0.5);

    Matrix RotationxRotation = outer_prod( rRotationVector, rRotationVector );

    RotationxRotation *= 0.25;

    Matrix DiagonalMatrix(3,3);
    noalias(DiagonalMatrix) = IdentityMatrix(3);

    rRotationTensor += RotationxRotation + DiagonalMatrix;

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
  int LargeDisplacementBeamEMCElement::Check(const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    // Perform base element checks
    int ErrorCode = 0;
    ErrorCode = LargeDisplacementBeamElement::Check(rCurrentProcessInfo);

    return ErrorCode;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************


  void LargeDisplacementBeamEMCElement::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementBeamElement )
    rSerializer.save("CurrentStrainResultantsVector",mCurrentStrainResultantsVector);
    rSerializer.save("PreviousStrainResultantsVector",mPreviousStrainResultantsVector);
  }

  void LargeDisplacementBeamEMCElement::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementBeamElement )
    rSerializer.load("CurrentStrainResultantsVector",mCurrentStrainResultantsVector);
    rSerializer.load("PreviousStrainResultantsVector",mPreviousStrainResultantsVector);

  }


} // Namespace Kratos


