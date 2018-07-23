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
#include "custom_elements/beam_elements/geometrically_exact_rod_element.hpp"

#include "solid_mechanics_application_variables.h"


namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  GeometricallyExactRodElement::GeometricallyExactRodElement(IndexType NewId,GeometryType::Pointer pGeometry)
    : LargeDisplacementBeamEMCElement(NewId, pGeometry)
  {

  }

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  GeometricallyExactRodElement::GeometricallyExactRodElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : LargeDisplacementBeamEMCElement(NewId, pGeometry, pProperties)
  {
    KRATOS_TRY

    KRATOS_CATCH( "" )
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  GeometricallyExactRodElement::GeometricallyExactRodElement(GeometricallyExactRodElement const& rOther)
    :LargeDisplacementBeamEMCElement(rOther)
    ,mInitialLocalDirectors(rOther.mInitialLocalDirectors)
    ,mCurrentLocalDirectors(rOther.mCurrentLocalDirectors)
    ,mPreviousLocalDirectors(rOther.mPreviousLocalDirectors)
    ,mInitialLocalDirectorsVelocities(rOther.mInitialLocalDirectorsVelocities)
    ,mCurrentLocalDirectorsVelocities(rOther.mCurrentLocalDirectorsVelocities)
    ,mPreviousLocalDirectorsVelocities(rOther.mPreviousLocalDirectorsVelocities)
  {
  }

  //*********************************OPERATIONS*****************************************
  //************************************************************************************

  Element::Pointer GeometricallyExactRodElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
  {
    return Kratos::make_shared<GeometricallyExactRodElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  GeometricallyExactRodElement::~GeometricallyExactRodElement()
  {
  }

  //************************************************************************************
  //************************************************************************************

  void GeometricallyExactRodElement::Initialize()
  {
    KRATOS_TRY

    LargeDisplacementBeamEMCElement::Initialize();

    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    const SizeType number_of_nodes  = GetGeometry().size();

    //Initial Local Directors initialization
    if ( mInitialLocalDirectors.size() != number_of_nodes )
      {
	mInitialLocalDirectors.resize(number_of_nodes);
      }

    //Current Local Directors initialization
    if ( mCurrentLocalDirectors.size() != number_of_nodes )
      {
	mCurrentLocalDirectors.resize(number_of_nodes);
      }

    //Previous Local Directors initialization
    if ( mPreviousLocalDirectors.size() != number_of_nodes )
      {
	mPreviousLocalDirectors.resize(number_of_nodes);
      }


    Matrix Identity = IdentityMatrix(dimension);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	mInitialLocalDirectors[i]   = Identity;
	mCurrentLocalDirectors[i]   = Identity;
	mPreviousLocalDirectors[i]  = Identity;
      }

    //*************//

    //Initial Local Directors Derivatives initialization
    if ( mInitialLocalDirectorsVelocities.size() != number_of_nodes )
      {
	mInitialLocalDirectorsVelocities.resize(number_of_nodes);
      }

    //Current Local Directors Derivatives initialization
    if ( mCurrentLocalDirectorsVelocities.size() != number_of_nodes )
      {
	mCurrentLocalDirectorsVelocities.resize(number_of_nodes);
      }

    //Previous Local Directors Derivatives initialization
    if ( mPreviousLocalDirectorsVelocities.size() != number_of_nodes )
      {
	mPreviousLocalDirectorsVelocities.resize(number_of_nodes);
      }

    Matrix MatrixZero(dimension,dimension);
    noalias(MatrixZero) = ZeroMatrix(dimension,dimension);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	mInitialLocalDirectorsVelocities[i]  = MatrixZero;
	mCurrentLocalDirectorsVelocities[i]  = MatrixZero;
	mPreviousLocalDirectorsVelocities[i] = MatrixZero;
      }

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void GeometricallyExactRodElement::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    LargeDisplacementBeamEMCElement::InitializeSolutionStep(rCurrentProcessInfo);

    for ( unsigned int i = 0; i < mCurrentLocalDirectors.size(); i++ )
      {
	mPreviousLocalDirectors[i] = mCurrentLocalDirectors[i] ;
      }

    for ( unsigned int i = 0; i < mCurrentLocalDirectorsVelocities.size(); i++ )
      {
	mPreviousLocalDirectorsVelocities[i] = mCurrentLocalDirectorsVelocities[i] ;
      }

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void GeometricallyExactRodElement::InitializeElementData(ElementDataType& rVariables, const ProcessInfo& rCurrentProcessInfo)
  {

    KRATOS_TRY

    LargeDisplacementBeamEMCElement::InitializeElementData(rVariables,rCurrentProcessInfo);

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();


    DirectorsVariables& Directors = rVariables.GetDirectors();

    Directors.Initial.resize(dimension);
    Directors.InitialDerivatives.resize(dimension);

    Directors.Current.resize(dimension);
    Directors.CurrentDerivatives.resize(dimension);

    Directors.Previous.resize(dimension);
    Directors.PreviousDerivatives.resize(dimension);

    Vector VectorZero(3);
    noalias(VectorZero) = ZeroVector(3);

    for ( unsigned int j = 0; j < 3; j++ )
      {
	Directors.Initial[j]               = VectorZero;
	Directors.InitialDerivatives[j]    = VectorZero;

	Directors.Current[j]               = VectorZero;
	Directors.CurrentDerivatives[j]    = VectorZero;

	Directors.Previous[j]              = VectorZero;
	Directors.PreviousDerivatives[j]   = VectorZero;
      }

    Directors.CurrentNode.resize(number_of_nodes);
    Directors.PreviousNode.resize(number_of_nodes);

    Directors.CurrentNodeVelocities.resize(number_of_nodes);
    Directors.PreviousNodeVelocities.resize(number_of_nodes);

    Matrix MatrixZero(dimension,dimension);
    noalias(MatrixZero) = ZeroMatrix(dimension,dimension);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {

	Directors.CurrentNode[i]            = MatrixZero;
	Directors.PreviousNode[i]           = MatrixZero;

	Directors.CurrentNodeVelocities[i]  = MatrixZero;
	Directors.PreviousNodeVelocities[i] = MatrixZero;

      }

    KRATOS_CATCH( "" )
  }


  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************

  void GeometricallyExactRodElement::CalculateKinematics(ElementDataType& rVariables, const unsigned int& rPointNumber)
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

    //TOTAL LAGRANGIAN (spatial coordinates)
    //Compute cartesian derivatives [dN/dx_0]
    if( mThisIntegrationMethod == mReducedIntegrationMethod ){
      rVariables.DN_DX = mInvJ0 * DN_De[rPointNumber];
      rVariables.detJ  = 1.0/mInvJ0;
    }

    if( mThisIntegrationMethod == mFullIntegrationMethod ){
      rVariables.DN_DX = mInvJ0 * DN_De[rPointNumber];
      rVariables.detJ  = 1.0/mInvJ0;
    }


    //********************************************************//

    //Compute centroid displacement derivatives:
    noalias(rVariables.InitialAxisPositionDerivatives)  = ZeroVector(dimension);
    noalias(rVariables.CurrentAxisPositionDerivatives)  = ZeroVector(dimension);
    noalias(rVariables.PreviousAxisPositionDerivatives) = ZeroVector(dimension);

    QuaternionType QuaternionValue;
    Matrix CurrentValueMatrix(3,3);
    noalias(CurrentValueMatrix) = ZeroMatrix(3,3);
    Matrix CurrentValueDirectors(3,3);
    noalias(CurrentValueDirectors) = ZeroMatrix(3,3);

    Vector InitialValueVector(3);
    noalias(InitialValueVector) = ZeroVector(3);
    Vector PreviousValueVector(3);
    noalias(PreviousValueVector) = ZeroVector(3);
    Vector CurrentValueVector(3);
    noalias(CurrentValueVector) = ZeroVector(3);

    DirectorsVariables& Directors = rVariables.GetDirectors();

    //strains due to displacements and rotations
    for ( SizeType i = 0; i < number_of_nodes; i++ )
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

	//C: Initial Nodes Position
     	CurrentValueVector = GetGeometry()[i].GetInitialPosition();
	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rPointNumber );

	//Initial Frame Axis Position derivative
	rVariables.InitialAxisPositionDerivatives +=  rVariables.DN_DX(i,0) * ( CurrentValueVector );


	//********************************************************//
	//UPDATE CURRENT NODE DIRECTORS:

	Directors.PreviousNode[i] = mPreviousLocalDirectors[i];

	//current directors (total)
	// CurrentValueVector = GetNodalCurrentValue( ROTATION, CurrentValueVector, i );
	// CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rPointNumber );

	// QuaternionValue = QuaternionType::FromRotationVector(CurrentValueVector);
	// QuaternionValue.ToRotationMatrix(CurrentValueMatrix);

	// noalias(CurrentValueDirectors) = prod(CurrentValueMatrix, mInitialLocalDirectors[i]);

	// Directors.CurrentNode[i] = CurrentValueDirectors;


	//current directors (incremental)
	CurrentValueVector = GetNodalCurrentValue( STEP_ROTATION, CurrentValueVector, i );
	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rPointNumber );

	QuaternionValue = QuaternionType::FromRotationVector(CurrentValueVector);
	QuaternionValue.ToRotationMatrix(CurrentValueMatrix);

	noalias(CurrentValueDirectors) = prod(CurrentValueMatrix, mPreviousLocalDirectors[i]);

	Directors.CurrentNode[i] = CurrentValueDirectors;

	//********************************************************//

	//Current Directors and Derivatives
	for ( unsigned int j = 0; j < 3; j++ )
	  {
	    for ( unsigned int k = 0; k < 3; k++ )
	      {
		InitialValueVector[k]  = mInitialLocalDirectors[i](k,j);
		PreviousValueVector[k] = mPreviousLocalDirectors[i](k,j);
	        CurrentValueVector[k]  = Directors.CurrentNode[i](k,j);
	      }

	    Directors.Initial[j]  += rVariables.N[i] * InitialValueVector;
	    Directors.Previous[j] += rVariables.N[i] * PreviousValueVector;
	    Directors.Current[j]  += rVariables.N[i] * CurrentValueVector;

	    Directors.InitialDerivatives[j]   +=  rVariables.DN_DX(i,0) * InitialValueVector;
	    Directors.PreviousDerivatives[j]  +=  rVariables.DN_DX(i,0) * PreviousValueVector;
	    Directors.CurrentDerivatives[j]   +=  rVariables.DN_DX(i,0) * CurrentValueVector;
	  }

      }


    //********************************************************//

    Matrix TensorAngularVelocity(3,3);
    noalias(TensorAngularVelocity) = ZeroMatrix(3,3);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {

	//UPDATE CURRENT NODE DIRECTORS VELOCITIES:

	//previous local directors velocities
	Directors.PreviousNodeVelocities[i] = mPreviousLocalDirectorsVelocities[i];

	//current local directors velocities
	CurrentValueVector  = GetNodalCurrentValue( ANGULAR_VELOCITY, CurrentValueVector, i );
	CurrentValueVector  = MapToInitialLocalFrame( CurrentValueVector, rPointNumber );

	noalias(TensorAngularVelocity) = ZeroMatrix(3,3);
	BeamMathUtilsType::VectorToSkewSymmetricTensor( CurrentValueVector, TensorAngularVelocity );

	noalias(CurrentValueMatrix) = ZeroMatrix(3,3);
	this->CalculateAlphaDirectors( CurrentValueMatrix, rVariables, i, rVariables.Alpha);

	Directors.CurrentNodeVelocities[i] = (-1) * mPreviousLocalDirectorsVelocities[i] +  2 * prod( TensorAngularVelocity, CurrentValueMatrix );
      }


    //compute local to global frame
    this->CalculateFrameMapping( rVariables, rPointNumber );


    //*************************************//

    //set current STRAIN RESULTANTS
    rVariables.CurrentStrainResultantsVector  = mPreviousStrainResultantsVector[rPointNumber];
    rVariables.PreviousStrainResultantsVector = mPreviousStrainResultantsVector[rPointNumber];

    //set current CURVATURES
    rVariables.CurrentCurvatureVector  = mPreviousCurvatureVectors[rPointNumber];
    rVariables.PreviousCurvatureVector = mPreviousCurvatureVectors[rPointNumber];


    KRATOS_CATCH( "" )
  }

  //*************************COMPUTE FRAME MAPPING*************************************
  //************************************************************************************

  void GeometricallyExactRodElement::CalculateFrameMapping(ElementDataType& rVariables,const unsigned int& rPointNumber)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();

    //*------------------------------*//
    DirectorsVariables& Directors = rVariables.GetDirectors();

    //directors previous
    rVariables.PreviousRotationMatrix.resize(3,3,false);
    noalias(rVariables.PreviousRotationMatrix) = ZeroMatrix(3,3);
    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
    	rVariables.PreviousRotationMatrix  += rVariables.N[i] * ( Directors.PreviousNode[i]);
      }


    //directors current
    rVariables.CurrentRotationMatrix.resize(3,3,false);
    noalias(rVariables.CurrentRotationMatrix) = ZeroMatrix(3,3);
    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
    	rVariables.CurrentRotationMatrix  += rVariables.N[i] * ( Directors.CurrentNode[i]);
      }


    //*------------------------------*//

    CalculateAlphaRotationMatrix( rVariables.PreviousRotationMatrix, rVariables.CurrentRotationMatrix, rVariables.AlphaRotationMatrix, rVariables.AlphaRotationMatrixAsterisk, rVariables.Alpha);

    //*------------------------------*//

    //set variables for the initialization update
    UpdateRotationVariables(rVariables,rPointNumber);

    KRATOS_CATCH( "" )
  }

  //*********************************SET ROTATION VARIABLES*****************************
  //************************************************************************************

  void GeometricallyExactRodElement::UpdateRotationVariables(ElementDataType& rVariables, const unsigned int& rPointNumber)
  {
    KRATOS_TRY

    DirectorsVariables& Directors = rVariables.GetDirectors();

    for ( unsigned int i = 0; i < mCurrentLocalDirectors.size(); i++ )
      {
	mCurrentLocalDirectors[i] = Directors.CurrentNode[i];
      }

    for ( unsigned int i = 0; i < mCurrentLocalDirectorsVelocities.size(); i++ )
      {
	mCurrentLocalDirectorsVelocities[i] = Directors.CurrentNodeVelocities[i];
      }

    KRATOS_CATCH( "" )
  }



  //*************************COMPUTE MAPPING TENSOR*************************************
  //************************************************************************************

  void GeometricallyExactRodElement::CalculateDirectorsMappingTensor(Matrix& rMappingTensor, ElementDataType& rVariables, const int& rNode, double alpha)
  {

    Matrix SkewSymDirector(3,3);
    noalias(SkewSymDirector) = ZeroMatrix(3,3);
    std::vector<Vector> DirectorVectorsAlpha;
    DirectorVectorsAlpha.resize(3);

    for( unsigned int i = 0; i < 3; i++ )
      {
	DirectorVectorsAlpha[i].resize(3,false);
	noalias(DirectorVectorsAlpha[i]) = ZeroVector(3);
      }

    Matrix DiagonalMatrix = IdentityMatrix(3);

    //nodal directors
    DirectorsVariables& Directors = rVariables.GetDirectors();
    Matrix& CurrentDirectors      = Directors.CurrentNode[rNode];
    Matrix& PreviousDirectors     = Directors.PreviousNode[rNode];

    for( unsigned int i = 0; i < 3; i++ )
      {
    	DirectorVectorsAlpha[0][i] = (1-alpha) * PreviousDirectors(i,0) + alpha * CurrentDirectors(i,0);
    	DirectorVectorsAlpha[1][i] = (1-alpha) * PreviousDirectors(i,1) + alpha * CurrentDirectors(i,1);
    	DirectorVectorsAlpha[2][i] = (1-alpha) * PreviousDirectors(i,2) + alpha * CurrentDirectors(i,2);
      }

    Matrix CurrentRotationMatrixAlpha =  (1-alpha) * PreviousDirectors + alpha * CurrentDirectors;

    //elemental directors
    // for( unsigned int i = 0; i < 3; i++ )
    //   {
    // 	DirectorVectorsAlpha[0][i] = (1-alpha) * Directors.Previous[0][i] + alpha * Directors.Current[0][i];
    // 	DirectorVectorsAlpha[1][i] = (1-alpha) * Directors.Previous[1][i] + alpha * Directors.Current[1][i];
    // 	DirectorVectorsAlpha[2][i] = (1-alpha) * Directors.Previous[2][i] + alpha * Directors.Current[2][i];
    //   }

    // std::cout<<" Director1 "<<DirectorVectorsAlpha[0]<<std::endl;
    // std::cout<<" Director2 "<<DirectorVectorsAlpha[1]<<std::endl;
    // std::cout<<" Director3 "<<DirectorVectorsAlpha[2]<<std::endl;


    //Mapping Matrix
    if( rMappingTensor.size1() !=12 || rMappingTensor.size2() != 6 )
      rMappingTensor.resize(12,6);

    noalias(rMappingTensor) = ZeroMatrix(12,6);

    BeamMathUtilsType::AddMatrix( rMappingTensor, DiagonalMatrix, 0, 0 );

    BeamMathUtilsType::VectorToSkewSymmetricTensor(DirectorVectorsAlpha[0], SkewSymDirector);
    SkewSymDirector *= (-1);
    BeamMathUtilsType::AddMatrix( rMappingTensor, SkewSymDirector, 3, 3 );

    BeamMathUtilsType::VectorToSkewSymmetricTensor(DirectorVectorsAlpha[1], SkewSymDirector);
    SkewSymDirector *= (-1);
    BeamMathUtilsType::AddMatrix( rMappingTensor, SkewSymDirector, 6, 3 );

    BeamMathUtilsType::VectorToSkewSymmetricTensor(DirectorVectorsAlpha[2], SkewSymDirector);
    SkewSymDirector *= (-1);
    BeamMathUtilsType::AddMatrix( rMappingTensor, SkewSymDirector, 9, 3 );

    //std::cout<<" Mapping Tensor "<<rMappingTensor<<std::endl;

  }


  //************************************************************************************
  //************************************************************************************


  void GeometricallyExactRodElement::CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
							       ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    //create and initialize element variables:
    ElementDataType Variables;
    this->InitializeElementData(Variables,rCurrentProcessInfo);

    DirectorsVariables Directors;
    Variables.SetDirectors(Directors);

    this->InitializeElementData(Variables,rCurrentProcessInfo);

    //reading integration points (in fact is the two nodes beam element, only one integration point)
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //auxiliary terms
    Vector VolumeForce(3);
    noalias(VolumeForce) = ZeroVector(3);

    double IntegrationWeight = 1.0;

    //(in fact is the two nodes beam element, only one integration point)
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {

        //compute element kinematics  ...
        this->CalculateKinematics(Variables,PointNumber);

	// std::cout<<" ID "<<this->Id()<<std::endl;
	// std::cout<<" Delta Position "<<Variables.DeltaPosition<<std::endl;
	// std::cout<<" CurrentRotationMatrix "<<Variables.CurrentRotationMatrix<<std::endl;

	//compute element ConstitutiveTensor
	this->CalculateConstitutiveMatrix(Variables);

	//compute element Strain and Stress Resultants and Couples
	this->CalculateStressResultants(Variables, PointNumber);

	// std::cout<<" StrainResultants "<<Variables.StrainVector<<std::endl;
	// std::cout<<" StressResultants "<<Variables.StressVector<<std::endl;

	IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
	IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );


	if ( rLocalSystem.CalculationFlags.Is(BeamElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
	  {
	    this->CalculateAndAddLHS( rLocalSystem, Variables, IntegrationWeight );
	  }

	if ( rLocalSystem.CalculationFlags.Is(BeamElement::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
	  {
	    //contribution to external forces
	    VolumeForce  = this->CalculateVolumeForce( VolumeForce, Variables.N );

	    this->CalculateAndAddRHS( rLocalSystem , Variables, VolumeForce, IntegrationWeight );
	  }

      }

    KRATOS_CATCH( "" )
  }



  //************************************************************************************
  //************************************************************************************


  void GeometricallyExactRodElement::CalculateDynamicSystem( LocalSystemComponents& rLocalSystem,
							     ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    //create and initialize element variables:
    ElementDataType Variables;
    this->InitializeElementData(Variables,rCurrentProcessInfo);

    DirectorsVariables Directors;
    Variables.SetDirectors(Directors);


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


    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();


    //reading shape functions
    Variables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod ));

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


	Directors.CurrentNode.resize(number_of_nodes);
	Directors.PreviousNode.resize(number_of_nodes);

	Directors.CurrentNodeVelocities.resize(number_of_nodes);
	Directors.PreviousNodeVelocities.resize(number_of_nodes);

	Matrix MatrixZero(dimension,dimension);
	noalias(MatrixZero) = ZeroMatrix(dimension,dimension);

	for ( SizeType i = 0; i < number_of_nodes; i++ )
	  {
	    Directors.CurrentNode[i].resize(dimension,dimension,false);
	    Directors.PreviousNode[i].resize(dimension,dimension,false);

	    Directors.CurrentNodeVelocities[i].resize(dimension,dimension,false);
	    Directors.PreviousNodeVelocities[i].resize(dimension,dimension,false);

	    noalias(Directors.CurrentNode[i])            = MatrixZero;
	    noalias(Directors.PreviousNode[i])           = MatrixZero;

	    noalias(Directors.CurrentNodeVelocities[i])  = MatrixZero;
	    noalias(Directors.PreviousNodeVelocities[i]) = MatrixZero;

	  }

	QuaternionType QuaternionValue;
	Matrix CurrentValueMatrix(3,3);
	noalias(CurrentValueMatrix) = ZeroMatrix(3,3);
	Matrix CurrentValueDirectors(3,3);
	noalias(CurrentValueDirectors) = ZeroMatrix(3,3);

	Vector CurrentValueVector(3);
	noalias(CurrentValueVector) = ZeroVector(3);
	Matrix TensorAngularVelocity(3,3);
	noalias(TensorAngularVelocity) = ZeroMatrix(3,3);

	//strains due to displacements and rotations
	for ( SizeType i = 0; i < number_of_nodes; i++ )
	  {

	    //********************************************************//
	    //UPDATE CURRENT NODE DIRECTORS:

	    Directors.PreviousNode[i] = mPreviousLocalDirectors[i];


	    //current directors (total)
	    // CurrentValueVector = GetNodalCurrentValue( ROTATION, CurrentValueVector, i );
	    // CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rPointNumber );

	    // QuaternionValue = QuaternionType::FromRotationVector(CurrentValueVector);
	    // QuaternionValue.ToRotationMatrix(CurrentValueMatrix);

	    // noalias(CurrentValueDirectors) = prod(CurrentValueMatrix, mInitialLocalDirectors[i]);

	    // Directors.CurrentNode[i] = CurrentValueDirectors;


	    //current directors (incremental)
	    CurrentValueVector = GetNodalCurrentValue( STEP_ROTATION, CurrentValueVector, i );
	    CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, PointNumber );

	    QuaternionValue = QuaternionType::FromRotationVector(CurrentValueVector);
	    QuaternionValue.ToRotationMatrix(CurrentValueMatrix);

	    noalias(CurrentValueDirectors) = prod(CurrentValueMatrix, mPreviousLocalDirectors[i]);

	    Directors.CurrentNode[i] = CurrentValueDirectors;


	    //********************************************************//
	    //UPDATE CURRENT NODE DIRECTORS VELOCITIES:

	    //previous local directors velocities
	    Directors.PreviousNodeVelocities[i] = mPreviousLocalDirectorsVelocities[i];

	    //current local directors velocities
	    CurrentValueVector  = GetNodalCurrentValue( ANGULAR_VELOCITY, CurrentValueVector, i );
	    CurrentValueVector  = MapToInitialLocalFrame( CurrentValueVector, PointNumber );

	    noalias(TensorAngularVelocity) = ZeroMatrix(3,3);
	    BeamMathUtilsType::VectorToSkewSymmetricTensor( CurrentValueVector, TensorAngularVelocity );

	    noalias(CurrentValueMatrix) = ZeroMatrix(3,3);
	    this->CalculateAlphaDirectors( CurrentValueMatrix, Variables, i, Variables.Alpha);

	    Directors.CurrentNodeVelocities[i] = (-1) * mPreviousLocalDirectorsVelocities[i] +  2 * prod( TensorAngularVelocity, CurrentValueMatrix );

	  }

	//compute local to global frame
	this->CalculateFrameMapping( Variables, PointNumber );


	//TOTAL LAGRANGIAN (spatial coordinates)
	Variables.detJ = 1.0/mInvJ0;

	double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
	IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );

	if ( rLocalSystem.CalculationFlags.Is(GeometricallyExactRodElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
	  {
	    LocalLeftHandSideMatrix.clear();

    	    this->CalculateAndAddInertiaLHS( LocalLeftHandSideMatrix, Variables, rCurrentProcessInfo, IntegrationWeight ); // (R_N+1, R_N)

	    BeamMathUtilsType::MapLocalToGlobal3D(mInitialLocalQuaternion,LocalLeftHandSideMatrix);

	    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();
	    rLeftHandSideMatrix += LocalLeftHandSideMatrix;
	  }

	if ( rLocalSystem.CalculationFlags.Is(GeometricallyExactRodElement::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
	  {
	    LocalRightHandSideVector.clear();

	    this->CalculateAndAddInertiaRHS( LocalRightHandSideVector, Variables, rCurrentProcessInfo, IntegrationWeight );

	    BeamMathUtilsType::MapLocalToGlobal3D(mInitialLocalQuaternion,LocalRightHandSideVector);

	    VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
	    rRightHandSideVector += LocalRightHandSideVector;
	  }

      }

    mThisIntegrationMethod = ThisIntegrationMethod;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void GeometricallyExactRodElement::CalculateCurrentStrainResultantsVector(ElementDataType& rVariables,
									    Vector& rCurrentStrainResultantsVector,
									    double Alpha)
  {
    KRATOS_TRY

      // current strain resultants
      if( rCurrentStrainResultantsVector.size() != 3 )
	rCurrentStrainResultantsVector.resize(3, false);

    noalias(rCurrentStrainResultantsVector) = ZeroVector(3);

    DirectorsVariables& Directors          = rVariables.GetDirectors();
    std::vector<Vector>& InitialDirectors  = Directors.Initial;
    std::vector<Vector>& CurrentDirectors  = Directors.Current;

    rCurrentStrainResultantsVector[0]  = inner_prod( CurrentDirectors[0], rVariables.CurrentAxisPositionDerivatives );
    rCurrentStrainResultantsVector[1]  = inner_prod( CurrentDirectors[1], rVariables.CurrentAxisPositionDerivatives );
    rCurrentStrainResultantsVector[2]  = inner_prod( CurrentDirectors[2], rVariables.CurrentAxisPositionDerivatives );


    //equivalent:
    //rCurrentStrainResultantsVector = prod( trans(rVariables.CurrentRotationMatrix), rVariables.CurrentAxisPositionDerivatives );

    //----------------------

    // Vector E1(3);
    // noalias(E1) = ZeroVector(3);
    // E1[0] = 1.0;
    // rCurrentStrainResultantsVector -= E1;  //next is equivalent

    rCurrentStrainResultantsVector[0] -= inner_prod( InitialDirectors[0], rVariables.InitialAxisPositionDerivatives );
    rCurrentStrainResultantsVector[1] -= inner_prod( InitialDirectors[1], rVariables.InitialAxisPositionDerivatives );
    rCurrentStrainResultantsVector[2] -= inner_prod( InitialDirectors[2], rVariables.InitialAxisPositionDerivatives );


    // Vector InitialStrainResultantsVector = ZeroVector(3);

    // InitialStrainResultantsVector[0] -= inner_prod( InitialDirectors[0], rVariables.InitialAxisPositionDerivatives );
    // InitialStrainResultantsVector[1] -= inner_prod( InitialDirectors[1], rVariables.InitialAxisPositionDerivatives );
    // InitialStrainResultantsVector[2] -= inner_prod( InitialDirectors[2], rVariables.InitialAxisPositionDerivatives );

    // std::cout<<" InitialStrainResultantsVector "<<InitialStrainResultantsVector<<std::endl;


    //std::cout<<" CurrentStrainResultants "<<rCurrentStrainResultantsVector<<std::endl;

    KRATOS_CATCH( "" )

  }

  //************************************************************************************
  //************************************************************************************

  void GeometricallyExactRodElement::CalculateCurrentCurvatureVector(ElementDataType& rVariables,
								     Vector& rCurrentCurvatureVector,
								     double Alpha)
  {
    KRATOS_TRY

    // current strain resultants
    if( rCurrentCurvatureVector.size() != 3 )
      rCurrentCurvatureVector.resize(3, false);

    noalias(rCurrentCurvatureVector) = ZeroVector(3);

    DirectorsVariables& Directors          = rVariables.GetDirectors();
    std::vector<Vector>& InitialDirectors  = Directors.Initial;
    std::vector<Vector>& CurrentDirectors  = Directors.Current;
    std::vector<Vector>& InitialDirectorsDerivatives  = Directors.InitialDerivatives;
    std::vector<Vector>& CurrentDirectorsDerivatives  = Directors.CurrentDerivatives;


    rCurrentCurvatureVector[0]  = 0.5 * inner_prod( CurrentDirectors[2], CurrentDirectorsDerivatives[1] );
    rCurrentCurvatureVector[1]  = 0.5 * inner_prod( CurrentDirectors[0], CurrentDirectorsDerivatives[2] );
    rCurrentCurvatureVector[2]  = 0.5 * inner_prod( CurrentDirectors[1], CurrentDirectorsDerivatives[0] );

    rCurrentCurvatureVector[0] -= 0.5 * inner_prod( CurrentDirectors[1], CurrentDirectorsDerivatives[2] );
    rCurrentCurvatureVector[1] -= 0.5 * inner_prod( CurrentDirectors[2], CurrentDirectorsDerivatives[0] );
    rCurrentCurvatureVector[2] -= 0.5 * inner_prod( CurrentDirectors[0], CurrentDirectorsDerivatives[1] );

    //----------------------

    //next is equivalent to zero

    rCurrentCurvatureVector[0] -= 0.5 * inner_prod( InitialDirectors[2], InitialDirectorsDerivatives[1] );
    rCurrentCurvatureVector[1] -= 0.5 * inner_prod( InitialDirectors[0], InitialDirectorsDerivatives[2] );
    rCurrentCurvatureVector[2] -= 0.5 * inner_prod( InitialDirectors[1], InitialDirectorsDerivatives[0] );

    rCurrentCurvatureVector[0] += 0.5 * inner_prod( InitialDirectors[1], InitialDirectorsDerivatives[2] );
    rCurrentCurvatureVector[1] += 0.5 * inner_prod( InitialDirectors[2], InitialDirectorsDerivatives[0] );
    rCurrentCurvatureVector[2] += 0.5 * inner_prod( InitialDirectors[0], InitialDirectorsDerivatives[1] );

    // Vector InitialCurvatureVector = ZeroVector(3);
    // InitialCurvatureVector[0] -= 0.5 * inner_prod( InitialDirectors[2], InitialDirectorsDerivatives[1] );
    // InitialCurvatureVector[1] -= 0.5 * inner_prod( InitialDirectors[0], InitialDirectorsDerivatives[2] );
    // InitialCurvatureVector[2] -= 0.5 * inner_prod( InitialDirectors[1], InitialDirectorsDerivatives[0] );

    // InitialCurvatureVector[0] += 0.5 * inner_prod( InitialDirectors[1], InitialDirectorsDerivatives[2] );
    // InitialCurvatureVector[1] += 0.5 * inner_prod( InitialDirectors[2], InitialDirectorsDerivatives[0] );
    // InitialCurvatureVector[2] += 0.5 * inner_prod( InitialDirectors[0], InitialDirectorsDerivatives[1] );

    // std::cout<<" InitialCurvatureVector "<<InitialCurvatureVector<<std::endl;

    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************

  void GeometricallyExactRodElement::CalculateConstitutiveMatrix(ElementDataType& rVariables)
  {
    KRATOS_TRY

    // Material Elastic constitutive matrix
    this->CalculateMaterialConstitutiveMatrix(rVariables.ConstitutiveMatrix, rVariables);

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void GeometricallyExactRodElement::CalculateStressResultants(ElementDataType& rVariables, const unsigned int& rPointNumber)
  {
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    //compute Strain Resultants and Couples
    Vector StrainResultants(3);
    noalias(StrainResultants) = ZeroVector(3);
    Vector StrainCouples(3);
    noalias(StrainCouples) = ZeroVector(3);

    this->CalculateStrainResultants(StrainResultants, rVariables, rVariables.Alpha);
    this->CalculateStrainCouples(StrainCouples, rVariables, rVariables.Alpha);

    //std::cout<<" CurrentAxisPositionDerivatives "<<rVariables.CurrentAxisPositionDerivatives<<std::endl;

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

    //******************************************//

    //Set current curvatures and strain resultans:
    CalculateCurrentStrainResultantsVector(rVariables, rVariables.CurrentStrainResultantsVector, rVariables.Alpha);
    CalculateCurrentCurvatureVector(rVariables, rVariables.CurrentCurvatureVector, rVariables.Alpha);

    //set variables for the initialization update
    this->UpdateStrainVariables(rVariables,rPointNumber);

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void GeometricallyExactRodElement::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
								   ElementDataType& rVariables,
								   Vector& rVolumeForce,
								   double& rIntegrationWeight)
  {
    KRATOS_TRY

    SizeType number_of_nodes  = GetGeometry().PointsNumber();
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    double DomainSize = rVariables.Section.Area;

    //gravity load
    Vector GravityLoad(dimension);
    noalias(GravityLoad) = ZeroVector(dimension);

    VectorType Fe(12);
    noalias(Fe) = ZeroVector(12);
    VectorType Fext(6);
    noalias(Fext) = ZeroVector(6);
    MatrixType MappingTensor;
    noalias(MappingTensor) = ZeroMatrix(12,6);

    unsigned int RowIndex = 0;
    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
      	RowIndex = i * (dimension * 2);
        GravityLoad =  rIntegrationWeight * rVariables.N[i] * rVolumeForce * DomainSize;

	BeamMathUtilsType::AddVector(GravityLoad, Fe, 0);

	this->CalculateDirectorsMappingTensor(MappingTensor, rVariables, i , rVariables.Alpha);

	Fext = prod( trans(MappingTensor), Fe );

	//Fext *= rVariables.Alpha;

	BeamMathUtilsType::AddVector(Fext, rRightHandSideVector, RowIndex);

	//std::cout<<" Fext "<<Fe<<std::endl;
	//std::cout<<" rRightHandSideVector "<<rRightHandSideVector<<std::endl;
      }

    //follower load forces (to implement)
    this->CalculateAndAddFollowerForces( rRightHandSideVector, rVariables, rIntegrationWeight );

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void GeometricallyExactRodElement::CalculateAndAddFollowerForces(VectorType& rRightHandSideVector,
								   ElementDataType & rVariables,
								   double& rIntegrationWeight)
  {
    KRATOS_TRY

    KRATOS_CATCH( "" )
  }



  //************************************************************************************
  //************************************************************************************

  void GeometricallyExactRodElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
								   ElementDataType & rVariables,
								   double& rIntegrationWeight)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    VectorType Fi(12);
    noalias(Fi) = ZeroVector(12);
    VectorType Fint(6);
    noalias(Fint) = ZeroVector(6);

    //Initialize Local Matrices
    MatrixType DifferentialOperatorI(12,6);
    noalias(DifferentialOperatorI) = ZeroMatrix(12,6);

    unsigned int RowIndex = 0;

    MatrixType MappingTensor(12,6);
    noalias(MappingTensor) = ZeroMatrix(12,6);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {

	RowIndex = i * (dimension * 2);

	noalias(Fi)   = ZeroVector(12);
	noalias(Fint) = ZeroVector(6);

 	this->CalculateDifferentialOperator(DifferentialOperatorI, rVariables, i , rVariables.Alpha);

	//std::cout<<" Differential B "<<DifferentialOperatorI<<std::endl;

	//nodal force vector
	Fi = prod( DifferentialOperatorI, rVariables.StressVector );

	Fi *= rIntegrationWeight;

	this->CalculateDirectorsMappingTensor(MappingTensor, rVariables, i , rVariables.Alpha);

	//std::cout<<" Btensor "<<prod(trans(MappingTensor), DifferentialOperatorI)<<std::endl;

	Fint = prod( trans(MappingTensor), Fi );

	//std::cout<<" Fint "<<Fint<<" Fi "<<Fi<<std::endl;
	BeamMathUtilsType::SubstractVector(Fint, rRightHandSideVector, RowIndex);

	//std::cout<<" Fi "<<Fint<<std::endl;

	//std::cout<<" rRightHandSideVector "<<rRightHandSideVector<<std::endl;

      }


    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************

  void GeometricallyExactRodElement::CalculateDifferentialOperator(MatrixType& rDifferentialOperator,
								   ElementDataType& rVariables,
								   const int& rNode,
								   double alpha)
  {
    KRATOS_TRY

    //Differencial operator transposed

    //Initialize Local Matrices
    if( rDifferentialOperator.size1() != 12 )
      rDifferentialOperator.resize(12, 6, false);

    noalias(rDifferentialOperator) = ZeroMatrix(12,6);

    Matrix OperatorTau(6,12);
    noalias(OperatorTau) = ZeroMatrix(6,12);
    Matrix OperatorOmegaI(6,12);
    noalias(OperatorOmegaI) = ZeroMatrix(6,12);
    Matrix OperatorOmegaII(6,12);
    noalias(OperatorOmegaII) = ZeroMatrix(6,12);

    std::vector<Vector> PreviousDirectors;
    PreviousDirectors.resize(3);
    std::vector<Vector> CurrentDirectors;
    CurrentDirectors.resize(3);

    DirectorsVariables& Directors  = rVariables.GetDirectors();
    Matrix& CurrentNodeDirectors   = Directors.CurrentNode[rNode];
    Matrix& PreviousNodeDirectors  = Directors.PreviousNode[rNode];

    //Current Directors and Derivatives
    for ( unsigned int j = 0; j < 3; j++ )
      {
	PreviousDirectors[j].resize(3,false);
	noalias(PreviousDirectors[j]) = ZeroVector(3);
    	CurrentDirectors[j].resize(3,false);
	noalias(CurrentDirectors[j]) = ZeroVector(3);

    	for ( unsigned int k = 0; k < 3; k++ )
    	  {
    	    PreviousDirectors[j][k] =  PreviousNodeDirectors(k,j);
    	    CurrentDirectors[j][k]  =  CurrentNodeDirectors(k,j);
    	  }
      }

    Matrix AlphaNodeRotationMatrix(3,3);
    noalias(AlphaNodeRotationMatrix) = ZeroMatrix(3,3);
    AlphaNodeRotationMatrix = (1-alpha) * PreviousNodeDirectors + alpha * CurrentNodeDirectors;

    Matrix AuxiliarRotationMatrix(6,6);
    noalias(AuxiliarRotationMatrix) = ZeroMatrix(6,6);

    //Building the rotation matrix for the local element matrix
    for (unsigned int kk=0; kk < 6; kk += 3)
      {
        for (unsigned int i=0; i<3; i++)
	  {
            for(unsigned int j=0; j<3; j++)
	      {
		AuxiliarRotationMatrix(i+kk,j+kk) = AlphaNodeRotationMatrix(i,j);
	      }
	  }
      }

    std::vector<Vector> DirectorsAlpha(3);
    std::vector<Vector> DirectorsDerivativesAlpha(3);

    std::vector<Vector>& CurrentElementDirectors  = Directors.Current;
    std::vector<Vector>& PreviousElementDirectors = Directors.Previous;

    std::vector<Vector>& CurrentElementDirectorsDerivatives  = Directors.CurrentDerivatives;
    std::vector<Vector>& PreviousElementDirectorsDerivatives = Directors.PreviousDerivatives;

    for ( unsigned int j = 0; j < 3; j++ )
      {
	DirectorsAlpha[j] = (1-alpha) * PreviousElementDirectors[j] + alpha * CurrentElementDirectors[j]; //elemental directors
	//DirectorsAlpha[j] = (1-alpha) * PreviousDirectors[j] + alpha * CurrentDirectors[j]; //nodal directors
	DirectorsDerivativesAlpha[j] = (1-alpha) * PreviousElementDirectorsDerivatives[j] + alpha * CurrentElementDirectorsDerivatives[j];
      }

    Vector AxisPositionDerivativesAlpha(3);
    AxisPositionDerivativesAlpha = (1-alpha) * rVariables.PreviousAxisPositionDerivatives + alpha * rVariables.CurrentAxisPositionDerivatives;

    //BtauT
    OperatorTau( 0, 0 )  =  rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[0][0] );
    OperatorTau( 0, 1 )  =  rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[0][1] );
    OperatorTau( 0, 2 )  =  rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[0][2] );

    OperatorTau( 1, 0 )  =  rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[1][0] );
    OperatorTau( 1, 1 )  =  rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[1][1] );
    OperatorTau( 1, 2 )  =  rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[1][2] );

    OperatorTau( 2, 0 )  =  rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[2][0] );
    OperatorTau( 2, 1 )  =  rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[2][1] );
    OperatorTau( 2, 2 )  =  rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[2][2] );


    OperatorTau( 0, 3 )  =  rVariables.N[rNode] * ( AxisPositionDerivativesAlpha[0] );
    OperatorTau( 0, 4 )  =  rVariables.N[rNode] * ( AxisPositionDerivativesAlpha[1] );
    OperatorTau( 0, 5 )  =  rVariables.N[rNode] * ( AxisPositionDerivativesAlpha[2] );

    OperatorTau( 1, 6 )  =  rVariables.N[rNode] * ( AxisPositionDerivativesAlpha[0] );
    OperatorTau( 1, 7 )  =  rVariables.N[rNode] * ( AxisPositionDerivativesAlpha[1] );
    OperatorTau( 1, 8 )  =  rVariables.N[rNode] * ( AxisPositionDerivativesAlpha[2] );

    OperatorTau( 2, 9  ) =  rVariables.N[rNode] * ( AxisPositionDerivativesAlpha[0] );
    OperatorTau( 2, 10 ) =  rVariables.N[rNode] * ( AxisPositionDerivativesAlpha[1] );
    OperatorTau( 2, 11 ) =  rVariables.N[rNode] * ( AxisPositionDerivativesAlpha[2] );

    //BomegaT

    double sign = (-1);

    //part 1
    OperatorOmegaI( 3, 6 )  =  rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[2][0] );
    OperatorOmegaI( 3, 7 )  =  rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[2][1] );
    OperatorOmegaI( 3, 8 )  =  rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[2][2] );

    OperatorOmegaI( 3, 9 )  =  sign * rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[1][0] );
    OperatorOmegaI( 3, 10 ) =  sign * rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[1][1] );
    OperatorOmegaI( 3, 11 ) =  sign * rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[1][2] );

    OperatorOmegaI( 4, 3 )  =  sign * rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[2][0] );
    OperatorOmegaI( 4, 4 )  =  sign * rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[2][1] );
    OperatorOmegaI( 4, 5 )  =  sign * rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[2][2] );

    OperatorOmegaI( 4, 9 )  =  rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[0][0] );
    OperatorOmegaI( 4, 10 ) =  rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[0][1] );
    OperatorOmegaI( 4, 11 ) =  rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[0][2] );

    OperatorOmegaI( 5, 3 )  =  rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[1][0] );
    OperatorOmegaI( 5, 4 )  =  rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[1][1] );
    OperatorOmegaI( 5, 5 )  =  rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[1][2] );

    OperatorOmegaI( 5, 6 )  =  sign * rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[0][0] );
    OperatorOmegaI( 5, 7 )  =  sign * rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[0][1] );
    OperatorOmegaI( 5, 8 )  =  sign * rVariables.DN_DX( rNode, 0 ) * ( DirectorsAlpha[0][2] );


    //part 2
    OperatorOmegaII( 3, 6 )  =  sign * rVariables.N[rNode] * ( DirectorsDerivativesAlpha[2][0] );
    OperatorOmegaII( 3, 7 )  =  sign * rVariables.N[rNode] * ( DirectorsDerivativesAlpha[2][1] );
    OperatorOmegaII( 3, 8 )  =  sign * rVariables.N[rNode] * ( DirectorsDerivativesAlpha[2][2] );

    OperatorOmegaII( 3, 9 )  =  rVariables.N[rNode] * ( DirectorsDerivativesAlpha[1][0] );
    OperatorOmegaII( 3, 10 ) =  rVariables.N[rNode] * ( DirectorsDerivativesAlpha[1][1] );
    OperatorOmegaII( 3, 11 ) =  rVariables.N[rNode] * ( DirectorsDerivativesAlpha[1][2] );

    OperatorOmegaII( 4, 3 )  =  rVariables.N[rNode] * ( DirectorsDerivativesAlpha[2][0] );
    OperatorOmegaII( 4, 4 )  =  rVariables.N[rNode] * ( DirectorsDerivativesAlpha[2][1] );
    OperatorOmegaII( 4, 5 )  =  rVariables.N[rNode] * ( DirectorsDerivativesAlpha[2][2] );

    OperatorOmegaII( 4, 9 )  =  sign * rVariables.N[rNode] * ( DirectorsDerivativesAlpha[0][0] );
    OperatorOmegaII( 4, 10 ) =  sign * rVariables.N[rNode] * ( DirectorsDerivativesAlpha[0][1] );
    OperatorOmegaII( 4, 11 ) =  sign * rVariables.N[rNode] * ( DirectorsDerivativesAlpha[0][2] );

    OperatorOmegaII( 5, 3 )  =  sign * rVariables.N[rNode] * ( DirectorsDerivativesAlpha[1][0] );
    OperatorOmegaII( 5, 4 )  =  sign * rVariables.N[rNode] * ( DirectorsDerivativesAlpha[1][1] );
    OperatorOmegaII( 5, 5 )  =  sign * rVariables.N[rNode] * ( DirectorsDerivativesAlpha[1][2] );

    OperatorOmegaII( 5, 6 )  =  rVariables.N[rNode] * ( DirectorsDerivativesAlpha[0][0] );
    OperatorOmegaII( 5, 7 )  =  rVariables.N[rNode] * ( DirectorsDerivativesAlpha[0][1] );
    OperatorOmegaII( 5, 8 )  =  rVariables.N[rNode] * ( DirectorsDerivativesAlpha[0][2] );

    OperatorTau += 0.5 * (OperatorOmegaI + OperatorOmegaII);

    //OperatorTau = prod( trans(AuxiliarRotationMatrix), OperatorTau);

    rDifferentialOperator = trans(OperatorTau);

    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************

  void GeometricallyExactRodElement::CalculateDiscreteOperatorN(MatrixType& rDiscreteOperator,
								ElementDataType& rVariables,
								const int& rNodeI,
								const int& rNodeJ,
								const int& rComponent)
  {
    KRATOS_TRY

    //Initialize Local Matrices
    if( rDiscreteOperator.size1() != 3 || rDiscreteOperator.size2() != 3)
      rDiscreteOperator.resize(3, 3, false);

    noalias(rDiscreteOperator) = ZeroMatrix(3,3);

    Matrix DiagonalMatrix(3,3);
    noalias(DiagonalMatrix) = IdentityMatrix(3);

    Vector StressResultants(3);
    noalias(StressResultants) = ZeroVector(3);

    for ( unsigned int i = 0; i < 3; i++ )
      {
	StressResultants[i] = rVariables.StressVector[i];
      }

    //permutation tensor
    double kronecker = 1;


    Vector DirectorKNodeI(3);
    noalias(DirectorKNodeI) = ZeroVector(3);

    this->CalculateAlphaDirectorVector( DirectorKNodeI, rVariables, rNodeI, rComponent, 1 );

    Vector AxisPositionDerivatives(3);
    noalias(AxisPositionDerivatives) = ZeroVector(3);

    AxisPositionDerivatives = ( 1-rVariables.Alpha ) * rVariables.PreviousAxisPositionDerivatives  + rVariables.Alpha * rVariables.CurrentAxisPositionDerivatives;


    rDiscreteOperator += outer_prod( DirectorKNodeI, AxisPositionDerivatives );


    rDiscreteOperator -= inner_prod( DirectorKNodeI, AxisPositionDerivatives ) * DiagonalMatrix;

    kronecker = BeamMathUtilsType::KroneckerDelta(rNodeI, rNodeJ);

    rDiscreteOperator *= kronecker * rVariables.N[rNodeI] * StressResultants[rComponent];


    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************

  void GeometricallyExactRodElement::CalculateDiscreteOperatorM(MatrixType& rDiscreteOperator,
								ElementDataType& rVariables,
								const int& rNodeI,
								const int& rNodeJ,
								const int& rComponent)
  {
    KRATOS_TRY

    //Initialize Local Matrices
    if( rDiscreteOperator.size1() != 3 || rDiscreteOperator.size2() != 3)
      rDiscreteOperator.resize(3, 3, false);

    noalias(rDiscreteOperator) = ZeroMatrix(3,3);

    Vector StressCouples(3);
    noalias(StressCouples) = ZeroVector(3);
    for ( unsigned int i = 0; i < 3; i++ )
      {
	StressCouples[i] = rVariables.StressVector[i+3];
      }

    //permutation tensor
    double epsilon   = 1;
    double kronecker = 1;


    for ( unsigned int j = 0; j < 3; j++ )
      {

        for ( unsigned int k = 0; k < 3; k++ )
	  {

	    Matrix SkewSymDirectorKNodeI(3,3);
	    noalias(SkewSymDirectorKNodeI) = ZeroMatrix(3,3);
	    this->CalculateAlphaDirectorSkewSymTensor( SkewSymDirectorKNodeI, rVariables, rNodeI, k, rVariables.Alpha );

	    rDiscreteOperator += ( (rVariables.DN_DX(rNodeI, 0) * rVariables.N[rNodeJ]) -  (rVariables.N[rNodeI] * rVariables.DN_DX(rNodeJ, 0)) ) *  SkewSymDirectorKNodeI;

	    Matrix SkewSymDirectorK(3,3);
	    noalias(SkewSymDirectorK) = ZeroMatrix(3,3);
	    this->CalculateAlphaDirectorSkewSymTensor ( SkewSymDirectorK, rVariables, k, rVariables.Alpha );

	    kronecker = BeamMathUtilsType::KroneckerDelta(rNodeI, rNodeJ);

	    rDiscreteOperator += kronecker * ( rVariables.DN_DX(rNodeI, 0) )* SkewSymDirectorK;
	    //rDiscreteOperator += kronecker * ( rVariables.DN_DX(rNodeI, 0) )* SkewSymDirectorKNodeI; //incorrect

	    Matrix SkewSymDirectorKDerivatives(3,3);
	    noalias(SkewSymDirectorKDerivatives) = ZeroMatrix(3,3);
	    this->CalculateDirectorDerivativesSkewSymTensor ( SkewSymDirectorKDerivatives, rVariables, k, rVariables.Alpha );

	    rDiscreteOperator -= kronecker * ( rVariables.N[rNodeI] ) * SkewSymDirectorKDerivatives;

	    epsilon = BeamMathUtilsType::LeviCivitaEpsilon(j, rComponent, k);

	    rDiscreteOperator *= epsilon;

	  }
      }

    rDiscreteOperator *= 0.5 * StressCouples[rComponent];

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************


  void GeometricallyExactRodElement::CalculateAlphaDirectors(Matrix& rDirectors, ElementDataType& rVariables, const int& rNode, double alpha)
  {
    KRATOS_TRY

    DirectorsVariables& Directors = rVariables.GetDirectors();

    rDirectors = (1 - alpha) * Directors.PreviousNode[rNode] + (alpha) * Directors.CurrentNode[rNode];

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************


  void GeometricallyExactRodElement::CalculateAlphaDirectorVector(Vector& rDirectorVector, ElementDataType& rVariables, const int& rNode, const int& rDirection, double alpha)
  {
    KRATOS_TRY

    // nodal directors
    Matrix Directors(3,3);
    noalias(Directors) = ZeroMatrix(3,3);
    this->CalculateAlphaDirectors( Directors, rVariables, rNode, alpha );

    for ( unsigned int i = 0; i < 3; i++ )
      {
    	rDirectorVector[i] = Directors(i,rDirection);
      }

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************


  void GeometricallyExactRodElement::CalculateAlphaDirectorSkewSymTensor(Matrix& rDirectorSkewSymTensor, ElementDataType& rVariables, const int& rNode, const int& rDirection, double alpha)
  {
    KRATOS_TRY

    Vector DirectorVector(3);
    noalias(DirectorVector) = ZeroVector(3);

    this->CalculateAlphaDirectorVector( DirectorVector, rVariables, rNode, rDirection, alpha );

    BeamMathUtilsType::VectorToSkewSymmetricTensor(DirectorVector, rDirectorSkewSymTensor);

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************


  void GeometricallyExactRodElement::CalculateAlphaDirectorVector(Vector& rDirectorVector, ElementDataType& rVariables, const int& rDirection, double alpha)
  {
    KRATOS_TRY

    DirectorsVariables& Directors = rVariables.GetDirectors();

    rDirectorVector.resize(3,false);
    noalias(rDirectorVector) = ZeroVector(3);

    rDirectorVector = (1 - alpha) * Directors.Previous[rDirection] + (alpha) * Directors.Current[rDirection] ;

    KRATOS_CATCH( "" )
  }



  //************************************************************************************
  //************************************************************************************


  void GeometricallyExactRodElement::CalculateAlphaDirectorSkewSymTensor(Matrix& rDirectorSkewSymTensor, ElementDataType& rVariables, const int& rDirection, double alpha)
  {

    KRATOS_TRY

    Vector DirectorVector(3);
    noalias(DirectorVector) = ZeroVector(3);

    this->CalculateAlphaDirectorVector( DirectorVector, rVariables, rDirection, alpha );

    BeamMathUtilsType::VectorToSkewSymmetricTensor(DirectorVector, rDirectorSkewSymTensor);


    KRATOS_CATCH( "" )

      }


  //************************************************************************************
  //************************************************************************************

  void GeometricallyExactRodElement::CalculateDirectorDerivativesVector(Vector& rDirectorDerivativesVector, ElementDataType& rVariables, const int& rDirection, double alpha)
  {
    KRATOS_TRY

    DirectorsVariables& Directors = rVariables.GetDirectors();

    rDirectorDerivativesVector.resize(3,false);
    noalias(rDirectorDerivativesVector) = ZeroVector(3);

    rDirectorDerivativesVector = (1 - alpha) * Directors.PreviousDerivatives[rDirection] + (alpha) * Directors.CurrentDerivatives[rDirection] ;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************


  void GeometricallyExactRodElement::CalculateDirectorDerivativesSkewSymTensor(Matrix& rDirectorDerivativesSkewSymTensor, ElementDataType& rVariables, const int& rDirection, double alpha)
  {
    KRATOS_TRY

    Vector DirectorDerivativesVector(3);
    noalias(DirectorDerivativesVector) = ZeroVector(3);

    this->CalculateDirectorDerivativesVector( DirectorDerivativesVector, rVariables, rDirection, alpha );

    BeamMathUtilsType::VectorToSkewSymmetricTensor(DirectorDerivativesVector, rDirectorDerivativesSkewSymTensor);

    KRATOS_CATCH( "" )
  }



  //************************************************************************************
  //************************************************************************************

  void GeometricallyExactRodElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
							 ElementDataType& rVariables,
							 double& rIntegrationWeight)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    //Initialize Local Matrices
    MatrixType DifferentialOperatorI(12,6);
    noalias(DifferentialOperatorI) = ZeroMatrix(12,6);
    MatrixType DifferentialOperatorJ(12,6);
    noalias(DifferentialOperatorJ) = ZeroMatrix(12,6);

    MatrixType Kij(12,12);
    noalias(Kij) = ZeroMatrix(12,12);
    MatrixType Kmij(6,6);
    noalias(Kmij) = ZeroMatrix(6,6);

    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    MatrixType MappingTensorI(12,6);
    noalias(MappingTensorI) = ZeroMatrix(12,6);
    MatrixType MappingTensorJ(12,6);
    noalias(MappingTensorJ) = ZeroMatrix(12,6);


    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	RowIndex = i * (dimension * 2);

	this->CalculateDifferentialOperator( DifferentialOperatorI, rVariables, i , rVariables.Alpha );

	this->CalculateDirectorsMappingTensor( MappingTensorI, rVariables, i , rVariables.Alpha );

	//std::cout<<" DifferentialOperatorI "<<DifferentialOperatorI<<std::endl;

	for ( SizeType j = 0; j < number_of_nodes; j++ )
	  {

	    noalias(Kij) = ZeroMatrix(12,12);

	    ColIndex = j * (dimension * 2);

	    this->CalculateDifferentialOperator( DifferentialOperatorJ, rVariables, j , 1 );

	    //std::cout<<" DifferentialOperatorJ "<<DifferentialOperatorJ<<std::endl;

	    noalias(Kij) = prod( DifferentialOperatorI, Matrix(prod( rVariables.ConstitutiveMatrix, trans(DifferentialOperatorJ) )) );

	    Kij *= rIntegrationWeight;

	    this->CalculateDirectorsMappingTensor(MappingTensorJ, rVariables, j , 1 );

	    noalias(Kmij) = ZeroMatrix(6,6);

	    noalias(Kmij) = prod( trans(MappingTensorI), Matrix(prod( Kij, MappingTensorJ )) );

	    Kmij *= rVariables.Alpha; //0.5;

	    //Building the Local Stiffness Matrix
	    BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, Kmij, RowIndex, ColIndex );
	  }
      }

    //std::cout<<" Kuum "<<rLeftHandSideMatrix<<std::endl;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void GeometricallyExactRodElement::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
							 ElementDataType& rVariables,
							 double& rIntegrationWeight)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();


    //Initialize Local Matrices
    MatrixType Kij(6,6);
    noalias(Kij) = ZeroMatrix(6,6);

    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    Matrix GabK(3,3);
    noalias(GabK) = ZeroMatrix(3,3);
    Matrix SkewSymDirectorKNodeI(3,3);
    noalias(SkewSymDirectorKNodeI) = ZeroMatrix(3,3);
    Matrix SkewSymDirectorKNodeJ(3,3);
    noalias(SkewSymDirectorKNodeJ) = ZeroMatrix(3,3);

    // compute KN
    Vector StressResultants(3);
    noalias(StressResultants) = ZeroVector(3);
    for ( unsigned int i = 0; i < 3; i++ )
      {
	StressResultants[i] = rVariables.StressVector[i];
      }

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	RowIndex = i * (dimension * 2);


	for ( SizeType j = 0; j < number_of_nodes; j++ )
	  {

	    noalias(Kij) = ZeroMatrix(6,6);

	    ColIndex = j * (dimension * 2);

	    for( unsigned int k = 0; k < 3; k++)
	      {

		//term 11 -> 0
		//term 12
		noalias(GabK) = ZeroMatrix(3,3);
		noalias(SkewSymDirectorKNodeJ) = ZeroMatrix(3,3);
		this->CalculateAlphaDirectorSkewSymTensor( SkewSymDirectorKNodeJ, rVariables, j, k, 1 );
		GabK = (-1) * (rVariables.DN_DX(i, 0) *  rVariables.N[j] * StressResultants[k] ) * SkewSymDirectorKNodeJ;
		//Building the Local Stiffness Matrix
		BeamMathUtilsType::AddMatrix( Kij, GabK, 0, 3 );

		//term 21
		noalias(GabK) = ZeroMatrix(3,3);
		noalias(SkewSymDirectorKNodeI) = ZeroMatrix(3,3);
		this->CalculateAlphaDirectorSkewSymTensor( SkewSymDirectorKNodeI, rVariables, i, k, rVariables.Alpha );
		GabK = (rVariables.N[i] * rVariables.DN_DX(j, 0) * StressResultants[k] ) * SkewSymDirectorKNodeI;
		//Building the Local Stiffness Matrix
		BeamMathUtilsType::AddMatrix( Kij, GabK, 3, 0 );

		//term 22
		noalias(GabK) = ZeroMatrix(3,3);
		this->CalculateDiscreteOperatorN(GabK,rVariables,i,j,k);

		//Building the Local Stiffness Matrix
		BeamMathUtilsType::AddMatrix( Kij, GabK, 3, 3 );
	      }

	    Kij *= rIntegrationWeight * rVariables.Alpha; //0.5;

	    //Building the Local Stiffness Matrix
	    BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, Kij, RowIndex, ColIndex );

	  }

      }


    // compute KM
    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {

	RowIndex = i * (dimension * 2);


	for ( SizeType j = 0; j < number_of_nodes; j++ )
	  {

	    noalias(Kij) = ZeroMatrix(6,6);

	    ColIndex = j * (dimension * 2);

	    for( unsigned int k = 0; k < 3; k++)
	      {
		//term 11 -> 0
		//term 12 -> 0
		//term 21 -> 0
		//term 22
		noalias(GabK) = ZeroMatrix(3,3);
		noalias(SkewSymDirectorKNodeJ) = ZeroMatrix(3,3);

		this->CalculateDiscreteOperatorM(GabK,rVariables,i,j,k);
		this->CalculateAlphaDirectorSkewSymTensor( SkewSymDirectorKNodeJ, rVariables, j, k, 1 );

		GabK = (-1) * prod(GabK, SkewSymDirectorKNodeJ);

		//Building the Local Stiffness Matrix
		BeamMathUtilsType::AddMatrix( Kij, GabK, 3, 3 );
	      }

	    Kij *= rIntegrationWeight * rVariables.Alpha; //0.5;

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

  void GeometricallyExactRodElement::CalculateAndAddKuuf(MatrixType& rLeftHandSideMatrix,
							 ElementDataType& rVariables,
							 double& rIntegrationWeight)
  {
    KRATOS_TRY

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void GeometricallyExactRodElement::CalculateAndAddInertiaLHS(MatrixType& rLeftHandSideMatrix, ElementDataType& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight )
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


    //block m(1,1) of the mass matrix

    MatrixType m11(3,3);
    noalias(m11)= ZeroMatrix(3,3);

    double TotalMass = 0;
    TotalMass  = this->CalculateTotalMass( Section, TotalMass );

    //block m(2,2) of the mass matrix

    MatrixType m22(3,3);
    noalias(m22) = ZeroMatrix(3,3);


    //2.-Get inertia dyadic
    Matrix InertiaDyadic(3,3);
    noalias(InertiaDyadic) = ZeroMatrix(3,3);
    this->CalculateInertiaDyadic( Section, InertiaDyadic );

    // InertiaDyadic = prod(rVariables.CurrentRotationMatrix,InertiaDyadic);
    // InertiaDyadic = prod(InertiaDyadic,trans(rVariables.CurrentRotationMatrix));

    Matrix AlgorithmicInertia(3,3);
    noalias(AlgorithmicInertia) = ZeroMatrix(3,3);

    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    Matrix DiagonalMatrix = identity_matrix<double> (3);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	noalias(m11) = ZeroMatrix(3,3);
	noalias(m22) = ZeroMatrix(3,3);

	RowIndex = i * (dimension * 2);

	for ( SizeType j = 0; j < number_of_nodes; j++ )
	  {

	    ColIndex = j * (dimension * 2);

	    m11 = TotalMass * rVariables.N[i] * rVariables.N[j] * rIntegrationWeight * DiagonalMatrix;

	    noalias(AlgorithmicInertia) = ZeroMatrix(3,3);

	    this->CalculateAlgorithmicInertia( AlgorithmicInertia, InertiaDyadic, rVariables, j, i , rVariables.Alpha );

	    m22 = rVariables.N[i] * rVariables.N[j] * rIntegrationWeight * AlgorithmicInertia;

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

  void GeometricallyExactRodElement::CalculateAndAddInertiaRHS(VectorType& rRightHandSideVector, ElementDataType& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize            = number_of_nodes * ( dimension * 2 );

    if(rRightHandSideVector.size() != MatSize)
      rRightHandSideVector.resize(MatSize, false);

    noalias(rRightHandSideVector) = ZeroVector(MatSize);


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

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	//Current Linear Velocity Vector
	CurrentValueVector = GetNodalCurrentValue( VELOCITY, CurrentValueVector, i );
	CurrentLinearVelocityVector   +=  rVariables.N[i] * CurrentValueVector;

	//Previous Linear Velocity Vector
	CurrentValueVector = GetNodalPreviousValue( VELOCITY, CurrentValueVector, i );
	PreviousLinearVelocityVector  +=  rVariables.N[i] * CurrentValueVector;
      }

    CurrentLinearVelocityVector  = MapToInitialLocalFrame( CurrentLinearVelocityVector, rVariables.PointNumber );
    PreviousLinearVelocityVector = MapToInitialLocalFrame( PreviousLinearVelocityVector, rVariables.PointNumber );


    //Compute Linear Term:

    Vector LinearInertialForceVector(3);
    noalias(LinearInertialForceVector) = ZeroVector(3);
    LinearInertialForceVector  = TotalMass * (CurrentLinearVelocityVector - PreviousLinearVelocityVector);

    //std::cout<<" Linear Velocity "<<CurrentLinearVelocityVector<<std::endl;

    //-----------------

    std::vector<Vector> CurrentDirectorsVelocity;
    CurrentDirectorsVelocity.resize(3);
    std::vector<Vector> PreviousDirectorsVelocity;
    PreviousDirectorsVelocity.resize(3);

    DirectorsVariables& Directors = rVariables.GetDirectors();

    for ( unsigned int i = 0; i < 3; i++ )
      {
	CurrentDirectorsVelocity[i].resize(3,false);
	PreviousDirectorsVelocity[i].resize(3,false);

	noalias(CurrentDirectorsVelocity[i]) = ZeroVector(3);
	noalias(PreviousDirectorsVelocity[i]) = ZeroVector(3);
      }

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	Matrix& CurrentVelocity  = Directors.CurrentNodeVelocities[i];
	Matrix& PreviousVelocity = Directors.PreviousNodeVelocities[i];

	for ( unsigned int j = 0; j < 3; j++ )
	  {
	    for ( unsigned int k = 0; k < 3; k++ )
	      {
		CurrentDirectorsVelocity[j][k]  +=  rVariables.N[i] * ( CurrentVelocity(k,j) );
		PreviousDirectorsVelocity[j][k] +=  rVariables.N[i] * ( PreviousVelocity(k,j) );
	      }
	  }
      }


    //Compute Angular Term:

    //Get inertia dyadic
    Matrix InertiaDyadic(3,3);
    noalias(InertiaDyadic) = ZeroMatrix(3,3);
    this->CalculateInertiaDyadic( Section, InertiaDyadic );

    //director I
    Vector AngularInertialForceVectorI(3);
    noalias(AngularInertialForceVectorI) = ZeroVector(3);

    AngularInertialForceVectorI = prod( InertiaDyadic,CurrentDirectorsVelocity[0] ) - prod( InertiaDyadic,PreviousDirectorsVelocity[0] );

    //director II
    Vector AngularInertialForceVectorII(3);
    noalias(AngularInertialForceVectorII) = ZeroVector(3);

    AngularInertialForceVectorII = prod( InertiaDyadic,CurrentDirectorsVelocity[1] ) - prod( InertiaDyadic,PreviousDirectorsVelocity[1] );


    // Build incremental momentum vector
    Vector TotalInertialForceVector(12);
    noalias(TotalInertialForceVector) = ZeroVector(12);

    BeamMathUtilsType::AddVector(LinearInertialForceVector, TotalInertialForceVector, 0);
    BeamMathUtilsType::AddVector(AngularInertialForceVectorI, TotalInertialForceVector, 3);
    BeamMathUtilsType::AddVector(AngularInertialForceVectorII, TotalInertialForceVector, 6);


    //-----------------
    VectorType Fi(6);
    noalias(Fi) = ZeroVector(6);
    VectorType Fmi(12);
    noalias(Fmi) = ZeroVector(12);

    //Initialize Local Matrices
    unsigned int RowIndex = 0;

    MatrixType MappingTensor(12,6);
    noalias(MappingTensor) = ZeroMatrix(12,6);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {

	RowIndex = i * (dimension * 2);

	noalias(Fi) = ZeroVector(6);
	noalias(Fmi) = ZeroVector(12);

	Fmi = TotalInertialForceVector * rVariables.N[i] * rIntegrationWeight / DeltaTime;

	this->CalculateDirectorsMappingTensor(MappingTensor, rVariables, i, rVariables.Alpha);

	// std::cout<<" MappingTensor "<<trans(MappingTensor)<<std::endl;
	// std::cout<<" TotalInertialForceVector "<<Fmi<<std::endl;

	Fi = prod( trans(MappingTensor), Fmi );

	BeamMathUtilsType::AddVector(Fi, rRightHandSideVector, RowIndex);

	//std::cout<<" Fi "<<Fi<<std::endl;
	//std::cout<<" rRightHandSideVector "<<rRightHandSideVector<<std::endl;

      }

    //std::cout<<" rRightHandSideDynamic "<<rRightHandSideVector<<std::endl;

    KRATOS_CATCH( "" )
  }



  //************************************************************************************
  //************************************************************************************

  void GeometricallyExactRodElement::CalculateAlgorithmicInertia(Matrix & rAlgorithmicInertia, const Matrix& rInertiaDyadic, ElementDataType& rVariables, const int& rNodeJ, const int& rNodeI, double alpha)
  {
    KRATOS_TRY

    if( rAlgorithmicInertia.size1() != 3 )
      rAlgorithmicInertia.resize(3, 3, false);

    noalias(rAlgorithmicInertia) = ZeroMatrix(3,3);

    Matrix DiagonalMatrix(3,3);
    noalias(DiagonalMatrix) = IdentityMatrix(3);

    Vector DirectorKNodeI(3);
    noalias(DirectorKNodeI) = ZeroVector(3);
    Vector DirectorKNodeJ(3);
    noalias(DirectorKNodeJ) = ZeroVector(3);

    // std::cout<<" Inertia "<<rInertiaDyadic<<std::endl;

    for( unsigned int k = 0; k < 2; k++)
      {
	this->CalculateAlphaDirectorVector( DirectorKNodeI, rVariables, rNodeI, k, alpha );
	this->CalculateAlphaDirectorVector( DirectorKNodeJ, rVariables, rNodeJ, k, 1 );

	// std::cout<<" Director NodeI "<<DirectorKNodeI<<std::endl;
	// std::cout<<" Director NodeJ "<<DirectorKNodeJ<<std::endl;

	rAlgorithmicInertia += rInertiaDyadic(k,k) * inner_prod( DirectorKNodeJ, DirectorKNodeI ) * DiagonalMatrix;

	rAlgorithmicInertia -= rInertiaDyadic(k,k) * outer_prod( DirectorKNodeJ, DirectorKNodeI );
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
  int  GeometricallyExactRodElement::Check(const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    // Perform base element checks
    int ErrorCode = 0;
    ErrorCode = LargeDisplacementBeamEMCElement::Check(rCurrentProcessInfo);

    return ErrorCode;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************


  void GeometricallyExactRodElement::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementBeamEMCElement )
    rSerializer.save("InitialLocalDirectors",mInitialLocalDirectors);
    rSerializer.save("mCurrentLocalDirectors",mCurrentLocalDirectors);
    rSerializer.save("PreviousLocalDirectors",mPreviousLocalDirectors);
    rSerializer.save("InitialLocalDirectorsVelocities",mInitialLocalDirectorsVelocities);
    rSerializer.save("CurrentLocalDirectorsVelocities",mCurrentLocalDirectorsVelocities);
    rSerializer.save("PreviousLocalDirectorsVelocities",mPreviousLocalDirectorsVelocities);
  }

  void GeometricallyExactRodElement::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementBeamEMCElement )
    rSerializer.load("InitialLocalDirectors",mInitialLocalDirectors);
    rSerializer.load("mCurrentLocalDirectors",mCurrentLocalDirectors);
    rSerializer.load("PreviousLocalDirectors",mPreviousLocalDirectors);
    rSerializer.load("InitialLocalDirectorsVelocities",mInitialLocalDirectorsVelocities);
    rSerializer.load("CurrentLocalDirectorsVelocities",mCurrentLocalDirectorsVelocities);
    rSerializer.load("PreviousLocalDirectorsVelocities",mPreviousLocalDirectorsVelocities);

  }



} // Namespace Kratos


