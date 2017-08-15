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
#include "custom_elements/large_displacement_beam_element.hpp"

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

/**
 * Flags related to the element computation
 */
KRATOS_CREATE_LOCAL_FLAG( LargeDisplacementBeamElement, COMPUTE_RHS_VECTOR,                 0 );
KRATOS_CREATE_LOCAL_FLAG( LargeDisplacementBeamElement, COMPUTE_LHS_MATRIX,                 1 );
KRATOS_CREATE_LOCAL_FLAG( LargeDisplacementBeamElement, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
KRATOS_CREATE_LOCAL_FLAG( LargeDisplacementBeamElement, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementBeamElement::LargeDisplacementBeamElement(IndexType NewId,GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!    
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************


LargeDisplacementBeamElement::LargeDisplacementBeamElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
    KRATOS_TRY

    mThisIntegrationMethod = GeometryData::GI_GAUSS_1;

    //DO NOT ADD DOFS HERE!!!

    KRATOS_CATCH( "" )

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LargeDisplacementBeamElement::LargeDisplacementBeamElement(LargeDisplacementBeamElement const& rOther)
    :Element(rOther)
    ,mThisIntegrationMethod(rOther.mThisIntegrationMethod)
    ,mCurrentCurvatureVectors(rOther.mCurrentCurvatureVectors)
    ,mInitialLocalQuaternionsReduced(rOther.mInitialLocalQuaternionsReduced)
    ,mCurrentLocalQuaternionsReduced(rOther.mCurrentLocalQuaternionsReduced)
    ,mPreviousLocalQuaternionsReduced(rOther.mPreviousLocalQuaternionsReduced)
    ,mInitialLocalQuaternionsFull(rOther.mInitialLocalQuaternionsFull)
    ,mCurrentLocalQuaternionsFull(rOther.mCurrentLocalQuaternionsFull)
    ,mPreviousLocalQuaternionsFull(rOther.mPreviousLocalQuaternionsFull)
    ,mIterationCounter(rOther.mIterationCounter)
    ,mInvJ0Reduced(rOther.mInvJ0Reduced)
    ,mDetJ0Reduced(rOther.mDetJ0Reduced)
    ,mInvJ0Full(rOther.mInvJ0Full)
    ,mDetJ0Full(rOther.mDetJ0Full)
    ,mEnergy(rOther.mEnergy)
{
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer LargeDisplacementBeamElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new LargeDisplacementBeamElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementBeamElement::~LargeDisplacementBeamElement()
{
}


//************* GETTING METHODS
//************************************************************************************
//************************************************************************************

LargeDisplacementBeamElement::IntegrationMethod  LargeDisplacementBeamElement::GetIntegrationMethod() const
{
  return mThisIntegrationMethod;
  //return GeometryData::GI_GAUSS_3; //writing results with 3 integration points....
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
{
    ElementalDofList.resize(0);

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
      {
	ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
	ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
	ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));

	ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
	ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
	ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));
      }
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int element_size          = number_of_nodes * ( dimension * 2 );

    if ( rResult.size() != element_size )
        rResult.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      int index = i * ( dimension * 2 );
      rResult[index]   = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
      rResult[index+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
      rResult[index+2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();

      rResult[index+3] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
      rResult[index+4] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
      rResult[index+5] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
    }
 
}


//*********************************DISPLACEMENT***************************************
//************************************************************************************

void LargeDisplacementBeamElement::GetValuesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * ( dimension * 2 );

    if ( rValues.size() != element_size ) rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      int index = i * ( dimension * 2 );
      rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
      rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );
      rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );

      rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( ROTATION_X, Step );
      rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Y, Step );
      rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Z, Step );
    }

    KRATOS_CATCH( "" )
}

//************************************VELOCITY****************************************
//************************************************************************************

//************************************************************************************
//************************************************************************************
void LargeDisplacementBeamElement::GetFirstDerivativesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * ( dimension * 2 );

    if ( rValues.size() != element_size ) rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      int index = i * ( dimension * 2 );
      rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
      rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );
      rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );

      rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_X, Step );
      rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_Y, Step );
      rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_Z, Step );
    }


    KRATOS_CATCH( "" )
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void LargeDisplacementBeamElement::GetSecondDerivativesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * ( dimension * 2 );

    if ( rValues.size() != element_size ) rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      int index = i * ( dimension * 2 );
      rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
      rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );
      rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );

      rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_X, Step );
      rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_Y, Step );
      rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_Z, Step );
    }


    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

//*********************************GET DOUBLE VALUE***********************************
//************************************************************************************

void  LargeDisplacementBeamElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
							    std::vector<double>& rValues,
							    const ProcessInfo& rCurrentProcessInfo )
{ 
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}
//**********************************GET VECTOR VALUE**********************************
//************************************************************************************

void LargeDisplacementBeamElement::GetValueOnIntegrationPoints( const Variable<array_1d<double, 3 > >& rVariable,
								std::vector< array_1d<double, 3 > >& rValues,
								const ProcessInfo& rCurrentProcessInfo )
{

    this->CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);

}


//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::Initialize()
{
    KRATOS_TRY

    //------------- REDUCED QUADRATURE INTEGRATION 
      
    IntegrationMethod ReducedIntegrationMethod = GeometryData::GI_GAUSS_1;

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( ReducedIntegrationMethod );

    //Curvature Initialization
    if ( mCurrentCurvatureVectors.size() != integration_points.size() )
      {
        mCurrentCurvatureVectors.resize( integration_points.size() );
      }
    
    for ( unsigned int i = 0; i < mCurrentCurvatureVectors.size(); i++ )
      {
	mCurrentCurvatureVectors[i] = ZeroVector(3);
      }


    //-------------

    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    Matrix Identity = IdentityMatrix(dimension);
   
    //Local Quaternions Initialization
    if ( mInitialLocalQuaternionsReduced.size() != integration_points.size() )
      {
        mInitialLocalQuaternionsReduced.resize( integration_points.size() );
      }

    //Local Quaternions Initialization
    if ( mCurrentLocalQuaternionsReduced.size() != integration_points.size() )
      {
        mCurrentLocalQuaternionsReduced.resize( integration_points.size() );
      }

    //Local Quaternions Initialization
    if ( mPreviousLocalQuaternionsReduced.size() != integration_points.size() )
      {
        mPreviousLocalQuaternionsReduced.resize( integration_points.size() );
      }


    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {

	mInitialLocalQuaternionsReduced[PointNumber]  = QuaternionType::FromRotationMatrix( Identity );

	mCurrentLocalQuaternionsReduced[PointNumber]  = QuaternionType::FromRotationMatrix( Identity );

	mPreviousLocalQuaternionsReduced[PointNumber] = QuaternionType::FromRotationMatrix( Identity ); 
	//---

	Matrix LocalTransformationMatrix = ZeroMatrix(dimension,dimension);
	this->CalculateLocalAxesMatrix( LocalTransformationMatrix );
    
	//Quaternions initialization
	mInitialLocalQuaternionsReduced[PointNumber] = QuaternionType::FromRotationMatrix( LocalTransformationMatrix );
      }



    //-------------

    //Resizing jacobian inverses container
    mInvJ0Reduced.resize( integration_points.size() );
    mDetJ0Reduced.resize( integration_points.size(), false );

    //Compute jacobian inverses and set the domain initial size:
    GeometryType::JacobiansType J0;
    J0 = GetGeometry().Jacobian( J0, ReducedIntegrationMethod );

    //calculating the inverse J0
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {
	Vector Jacobian = ZeroVector(3);
    
	//Calculating the jacobian [dx_0/d£]
	for( unsigned int i=0; i<dimension; i++)
	  {
	    Jacobian[i] = J0[PointNumber](i,0);
	  }


	mDetJ0Reduced[PointNumber] = norm_2(Jacobian); 

	//Calculating the inverse of the jacobian and the parameters needed [d£/dx_0]
	double InverseJacobian = 1.0/ mDetJ0Reduced[PointNumber]; 
      
	mInvJ0Reduced[PointNumber] = ZeroVector(dimension);
      
	for( unsigned int i=0; i<dimension; i++)
	  {
	    mInvJ0Reduced[PointNumber][i] = InverseJacobian;
	  }

				     
      }


    //------------- FULL QUADRATURE INTEGRATION 

    IntegrationMethod FullIntegrationMethod = GeometryData::GI_GAUSS_2;

    const GeometryType::IntegrationPointsArrayType& integration_points_full = GetGeometry().IntegrationPoints( FullIntegrationMethod );

   
    //Local Quaternions Initialization
    if ( mInitialLocalQuaternionsFull.size() != integration_points_full.size() )
      {
        mInitialLocalQuaternionsFull.resize( integration_points_full.size() );
      }

    //Local Quaternions Initialization
    if ( mCurrentLocalQuaternionsFull.size() != integration_points_full.size() )
      {
        mCurrentLocalQuaternionsFull.resize( integration_points_full.size() );
      }
   
    //Local Quaternions Initialization
    if ( mPreviousLocalQuaternionsFull.size() != integration_points_full.size() )
      {
        mPreviousLocalQuaternionsFull.resize( integration_points_full.size() );
      }

    for ( unsigned int PointNumber = 0; PointNumber < integration_points_full.size(); PointNumber++ )
      {

	mInitialLocalQuaternionsFull[PointNumber]  = QuaternionType::FromRotationMatrix( Identity );

       	mCurrentLocalQuaternionsFull[PointNumber]  = QuaternionType::FromRotationMatrix( Identity );

	mPreviousLocalQuaternionsFull[PointNumber] = QuaternionType::FromRotationMatrix( Identity );

	//---

	Matrix LocalTransformationMatrix = ZeroMatrix(dimension,dimension);
	this->CalculateLocalAxesMatrix( LocalTransformationMatrix );
    
	//Quaternions initialization
	mInitialLocalQuaternionsFull[PointNumber] = QuaternionType::FromRotationMatrix( LocalTransformationMatrix );

      }


    //-------------

    //Resizing jacobian inverses container
    mInvJ0Full.resize( integration_points_full.size() );
    mDetJ0Full.resize( integration_points_full.size(), false );

    //Compute jacobian inverses and set the domain initial size:
    J0 = GetGeometry().Jacobian( J0, FullIntegrationMethod );

    //calculating the inverse J0
    for ( unsigned int PointNumber = 0; PointNumber < integration_points_full.size(); PointNumber++ )
      {
	Vector Jacobian = ZeroVector(3);
    
	//Calculating the jacobian [dx_0/d£]
	for( unsigned int i=0; i<dimension; i++)
	  {
	    Jacobian[i] = J0[PointNumber](i,0);
	  }


	mDetJ0Full[PointNumber] = norm_2(Jacobian); 

	//Calculating the inverse of the jacobian and the parameters needed [d£/dx_0]
	double InverseJacobian = 1.0/ mDetJ0Full[PointNumber]; 
      
	mInvJ0Full[PointNumber] = ZeroVector(dimension);
      
	for( unsigned int i=0; i<dimension; i++)
	  {
	    mInvJ0Full[PointNumber][i] = InverseJacobian;
	  }
			     
      }


    mEnergy.clear();
    

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    mIterationCounter = 0;

    for ( unsigned int i = 0; i < mPreviousLocalQuaternionsReduced.size(); i++ )
      {
	mPreviousLocalQuaternionsReduced[i] =  mCurrentLocalQuaternionsReduced[i];
      }

    for ( unsigned int i = 0; i < mPreviousLocalQuaternionsFull.size(); i++ )
      {
	mPreviousLocalQuaternionsFull[i] =  mCurrentLocalQuaternionsFull[i];
      }

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //Curvature Update
    if( mIterationCounter == 0 ){ //predict is done after initialize solution step
      this->CalculateCurvatureUpdate( mCurrentCurvatureVectors, STEP_ROTATION, rCurrentProcessInfo ); 
    }

    mEnergy.clear();

    KRATOS_CATCH( "" )
}

////************************************************************************************
////************************************************************************************

void LargeDisplacementBeamElement::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
     KRATOS_TRY
    
     //Curvature Update (on 14/07/2016 the FinalizeNonLinearIteration changed in solver, it is done after movemesh...problems¿?)
     this->CalculateCurvatureUpdate( mCurrentCurvatureVectors, DELTA_ROTATION, rCurrentProcessInfo ); 

     mIterationCounter++;

     KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
							    VectorType& rRightHandSideVector,
							    Flags& rCalculationFlags)

{

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * ( dimension * 2 );

    if ( rCalculationFlags.Is(LargeDisplacementBeamElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( rCalculationFlags.Is(LargeDisplacementBeamElement::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
	    rRightHandSideVector.resize( MatSize, false );
      
	rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
	  
    }
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//************************************************************************************
//************************************************************************************

Vector& LargeDisplacementBeamElement::MapToInitialLocalFrame(Vector& rVariable, unsigned int PointNumber)
{
    KRATOS_TRY
      
    // Global System coordinate to local beam system coordinates (Total Lagrangian Local Reference Configuration)
    if( GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod) == 1 ){
      BeamMathUtilsType::MapToCurrentLocalFrame(mInitialLocalQuaternionsReduced[PointNumber], rVariable);    
    }
    if( GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod) == 2 ){      
      BeamMathUtilsType::MapToCurrentLocalFrame(mInitialLocalQuaternionsFull[PointNumber], rVariable);
    }   
    
    return rVariable;

    KRATOS_CATCH( "" )
}


///************************************************************************************
//************************************************************************************


void LargeDisplacementBeamElement::CalculateCurvatureUpdate(std::vector<Vector>& rCurvatureVectors, const Variable<array_1d<double, 3 > >& rVariable, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    //(in fact is the two nodes beam element, only one integration point)
    for ( unsigned int PointNumber = 0; PointNumber < rCurvatureVectors.size(); PointNumber++ )
      {
	//compute element kinematics  ...

	//std::cout<<" Curvature Update "<<std::endl;
	this->CalculateKinematics(Variables,PointNumber);
	
	Variables.PreviousCurvatureVector = rCurvatureVectors[PointNumber];

	this->CalculateCurrentCurvature(Variables,rVariable);

	rCurvatureVectors[PointNumber] = Variables.CurrentCurvatureVector;
      }

    
    KRATOS_CATCH( "" )
}


///************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateCurrentCurvature(GeneralVariables& rVariables, const Variable<array_1d<double, 3 > >& rVariable)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();

    Vector CurrentStepRotationDerivativesVector = ZeroVector(3);
    Vector CurrentStepRotationVector            = ZeroVector(3);
    Vector CurrentValueVector                   = ZeroVector(3);

    //jelenic:crisfield [j-c]:start
    QuaternionType QuaternionValue;
    std::vector<QuaternionType> CurrentStepRotationQuaternions;
    //------------------[j-c]:end

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
    	//Current Rotation Derivatives	
    	CurrentValueVector = GetNodalCurrentValue( rVariable, CurrentValueVector, i );
    	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rVariables.PointNumber );

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
    // Vector EquationVectorPartA = ZeroVector(3);
    // Vector EquationVectorPartB = ZeroVector(3);
    // Vector EquationVectorPartC = ZeroVector(3);
    //------------------[option 0]:end

    Matrix EquationMatrixPartA = ZeroMatrix(3,3);
    Matrix EquationMatrixPartB = ZeroMatrix(3,3);
    Matrix EquationMatrixPartC = ZeroMatrix(3,3);
    
    Matrix DiagonalMatrix = IdentityMatrix(3);
	
    Matrix SkewSymCurrentStepRotation = ZeroMatrix(3,3);
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

void LargeDisplacementBeamElement::MapToSpatialFrame(const GeneralVariables& rVariables, Matrix& rVariable)
{

    KRATOS_TRY
     
    BeamMathUtilsType::MapLocalToGlobal3D(rVariables.CurrentRotationMatrix, rVariable);

    KRATOS_CATCH( "" )

}
//************************************************************************************
//************************************************************************************

Vector& LargeDisplacementBeamElement::GetLocalCurrentValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const Vector& rN)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    Vector CurrentValueVector = ZeroVector(dimension);

    //strains due to displacements and rotations
    for ( unsigned int i = 0; i <number_of_nodes ; i++ )
      {
	CurrentValueVector = GetNodalCurrentValue( rVariable, CurrentValueVector, i );
	rValue += rN[i] * CurrentValueVector;
      }
    
    //Current Frame is the Local Frame
    rValue = MapToInitialLocalFrame( rValue );

    return rValue;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

Vector& LargeDisplacementBeamElement::GetLocalPreviousValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const Vector& rN)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    Vector PreviousValueVector = ZeroVector(dimension);

    //strains due to displacements and rotations
    for ( unsigned int i = 0; i <number_of_nodes ; i++ )
      {
	PreviousValueVector = GetNodalPreviousValue( rVariable, PreviousValueVector, i );
	rValue += rN[i] * PreviousValueVector;
      }
    
    //Current Frame is the Local Frame
    rValue = MapToInitialLocalFrame( rValue );

    return rValue;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

Vector& LargeDisplacementBeamElement::GetNodalCurrentValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const unsigned int& rNode)
{
    KRATOS_TRY
      
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    
    if( rValue.size() != dimension )
      rValue.resize(dimension, false);

    rValue = GetGeometry()[rNode].FastGetSolutionStepValue( rVariable );
    
    return rValue;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

Vector& LargeDisplacementBeamElement::GetNodalPreviousValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const unsigned int& rNode)
{
    KRATOS_TRY

    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if( rValue.size() != dimension )
      rValue.resize(dimension, false);

    rValue = GetGeometry()[rNode].FastGetSolutionStepValue( rVariable, 1 );
    
    return rValue;

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::InitializeGeneralVariables(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rVariables.Alpha = 1;

    //Compute Section Properties:
    this->CalculateSectionProperties(rVariables.Section);

    rVariables.Length = GetGeometry().Length();

    if(rVariables.Length == 0.00)
      KRATOS_THROW_ERROR(std::invalid_argument, "Zero length found in element #", this->Id())


    unsigned int voigtsize = 3;
    if( dimension == 3 )
    {
        voigtsize  = 6;
    }

    rVariables.ConstitutiveMatrix.resize( voigtsize, voigtsize );

    rVariables.CurrentRotationMatrix.resize( dimension, dimension );
 
    rVariables.CurrentAxisPositionDerivatives.resize( dimension );
   
    rVariables.StrainVector.resize( voigtsize );
    rVariables.StressVector.resize( voigtsize );

    rVariables.DN_DX.resize( number_of_nodes, 1 ); //they are the same in all dimensions local 1D

    //set variables including all integration points values

    //reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    //reading shape functions local gradients
    rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

    //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );

    //Calculate Delta Position
    rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateSectionProperties(SectionProperties & rSection)
{
    KRATOS_TRY
   
    if( GetProperties().Has(CROSS_SECTION_AREA) ){
        rSection.Area = GetProperties()[CROSS_SECTION_AREA];
    }
    else{
        rSection.Area = GetValue(CROSS_SECTION_AREA);
    }


    if( GetProperties().Has(LOCAL_INERTIA_TENSOR) )
    {
        Matrix& inertia = GetProperties()[LOCAL_INERTIA_TENSOR];
        rSection.Inertia_z = inertia(0,0);
        rSection.Inertia_y = inertia(1,1);
        rSection.Polar_Inertia = inertia(0,1);
    }
    else
    {
        Matrix& inertia = GetValue(LOCAL_INERTIA_TENSOR);
        rSection.Inertia_z = inertia(0,0);
        rSection.Inertia_y = inertia(1,1);
        rSection.Polar_Inertia = inertia(0,1);
    }

    rSection.Rotational_Inertia = rSection.Polar_Inertia;

    KRATOS_CATCH( "" )

}


//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateKinematics(GeneralVariables& rVariables, const unsigned int& rPointNumber)
{
    KRATOS_TRY


    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().size();

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //Parent to reference configuration
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    //point number
    rVariables.PointNumber = rPointNumber;

    //Set Shape Functions Values for this integration point
    rVariables.N=row( Ncontainer, rPointNumber);
    
    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();

    //TOTAL LAGRANGIAN
    //Compute cartesian derivatives [dN/dx_0]    
    if( GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod) == 1 ){
      rVariables.DN_DX = mInvJ0Reduced[rPointNumber][0] * DN_De[rPointNumber]; 
      rVariables.detJ  = mDetJ0Reduced[rPointNumber];
    }
    
    if( GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod) == 2 ){
      rVariables.DN_DX = mInvJ0Full[rPointNumber][0] * DN_De[rPointNumber]; 
      rVariables.detJ  = mDetJ0Full[rPointNumber];
    }

    //Compute centroid displacement derivatives:
    Vector CurrentDisplacements  = ZeroVector(dimension);
    Vector DeltaDisplacements    = ZeroVector(dimension);

    rVariables.CurrentAxisPositionDerivatives = ZeroVector(dimension);
   
    //compute local to global frame
    this->CalculateFrameMapping(rVariables,rPointNumber);

    Vector CurrentValueVector = ZeroVector(3);
    Matrix CurrentValueMatrix = ZeroMatrix(3,3);
    QuaternionType QuaternionValue;

    //strains due to displacements and rotations
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {

	//A: Current Nodes Position
	CurrentValueVector = GetGeometry()[i].Coordinates();
	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rPointNumber );

	//Current Frame Axis Position derivative
	rVariables.CurrentAxisPositionDerivatives +=  rVariables.DN_DX(i,0) * ( CurrentValueVector );

      }
    

    //set current CURVATURES
    rVariables.CurrentCurvatureVector = mCurrentCurvatureVectors[rPointNumber];

    KRATOS_CATCH( "" )
}



//*************************COMPUTE DELTA POSITION*************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateFrameMapping(GeneralVariables& rVariables,const unsigned int& rPointNumber)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();

    Vector CurrentValueVector  = ZeroVector(3);

    std::vector<QuaternionType> NodeQuaternions;
    QuaternionType QuaternionValue;
    
    //strains due to displacements and rotations
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {

	//Rotations	
	CurrentValueVector = GetNodalCurrentValue( ROTATION, CurrentValueVector, i );

	//Current Frame is the Local Frame
	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rPointNumber );

	QuaternionValue = QuaternionType::FromRotationVector(CurrentValueVector);
	
	NodeQuaternions.push_back(QuaternionValue);

      }


    QuaternionType PreviousLocalQuaternion;
    if( GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod) == 1 ){
	PreviousLocalQuaternion = mPreviousLocalQuaternionsReduced[rPointNumber];
	mPreviousLocalQuaternionsReduced[rPointNumber].ToRotationMatrix(rVariables.PreviousRotationMatrix);
    }

    if( GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod) == 2 ){
      PreviousLocalQuaternion = mPreviousLocalQuaternionsFull[rPointNumber];
      mPreviousLocalQuaternionsFull[rPointNumber].ToRotationMatrix(rVariables.PreviousRotationMatrix);
    }

    //option frame 1:
    Vector CurrentStepRotationVector  = ZeroVector(3);
    this->GetLocalCurrentValue(STEP_ROTATION, CurrentStepRotationVector, rVariables.N);

    Matrix ExponentialRotationMatrix = ZeroMatrix(3,3);

    QuaternionValue = QuaternionType::FromRotationVector(CurrentStepRotationVector);

    QuaternionValue.ToRotationMatrix(ExponentialRotationMatrix);

    rVariables.CurrentRotationMatrix = prod(ExponentialRotationMatrix, rVariables.PreviousRotationMatrix);


    //option frame 2:
    // QuaternionValue = (NodeQuaternions.front()).conjugate() * (NodeQuaternions.back());

    // QuaternionValue.ToRotationVector(CurrentValueVector);

    // CurrentValueVector *= rVariables.N[0];
    
    // QuaternionValue = QuaternionType::FromRotationVector(CurrentValueVector);

    // QuaternionValue = (NodeQuaternions.front()) * QuaternionValue;
    
    // // Current Rotation Matrix
    // QuaternionValue.ToRotationMatrix(rVariables.CurrentRotationMatrix);

    //*------------------------------*//

    //update the global to local transformation matrix with the frame rotation matrix
    rVariables.LocalTransformationMatrix = ZeroMatrix(3,3);

    QuaternionType SpatialLocalQuaternion;
    if( GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod) == 1 ){
      SpatialLocalQuaternion = mInitialLocalQuaternionsReduced[rPointNumber];
    }
     
    if( GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod) == 2 ){
      SpatialLocalQuaternion = mInitialLocalQuaternionsFull[rPointNumber];
    }

    SpatialLocalQuaternion.ToRotationMatrix(rVariables.LocalTransformationMatrix);

    //set variables for the initialization update
    this->UpdateRotationVariables(rVariables,rPointNumber);

    KRATOS_CATCH( "" )
}

//*********************************SET ROTATION VARIABLES*****************************
//************************************************************************************

void LargeDisplacementBeamElement::UpdateRotationVariables(GeneralVariables& rVariables, const unsigned int& rPointNumber)
{
    KRATOS_TRY

    if( GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod) == 1 ){
      mCurrentLocalQuaternionsReduced[rPointNumber] = QuaternionType::FromRotationMatrix(rVariables.CurrentRotationMatrix);
    }

    if( GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod) == 2 ){
      mCurrentLocalQuaternionsFull[rPointNumber] = QuaternionType::FromRotationMatrix(rVariables.CurrentRotationMatrix);
    }

    KRATOS_CATCH( "" )
}

//*********************************SET STRAIN VARIABLES*******************************
//************************************************************************************

void LargeDisplacementBeamElement::UpdateStrainVariables(GeneralVariables& rVariables, const unsigned int& rPointNumber)
{
    KRATOS_TRY


    KRATOS_CATCH( "" )
}

//*************************COMPUTE DELTA POSITION*************************************
//************************************************************************************


Matrix& LargeDisplacementBeamElement::CalculateDeltaPosition(Matrix & rDeltaPosition)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rDeltaPosition = zero_matrix<double>( number_of_nodes , dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rDeltaPosition(i,j) = CurrentDisplacement[j]-PreviousDisplacement[j];
        }

    }

    return rDeltaPosition;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************


void LargeDisplacementBeamElement::CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
							     ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //create and initialize element variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);


    IntegrationMethod ThisIntegrationMethod = mThisIntegrationMethod;
    mThisIntegrationMethod = GeometryData::GI_GAUSS_1;

    //reading integration points (in fact is the two nodes beam element, only one integration point)
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //auxiliary terms
    Vector VolumeForce;
    
    //strain energy clear
    mEnergy.Deformation = 0;

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

	// std::cout<<" ConstitutiveMatrix "<<Variables.ConstitutiveMatrix<<std::endl;

	//compute element Strain and Stress Resultants and Couples
	this->CalculateStressResultants(Variables, PointNumber, Variables.Alpha);

	// std::cout<<" StrainResultants "<<Variables.StrainVector<<std::endl;
	// std::cout<<" StressResultants "<<Variables.StressVector<<std::endl;

	
	double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
	IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );


	if ( rLocalSystem.CalculationFlags.Is(LargeDisplacementBeamElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
	  {
	    this->CalculateAndAddLHS( rLocalSystem, Variables, IntegrationWeight );
	  }

	if ( rLocalSystem.CalculationFlags.Is(LargeDisplacementBeamElement::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
	  {
	    //contribution to external forces
	    VolumeForce  = this->CalculateVolumeForce( VolumeForce, Variables.N );

	    this->CalculateAndAddRHS( rLocalSystem , Variables, VolumeForce, IntegrationWeight );
	  }

	//Strain energy calculation
	this->CalculateStrainEnergy(Variables, rCurrentProcessInfo, IntegrationWeight);
      }


    // VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
    // std::cout<<"["<<this->Id()<<"] Beam RotatedStatic rRightHandSideVector "<<rRightHandSideVector<<std::endl;
    // MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();
    // std::cout<<"["<<this->Id()<<"] Beam RotatedStatic rLeftHandSideMatrix "<<rLeftHandSideMatrix<<std::endl;

    mThisIntegrationMethod = ThisIntegrationMethod;
    
    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************


void LargeDisplacementBeamElement::CalculateDynamicSystem( LocalSystemComponents& rLocalSystem,
							   ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //create and initialize element variables:
    GeneralVariables Variables;

    IntegrationMethod ThisIntegrationMethod = mThisIntegrationMethod;
    //full quadrature integration:
    mThisIntegrationMethod = GeometryData::GI_GAUSS_2;

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    // initialize variables short version:

    //std::cout<<" ID "<<this->Id()<<std::endl;

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
	Variables.N=row( Ncontainer, PointNumber);

	//compute local to global frame
	this->CalculateFrameMapping( Variables, PointNumber );
	
	Variables.detJ = mDetJ0Full[PointNumber];
 
	double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
	IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );
	
	if ( rLocalSystem.CalculationFlags.Is(LargeDisplacementBeamElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
	  {
	    LocalLeftHandSideMatrix.clear();
	    
    	    this->CalculateAndAddInertiaLHS( LocalLeftHandSideMatrix, Variables, rCurrentProcessInfo, IntegrationWeight ); // (R_N+1, R_N)

	    BeamMathUtilsType::MapLocalToGlobal3D(Variables.LocalTransformationMatrix,LocalLeftHandSideMatrix);

	    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();
	    rLeftHandSideMatrix += LocalLeftHandSideMatrix;

	    //std::cout<<"["<<this->Id()<<"] Beam RotatedDynamic rLeftHandSideMatrix "<<rLeftHandSideMatrix<<std::endl;

	  }

	if ( rLocalSystem.CalculationFlags.Is(LargeDisplacementBeamElement::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
	  {
	    LocalRightHandSideVector.clear();
	    
	    this->CalculateAndAddInertiaRHS( LocalRightHandSideVector, Variables, rCurrentProcessInfo, IntegrationWeight );
	    
	    BeamMathUtilsType::MapLocalToGlobal3D(Variables.LocalTransformationMatrix,LocalRightHandSideVector);

	    VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
	    rRightHandSideVector += LocalRightHandSideVector;

	    //std::cout<<"["<<this->Id()<<"] Beam RotatedDynamic rRightHandSideVector "<<rRightHandSideVector<<std::endl;
	  }

	//Kinetic energy calculation
	this->CalculateKineticEnergy(Variables, rCurrentProcessInfo, IntegrationWeight);

	//Momentum calculation
	this->CalculateMomentumRelations(Variables, rCurrentProcessInfo, IntegrationWeight);

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

void LargeDisplacementBeamElement::GetMaterialConstitutiveMatrix(Matrix& rConstitutiveMatrix, GeneralVariables& rVariables)
{
    KRATOS_TRY

    if( rConstitutiveMatrix.size1() != 6 || rConstitutiveMatrix.size2() != 6)
      rConstitutiveMatrix.resize(6,6, false);

    rConstitutiveMatrix = ZeroMatrix(6,6);

    if(GetProperties().Has(LOCAL_CONSTITUTIVE_MATRIX)){

      //Axis local E3
      rConstitutiveMatrix = GetProperties()[LOCAL_CONSTITUTIVE_MATRIX];

    }
    else{
      
      const double PoissonCoefficient = GetProperties()[POISSON_RATIO];
      const double YoungModulus       = GetProperties()[YOUNG_MODULUS];
      const double ShearModulus       = YoungModulus*0.5/(1.0 + PoissonCoefficient);


      //Axis local E3
      rConstitutiveMatrix( 0, 0 ) = ShearModulus * rVariables.Section.Area;
      rConstitutiveMatrix( 1, 1 ) = ShearModulus * rVariables.Section.Area;
      rConstitutiveMatrix( 2, 2 ) = YoungModulus * rVariables.Section.Area;

      rConstitutiveMatrix( 3, 3 ) = YoungModulus * rVariables.Section.Inertia_y;
      rConstitutiveMatrix( 4, 4 ) = YoungModulus * rVariables.Section.Inertia_z;
      rConstitutiveMatrix( 5, 5 ) = ShearModulus * rVariables.Section.Polar_Inertia;

    }
    //std::cout<<" ConstitutiveMatrix "<<rConstitutiveMatrix<<std::endl;

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateConstitutiveMatrix(GeneralVariables& rVariables)
{
    KRATOS_TRY
      
    //Material Elastic constitutive matrix
    this->GetMaterialConstitutiveMatrix(rVariables.ConstitutiveMatrix, rVariables);
    

    //Spatial Elastic constitutive matrix
    this->MapToSpatialFrame( rVariables, rVariables.ConstitutiveMatrix);


    KRATOS_CATCH( "" )
}




//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateStressResultants(GeneralVariables& rVariables, const unsigned int& rPointNumber, double alpha)
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    Vector StrainResultants = rVariables.CurrentAxisPositionDerivatives;
    Vector StrainCouples    = rVariables.CurrentCurvatureVector;

    //Reference frame given by the Frame Rotation
    StrainResultants = prod( trans(rVariables.CurrentRotationMatrix), StrainResultants );
    StrainCouples    = prod( trans(rVariables.CurrentRotationMatrix), StrainCouples );

    Vector E3 = ZeroVector(3);
    E3[2] = 1.0;
    StrainResultants -= E3;


    for ( unsigned int i = 0; i < dimension; i++ )
      {
	rVariables.StrainVector[i]   = StrainResultants[i];
	rVariables.StrainVector[i+3] = StrainCouples[i];
      }

    
    Matrix ConstitutiveMatrix  = ZeroMatrix(6,6);
    this->GetMaterialConstitutiveMatrix(ConstitutiveMatrix, rVariables);

    //Reference Stress Vector
    rVariables.StressVector = prod( ConstitutiveMatrix, rVariables.StrainVector );

    Vector StressResultants = ZeroVector(dimension);
    Vector StressCouples    = ZeroVector(dimension);

    for ( unsigned int i = 0; i < dimension; i++ )
      {
	StressResultants[i] = rVariables.StressVector[i];
	StressCouples[i]    = rVariables.StressVector[i+3];
      }

    //Current frame given by the Frame Rotation
    StressResultants = prod( rVariables.CurrentRotationMatrix, StressResultants );
    StressCouples    = prod( rVariables.CurrentRotationMatrix, StressCouples );

    for ( unsigned int i = 0; i < dimension; i++ )
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
void LargeDisplacementBeamElement::CalculateStrainEnergy(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight)
{
    KRATOS_TRY

    //Internal Energy Calculation: alpha = 1     
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //compute Strain Resultants and Couples
    Vector StrainResultants = rVariables.CurrentAxisPositionDerivatives;
    Vector StrainCouples    = rVariables.CurrentCurvatureVector;

    //Reference frame given by the Frame Rotation
    StrainResultants = prod( trans(rVariables.CurrentRotationMatrix), StrainResultants );
    StrainCouples    = prod( trans(rVariables.CurrentRotationMatrix), StrainCouples );

    Vector E3 = ZeroVector(3);
    E3[2] = 1.0;
    StrainResultants -= E3;

    //----------------

    Matrix ConstitutiveMatrix  = ZeroMatrix(6,6);
    this->GetMaterialConstitutiveMatrix(ConstitutiveMatrix, rVariables);

    //----------------

    Vector StrainVector = ZeroVector(6);
    for ( unsigned int i = 0; i < dimension; i++ )
      {
    	StrainVector[i]   = StrainResultants[i];
    	StrainVector[i+3] = StrainCouples[i];
      }

    //Reference Stress Vector
    Vector StressVector = prod( ConstitutiveMatrix, StrainVector );

    mEnergy.Deformation += 0.5 * (inner_prod(StressVector, StrainVector)) * rIntegrationWeight ;

    
    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

double& LargeDisplacementBeamElement::CalculateIntegrationWeight(double& rIntegrationWeight)
{
    KRATOS_TRY

    return rIntegrationWeight;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
{
    KRATOS_TRY
 
    Flags       LocalFlags;
    MatrixType  LocalLeftHandSideMatrix;
    VectorType  LocalRightHandSideVector;

    //Initialize sizes for the local system components:
    LocalFlags.Set(LargeDisplacementBeamElement::COMPUTE_LHS_MATRIX);
    this->InitializeSystemMatrices( LocalLeftHandSideMatrix, LocalRightHandSideVector, LocalFlags );
    

    // Local material stiffness
    this->CalculateAndAddKuum( LocalLeftHandSideMatrix, rVariables, rIntegrationWeight );

    //std::cout<<" rLeftHandSideMatrix Kuum "<<LocalLeftHandSideMatrix<<std::endl;

    // Local geometrical stiffness (not used in quasi-static additive rotations it degradates convergence)
    if( mIterationCounter < 2 || mIterationCounter > 5 ) //TEST
      this->CalculateAndAddKuug( LocalLeftHandSideMatrix, rVariables, rIntegrationWeight );
    

    // Local follower load stiffness
    this->CalculateAndAddKuuf( LocalLeftHandSideMatrix, rVariables, rIntegrationWeight );


    // LocalToGlobalSystem for the correct assembly
    BeamMathUtilsType::MapLocalToGlobal3D(rVariables.LocalTransformationMatrix,LocalLeftHandSideMatrix);
    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();
    rLeftHandSideMatrix += LocalLeftHandSideMatrix;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
{
    KRATOS_TRY

    Flags       LocalFlags;
    MatrixType  LocalLeftHandSideMatrix;
    VectorType  LocalRightHandSideVector;

    //Initialize sizes for the local system components:
    LocalFlags.Set(LargeDisplacementBeamElement::COMPUTE_RHS_VECTOR);
    this->InitializeSystemMatrices( LocalLeftHandSideMatrix, LocalRightHandSideVector, LocalFlags );
    
    // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
    this->CalculateAndAddExternalForces( LocalRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

    // operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
    this->CalculateAndAddInternalForces( LocalRightHandSideVector, rVariables, rIntegrationWeight );

    // follower load forces (to implement)
    this->CalculateAndAddFollowerForces( LocalRightHandSideVector, rVariables, rIntegrationWeight );
    

    BeamMathUtilsType::MapLocalToGlobal3D(rVariables.LocalTransformationMatrix,LocalRightHandSideVector);
    VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
    rRightHandSideVector += LocalRightHandSideVector;

    KRATOS_CATCH( "" )

}


//************************************************************************************
//************************************************************************************

void  LargeDisplacementBeamElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
							   ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(LargeDisplacementBeamElement::COMPUTE_RHS_VECTOR);

    MatrixType LeftHandSideMatrix = Matrix();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    //Calculate elemental system
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH( "" )

}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags   
    LocalSystem.CalculationFlags.Set(LargeDisplacementBeamElement::COMPUTE_LHS_MATRIX);

    VectorType RightHandSideVector = Vector();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector,  LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(RightHandSideVector);

    //Calculate elemental system
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH( "" )

}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags 
    LocalSystem.CalculationFlags.Set(LargeDisplacementBeamElement::COMPUTE_RHS_VECTOR);
    LocalSystem.CalculationFlags.Set(LargeDisplacementBeamElement::COMPUTE_LHS_MATRIX);

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    //Calculate elemental system
    CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH( "" )

}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
								 GeneralVariables& rVariables,
								 Vector& rVolumeForce,
								 double& rIntegrationWeight)
{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    double DomainSize = rVariables.Section.Area; 

    //gravity load
    Vector GravityLoad = ZeroVector(dimension);
    Vector Fe = ZeroVector(6);

    unsigned int RowIndex = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
      	RowIndex = i * (dimension * 2);

        GravityLoad = rIntegrationWeight * rVariables.N[i] * rVolumeForce * DomainSize; 	
	BeamMathUtilsType::AddVector(GravityLoad, rRightHandSideVector, RowIndex);
	
      }
    
    //---------------------

    // External Energy Calculation:
    Vector CurrentValueVector           = ZeroVector(3);
    Vector CurrentDisplacementVector    = ZeroVector(3);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	//Current Displacement Vector
	CurrentValueVector = GetNodalCurrentValue( DISPLACEMENT, CurrentValueVector, i );
	CurrentDisplacementVector += rVariables.N[i] * CurrentValueVector;
      }

    CurrentDisplacementVector     = MapToInitialLocalFrame( CurrentDisplacementVector, rVariables.PointNumber  );

    //for more than one integration point
    Vector PotencialForceVector = ZeroVector(3);
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
    	PotencialForceVector  += rIntegrationWeight * rVariables.N[i] * rVolumeForce * DomainSize;	  
      }
    
	
    mEnergy.External += inner_prod(PotencialForceVector,CurrentDisplacementVector);
 
    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateAndAddFollowerForces(VectorType& rRightHandSideVector,
								 GeneralVariables & rVariables,
								 double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //Initialize Local Matrices
    VectorType Fi = ZeroVector(6);
    Vector ResultantsVector = ZeroVector(6);

    //locate follower load resultants in skew-symmetric form
    Vector FollowerLoad   = ZeroVector(dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      if( GetGeometry()[i].SolutionStepsDataHas( FOLLOWER_LINE_LOAD ) ) 
    	FollowerLoad += rVariables.N[i] * GetGeometry()[i].FastGetSolutionStepValue( FOLLOWER_LINE_LOAD );
    }
    
    //Current Frame given by the frame rotation
    FollowerLoad = prod( rVariables.CurrentRotationMatrix, FollowerLoad );
    
    BeamMathUtilsType::AddVector(FollowerLoad, ResultantsVector, 3);

    //Initialize Local Matrices
    MatrixType OperatorI = ZeroMatrix(6,6);

    unsigned int RowIndex = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
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
								 GeneralVariables & rVariables,
								 double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = rRightHandSideVector.size();

    VectorType Fi = ZeroVector(6);

    //Initialize Local Matrices
    MatrixType DifferentialOperatorI = ZeroMatrix(MatSize,MatSize);

    unsigned int RowIndex = 0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {

	RowIndex = i * (dimension * 2);

	Fi = ZeroVector(6);
	
 	this->CalculateDifferentialOperator(DifferentialOperatorI, rVariables, i, rVariables.Alpha );

	//nodal force vector
	Fi = prod( DifferentialOperatorI, rVariables.StressVector );

	Fi *= rIntegrationWeight;
	
	BeamMathUtilsType::SubstractVector( Fi, rRightHandSideVector, RowIndex );

	// bool write = true;
	// if( write ){

	//   Vector CurrentStepRotationVector  = ZeroVector(3);
	//   this->GetLocalCurrentValue(STEP_ROTATION, CurrentStepRotationVector, rVariables.N);

	//   Vector CurrentDisplacementVector  = ZeroVector(3);
	//   this->GetLocalCurrentValue(DISPLACEMENT, CurrentDisplacementVector, rVariables.N);

	//   std::cout<<" STEP_ROTATION  ("<<this->Id()<<"): "<<CurrentStepRotationVector<<std::endl;
	//   std::cout<<" DISPLACEMENT   ("<<this->Id()<<"): "<<CurrentDisplacementVector<<std::endl;
	//   std::cout<<" INTERNAL_FORCE ("<<this->Id()<<"): "<<Fi<<" StressVector "<<rVariables.StressVector<<std::endl;
	  
	//   //std::cout<<" Btensor "<<DifferentialOperatorI<<std::endl;

	//   // std::cout<<" rRightHandSideVector "<<rRightHandSideVector<<std::endl;
	// }

	//---------------------

	// Internal Energy Calculation:
	Vector Movements;
	this->GetCurrentNodalMovements( Movements, i, rVariables.PointNumber);
	
	mEnergy.Internal += inner_prod( Fi, Movements );
      }
 

    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::AddExplicitContribution(const VectorType& rRHSVector, 
							   const Variable<VectorType>& rRHSVariable, 
							   Variable<array_1d<double,3> >& rDestinationVariable, 
							   const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
      
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if( (rRHSVariable == RESIDUAL_VECTOR) ){

      if ( rDestinationVariable == FORCE_RESIDUAL )
	{

	  for(unsigned int i=0; i< number_of_nodes; i++)
	    {
	      int index = (dimension * 2) * i;

	      GetGeometry()[i].SetLock();

	      array_1d<double, 3 > &ForceResidual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

	      for(unsigned int j=0; j<dimension; j++)
		{
		  ForceResidual[j] += rRHSVector[index + j];
		}

	      GetGeometry()[i].UnSetLock();
	    }
	}
      else if( rDestinationVariable == MOMENT_RESIDUAL )
	{

	  for(unsigned int i=0; i< number_of_nodes; i++)
	    {
	      int index = dimension + (dimension * 2) * i;

	      GetGeometry()[i].SetLock();

	      array_1d<double, 3 > &MomentResidual = GetGeometry()[i].FastGetSolutionStepValue(MOMENT_RESIDUAL);

	      for(unsigned int j=0; j<dimension; j++)
		{
		  MomentResidual[j] += rRHSVector[index + j];
		}

	      GetGeometry()[i].UnSetLock();
	    }

	}

      
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
    
    rOperator = ZeroMatrix(6,6);

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
								 GeneralVariables& rVariables,
								 const int& rNode,
								 double alpha)
{
    KRATOS_TRY

    //Differencial operator transposed

    //Initialize Local Matrices
    if( rDifferentialOperator.size1() != 6 )
      rDifferentialOperator.resize(6, 6, false);
    
    rDifferentialOperator = ZeroMatrix(6,6);

    rDifferentialOperator( 0, 0 ) =  rVariables.DN_DX( rNode, 0 );
    rDifferentialOperator( 1, 1 ) =  rVariables.DN_DX( rNode, 0 );
    rDifferentialOperator( 2, 2 ) =  rVariables.DN_DX( rNode, 0 );
    rDifferentialOperator( 3, 3 ) =  rVariables.DN_DX( rNode, 0 );
    rDifferentialOperator( 4, 4 ) =  rVariables.DN_DX( rNode, 0 );
    rDifferentialOperator( 5, 5 ) =  rVariables.DN_DX( rNode, 0 );


    //locate stress resultants in skew-symmetric "transposed" form 
    Matrix SkewSymResultants = ZeroMatrix(3,3);

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
							     GeneralVariables& rVariables,
							     const int& rNode)
{
    KRATOS_TRY

    //Initialize Local Matrices
    if( rDiscreteOperator.size1() != 6 || rDiscreteOperator.size2() != 9)
      rDiscreteOperator.resize(6, 9, false);
    
    rDiscreteOperator = ZeroMatrix(6,9);

    rDiscreteOperator( 0, 0 ) =  rVariables.DN_DX( rNode, 0 );
    rDiscreteOperator( 1, 1 ) =  rVariables.DN_DX( rNode, 0 );
    rDiscreteOperator( 2, 2 ) =  rVariables.DN_DX( rNode, 0 );
    rDiscreteOperator( 3, 3 ) =  rVariables.DN_DX( rNode, 0 );
    rDiscreteOperator( 4, 4 ) =  rVariables.DN_DX( rNode, 0 );
    rDiscreteOperator( 5, 5 ) =  rVariables.DN_DX( rNode, 0 );

    rDiscreteOperator( 3, 6 ) =  rVariables.N[ rNode ];
    rDiscreteOperator( 4, 7 ) =  rVariables.N[ rNode ];
    rDiscreteOperator( 5, 8 ) =  rVariables.N[ rNode ];


    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

//Geometric operator
void LargeDisplacementBeamElement::CalculateBmatrix(MatrixType& rBmatrix,
						    GeneralVariables& rVariables)
{
    KRATOS_TRY

    //Initialize Local Matrices
    if( rBmatrix.size1() != 9 )
      rBmatrix.resize(9, 9, false);
    
    rBmatrix = ZeroMatrix(9,9);

    Vector StressResultants = ZeroVector(3);
    Vector StressCouples    = ZeroVector(3);


    StressResultants[0] = rVariables.StressVector[0];
    StressResultants[1] = rVariables.StressVector[1];
    StressResultants[2] = rVariables.StressVector[2];

    StressCouples[0] = rVariables.StressVector[3];
    StressCouples[1] = rVariables.StressVector[4];
    StressCouples[2] = rVariables.StressVector[5];


    //locate stress resultants in skew-symmetric form
    Matrix SkewSymStressResultants = ZeroMatrix(3,3);

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
    Matrix SkewSymStressCouples = ZeroMatrix(3,3);
    
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


    Matrix AxialStressMatrix = ZeroMatrix(3,3);

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

    
    // std::cout<<" AxialStressMatrix "<<AxialStressMatrix<<std::endl;

    //axial skew-symmetric matrix (Crisfield)
    Matrix AxialSkewSymMatrix = ZeroMatrix(3,3);
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
						       GeneralVariables& rVariables,
						       double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = dimension * 2;
    
    //Initialize Local Matrices
    MatrixType DifferentialOperatorI = ZeroMatrix(MatSize,MatSize);
    MatrixType DifferentialOperatorJ = ZeroMatrix(MatSize,MatSize);

    MatrixType Kij = ZeroMatrix(MatSize,MatSize);

    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	RowIndex = i * (dimension * 2);

	this->CalculateDifferentialOperator( DifferentialOperatorI, rVariables, i, rVariables.Alpha );
	
	for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {

	    Kij = ZeroMatrix(6,6);

	    ColIndex = j * (dimension * 2);

	    this->CalculateDifferentialOperator( DifferentialOperatorJ, rVariables, j, 1 );

	    noalias(Kij) = prod( DifferentialOperatorI, Matrix(prod( rVariables.ConstitutiveMatrix, trans(DifferentialOperatorJ) )) );

	    Kij *= rIntegrationWeight;

	    //std::cout<<" Kij "<<Kij<<std::endl;

	    //Building the Local Stiffness Matrix
	    MathUtils<double>::AddMatrix( rLeftHandSideMatrix, Kij, RowIndex, ColIndex );
	  }
      }

    //std::cout<<" Kuum "<<rLeftHandSideMatrix<<std::endl;

   
    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
						       GeneralVariables& rVariables,
						       double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();


    //MatrixType Kuum = rLeftHandSideMatrix;

    //Initialize Local Matrices
    MatrixType Bmatrix = ZeroMatrix(9, 9);
    this->CalculateBmatrix( Bmatrix, rVariables );

    MatrixType DiscreteOperatorI = ZeroMatrix(6,9);
    MatrixType DiscreteOperatorJ = ZeroMatrix(6,9);

    MatrixType Kij = ZeroMatrix(6,6);

    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	
	RowIndex = i * (dimension * 2);

 	this->CalculateDiscreteOperator( DiscreteOperatorI, rVariables, i );

	for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {

	    Kij = ZeroMatrix(6,6);

	    ColIndex = j * (dimension * 2);

	    this->CalculateDiscreteOperator( DiscreteOperatorJ, rVariables, j );

	    noalias(Kij) = prod( DiscreteOperatorI, Matrix( prod( Bmatrix, trans(DiscreteOperatorJ) ) ) );
	    
	    Kij *= rIntegrationWeight;

	    //Building the Local Stiffness Matrix
	    MathUtils<double>::AddMatrix( rLeftHandSideMatrix, Kij, RowIndex, ColIndex );
	    
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

    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************


void LargeDisplacementBeamElement::CalculateAndAddKuug2(MatrixType& rLeftHandSideMatrix,
							GeneralVariables& rVariables,
							double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();


    // MatrixType Kuum = rLeftHandSideMatrix;

    //Initialize Local Matrices
    MatrixType Kij        = ZeroMatrix(6,6);
    Matrix GabK           = ZeroMatrix(3,3);
    Matrix DiagonalMatrix = IdentityMatrix(3);   

    Vector StressResultants = ZeroVector(3);
    for ( unsigned int i = 0; i < 3; i++ )
      {
	StressResultants[i] = rVariables.StressVector[i];
      }

    Vector StressCouples = ZeroVector(3);
    for ( unsigned int i = 0; i < 3; i++ )
      {
	StressCouples[i] = rVariables.StressVector[i+3];
      }

    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	
	RowIndex = i * (dimension * 2);

	for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {

	    Kij = ZeroMatrix(6,6);

	    ColIndex = j * (dimension * 2);

	    
	    //term 11 -> 0
	    //term 12
	    GabK = ZeroMatrix(3,3);
	    Matrix SkewSymStressResultants = ZeroMatrix(3,3);
	    BeamMathUtilsType::VectorToSkewSymmetricTensor(StressResultants, SkewSymStressResultants);
	    GabK = (-1) * (rVariables.DN_DX(i, 0) * rVariables.N[j]) * SkewSymStressResultants; 
	    //Building the Local Stiffness Matrix
	    MathUtils<double>::AddMatrix( Kij, GabK, 0, 3 );
		
	    //term 21
	    GabK = ZeroMatrix(3,3);
	    GabK = (rVariables.N[i] * rVariables.DN_DX(j, 0) ) * SkewSymStressResultants;
	    //Building the Local Stiffness Matrix
	    MathUtils<double>::AddMatrix( Kij, GabK, 3, 0 );
	

	    //term 22
	    GabK = ZeroMatrix(3,3);

	    Matrix SkewSymStressCouples = ZeroMatrix(3,3);
	    BeamMathUtilsType::VectorToSkewSymmetricTensor(StressCouples, SkewSymStressCouples);

	    GabK  = (-1) * (rVariables.DN_DX(i, 0) * rVariables.N[j]) * SkewSymStressCouples; 
	    GabK += ( rVariables.N[i] * rVariables.N[j]) * outer_prod( StressResultants, rVariables.CurrentAxisPositionDerivatives );
		
	    GabK -= ( rVariables.N[i] * rVariables.N[j]) * inner_prod( StressResultants, rVariables.CurrentAxisPositionDerivatives ) * DiagonalMatrix;

	    //Building the Local Stiffness Matrix
	    MathUtils<double>::AddMatrix( Kij, GabK, 3, 3 );

	    Kij *= rIntegrationWeight;

	    //Building the Local Stiffness Matrix
	    MathUtils<double>::AddMatrix( rLeftHandSideMatrix, Kij, RowIndex, ColIndex );
	    
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

    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateAndAddKuuf(MatrixType& rLeftHandSideMatrix,
						       GeneralVariables& rVariables,
						       double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //Initialize Local Matrices
    MatrixType Kij = ZeroMatrix(6,6);
    MatrixType Lij = ZeroMatrix(3,3);

    //locate follower load resultants in skew-symmetric form
    Vector FollowerLoad   = ZeroVector(dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      if( GetGeometry()[i].SolutionStepsDataHas( FOLLOWER_LINE_LOAD ) ) 
	FollowerLoad += rVariables.N[i] *  GetGeometry()[i].FastGetSolutionStepValue( FOLLOWER_LINE_LOAD );
    }

    //Current Frame given by the frame rotation
    FollowerLoad = prod( rVariables.CurrentRotationMatrix, FollowerLoad );

    Matrix SkewSymResultants = ZeroMatrix(dimension,dimension);

    BeamMathUtilsType::VectorToSkewSymmetricTensor(FollowerLoad, SkewSymResultants);

    SkewSymResultants *= (-1);

    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {

	Kij = ZeroMatrix(6,6);

	RowIndex = i * (dimension * 2);


	for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {

	    ColIndex = j * (dimension * 2);

	    Lij =  rIntegrationWeight * rVariables.N[i] * rVariables.N[j] * SkewSymResultants;
	    
	    //Building the Local Stiffness Matrix
	    MathUtils<double>::AddMatrix( Kij, Lij, 0, 3 );
	    	    
	    //Building the Local Stiffness Matrix
	    MathUtils<double>::AddMatrix( rLeftHandSideMatrix, Kij, RowIndex, ColIndex );
	    
	  }
      }

   
    KRATOS_CATCH( "" )

}


//****************GID DEFINITION OF THE AUTOMATIC LOCAL AXES******************
//*****************************************************************************

//Local E3 beam axis
void LargeDisplacementBeamElement::CalculateLocalAxesMatrix(Matrix& rRotationMatrix)
{

    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       size            = number_of_nodes * dimension;

    Vector LocalZ               = ZeroVector(3);
    Vector ReferenceCoordinates = ZeroVector(size);

    ReferenceCoordinates[0] = GetGeometry()[0].X();
    ReferenceCoordinates[1] = GetGeometry()[0].Y();
    ReferenceCoordinates[2] = GetGeometry()[0].Z();     

    int k = number_of_nodes - 1 ;

    ReferenceCoordinates[3] = GetGeometry()[k].X();
    ReferenceCoordinates[4] = GetGeometry()[k].Y();
    ReferenceCoordinates[5] = GetGeometry()[k].Z();     
   
    for( unsigned int i = 0; i < dimension; i++ )
    {
      LocalZ[i]  = (ReferenceCoordinates[i+3] - ReferenceCoordinates[i]);
    }

    BeamMathUtilsType::CalculateLocalAxesMatrix(LocalZ,rRotationMatrix);

    KRATOS_CATCH( "" )

}


//************************************CALCULATE VOLUME ACCELERATION*******************
//************************************************************************************

Vector&  LargeDisplacementBeamElement::CalculateVolumeForce( Vector& rVolumeForce, const Vector &rN)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rVolumeForce = ZeroVector(dimension);
    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
      //temporary, will be checked once at the beginning only
      if( GetGeometry()[j].SolutionStepsDataHas(VOLUME_ACCELERATION) ){
	rVolumeForce += rN[j] * GetGeometry()[j].FastGetSolutionStepValue(VOLUME_ACCELERATION);
      }
    }

    rVolumeForce *= GetProperties()[DENSITY]; 

    //Current Frame is the local frame
    rVolumeForce = MapToInitialLocalFrame( rVolumeForce );

    return rVolumeForce;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY


    bool ComputeDynamicTangent = false;
    if( rCurrentProcessInfo.Has(COMPUTE_DYNAMIC_TANGENT) ){

      if(rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true){
	ComputeDynamicTangent = true;
      }
    }

    if( ComputeDynamicTangent == true ){

      //create local system components
      LocalSystemComponents LocalSystem;

      //calculation flags 
      LocalSystem.CalculationFlags.Set(LargeDisplacementBeamElement::COMPUTE_RHS_VECTOR);
      LocalSystem.CalculationFlags.Set(LargeDisplacementBeamElement::COMPUTE_LHS_MATRIX);
    
      //Initialize sizes for the system components:
      this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

      //Set Variables to Local system components
      LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
      LocalSystem.SetRightHandSideVector(rRightHandSideVector);
	
      //Calculate elemental system
      CalculateDynamicSystem( LocalSystem, rCurrentProcessInfo );

    }
    else{

      //1.-Calculate Tangent Inertia Matrix:
      this->CalculateMassMatrix( rLeftHandSideMatrix, rCurrentProcessInfo );

      double MatSize = rLeftHandSideMatrix.size1();

      //2.-Calculate Inertial Forces:
      if ( rRightHandSideVector.size() != MatSize )
	rRightHandSideVector.resize( MatSize, false );
      
      rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
    }


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
      
    bool ComputeDynamicTangent = false;
    if( rCurrentProcessInfo.Has(COMPUTE_DYNAMIC_TANGENT) )
      if(rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true)
	ComputeDynamicTangent = true;
    
    if( ComputeDynamicTangent == true ){

      //create local system components
      LocalSystemComponents LocalSystem;

      //calculation flags   
      LocalSystem.CalculationFlags.Set(LargeDisplacementBeamElement::COMPUTE_LHS_MATRIX);

      VectorType RightHandSideVector = Vector();

      //Initialize sizes for the system components:
      this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector,  LocalSystem.CalculationFlags );
	
      //Set Variables to Local system components
      LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
      LocalSystem.SetRightHandSideVector(RightHandSideVector);
	
      //Calculate elemental system
      CalculateDynamicSystem( LocalSystem, rCurrentProcessInfo );    

    }
    else{

      //1.-Calculate Tangent Inertia Matrix:
      this->CalculateMassMatrix( rLeftHandSideMatrix, rCurrentProcessInfo );

    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    bool ComputeDynamicTangent = false;
    if( rCurrentProcessInfo.Has(COMPUTE_DYNAMIC_TANGENT) )
      if(rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true)
	ComputeDynamicTangent = true;
    
    if( ComputeDynamicTangent == true ){

      //create local system components
      LocalSystemComponents LocalSystem;

      //calculation flags
      LocalSystem.CalculationFlags.Set(LargeDisplacementBeamElement::COMPUTE_RHS_VECTOR);

      MatrixType LeftHandSideMatrix = Matrix();

      //Initialize sizes for the system components:
      this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

      //Set Variables to Local system components
      LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
      LocalSystem.SetRightHandSideVector(rRightHandSideVector);

      //Calculate elemental system
      CalculateDynamicSystem( LocalSystem, rCurrentProcessInfo );
 
    }
    else{

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
      unsigned int MatSize               = number_of_nodes * ( dimension * 2 );

      //2.-Calculate Inertial Forces:
      if ( rRightHandSideVector.size() != MatSize )
	rRightHandSideVector.resize( MatSize, false );
      
      rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
    }


    KRATOS_CATCH( "" )   
}

//************************************************************************************
//************************************************************************************

//Inertia in the SPATIAL configuration 
void LargeDisplacementBeamElement::CalculateAndAddInertiaLHS(MatrixType& rLeftHandSideMatrix, GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight )
{

    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize               = number_of_nodes * ( dimension * 2 );

    if(rLeftHandSideMatrix.size1() != MatSize)
        rLeftHandSideMatrix.resize (MatSize, MatSize, false);

    rLeftHandSideMatrix = ZeroMatrix( MatSize, MatSize );

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

    MatrixType m11 = ZeroMatrix(3,3);

    double TotalMass = 0;
    TotalMass  = this->CalculateTotalMass( Section, TotalMass );
    
    //block m(2,2) of the mass matrix
     
    MatrixType m22 = ZeroMatrix(3,3);

    Vector CurrentCompoundRotationVector   = ZeroVector(3);
    Vector ReferenceCompoundRotationVector = ZeroVector(3);

    Vector CurrentStepRotationVector = ZeroVector(3);
    Vector PreviousStepRotationVector = ZeroVector(3);

    Vector AngularVelocityVector     = ZeroVector(3);
    Vector AngularAccelerationVector = ZeroVector(3);
    Vector CurrentAngularAccelerationVector = ZeroVector(3);
    Vector PreviousAngularAccelerationVector = ZeroVector(3);

    Vector CurrentValueVector   = ZeroVector(3);
    Vector PreviousValueVector = ZeroVector(3);

    std::vector<QuaternionType> ReferenceNodeQuaternions;
    QuaternionType QuaternionValue;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {

	//Current Compound Rotation Vector 
	CurrentValueVector = GetNodalCurrentValue( ROTATION, CurrentValueVector, i );
	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector );

	CurrentCompoundRotationVector   += rVariables.N[i] * CurrentValueVector;
	
	//Reference Compound Rotation Vector 
	PreviousValueVector = GetNodalPreviousValue( ROTATION, PreviousValueVector, i );
	PreviousValueVector = MapToInitialLocalFrame( PreviousValueVector );

	ReferenceCompoundRotationVector += rVariables.N[i] * PreviousValueVector;	

	//Reference Frame is the Local Frame
	QuaternionValue = QuaternionType::FromRotationVector(PreviousValueVector);
	
	ReferenceNodeQuaternions.push_back(QuaternionValue);


	//Current Step Rotation Vector 
	CurrentValueVector = GetNodalCurrentValue( STEP_ROTATION, CurrentValueVector, i );
	CurrentStepRotationVector       += rVariables.N[i] * CurrentValueVector;

	//Previous Step Rotation Vector 
	PreviousValueVector = GetNodalPreviousValue( STEP_ROTATION, PreviousValueVector, i );
	PreviousStepRotationVector      += rVariables.N[i] * PreviousValueVector;

	//Angular Velocity Vector
	CurrentValueVector = GetNodalCurrentValue( ANGULAR_VELOCITY, CurrentValueVector, i );
	AngularVelocityVector           += rVariables.N[i] * CurrentValueVector;

	//Current Angular Acceleration Vector
	CurrentValueVector = GetNodalCurrentValue( ANGULAR_ACCELERATION, CurrentValueVector, i );
	CurrentAngularAccelerationVector       += rVariables.N[i] * CurrentValueVector;
        
        //Previous Angular Acceleration Vector
	CurrentValueVector = GetNodalPreviousValue( ANGULAR_ACCELERATION, CurrentValueVector, i );
	PreviousAngularAccelerationVector       += rVariables.N[i] * CurrentValueVector;
	
      }
     

    //Set step variables to local frame (current Frame is the local frame)
    PreviousStepRotationVector        = MapToInitialLocalFrame( PreviousStepRotationVector );
    CurrentStepRotationVector         = MapToInitialLocalFrame( CurrentStepRotationVector );
    AngularVelocityVector             = MapToInitialLocalFrame( AngularVelocityVector );
    CurrentAngularAccelerationVector  = MapToInitialLocalFrame( CurrentAngularAccelerationVector );
    PreviousAngularAccelerationVector = MapToInitialLocalFrame( PreviousAngularAccelerationVector );
    
    double AlphaM = 0;
    if( rCurrentProcessInfo.Has(BOSSAK_ALPHA) ){
      AlphaM = rCurrentProcessInfo[BOSSAK_ALPHA];
    }

    AngularAccelerationVector = (1.0-AlphaM)*CurrentAngularAccelerationVector + AlphaM*(PreviousAngularAccelerationVector);
      
    //1.-Get inertia dyadic
    Matrix InertiaDyadic = ZeroMatrix(3,3);
    this->CalculateInertiaDyadic( Section, InertiaDyadic );
    InertiaDyadic = prod(rVariables.CurrentRotationMatrix,InertiaDyadic);
    InertiaDyadic = prod(InertiaDyadic,trans(rVariables.CurrentRotationMatrix));

    //2.- Compute Term 1:
    
    Matrix MassTerm1 = ZeroMatrix(3,3);
    
    Vector InertiaxAngularVelocity     = prod( InertiaDyadic, AngularVelocityVector );
    Vector InertiaxAngularAcceleration = prod( InertiaDyadic, AngularAccelerationVector );
    
    Matrix TensorAngularVelocity = ZeroMatrix(3,3);

    BeamMathUtilsType::VectorToSkewSymmetricTensor( AngularVelocityVector, TensorAngularVelocity );

    Vector VectorTerm1 = ZeroVector(3);

    VectorTerm1  = prod( TensorAngularVelocity, InertiaxAngularVelocity );

    VectorTerm1 += InertiaxAngularAcceleration;
    
    BeamMathUtilsType::VectorToSkewSymmetricTensor(VectorTerm1, MassTerm1);
    
    
    //3.- Compute Term 2:
    
    Matrix MassTerm2 = ZeroMatrix(3,3);

    Matrix InertiaxAngularVelocityTensor = ZeroMatrix(3,3);
    BeamMathUtilsType::VectorToSkewSymmetricTensor( InertiaxAngularVelocity, InertiaxAngularVelocityTensor );

    Matrix TensorAngularVelocityxInertia = ZeroMatrix(3,3);
    noalias(TensorAngularVelocityxInertia) = prod( TensorAngularVelocity, InertiaDyadic );
    
    Matrix InertiaxTensorAngularVelocity = ZeroMatrix(3,3);
    noalias(InertiaxTensorAngularVelocity) = prod( InertiaDyadic, TensorAngularVelocity );
    
    //MassTerm2 = (1.0-AlphaM) * InertiaDyadic - Newmark2 * InertiaxAngularVelocityTensor + Newmark2 * TensorAngularVelocityxInertia + Newmark2 * InertiaxTensorAngularVelocity; //Cardona

    MassTerm2 = (1.0-AlphaM) * InertiaDyadic - Newmark2 * InertiaxAngularVelocityTensor + Newmark2 * TensorAngularVelocityxInertia; //Simo       

    MassTerm2 *= Newmark1;
    
    MatrixType MassMatrixBlock2 = ZeroMatrix(3,3);

    MassMatrixBlock2 = (-1) * MassTerm1 + MassTerm2;
    	    
    // Compute Linear Part of the Step Rotation
    Matrix LinearPartRotationTensor = ZeroMatrix(3,3);

    Vector StepRotationVector = PreviousStepRotationVector; //CurrentStepRotationVector;

    double NormStepRotation =  norm_2(StepRotationVector);
    if( NormStepRotation != 0 ){
      this->CalculateRotationLinearPartTensor( StepRotationVector , LinearPartRotationTensor );
    }
    else{
      LinearPartRotationTensor = IdentityMatrix(3);
    }

    
    MassMatrixBlock2 = prod( MassMatrixBlock2, LinearPartRotationTensor );
    
    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    Matrix DiagonalMatrix = IdentityMatrix(3);   
    
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	m11 = ZeroMatrix(3,3);
	m22 = ZeroMatrix(3,3);

	RowIndex = i * (dimension * 2);


	for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {

	    ColIndex = j * (dimension * 2);

	    m11 = (1.0-AlphaM) * Newmark1 * TotalMass * rVariables.N[i] * rVariables.N[j] * rIntegrationWeight * DiagonalMatrix;

	    m22 = rVariables.N[i] * rVariables.N[j] * rIntegrationWeight * MassMatrixBlock2;
	    

	    //Building the Local Tangent Inertia Matrix
	    MathUtils<double>::AddMatrix( rLeftHandSideMatrix, m11, RowIndex, ColIndex );
	    MathUtils<double>::AddMatrix( rLeftHandSideMatrix, m22, RowIndex+3, ColIndex+3 );
	    
	  }
      }
   
    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

//Inertia in the SPATIAL configuration 
void LargeDisplacementBeamElement::CalculateAndAddInertiaRHS(VectorType& rRightHandSideVector, GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize               = number_of_nodes * ( dimension * 2 );

    if(rRightHandSideVector.size() != MatSize)
        rRightHandSideVector.resize(MatSize, false);

    rRightHandSideVector = ZeroVector( MatSize );

    SectionProperties Section;
    this->CalculateSectionProperties(Section);

  
    double TotalMass = 0;
    TotalMass = this->CalculateTotalMass( Section, TotalMass );     
 
    //displacements and rotations vector

    Vector LinearAccelerationVector   = ZeroVector(3);
    Vector CurrentLinearAccelerationVector   = ZeroVector(3);
    Vector PreviousLinearAccelerationVector  = ZeroVector(3);

    Vector AngularVelocityVector     = ZeroVector(3);
    Vector AngularAccelerationVector = ZeroVector(3);
    Vector CurrentAngularAccelerationVector = ZeroVector(3);
    Vector PreviousAngularAccelerationVector = ZeroVector(3);
 

    Vector CurrentValueVector = ZeroVector(3);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
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
    CurrentLinearAccelerationVector   = MapToInitialLocalFrame( CurrentLinearAccelerationVector );
    PreviousLinearAccelerationVector  = MapToInitialLocalFrame( PreviousLinearAccelerationVector );

    AngularVelocityVector             = MapToInitialLocalFrame( AngularVelocityVector );

    CurrentAngularAccelerationVector  = MapToInitialLocalFrame( CurrentAngularAccelerationVector );
    PreviousAngularAccelerationVector = MapToInitialLocalFrame( PreviousAngularAccelerationVector );

    double AlphaM = 0;
    if( rCurrentProcessInfo.Has(BOSSAK_ALPHA) ){
      AlphaM = rCurrentProcessInfo[BOSSAK_ALPHA];
    }

    LinearAccelerationVector  = (1.0-AlphaM) * CurrentLinearAccelerationVector  + AlphaM * (PreviousLinearAccelerationVector);
    AngularAccelerationVector = (1.0-AlphaM) * CurrentAngularAccelerationVector + AlphaM * (PreviousAngularAccelerationVector);

    //-----------------
    //block m(1) of the inertial force vector

    //Compute Linear Term
    Vector LinearInertialForceVector = ZeroVector(3);
    LinearInertialForceVector = TotalMass * LinearAccelerationVector;


    //-----------------
    //block m(2,2) of the inertial force vector (rotations part::to be defined)

    //Get inertia dyadic
    Matrix InertiaDyadic = ZeroMatrix(3,3);
    this->CalculateInertiaDyadic( Section, InertiaDyadic );
    InertiaDyadic = prod(rVariables.CurrentRotationMatrix,InertiaDyadic);
    InertiaDyadic = prod(InertiaDyadic,trans(rVariables.CurrentRotationMatrix));


    //Compute Angular Term:  
    Vector InertiaxAngularVelocity     = prod( InertiaDyadic, AngularVelocityVector );
    Vector InertiaxAngularAcceleration = prod( InertiaDyadic, AngularAccelerationVector );
    
    Matrix TensorAngularVelocity = ZeroMatrix(3,3);
    BeamMathUtilsType::VectorToSkewSymmetricTensor( AngularVelocityVector, TensorAngularVelocity );

    Vector AngularInertialForceVector = ZeroVector(3);
    AngularInertialForceVector  = prod( TensorAngularVelocity, InertiaxAngularVelocity );

    // CROSS PRODUCT of AxB = prod( skewA, B )  where (skewA) = [Ax] = [A]^v =>  (skewA)^T = hat(A) (nomenclature)
    AngularInertialForceVector += InertiaxAngularAcceleration;
 

    //compose total acceleration integral function:
    Vector TotalInertialForceVector = ZeroVector(6);
    
    BeamMathUtilsType::AddVector(LinearInertialForceVector, TotalInertialForceVector, 0);
    BeamMathUtilsType::AddVector(AngularInertialForceVector, TotalInertialForceVector, 3);

    //Initialize Local Matrices
    VectorType Fi = ZeroVector(6);
    MatrixType OperatorI = ZeroMatrix(6,6);
    unsigned int RowIndex = 0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	RowIndex = i * (dimension * 2);

	Fi = ZeroVector(6);
	
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
void LargeDisplacementBeamElement::CalculateKineticEnergy(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();

    SectionProperties Section;
    this->CalculateSectionProperties(Section);

    //Get total mass
    double TotalMass = 0;
    TotalMass = this->CalculateTotalMass( Section, TotalMass );     
 
    //Get inertia dyadic
    Matrix InertiaDyadic = ZeroMatrix(3,3);
    this->CalculateInertiaDyadic( Section, InertiaDyadic );

    Matrix CurrentInertiaDyadic = prod(rVariables.CurrentRotationMatrix,InertiaDyadic);
    CurrentInertiaDyadic = prod(CurrentInertiaDyadic,trans(rVariables.CurrentRotationMatrix));

    // Kinetic Energy Calculation:
    Vector CurrentNodalVelocities    = ZeroVector(6);
    Vector CurrentVelocitiesVector   = ZeroVector(6);

    Matrix KineticMatrix;
    this->GetKineticMatrix( KineticMatrix, TotalMass, CurrentInertiaDyadic );    


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {

	CurrentNodalVelocities    = ZeroVector(6);
	this->GetCurrentNodalVelocities( CurrentNodalVelocities, i, rVariables.PointNumber );
	CurrentVelocitiesVector += rVariables.N[i] * CurrentNodalVelocities;

      }

 
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	mEnergy.Kinetic += 0.5 * rVariables.N[i] * inner_prod( CurrentVelocitiesVector, prod( KineticMatrix, CurrentVelocitiesVector) ) * rIntegrationWeight;
      }
  
      
    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

//Linear Momentum and Angular Momentum Calculation 
void LargeDisplacementBeamElement::CalculateMomentumRelations(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();

    SectionProperties Section;
    this->CalculateSectionProperties(Section);

    //Get total mass
    double TotalMass = 0;
    TotalMass = this->CalculateTotalMass( Section, TotalMass );     

    //Get inertia dyadic
    Matrix InertiaDyadic = ZeroMatrix(3,3);
    this->CalculateInertiaDyadic( Section, InertiaDyadic );

    //spatial
    Matrix CurrentInertiaDyadic = prod(rVariables.CurrentRotationMatrix,InertiaDyadic);
    CurrentInertiaDyadic = prod(CurrentInertiaDyadic,trans(rVariables.CurrentRotationMatrix));

    Matrix DiagonalMatrix = identity_matrix<double> (3);   

    Vector CurrentValueVector           = ZeroVector(3);
    Vector CurrentPositionVector        = ZeroVector(3);
    Vector CurrentDisplacementVector    = ZeroVector(3);
    Vector CurrentLinearVelocityVector  = ZeroVector(3);
    Vector CurrentAngularVelocityVector = ZeroVector(3);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
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

    CurrentPositionVector         = MapToInitialLocalFrame( CurrentPositionVector, rVariables.PointNumber );
    CurrentDisplacementVector     = MapToInitialLocalFrame( CurrentDisplacementVector, rVariables.PointNumber  );
    CurrentLinearVelocityVector   = MapToInitialLocalFrame( CurrentLinearVelocityVector, rVariables.PointNumber  );  
    CurrentAngularVelocityVector  = MapToInitialLocalFrame( CurrentAngularVelocityVector, rVariables.PointNumber );



    // Calculate Linear Momentum and Angular Momentum
    Vector LinearMomentumVector  = ZeroVector(3);
    Vector AngularMomentumVector = ZeroVector(3);


    //for more than one integration point
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
    	LinearMomentumVector  += TotalMass * rVariables.N[i] * prod( DiagonalMatrix, CurrentLinearVelocityVector ) * rIntegrationWeight;	  
      }
    
    AngularMomentumVector  =  MathUtils<double>::CrossProduct(CurrentPositionVector, LinearMomentumVector);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
    	AngularMomentumVector += prod( CurrentInertiaDyadic, CurrentAngularVelocityVector ) * rVariables.N[i] * rIntegrationWeight;
      }

    //for only one integration point
    //LinearMomentumVector   = TotalMass * prod( DiagonalMatrix, CurrentLinearVelocityVector ) * rIntegrationWeight;
    //AngularMomentumVector  = MathUtils<double>::CrossProduct(CurrentPositionVector, LinearMomentumVector);
    //AngularMomentumVector += prod( CurrentInertiaDyadic, CurrentAngularVelocityVector ) * rIntegrationWeight;



    LinearMomentumVector  = prod(LinearMomentumVector,  trans(rVariables.LocalTransformationMatrix) );
    AngularMomentumVector = prod(AngularMomentumVector, trans(rVariables.LocalTransformationMatrix) );

    for( unsigned int i=0; i<3; i++ )
      {
	mEnergy.LinearMomentum[i]  += LinearMomentumVector[i];
	mEnergy.AngularMomentum[i] += AngularMomentumVector[i];
      }


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::GetKineticMatrix(Matrix& rKineticMatrix, const double& rMass, const Matrix& rInertia)
{
  if( rKineticMatrix.size1() != 6 || rKineticMatrix.size2() != 6 )
    rKineticMatrix.resize(6,6);

  rKineticMatrix = ZeroMatrix(6,6);

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

}


//************************************************************************************
//************************************************************************************
void LargeDisplacementBeamElement::GetCurrentNodalVelocities(Vector& rValues, const int& rNode, unsigned int PointNumber)
{
  if( rValues.size() != 6 )
    rValues.resize(6);

  rValues = ZeroVector(6);

  Vector CurrentValueVector = ZeroVector(3);
  CurrentValueVector = GetNodalCurrentValue( VELOCITY, CurrentValueVector, rNode );
  CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, PointNumber );

  rValues[0] = CurrentValueVector[0];
  rValues[1] = CurrentValueVector[1];
  rValues[2] = CurrentValueVector[2];
	

  CurrentValueVector = ZeroVector(3);
  CurrentValueVector = GetNodalCurrentValue( ANGULAR_VELOCITY, CurrentValueVector, rNode );
  CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, PointNumber );

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

  rValues = ZeroVector(6);
  
  Vector CurrentValueVector = ZeroVector(3);
  CurrentValueVector = GetNodalCurrentValue( DISPLACEMENT, CurrentValueVector, rNode );
  CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, PointNumber );

  rValues[0] = CurrentValueVector[0];
  rValues[1] = CurrentValueVector[1];
  rValues[2] = CurrentValueVector[2];     

  CurrentValueVector = ZeroVector(3);
  CurrentValueVector = GetNodalCurrentValue( ROTATION, CurrentValueVector, rNode );
  CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, PointNumber );

  rValues[3] = CurrentValueVector[0];
  rValues[4] = CurrentValueVector[1];
  rValues[5] = CurrentValueVector[2];

}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize               = number_of_nodes * ( dimension * 2 );

    if(rMassMatrix.size1() != MatSize)
        rMassMatrix.resize (MatSize, MatSize, false);

    rMassMatrix = ZeroMatrix( MatSize, MatSize );

    SectionProperties Section;
    this->CalculateSectionProperties(Section);

    //block m(1,1) of the mass matrix

    MatrixType m11 = ZeroMatrix(3,3);

    double TotalMass = 0;
    TotalMass = this->CalculateTotalMass( Section, TotalMass );
    TotalMass *= GetGeometry().Length();

    Vector LumpFact = ZeroVector(number_of_nodes); 

    LumpFact = GetGeometry().LumpingFactors(LumpFact);


    //block m(2,2) of the mass matrix
    
    MatrixType m22 = ZeroMatrix(3,3);

    Matrix InertiaDyadic = ZeroMatrix(3,3);
    this->CalculateInertiaDyadic( Section, InertiaDyadic );


    for( unsigned int i=0; i < number_of_nodes; i++ )
      {

        double temp = LumpFact[i] * TotalMass;

	int RowIndex = i * (dimension * 2);

	for( unsigned int k=0; k < dimension; k++ )
	  {
	    m11(k,k) = temp;
	  }

	m22 = InertiaDyadic * temp;

	
	//Building the Local Tangent Inertia Matrix
	MathUtils<double>::AddMatrix( rMassMatrix, m11, RowIndex, RowIndex );
	MathUtils<double>::AddMatrix( rMassMatrix, m22, RowIndex+3, RowIndex+3 );
	
      }
 
    // Note: 
    // Variables.LocalTransformationMatrix is a rotation matrix with new base in columns
    // That means that the standard rotation K = Q·K'·QT and F = Q·F' is the correct transformation

    // initialize local transformation/rotation matrix
    Matrix LocalTransformationMatrix = ZeroMatrix(3,3);

    mInitialLocalQuaternionsReduced[0].ToRotationMatrix(LocalTransformationMatrix);

    // Transform Local to Global LHSMatrix:
    BeamMathUtilsType::MapLocalToGlobal3D(LocalTransformationMatrix,rMassMatrix);


    KRATOS_CATCH( "" )

}


//************************************CALCULATE TOTAL MASS****************************
//************************************************************************************

double& LargeDisplacementBeamElement::CalculateTotalMass( SectionProperties& Section, double& rTotalMass )
{
    KRATOS_TRY

    //const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

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
    
    rInertiaDyadic = ZeroMatrix(3,3);
       
    //if the local axes are the principal axes of the cross section

    //Axis Local E3
    rInertiaDyadic(0,0) = rSection.Inertia_y; //vertial axis
    rInertiaDyadic(1,1) = rSection.Inertia_z; //horizontal axis
    rInertiaDyadic(2,2) = rSection.Rotational_Inertia; //beam axis

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
    
    rRotationTensor = ZeroMatrix(3,3);

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


    if ( rVariable == KINETIC_ENERGY )
    {
      //reading integration points
      for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
        {
	  rOutput[PointNumber] = mEnergy.Kinetic;
	}
    }

    if ( rVariable == INTERNAL_ENERGY )
    {
      //reading integration points
      for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
        {
	  //rOutput[PointNumber] = mEnergy.Internal;
	  rOutput[PointNumber] = mEnergy.Deformation;
	}
    }

    if ( rVariable == EXTERNAL_ENERGY )
    {
      //reading integration points
      for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
        {
	  rOutput[PointNumber] = mEnergy.External;
	}
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


      if ( rVariable == LINEAR_MOMENTUM )
	{
	  //reading integration points
	  for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
	    {
	      rOutput[PointNumber] = mEnergy.LinearMomentum;
	    }

	}
      
      if ( rVariable == ANGULAR_MOMENTUM )
	{
	  //reading integration points
	  for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
	    {
	      rOutput[PointNumber] = mEnergy.AngularMomentum;
	    }

	}

      
      if(rVariable==MOMENT || rVariable==FORCE){ 

	//create and initialize element variables:
	GeneralVariables Variables;
	this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

	//reading integration points (in fact is the two nodes beam element, only one integration point)
	const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

	const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

	//(in fact is the two nodes beam element, only one integration point)
	for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
	  {
	  
	    //compute element kinematics  ...
	    //std::cout<<" Integration Points "<<std::endl;
	    this->CalculateKinematics(Variables,PointNumber);
	  
	    //compute element ConstitutiveTensor
	    this->CalculateConstitutiveMatrix(Variables);
	  
	    //compute element Strain and Stress Resultants and Couples
	    this->CalculateStressResultants(Variables, PointNumber, Variables.Alpha);
     
	    if(rVariable==MOMENT)
	      {
		Vector StressCouples    = ZeroVector(dimension);

		for ( unsigned int i = 0; i < dimension; i++ )
		  {
		    StressCouples[i]  = Variables.StressVector[i+3];
		  }

		StressCouples = BeamMathUtilsType::MapToReferenceLocalFrame( mInitialLocalQuaternionsReduced[PointNumber], StressCouples );

		for( unsigned int i=0; i<dimension; i++ )
		  {
		    rOutput[PointNumber][i] = StressCouples[i]; 
		  }
	      }

	    if(rVariable==FORCE)
	      {

		Vector StressResultants = ZeroVector(dimension);
		
		for ( unsigned int i = 0; i < dimension; i++ )
		  {
		    StressResultants[i] = Variables.StressVector[i];
		  }

		StressResultants = BeamMathUtilsType::MapToReferenceLocalFrame( mInitialLocalQuaternionsReduced[PointNumber], StressResultants );

		for( unsigned int i=0; i<dimension; i++ )
		  {
		    rOutput[PointNumber][i] = StressResultants[i]; 
		  }
	      }
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
int  LargeDisplacementBeamElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (GetGeometry().WorkingSpaceDimension() != 3 || GetGeometry().size()!=2 )
    {
      KRATOS_THROW_ERROR( std::invalid_argument, "This element works only in 3D and with 2 noded linear elements", "")
    }

    //verify that the variables are correctly initialized
    if(VELOCITY.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"VELOCITY has Key zero! (check if the application is correctly registered", "" )
    if(DISPLACEMENT.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered", "" )
    if(ACCELERATION.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"ACCELERATION has Key zero! (check if the application is correctly registered", "" )
    if(DENSITY.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"DENSITY has Key zero! (check if the application is correctly registered", "" )
     if(VOLUME_ACCELERATION.Key() == 0)
       //KRATOS_THROW_ERROR( std::invalid_argument,"VOLUME_ACCELERATION has Key zero! (check if the application is correctly registered", "" )
    if(CROSS_SECTION_AREA.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"CROSS_SECTION_AREA has Key zero! (check if the application is correctly registered", "" )
    if(LOCAL_CONSTITUTIVE_MATRIX.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"LOCAL_CONSTITUTIVE_MATRIX has Key zero! (check if the application is correctly registered", "" )
	  
    if(LOCAL_INERTIA_TENSOR.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"LOCAL_INERTIA_TENSOR has Key zero! (check if the application is correctly registered", "" )
    if(ROTATION.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"ROTATION has Key zero! (check if the application is correctly registered", "" )

    //verify that the dofs exist
    for(unsigned int i=0; i<this->GetGeometry().size(); i++)
    {
        if(this->GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false)
            KRATOS_THROW_ERROR( std::invalid_argument,"missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id() )
        if(this->GetGeometry()[i].HasDofFor(DISPLACEMENT_X) == false || this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false || this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Z) == false)
            KRATOS_THROW_ERROR( std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id() )
    }

    //verify that the area is given by properties
    if (this->GetProperties().Has(CROSS_SECTION_AREA)==false)
    {
        if( GetValue(CROSS_SECTION_AREA) == 0.0 )
            KRATOS_THROW_ERROR( std::logic_error,"CROSS_SECTION_AREA not provided for this element", this->Id() )
    }

    //verify that the inertia is given by properties
    if (this->GetProperties().Has(LOCAL_INERTIA_TENSOR)==false)
    {
        if( GetValue(LOCAL_INERTIA_TENSOR)(0,0) == 0.0 )
	  KRATOS_THROW_ERROR( std::logic_error,"LOCAL_INERTIA_TENSOR not provided for this element ", this->Id() )
    }


    return 0;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    int IntMethod = int(mThisIntegrationMethod);
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("CurrentCurvatureVectors",mCurrentCurvatureVectors);
    rSerializer.save("InitialLocalQuaternionsReduced",mInitialLocalQuaternionsReduced);
    rSerializer.save("CurrentLocalQuaternionsReduced",mCurrentLocalQuaternionsReduced);
    rSerializer.save("PreviousLocalQuaternionsReduced",mPreviousLocalQuaternionsReduced);
    rSerializer.save("InitialLocalQuaternionsFull",mInitialLocalQuaternionsFull);
    rSerializer.save("CurrentLocalQuaternionsFull",mCurrentLocalQuaternionsFull);
    rSerializer.save("PreviousLocalQuaternionsFull",mPreviousLocalQuaternionsFull);
    rSerializer.save("IterationCounter",mIterationCounter);
    rSerializer.save("InvJ0Reduced",mInvJ0Reduced);
    rSerializer.save("DetJ0Reduced",mDetJ0Reduced);
    rSerializer.save("InvJ0Full",mInvJ0Full);
    rSerializer.save("DetJ0Full",mDetJ0Full);
}

void LargeDisplacementBeamElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("CurrentCurvatureVectors",mCurrentCurvatureVectors);
    rSerializer.load("InitialLocalQuaternionsReduced",mInitialLocalQuaternionsReduced);
    rSerializer.load("CurrentLocalQuaternionsReduced",mCurrentLocalQuaternionsReduced);
    rSerializer.load("PreviousLocalQuaternionsReduced",mPreviousLocalQuaternionsReduced);
    rSerializer.load("InitialLocalQuaternionsFull",mInitialLocalQuaternionsFull);
    rSerializer.load("CurrentLocalQuaternionsFull",mCurrentLocalQuaternionsFull);
    rSerializer.load("PreviousLocalQuaternionsFull",mPreviousLocalQuaternionsFull);
    rSerializer.load("IterationCounter",mIterationCounter);
    rSerializer.load("InvJ0Reduced",mInvJ0Reduced);
    rSerializer.load("DetJ0Reduced",mDetJ0Reduced);
    rSerializer.load("InvJ0Full",mInvJ0Full);
    rSerializer.load("DetJ0Full",mDetJ0Full);
}


} // Namespace Kratos


