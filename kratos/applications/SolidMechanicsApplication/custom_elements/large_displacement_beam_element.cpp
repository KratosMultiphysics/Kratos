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
#include "custom_elements/small_displacement_beam_element_3D2N.hpp"
#include "utilities/math_utils.h"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

#include "solid_mechanics_application.h"


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
  //return mThisIntegrationMethod;
  return GeometryData::GI_GAUSS_3; //writing results with 3 integration points....
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

      rValues[index + 3] = 0.0;
      rValues[index + 4] = 0.0;
      rValues[index + 5] = 0.0;
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

      rValues[index + 3] = 0.0;
      rValues[index + 4] = 0.0;
      rValues[index + 5] = 0.0;
    }


    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

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

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    //Nodal Rotations
    if ( mTotalCompoundRotations.size() != number_of_nodes )
      {
        mTotalCompoundRotations.resize( number_of_nodes );
      }


    //Nodal Curvatures
    if ( mCurvatures.size() != number_of_nodes )
      {
        mCurvatures.resize( number_of_nodes );
      }

    
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	mTotalCompoundRotations[i].Quaternion = QuaternionType::Identity();
	mTotalCompoundRotations[i].Center     = GetGeometry()[i].Coordinates();
	
	mCurvatures[i] = ZeroVector(3);
      }


    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //Global to Local Transformation Matrix for the integration frame (the beam)
    
    Matrix GlobalToLocalRotationMatrix = ZeroMatrix(dimension);

    this->CalculateGlobalToLocalTransformationMatrix( GlobalToLocalRotationMatrix );
    
    mGlobalToLocalQuaternion = QuaternionType::FromRotationMatrix( GlobalToLocalRotationMatrix );


    //Set step variables to be kept and updated during the iteration
    this->InitializeStepVariables( mStepVariables );

    //Set step variables to local frame
    this->TransformVariablesToLocal( mStepVariables, mGlobalToLocalQuaternion );

    KRATOS_CATCH( "" )
}



////************************************************************************************
////************************************************************************************
void LargeDisplacementBeamElement::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //Global to Local Transformation Matrix for the integration frame (the beam)
    
    Matrix GlobalToLocalRotationMatrix = ZeroMatrix(dimension);

    this->CalculateGlobalToLocalTransformationMatrix( GlobalToLocalRotationMatrix );
    
    mGlobalToLocalQuaternion = QuaternionType::FromRotationMatrix( GlobalToLocalRotationMatrix );

    KRATOS_CATCH( "" )
   

}

////************************************************************************************
////************************************************************************************

void LargeDisplacementBeamElement::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
     KRATOS_TRY
    
     //Set step variables are updated in each iteration
     this->UpdateStepVariables( mStepVariables, rCurrentProcessInfo );
 
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

void LargeDisplacementBeamElement::InitializeStepVariables(IterationVariables& rGlobalStepVariables)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    // initialize variables
    rGlobalStepVariables.TotalCompoundRotationVectors.resize(number_of_nodes);
    rGlobalStepVariables.CurrentStepRotationVectors.resize(number_of_nodes);
    rGlobalStepVariables.PreviousIterationRotationVectors.resize(number_of_nodes);
    rGlobalStepVariables.AngularVelocityVectors.resize(number_of_nodes);
    rGlobalStepVariables.AngularAccelerationVectors.resize(number_of_nodes);
    rGlobalStepVariables.CurvatureVectors.resize(number_of_nodes);

    Vector VectorZero = ZeroVector(dimension);

    Vector ReferenceRotation = VectorZero;

    QuaternionType PreviouStepQuaternion;

    // set values from initial prediction
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {

	rGlobalStepVariables.TotalCompoundRotationVectors[i]     = VectorZero;
	rGlobalStepVariables.CurrentStepRotationVectors[i]       = VectorZero;
	rGlobalStepVariables.PreviousIterationRotationVectors[i] = VectorZero;
	rGlobalStepVariables.AngularVelocityVectors[i]           = VectorZero;
	rGlobalStepVariables.AngularAccelerationVectors[i]       = VectorZero;
	rGlobalStepVariables.CurvatureVectors[i]                 = VectorZero;

	
	//1.-Total Rotation Vector
	rGlobalStepVariables.TotalCompoundRotationVectors[i] = GetCurrentValue( ROTATION, rGlobalStepVariables.TotalCompoundRotationVectors[i], i );
	
	//2.-Step Rotation Vector 
	ReferenceRotation = GetReferenceValue( ROTATION, ReferenceRotation, i );

	PreviousStepQuaternion = QuaternionType::FromRotationVector( ReferenceRotation );

	PreviousStepQuaternion.RotateVector3( rGlobalStepVariables.TotalCompoundRotationVectors[i], rGlobalStepVariables.CurrentStepRotationVectors[i] );   
		  

	//4.-Angular Velocity Vector
	rGlobalStepVariables.AngularVelocityVectors[i] = GetCurrentValue( ANGULAR_VELOCITY, rGlobalStepVariables.AngularVelocityVectors[i], i );

	//5.-Angular Acceleration Vector
	rGlobalStepVariables.AngularAccelerationVectors[i] = GetCurrentValue( ANGULAR_ACCELERATION, rGlobalStepVariables.AngularAccelerationVectors[i], i );

	//6.-Curvature Vector
	rGlobalStepVariables.CurvatureVectors[i] = GetCurrentValue( CURVATURE, rGlobalStepVariables.CurvatureVectors[i], i );


      }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::UpdateStepVariables(IterationVariables& rGlobalStepVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    std::vector<Vector> IncrementRotationVectors;
    IncrementRotationVectors.resize( number_of_nodes );

    Vector IncrementStepRotation;

    QuaternionType TotalQuaternion;
    QuaternionType CurrentStepQuaternion;
    
    QuaternionType PreviouStepQuaternion;
    QuaternionType IncrementRotationQuaternion;    

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {

	//0.- Compute Increment Rotation Vector and Quaternion: increment of the rotation from iteration i to i+1
	PreviousStepQuaternion = QuaternionType::FromRotationVector( rGlobalStepVariables.PreviousIterationRotationVectors[i] );

	PreviousStepQuaternion.RotateVector3( rGlobalStepVariables.CurrentStepRotationVectors[i], IncrementRotationVectors[i] );

	IncrementRotationQuaternion = QuaternionType::FromRotationVector( IncrementRotationVectors[i] );
	
	//also linear rotation increment
	IncrementIterationStepRotation = rGlobalStepVariables.CurrentStepRotationVectors[i] - rGlobalStepVariables.PreviousIterationRotationVectors[i];

	//1.-Total Rotation Vector
	IncrementRotationQuaternion.RotateVector3( rGlobalStepVariables.TotalCompoundRotationVectors[i] );

	TotalQuaternion = QuaternionType::FromRotationVector( rGlobalStepVariables.TotalCompoundRotationVectors[i] );
	
	//2.-Step Rotation Vector
	IncrementRotationQuaternion.RotateVector3( rGlobalStepVariables.CurrentStepRotationVectors[i] );

	CurrentStepQuaternion = QuaternionType::FromRotationVector( rGlobalStepVariables.CurrentStepRotationVectors[i] );
		
	//4.-Angular Velocity Vector

	//rCurrentProcessInfo must give it:
	double beta  = 1;
	double gamma = 1;
	double DeltaTime = 1;

	TotalQuaternion = TotalQuaternion.conjugate() * CurrentStepQuaternion;

	TotalQuaternion.RotateVector3( IncrementIterationStepRotation );
	
	rGlobalStepVariables.AngularVelocityVectors[i] += (gamma/ ( beta * DeltaTime ) ) * IncrementIterationStepRotation;

	//5.-Angular Acceleration Vector
	rGlobalStepVariables.AngularAccelerationVectors[i] += (1.0/ ( beta * DeltaTime * DeltaTime ) ) * IncrementIterationStepRotation;

	//6.-Curvature Vector

	//6.1- Update Term of the axial vector for the determination of the curvature
	IncrementRotationQuaternion.RotateVector3( rGlobalStepVariables.CurvatureVectors[i] );
	
      }


    //6.2.-Calculate the Curvature rotation due to the iteration incremental rotation
    Vector DeltaCurvature = ZeroVector(3);
    
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {

	 DeltaCurvature = this->CalculteDeltaCurvature( IncrementRotationVectors[i], i );
		    
	 rGlobalStepVariables.CurvatureVectors[i] += DeltaCurvature;

      }
  

    
    //7.-Update Previous Iteration Rotation Vectors
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	rGlobalStepVariables.PreviousIterationRotationVectors[i] = rGlobalStepVariables.CurrentStepRotationVectors[i];
      }
	
    
    KRATOS_CATCH( "" )
}


///************************************************************************************
//************************************************************************************

Vector LargeDisplacementBeamElement::CalculateDeltaCurvature(std::vector<Vector>& rIncrementRotationVectors, const unsigned int& rNode)
{
    KRATOS_TRY

    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    //reading integration points (in fact is the two nodes beam element, only one integration point)
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
	

    Vector DeltaIncrementRotationVector = ZeroVector();
    //(in fact is the two nodes beam element, only one integration point)
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {
	
	//compute element kinematics  ...
	this->CalculateKinematics(Variables,PointNumber);
	    	    
	//Calculate the Curvature rotation due to the iteration incremental rotation:
	for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
	    DeltaIncrementRotationVector += Variables.DN_DX( rNode, 0 ) * IncrementRotationVectors[i];
	  }
      }
      

    double RotationModulus = IncrementRotationVectors[rNode].norm_2();
    double RotationSinus1  = std::sin(RotationModulus) / RotationModulus;
    double RotationSinus2  = std::sin(0.5 * RotationModulus) * 0.5 / RotationModulus ;
	
    
    Vector Term1 = (RotationSinus1) * DeltaIncrementRotationVector;
    Vector Term2 = ( 1.0 - RotationSinus1 ) * ( IncrementRotationVectors[rNode] * DeltaIncrementRotationVector / RotationModulus ) * ( 1.0 / RotationModulus ) * IncrementRotationVectors[rNode];
    Vector Term3 = 0.5 * RotationSinus2 * RotationSinus2 * MathUtils::CrossProduct(IncrementRotationVectors[rNode], DeltaIncrementRotationVector);
    
    Vector DeltaCurvature = Term1 + Term2 + Term3;

    return DeltaCurvature;

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::TransformVariablesToLocal(IterationVariables& rGlobalStepVariables, QuaternionType& rGlobalToLocalQuaternion)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	rGlobalToLocalQuaternion.RotateVector3(rGlobalStepVariables.TotalCompoundRotationVectors[i]);
	rGlobalToLocalQuaternion.RotateVector3(rGlobalStepVariables.CurrentStepRotationVectors[i]);
	rGlobalToLocalQuaternion.RotateVector3(rGlobalStepVariables.PreviousIterationRotationVectors[i]);
	rGlobalToLocalQuaternion.RotateVector3(rGlobalStepVariables.AngularVelocityVectors[i]);
	rGlobalToLocalQuaternion.RotateVector3(rGlobalStepVariables.AngularAccelerationVectors[i]);
	rGlobalToLocalQuaternion.RotateVector3(rGlobalStepVariables.CurvatureVectors[i]);
      }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::TransformVariablesToGlobal(IterationVariables& rGlobalStepVariables, QuaternionType& rGlobalToLocalQuaternion)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    QuaternionType LocalToGlobalQuaternion = rGlobalToLocalQuaternion.conjugate();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	LocalToGlobalQuaternion.RotateVector3(rGlobalStepVariables.TotalCompoundRotationVectors[i]);
	LocalToGlobalQuaternion.RotateVector3(rGlobalStepVariables.CurrentStepRotationVectors[i]);
	LocalToGlobalQuaternion.RotateVector3(rGlobalStepVariables.PreviousIterationRotationVectors[i]);
	LocalToGlobalQuaternion.RotateVector3(rGlobalStepVariables.AngularVelocityVectors[i]);
	LocalToGlobalQuaternion.RotateVector3(rGlobalStepVariables.AngularAccelerationVectors[i]);
	LocalToGlobalQuaternion.RotateVector3(rGlobalStepVariables.CurvatureVectors[i]);
      }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

Vector& LargeDisplacementBeamElement::GetCurrentValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const unsigned int& rNode)
{
    KRATOS_TRY

    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    array_1d<double,3> ArrayValue;
    ArrayValue = GetGeometry()[rNode].FastGetSolutionStepValue( rVariable );
    
    if( rValue.size() != dimension )
      rValue.resize(dimension, false);

    for( int i=0; i<dimension; i++ )
      {
	rValue[i] = ArrayValue[i];
      }
   

    return rValue;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

Vector& LargeDisplacementBeamElement::GetReferenceValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const unsigned int& rNode)
{
    KRATOS_TRY

    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    array_1d<double,3> ArrayValue;
    ArrayValue = GetGeometry()[rNode].FastGetSolutionStepValue( rVariable, 1 );
    
    if( rValue.size() != dimension )
      rValue.resize(dimension, false);

    for( int i=0; i<dimension; i++ )
      {
	rValue[i] = ArrayValue[i];
      }
   

    return rValue;

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::InitializeGeneralVariables(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //Compute Section Properties:
    this->CalculateSectionProperties(rVariables.Section);

    rVariables.Length = GetGeometry().Length();

    if(rVariables.Length == 0.00)
      KRATOS_ERROR(std::invalid_argument, "Zero length found in elemnet #", this->Id());


    unsigned int voigtsize = 3;
    if( dimension == 3 )
    {
        voigtsize  = 6;
    }

    rVariables.ConstitutiveMatrix.resize( voigtsize, voigtsize );

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

    
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateSectionProperties(SectionProperties & rSection)

{
    KRATOS_TRY

   
    if( GetProperties().Has(CROSS_AREA) ){
        rSection.Area = GetProperties()[CROSS_AREA];
    }
    else{
        rSection.Area = GetValue(AREA);
    }


    if( GetProperties().Has(LOCAL_INERTIA) )
    {
        Matrix& inertia = GetProperties()[LOCAL_INERTIA];
        rSection.Inertia_z = inertia(0,0);
        rSection.Inertia_y = inertia(1,1);
        rSection.Polar_Inertia = inertia(0,1);
    }
    else
    {
        Matrix& inertia = GetValue(LOCAL_INERTIA);
        rSection.Inertia_z = inertia(0,0);
        rSection.Inertia_y = inertia(1,1);
        rSection.Polar_Inertia = inertia(0,1);
    }

    rSection.Rotational_Inertia = rSection.Polar_Inertia;

    KRATOS_CATCH( "" )

}


//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************

void LargeDisplacementElement::CalculateKinematics(GeneralVariables& rVariables, const double& rPointNumber)
{
    KRATOS_TRY


    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //initialize local transformation/rotation matrix
    rVariables.LocalToGlobalRotationMatrix = ZeroMatrix(3);

    (mGlobalToLocalQuaternion.conjugate()).ToRotationMatrix(rVariables.LocalToGlobalRotationMatrix);

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //Parent to reference configuration
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    //Set Shape Functions Values for this integration point
    rVariables.N=row( Ncontainer, rPointNumber);
    
    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n+1]
    Vector Jacobian = ZeroVector(3);
    Vector InverseJacobian = ZeroVector(3);

    for(int n=0; n<dimension; n++)
      Jacobian[n] = rVariables.j[rPointNumber](n, 0);

    noalias(InverseJacobian) = prod( trans(rVariables.LocalToGlobalRotationMatrix), Jacobian );

    //Select the axis direction (lets take axis e1 as local beam axis)
    rVariables.detJ = InverseJacobian[0];
				     
    //Compute cartesian derivatives [dN/dx_n+1]
    rVariables.DN_DX = rVariables.detJ * DN_De[rPointNumber]; //overwrites DX now is the current position dx
   

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

    //reading integration points (in fact is the two nodes beam element, only one integration point)
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //auxiliary terms
    Vector VolumeForce;

    //(in fact is the two nodes beam element, only one integration point)
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {

        //compute element kinematics  ...
        this->CalculateKinematics(Variables,PointNumber);

	//compute element ConstitutiveTensor
	this->CalculateConstitutiveMatrix(Variables);

	//compute element Strain and Stress Resultants and Couples
	this->CalculateStressResultants(Variables);

	double IntegrationWeight = integration_points[PointNumber].Weight() * rVariables.detJ;

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

      }


    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateConstitutiveMatrix(GeneralVariables& rVariables)
{
    KRATOS_TRY
      
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int size = 2 * dimension;

    if( rVariables.ConstitutiveMatrix.size1() != size )
      rVariables.ConstitutiveMatrix.resize(6,6, false);

    Matrix ConstitutiveMatrix  = ZeroMatrix(6,6);
    rVariables.ConstitutiveMatrix = ZeroMatrix(6,6);

    const double PoissonCoefficient = GetProperties()[POISSON_RATIO];
    const double YoungModulus       = GetProperties()[YOUNG_MODULUS];
    const double ShearModulus       = YoungModulus*0.5/(1.0 + PoissonCoefficient);

    //Calculate Local Constitutive Matrix:

    
    const unsigned int number_of_nodes = GetGeometry().size();

    //displacements and rotations vector
    Matrix TotalRotation(3);
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	TotalRotation  =  rVariables.N[i] * rVariables.TotalCompoundRotations[i];
      }

    //build the global transformation matrix for the constitutive matrix
    Matrix RotationMatrix = ZeroMatrix(6);
    for ( unsigned int i = 0; i < 3; i++ )
      {
	for ( unsigned int j = 0; j < 3; j++ )
	  {
	    RotationMatrix(i,j) = TotalRotation(i,j);
	  }
      }

    for ( unsigned int i = 0; i < 3; i++ )
      {
	for ( unsigned int j = 0; j < 3; j++ )
	  {
	    RotationMatrix(i+3,j+3) = TotalRotation(i,j);
	  }
      }
    
    //reference elastic constitutive matrix
    ConstitutiveMatrix( 0, 0 ) = YoungModulus * rVariables.Section.Area;
    ConstitutiveMatrix( 1, 1 ) = ShearModulus * rVariables.Section.Area;
    ConstitutiveMatrix( 2, 2 ) = ShearModulus * rVariables.Section.Area;

    ConstitutiveMatrix( 3, 3 ) = ShearModulus * rVariables.Section.Inertia_zy;
    ConstitutiveMatrix( 4, 4 ) = YoundModulus * rVariables.Section.Inertia_y;
    ConstitutiveMatrix( 5, 5 ) = YoundModulus * rVariables.Section.Inertia_z;
    
    //current elastic constitutive matrix
    noalias(rVariables.ConstitutiveMatrix = prod( RotationMatrix, prod( rConstitutiveMatrix, trans(RotationMatrix) ) );


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateStressResultants(GeneralVariables& rVariables)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();     
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    Vector StrainResultants = ZeroVector(dimension);
    Vector StrainCouples    = ZeroVector(dimension);

    Matrix RotationVector   = ZeroVector(dimension);
    Vector Displacements    = ZeroVector(dimension);

    //strains due to displacements and rotations
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	Displacements = GetCurrentValue( DISPLACEMENT, Displacements, i );
    
	StrainResultants +=  rVariables.DN_DX(i,0) * Displacements;
	StrainCouples    +=  rVariables.N[i] * mStepVariables.CurvatureVectors[i];
	RotationVector   +=  rVariables.N[i] * mStepVariables.TotalCompoundRotationVectors[i];
      }

    //axis in e1 direction
    StrainResultants[0] = -1.0;

    RotationVector *=(-1.0);
    QuaternionType TotalQuaternion = QuaternionType::FromRotationVector( RotationVector );

    TotalQuaternion.RotateVector3( StrainResultants );
    TotalQuaternion.RotateVector3( StrainCouples );       

    for ( unsigned int i = 0; i < dimension; i++ )
      {
	rVariables.StrainVector[i]   = StrainResultants;
	rVariables.StrainVector[i+3] = StrainCouples;	
      }
    
    Matrix ConstitutiveMatrix  = ZeroMatrix(6,6);
    rVariables.ConstitutiveMatrix = ZeroMatrix(6,6);

    const double PoissonCoefficient = GetProperties()[POISSON_RATIO];
    const double YoungModulus       = GetProperties()[YOUNG_MODULUS];
    const double ShearModulus       = YoungModulus*0.5/(1.0 + PoissonCoefficient);

    //Calculate Local Constitutive Matrix:
    ConstitutiveMatrix( 0, 0 ) = YoungModulus * rVariables.Section.Area;
    ConstitutiveMatrix( 1, 1 ) = ShearModulus * rVariables.Section.Area;
    ConstitutiveMatrix( 2, 2 ) = ShearModulus * rVariables.Section.Area;

    ConstitutiveMatrix( 3, 3 ) = ShearModulus * rVariables.Section.Polar_Inertia;
    ConstitutiveMatrix( 4, 4 ) = YoundModulus * rVariables.Section.Inertia_y;
    ConstitutiveMatrix( 5, 5 ) = YoundModulus * rVariables.Section.Inertia_z;
    
    //current elastic constitutive matrix
    noalias(rVariables.StressVector) = prod( rConstitutiveMatrix, StrainVector );

    Vector StressResultants = ZeroVector(dimension);
    Vector StressCouples    = ZeroVector(dimension);

    for ( unsigned int i = 0; i < dimension; i++ )
      {
	StressResultants[i] = rVariables.StressVector[i];
	StressCouples[i]    = rVariables.StressVector[i+3];
      }


    TotalQuaternion = TotalQuaternion.conjugate();

    TotalQuaternion.RotateVector3( StressResultants );
    TotalQuaternion.RotateVector3( StressCouples );       

    for ( unsigned int i = 0; i < dimension; i++ )
      {
	rVariables.StressVector[i]   = StressResultants;
	rVariables.StressVector[i+3] = StressCouples;	
      }

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
 
    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

    this->CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    this->CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    //follower load stiffness (to implement)

    //std::cout<<" rLeftHandSideMatrix "<<rLeftHandSideMatrix<<std::endl;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& VolumeForce, double& rIntegrationWeight)
{
    KRATOS_TRY

    VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector(); 
  
    // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
    this->CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

    // operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
    this->CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight );

    //follower load forces (to implement)


    //inertial forces (to implement)


    //std::cout<<" rRightHandSideVector "<<rRightHandSideVector<<std::endl;

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

    double DomainSize = (rVariables.Section.Area / rVariables.detJ );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = (dimension * 2) * i;
        for ( unsigned int j = 0; j < dimension; j++ )
        {
	  rRightHandSideVector[index + j] += rIntegrationWeight * rVariables.N[i] * rVolumeForce[j] * DomainSize;
        }
    }

    KRATOS_CATCH( "" )

}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::VectorToSkewSymmetricTensor( const Vector& rVector, 
								Matrix& rSkewSymmetricTensor )
{
    KRATOS_TRY

    //Initialize Local Matrices
    if( rSkewSymmetricTensor.size1() != 3 )
      rSkewSymmetricTensor.resize(3, 3, false);
    
    rSkewSymmetricTensor = ZeroMatrix(3);

    SkewSymmetricTensor( 0, 1 ) =  rVector[2];
    SkewSymmetricTensor( 0, 2 ) = -rVector[1];
    SkewSymmetricTensor( 1, 2 ) =  rVector[0];

    SkewSymmetricTensor( 1, 0 ) = -rVector[2];
    SkewSymmetricTensor( 2, 0 ) =  rVector[1];
    SkewSymmetricTensor( 2, 1 ) = -rVector[0];

    
    KRATOS_CATCH( "" )

}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateOperator(MatrixType& rOperator,
						     GeneralVariables& rVariables,
						     const int& rNode)
{
    KRATOS_TRY

    //Initialize Local Matrices
    if( rOperator.size1() != 6 )
      rOperator.resize(6, 6, false);
    
    rOperator = ZeroMatrix(6);

    rOperator( 0, 0 ) =  rVariables.N[rNode];
    rOperator( 1, 1 ) =  rVariables.N[rNode];
    rOperator( 2, 2 ) =  rVariables.N[rNode];
    rOperator( 3, 3 ) =  rVariables.N[rNode];
    rOperator( 4, 4 ) =  rVariables.N[rNode];
    rOperator( 5, 5 ) =  rVariables.N[rNode];


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
    unsigned int MatSize = rRightHandSideVector.size();

    Vector NodalForceVector(MatSize);

    //Initialize Local Matrices
    MatrixType DifferentialOperatorI = ZeroMatrix(MatSize);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	NodalForveVector = ZeroVector(MatSize);
 	this->CalculateDifferentialOperator( DifferencialOperatorI, rVariables, i );

	//nodal force vector
	noalias(NodalForceVector) = prod( DifferencialOperatorI, rVariables.StressVector );

	noalias(rRightHandSideVector) -= prod(rVariables.LocalToGlobalRotationMatrix, NodalForceVector);
      }
 

    KRATOS_CATCH( "" )

}



//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateDifferentialOperator(MatrixType& rDifferentialOperator,
								 GeneralVariables& rVariables,
								 const int& rNode)
{
    KRATOS_TRY

    //Initialize Local Matrices
    if( rDifferentialOperator.size1() != 6 )
      rDifferentialOperator.resize(6, 6, false);
    
    rDifferentialOperator = ZeroMatrix(6);

    rDifferentialOperator( 0, 0 ) =  rVariables.DN_DX( rNode, 0 );
    rDifferentialOperator( 1, 1 ) =  rVariables.DN_DX( rNode, 0 );
    rDifferentialOperator( 2, 2 ) =  rVariables.DN_DX( rNode, 0 );
    rDifferentialOperator( 3, 3 ) =  rVariables.DN_DX( rNode, 0 );
    rDifferentialOperator( 4, 4 ) =  rVariables.DN_DX( rNode, 0 );
    rDifferentialOperator( 5, 5 ) =  rVariables.DN_DX( rNode, 0 );


    //locate stress resultants in skew-symmetric transpose form
    Matrix SkewSymTotalCompoundRotations = ZeroMatrix(3);
    
    SkewSymTotalCompoundRotations( 0, 1 ) =  mStepVariables.TotalCompoundRotationVectors[rNode][2] * rVariables.DN_DX( 2, 0 );
    SkewSymTotalCompoundRotations( 0, 2 ) = -mStepVariables.TotalCompoundRotationVectors[rNode][1] * rVariables.DN_DX( 1, 0 );
    SkewSymTotalCompoundRotations( 1, 2 ) =  mStepVariables.TotalCompoundRotationVectors[rNode][0] * rVariables.DN_DX( 0, 0 );

    SkewSymTotalCompoundRotations( 1, 0 ) = -mStepVariables.TotalCompoundRotationVectors[rNode][2] * rVariables.DN_DX( 2, 0 );
    SkewSymTotalCompoundRotations( 2, 0 ) =  mStepVariables.TotalCompoundRotationVectors[rNode][1] * rVariables.DN_DX( 1, 0 );
    SkewSymTotalCompoundRotations( 2, 1 ) = -mStepVariables.TotalCompoundRotationVectors[rNode][0] * rVariables.DN_DX( 0, 0 );


    for ( unsigned int i = 0; i < 3; i++ )
      {
	for ( unsigned int j = 0; j < 3; j++ )
	  {
	    rDifferentialOperator( i, j+3 ) = SkewSymTotalCompoundRotations(i,j);
   	  }
      }


    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
						       GeneralVariables& rVariables,
						       double& rIntegrationWeight)
{
    KRATOS_TRY

    unsigned int MatSize = rLeftHandSideMatrix.size1();

    MatrixType LocalStiffnessMatrix(MatSize);

    LocalStiffnessMatrix = ZeroMatrix(MatSize);

    //Initialize Local Matrices
    MatrixType DifferentialOperatorI = ZeroMatrix(MatSize);
    MatrixType DifferentialOperatorJ = ZeroMatrix(MatSize);

    const unsigned int number_of_nodes = GetGeometry().size();

    MatrixType NodalStiffnessMatrix = ZeroMatrix(MatSize);
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	NodalStiffnessMatrix = ZeroMatrix(MatSize);
 	this->CalculateDifferentialOperator( DifferencialOperatorI, rVariables, i );

	for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {
	    this->CalculateDifferentialOperator( DifferencialOperatorJ, rVariables, j );

	    noalias(NodalStiffnessMatrix) = prod( DifferentialOperatorI, prod( rVariables.ConstitutiveMatrix, trans(DifferentialOperatorJ) ) );

	    
	    //Building the Local Stiffness Matrix
	    for (unsigned int ii=0; ii<dimension; ii++)
	      {
		for(unsigned int jj=0; jj<dimension; jj++)
		  {
		    LocalStiffnessMatrix(ii+i,jj+j) = NodalStiffnessMatrix(ii,jj);
		  }
	      }
	    
	  }
      }


    Matrix RotationMatrix( MatSize );
    //Building the rotation matrix for the local element matrix
    for (unsigned int kk=0; kk < MatSize; kk += dimension)
    {
        for (unsigned int i=0; i<dimension; i++)
        {
            for(unsigned int j=0; j<dimension; j++)
            {
	      RotationMatrix(i+kk,j+kk) = rVariables.LocalToGlobalRotationMatrix(i,j);
            }
        }
    }

    //Rotate Local Stiffness Matrix
    Matrix aux_matrix   = ZeroMatrix(MatSize);
    noalias(aux_matrix) = prod(RotationMatrix, LocalStiffnessMatrix);

    //Stiffness Matrix
    LocalStiffnessMatrix = ZeroMatrix(MatSize);
    noalias(LocalStiffnessMatrix) = prod(aux_matrix,Matrix(trans(RotationMatrix)));

    noalias(rLeftHandSideMatrix) += LocalStiffnessMatrix;


    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateBmatrix(MatrixType& rBmatrix,
						    GeneralVariables& rVariables)
{
    KRATOS_TRY

    //Initialize Local Matrices
    if( rDifferenctialOperator.size1() != 9 )
      rDiscreteOperator.resize(9, 9, false);
    
    rBmatrix = ZeroMatrix(9);

    //locate stress resultants in skew-symmetric form
    Matrix SkewSymStressResultants = ZeroMatrix(3);
    
    SkewSymStressResultants( 0, 1 ) =  rVariables.StressVector[2];
    SkewSymStressResultants( 0, 2 ) = -rVariables.StressVector[1];
    SkewSymStressResultants( 1, 2 ) =  rVariables.StressVector[0];

    SkewSymStressResultants( 1, 0 ) = -rVariables.StressVector[2];
    SkewSymStressResultants( 2, 0 ) =  rVariables.StressVector[1];
    SkewSymStressResultants( 2, 1 ) = -rVariables.StressVector[0];

    for ( unsigned int i = 0; i < 3; i++ )
      {
	for ( unsigned int j = 0; j < 3; j++ )
	  {
	    rBmatrix( i+6, j ) =  SkewSymStressResultants(i,j);
	    rBmatrix( i, j+6 ) = -SkewSymStressResultants(i,j);  
   	  }
      }


    //locate stress couples in skew-symmetric form
    Matrix SkewSymStressCouples = ZeroMatrix(3);
    
    SkewSymStressCouples( 0, 1 ) =  rVariables.StressVector[5];
    SkewSymStressCouples( 0, 2 ) = -rVariables.StressVector[4];
    SkewSymStressCouples( 1, 2 ) =  rVariables.StressVector[3];

    SkewSymStressCouples( 1, 0 ) = -rVariables.StressVector[5];
    SkewSymStressCouples( 2, 0 ) =  rVariables.StressVector[4];
    SkewSymStressCouples( 2, 1 ) = -rVariables.StressVector[3];

 
    //locate stress couples in skew-symmetric form

    for ( unsigned int i = 0; i < 3; i++ )
      {
	for ( unsigned int j = 0; j < 3; j++ )
	  {
	    rBmatrix( i+3, j+6 ) = -SkewSymStressCouples(i,j);
   	  }
      }


    //calculate stress delta rotation matrix

    Vector StressResultants = ZeroVector(3);

    StressResultants[0] = rVariables.StressVector[0];
    StressResultants[1] = rVariables.StressVector[1];
    StressResultants[2] = rVariables.StressVector[2];

    Vector DeltaTotalCompoundRotations = ZeroVector(3);
    
    const unsigned int number_of_nodes = GetGeometry().size();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      DeltaTotalCompoundRotations[0] += mStepVariables.TotalCompoundRotationVectors[i][0] * rVariables.DN_DX( i, 0 );
      DeltaTotalCompoundRotations[1] += mStepVariables.TotalCompoundRotationVectors[i][1] * rVariables.DN_DX( i, 0 );
      DeltaTotalCompoundRotations[2] += mStepVariables.TotalCompoundRotationVectors[i][2] * rVariables.DN_DX( i, 0 );
    }

    Matrix StressDeltaRotationMatrix  = outer_prod( StressResultants, DeltaTotalCompoundRotation );
    
    double stress_delta_rotation = inner_prod( StressResultants, DeltaTotalCompoundRotation );
    
    GeometricTerm( 0, 0 ) -= stress_delta_rotation;
    GeometricTerm( 1, 1 ) -= stress_delta_rotation;
    GeometricTerm( 2, 2 ) -= stress_delta_rotation;


    for ( unsigned int i = 0; i < 3; i++ )
      {
	for ( unsigned int j = 0; j < 3; j++ )
	  {
	    rBmatrix( i+6, j+6 ) = StressDeltaRotationMatrix(i,j);
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
    if( rDifferenctialOperator.size1() != 6 )
      rDiscreteOperator.resize(6, 6, false);
    
    rDiscreteOperator = ZeroMatrix(6);

    rDiscreteOperator( 0, 0 ) =  rVariables.DN_DX( rNode, 0 );
    rDiscreteOperator( 1, 1 ) =  rVariables.DN_DX( rNode, 0 );
    rDiscreteOperator( 2, 2 ) =  rVariables.DN_DX( rNode, 0 );
    rDiscreteOperator( 3, 3 ) =  rVariables.DN_DX( rNode, 0 );
    rDiscreteOperator( 4, 4 ) =  rVariables.DN_DX( rNode, 0 );
    rDiscreteOperator( 5, 5 ) =  rVariables.DN_DX( rNode, 0 );

    rDiscreteOperator( 3, 6 ) =  1.0; 
    rDiscreteOperator( 4, 7 ) =  1.0;
    rDiscreteOperator( 5, 8 ) =  1.0;

    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
						       GeneralVariables& rVariables,
						       double& rIntegrationWeight)
{
    KRATOS_TRY

    unsigned int MatSize = rLeftHandSideMatrix.size1();

    MatrixType LocalStiffnessMatrix(MatSize);

    LocalStiffnessMatrix = ZeroMatrix(MatSize);

    //Initialize Local Matrices
    MatrixType Bmatrix = ZeroMatrix(MatSize+3, MatSize+3);
    MatrixType DiscreteOperatorI = ZeroMatrix(MatSize, MatSize+3);
    MatrixType DiscreteOperatorJ = ZeroMatrix(MatSize, MatSize+3);

    const unsigned int number_of_nodes = GetGeometry().size();

    MatrixType NodalStiffnessMatrix = ZeroMatrix(MatSize);
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	NodalStiffnessMatrix = ZeroMatrix(MatSize);

 	this->CalculateDiscreteOperator( DiscreteOperatorI, rVariables, i );

	for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {
	    this->CalculateDiscreteOperator( DiscreteOperatorJ, rVariables, j );

	    noalias(NodalStiffnessMatrix) = prod( DiscreteOperatorI, prod( Bmatrix, trans(DiscreteOperatorJ) ) );

	    //Building the Local Stiffness Matrix
	    for (unsigned int ii=0; ii<dimension; ii++)
	      {
		for(unsigned int jj=0; jj<dimension; jj++)
		  {
		    LocalStiffnessMatrix(ii+i,jj+j) = NodalStiffnessMatrix(ii,jj);
		  }
	      }
	    
	  }
      }


    Matrix RotationMatrix( MatSize );
    //Building the rotation matrix for the local element matrix
    for (unsigned int kk=0; kk < MatSize; kk += dimension)
    {
        for (unsigned int i=0; i<dimension; i++)
        {
            for(unsigned int j=0; j<dimension; j++)
            {
	      RotationMatrix(i+kk,j+kk) = rVariables.LocalToGlobalRotationMatrix(i,j);
            }
        }
    }

    //Rotate Local Stiffness Matrix
    Matrix aux_matrix   = ZeroMatrix(MatSize);
    noalias(aux_matrix) = prod(RotationMatrix, LocalStiffnessMatrix);

    //Stiffness Matrix
    LocalStiffnessMatrix = ZeroMatrix(MatSize);
    noalias(LocalStiffnessMatrix) = prod(aux_matrix,Matrix(trans(RotationMatrix)));

    noalias(rLeftHandSideMatrix) += LocalStiffnessMatrix;


    KRATOS_CATCH( "" )

}


//*****************************************************************************
//*****************************************************************************

void LargeDisplacementBeamElement::CalculateGlobalToLocalTransformationMatrix(Matrix& rRotationMatrix)

{

    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    Matrix AuxiliarRotationMatrix = ZeroMatrix(dimension);

    this->CalculateLocalToGlobalTransformationMatrix(AuxiliarRotationMatrix);

    rRotationMatrix = trans( RotationMatrix );
      
    KRATOS_CATCH( "" )

}

//*****************************************************************************
//*****************************************************************************

void LargeDisplacementBeamElement::CalculateLocalToGlobalTransformationMatrix(Matrix& rRotationMatrix)

{

    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       size            = number_of_nodes * dimension;
    unsigned int       MatSize         = 2 * size;

    Vector NormalVectors        = ZeroVector(9);    // vector containing director cosines
    Vector ReferenceCoordinates = ZeroVector(size);

    double nx, ny, nz, theta;
    
    ReferenceCoordinates[0] = GetGeometry()[0].X();
    ReferenceCoordinates[1] = GetGeometry()[0].Y();
    ReferenceCoordinates[2] = GetGeometry()[0].Z();     

    int k = number_of_nodes - 1 ;

    ReferenceCoordinates[3] = GetGeometry()[k].X();
    ReferenceCoordinates[4] = GetGeometry()[k].Y();
    ReferenceCoordinates[5] = GetGeometry()[k].Z();     

    double length_inverse = ( 1.00 / rVariables.Length );
    
    for( unsigned int i = 0; i < dimension; i++ )
    {
      NormalVectors[i]  = (ReferenceCoordinates[i+3] - ReferenceCoordinates[i]);
      NormalVectors[i] *= length_inverse;
    }

    // local x-axis (e1_local) is the beam axis 
    nx = NormalVectors[0];
    ny = NormalVectors[1];
    nz = NormalVectors[2];


    // local axis parallel to the global plane x-y
    if (nx == 0.0)
    {
	theta = 0;

	// local y-axis (e2_local) is the vector product of the local e1(beam axis) and the local e2(parallel to the global x-y)
	NormalVectors[6] = -nz*cos(theta);
	NormalVectors[7] = -nz*sin(theta);
	NormalVectors[8] =  nx*cos(theta) + ny*sin(theta);


	// local z-axis (e3_local) is always parallel to the global plane x-y
	NormalVectors[6] = -sin(theta);
	NormalVectors[7] =  cos(theta);
	NormalVectors[8] =  0.0;

    }
    else
    {
        theta = atan(ny/nx);

	if(nx < 0.0)
	  theta = theta + PI;

	// local y-axis (e2_local) is always parallel to the global plane x-y
	NormalVectors[3] = -sin(theta);
	NormalVectors[4] =  cos(theta);
	NormalVectors[5] =  0.0;
	
	// local z-axis (e3_local) is the vector product of the local e1(beam axis) and the local e2(parallel to the global x-y)
	NormalVectors[6] = -nz*cos(theta);
	NormalVectors[7] = -nz*sin(theta);
	NormalVectors[8] =  nx*cos(theta) + ny*sin(theta);

    }

    //Transformation matrix T = [e1_local, e2_local, e3_local] for each column

    //Building the rotation matrix
    for (unsigned int i=0; i<dimension; i++)
      {
	for(unsigned int j=0; j<dimension; j++)
	  {
	    rRotationMatrix(i,j) = NormalVectors[3*j+i];
	  }
      }


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
        if( GetGeometry()[j].SolutionStepsDataHas(VOLUME_ACCELERATION) ) //temporary, will be checked once at the beginning only
            rVolumeForce += rN[j] * GetGeometry()[j].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }

    rVolumeForce *= GetProperties()[DENSITY];

    return rVolumeForce;

    KRATOS_CATCH( "" )
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

    double TotalMass = 0;
    TotalMass = this->CalculateTotalMass( Section, TotalMass );

    Vector LumpFact = ZeroVector(number_of_nodes); 

    LumpFact = GetGeometry().LumpingFactors(LumpFact);

    for( unsigned int i=0; i < number_of_nodes; i++ )
      {
        double temp = LumpFact[i] * TotalMass;

        for( unsigned int j=0; j < dimension; j++ )
	  {
 	    unsigned int index = i * (dimension * 2) + j;

            rMassMatrix(index,index) = temp;
	  }

      }

    //block m(2,2) of the mass matrix
    
 
    Vector AngularVelocityVector(3);    
    AngularVelocityVector = ZeroVector(3);
    Vector AngularAccelerationVector(3);    
    AngularAccelerationVector = ZeroVector(3);
    Vector CurrentCompoundRotationVector(3);
    CurrentCompoundRotationVector = ZeroVector(3);
    Vector ReferenceCompoundRotationVector(3);
    ReferenceCompoundRotationVector = ZeroVector(3);
    Vector CurrentStepRotationVector(3);
    CurrentStepRotationVector = ZeroVector(3);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	AngularVelocityVector           += N[i] * mStepVariables.AngularVelocityVectors[i];
	AngularAccelerationVector       += N[i] * mStepVariables.AngularAccelerationVectors[i];
	CurrentCompoundRotationVector   += N[i] * mStepVariables.CurrentStepRotationVectors[i];

	CurrentStepRotationVector       += N[i] * mStepVariables.CurrentStepRotationVectors[i];

	//Reference Rotation Vector 
	Vector ReferenceRotation = GetReferenceValue( ROTATION, ReferenceRotation, i );
	
	ReferenceCompoundRotationVector += N[i] * ReferenceRotationVector;	
      }
    

    //rCurrentProcessInfo must give it:
    double beta  = 1;
    double gamma = 1;
    double DeltaTime = 1;
    

    Matrix CurrentRotationMatrix(3);
    Matrix ReferenceRotationMatrix(3);
    CurrentRotationMatrix = ZeroMatrix(3);
    ReferenceRotationMatrix = ZeroMatrix(3);
    

    Matrix MassTerm1(3);
    Matrix MassTerm2(3);
    MassTerm1 = ZeroMatrix(3);
    MassTerm2 = ZeroMatrix(3);
    
    //1.-Get rotation matrix
    QuaternionType TotalQuaternion;    

    //1.1.-Current Rotation Matrix
    TotalQuaternion = QuaternionType::FromRotationVector( CurrentCompoundRotationVector );		
    TotalQuaternion.ToRotationMatrix( CurrentRotationMatrix );
    
    
    //1.2.-Reference Rotation Matrix
    TotalQuaternion = QuaternionType::FromRotationVector( ReferenceCompoundRotationVector );
    TotalQuaternion.ToRotationMatrix( ReferenceRotationMatrix );
    
    
    //2.-Get inertia dyadic
    Matrix InertiaDyadic(3);
    this->CalculateInertiaDyadic( Section, InertiaDyadic );
    
    //2.- Compute Term 1:
    
    Vector InertiaxAngularVelocity     = prod( InertiaDyadic, AngularVelocityVector );
    Vector InertiaxAngularAcceleration = prod( InertiaDyadic, AngularAccelerationVector );
    
    Vector VectorTerm1 = MathUtils::CrossProduct( AngularVelocityVector, InertiaxAngularVelocity);
    VectorTerm1 += InertiaxAngularAcceleration;
    
    this->VectorToSkewSymmetricTensor(VectorTerm1, MassTerm1);
    
    MassTerm1 = prod( CurrentRotationMatrix, MassTerm1 );
    
    MassTerm1 *= (-1);
    
    //3.- Compute Term 2:
    
    Matrix TensorInertiaxAngularVelocity(3);
    this->VectorToSkewSymmetricTensor( InertiaxAngularVelocity, TensorInertiaxAngularVelocity );
    
    Matrix TensorAngularVelocityxInertia(3);
    Matrix TensorAngularVelocity(3);
    this->VectorToSkewSymmetricTensor(AngularVelocityVector, TensorAngularVelocity );
    noalias(TensorAngularVelocityxInertia) = prod( TensorAngularVelocity, InertiaDyadic );
    
    
    MassTerm2 = InertiaDyadic - (DeltaTime * gamma) * TensorInertiaxAngularVelocity + (DeltaTime * gamma) * TensorAngularVelocityxInertia ;
    
    
    MassTerm2 = prod( CurrentRotationMatrix, MassTerm1 );
    
    MassTerm2 *= (1.0/(DeltaTime * beta));
    
    
    Matrix MassMatrixBlock2 = MassTerm1 + MassTerm2;
    
	
    MassMatrixBlock2 =  prod( MassMatrixBlock2, trans(ReferenceRotationMatrix));
    
    
    // Compute Linear Part of the Step Rotation
    Matrix LinearPartRotationTensor(3);
    this->CalculateRotationLinearPartTensor( CurrentStepRotationVector , LinearPartRotationTensor );
    
	
    MassMatrixBlock2 = prod( MassMatrixBlock2, LinearPartRotationTensor );
    

    for( unsigned int i=0; i < number_of_nodes; i++ )
      {
	
	for( unsigned int j=0; j < number_of_nodes; j++ )
	  {
	
	    for( unsigned int k=0; k < dimension; k++ )
	      {
		unsigned int index_i = i * (dimension * 2) + k;

		for( unsigned int l=0; l < dimension; l++ )
		  {
		    unsigned int index_j = j * (dimension * 2) + l;

		    rMassMatrix(index_i,index_j) += MassMatrixBlock2( k, l ) * N[i] * N[j] * GetGeometry().Length() * 0.5;
		  }
	      }
	  }
      }

    

    //initialize local transformation/rotation matrix
    Matrix RotationMatrix( MatSize );

    Matrix LocalToGlobalRotationMatrix( 3 );   
    LocalToGlobalRotationMatrix = ZeroMatrix(3);

    (mGlobalToLocalQuaternion.conjugate()).ToRotationMatrix(LocalToGlobalRotationMatrix);

    //Building the rotation matrix for the local element matrix
    for (unsigned int kk=0; kk < MatSize; kk += dimension)
    {
        for (unsigned int i=0; i<dimension; i++)
        {
            for(unsigned int j=0; j<dimension; j++)
            {
	      RotationMatrix(i+kk,j+kk) = LocalToGlobalRotationMatrix(i,j);
            }
        }
    }

    //Rotate Local Stiffness Matrix
    Matrix aux_matrix   = ZeroMatrix(MatSize);
    noalias(aux_matrix) = prod(RotationMatrix, LocalStiffnessMatrix);

    //Stiffness Matrix
    LocalStiffnessMatrix = ZeroMatrix(MatSize);
    noalias(LocalStiffnessMatrix) = prod(aux_matrix,Matrix(trans(RotationMatrix)));

    noalias(rLeftHandSideMatrix) += LocalStiffnessMatrix;


    KRATOS_CATCH( "" )

}


//************************************CALCULATE TOTAL MASS****************************
//************************************************************************************

double& LargeDisplacementBeamElement::CalculateTotalMass( SectionProperties& Section, double& rTotalMass )
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    double Length = GetGeometry().Length();

    rTotalMass = ( Section.Area * Length ) * GetProperties()[DENSITY];

    return rTotalMass;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateInertiaDyadic(SectionProperties & rSection, Matrix & rInertiaDyadic)

{
    KRATOS_TRY

    if( rInertiaDyadic.size1() != 3 )
      rInertiaDyadic.resize(3, 3, false);
    
    rInertiaDyadic = ZeroMatrix(3);

       
    //if the local axis are the principal axis of the cross section

    rInertiaDyadic(0,0) = rSection.Rotational_Inertia; //beam axis
    rInertiaDyadic(1,1) = rSection.Inertia_y; //vertial axis
    rInertiaDyadic(2,2) = rSection.Inertia_z; //horizontal axis


    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamElement::CalculateRotationLinearPartTensor(Vector & rRotationVector, Matrix & rRotationTensor)

{
    KRATOS_TRY

    if( rRotationTensor.size1() != 3 )
      rRotationTensor.resize(3, 3, false);
    
    rRotationTensor = ZeroMatrix(3);

    // adding term 3
    this->VectorToSkewSymmetricTensor( rRotationVector, rRotationTensor );

    rRotationTensor *= -0.5;
    
    double NormRotation =  norm_2(rRotationVector);
    Matrix RotationxRotation = outer_prod( rRotationVector, rRotationVector );
    RotationxRotation *= ( 1.0/( NormRotation * NormRotation ) );

    // adding term 1
    rRotationTensor += RotationxRotation;
    
    double Coefficient = 0.5 * NormRotation / std::tan( 0.5 * NormRotation );
    
    for(unsigned int i=0; i<3; i++)
      RotationxRotation(i,i) -= 1.0;

    RotationxRotation *= (-1)*Coefficient;

    // adding term 2
    rRotationTensor += RotationxRotation;


    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************


void LargeDisplacementBeamElement::CalculateDistributedBodyForce(const int Direction, Vector& Load, Vector& rVolumeForce)
{

    KRATOS_TRY

    Vector Weight = rVariables.Section.Area * rVolumeForce ;

    array_1d<double, 6 > ReferenceCoordinates;
    ReferenceCoordinates(0)= GetGeometry()[0].X0();
    ReferenceCoordinates(1)= GetGeometry()[0].Y0();
    ReferenceCoordinates(2)= GetGeometry()[0].Z0();
    ReferenceCoordinates(3)= GetGeometry()[1].X0();
    ReferenceCoordinates(4)= GetGeometry()[1].Y0();
    ReferenceCoordinates(5)= GetGeometry()[1].Z0();

    Vector LocalBeamAxis;
    LocalBeamAxis.resize(3,false);

    for (unsigned int i=0; i<3; i++)
    {
        LocalBeamAxis[i] = ReferenceCoordinates[i+3] - ReferenceCoordinates[i];
    }

    if( Load.size() != 2 )
      Load.resize(2, false);

    double alpha  =  0.00;
    double sign  =  1.00;

    double sinus;
    double cosinus;

    Vector NormalPlaneProjection;
    NormalPlaneProjection.resize(3,false);

    if(Direction==0)  // 0->x, 1->y, 2->z
    {
	NormalPlaneProjection[0] = 0.00;
	NormalPlaneProjection[1] = LocalBeamAxis[1];
	NormalPlaneProjection[2] = LocalBeamAxis[2];

	if (LocalBeamAxis[0]<0)
	  {
	    sign =-1.00;
	  }
	if( norm_2(NormalPlaneProjection)==0 || norm_2(LocalBeamAxis)==0  )
	  {
	    alpha = sign * PI * 0.5;
	  }
	else
	  {
	    alpha = inner_prod(NormalPlaneProjection,LocalBeamAxis)/(norm_2(LocalBeamAxis)*norm_2(NormalPlaneProjection));
	    alpha = sign*acos(alpha);
	  }

	sinus = sin(alpha);
	cosinus = cos(alpha);

	if(fabs(sinus) < 1E-7) sinus = 0.00;
	if(fabs(cosinus) < 1E-7) cosinus = 0.00;

	// self weight only
	Load[0]= Weight[0]*sinus;         // Axial load
	Load[1]= Weight[0]*cosinus;       // Global Vertical load
    }

    if(Direction==1) // 0->x, 1->y, 2->z
    {
        NormalPlaneProjection = ZeroVector(3);
        NormalPlaneProjection[0] = LocalBeamAxis[0];
        NormalPlaneProjection[1] = 0.00 ;
        NormalPlaneProjection[2] = LocalBeamAxis[2];

        if (LocalBeamAxis[1]<0)
        {
            sign =-1.00;
        }
        if( norm_2(NormalPlaneProjection)==0 || norm_2(LocalBeamAxis)==0  )
        {
            alpha = sign * PI * 0.5;
        }
        else
        {
            alpha = inner_prod(NormalPlaneProjection,LocalBeamAxis)/(norm_2(LocalBeamAxis)*norm_2(NormalPlaneProjection));
            alpha = sign*acos(alpha);
        }

        sinus = sin(alpha);
        cosinus = cos(alpha);

        if(fabs(sinus) < 1E-7) sinus = 0.00;
        if(fabs(cosinus) < 1E-7) cosinus = 0.00;


        // self weight only
        Load[0]= Weight[1]*sinus;         // Axial load
        Load[1]= Weight[1]*cosinus;       // Global Vertical load
    }

    if(Direction==2) // 0->x, 1->y, 2->z
    {
        NormalPlaneProjection[0] = LocalBeamAxis[0] ;
        NormalPlaneProjection[1] = LocalBeamAxis[1] ;
        NormalPlaneProjection[2] = 0.00;

        if (LocalBeamAxis[2]<0)
        {
            sign =-1.00;
        }
        if( norm_2(NormalPlaneProjection)==0 || norm_2(LocalBeamAxis)==0  )
        {
            alpha = sign * PI * 0.5;
        }
        else
        {
            alpha = inner_prod(NormalPlaneProjection,LocalBeamAxis)/(norm_2(LocalBeamAxis)*norm_2(NormalPlaneProjection));
            alpha = sign*acos(alpha);
        }

        sinus = sin(alpha);
        cosinus = cos(alpha);

        if(fabs(sinus) < 1E-7) sinus = 0.00;
        if(fabs(cosinus) < 1E-7) cosinus = 0.00;

        // self weight only
        Load[0]= Weight[2]*sinus;         // Axial load
        Load[1]= Weight[2]*cosinus;       // Global Vertical load
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
   const unsigned int dimension                  = GetGeometry().WorkingSpaceDimension();
   
   if ( rOutput.size() != integration_points_number )
     rOutput.resize( integration_points_number );

    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );
    
    Vector Stress;
    std::vector<Vector> Load(dimension);

    //auxiliary terms
    int factor = 1;

    //(in fact is the two nodes beam element, only one integration point)
    for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
      {
	
	Vector N = row( Ncontainer, PointNumber);
	
	//contribution to external forces
	Vector VolumeForce;
	VolumeForce = this->CalculateVolumeForce( VolumeForce, N );

	this->CalculateLocalNodalStress(Stress, VolumeForce);

	//dangerous:
	for(unsigned int i = 0; i<Stress.size(); i++)
	  {
	    if( std::fabs(Stress[i])< 1E-6) Stress[i] = 0.00;
	  }
	
	
	double x_tolerance     = GetGeometry()[1].X0() - GetGeometry()[0].X0();
	double y_tolerance     = GetGeometry()[1].Y0() - GetGeometry()[0].Y0();
	const double tolerance = 1E-6;

	if(fabs(x_tolerance)>tolerance)
	  {
	    if(GetGeometry()[1].X0() > GetGeometry()[0].X0())
	      factor = 1;
	  }
	else if(fabs(y_tolerance)>tolerance)
	  {
	    if(GetGeometry()[1].Y0() > GetGeometry()[0].Y0())
	      factor = 1;
	  }
	else
	  {
	    factor = 1; //-1;
	  }
	
	for( unsigned int i= 0; i< dimension; i++ )
	  CalculateDistributedBodyForce(i, Load[i], VolumeForce);

	//std::cout<<" VolumeForce "<<VolumeForce<<std::endl;
      }

  
    //only moment in z axis (global ?)
    if(rVariable==MOMENT)
      {
	//internal moment in not using the given load

	/// Punto Inical
	rOutput[0][0] = factor * CalculateInternalMoment(Stress[3], Stress[9], Load[1][1], 1.00 * 0.25);  //Stress[3];
	rOutput[0][1] = factor * CalculateInternalMoment(Stress[4], Stress[10], Load[1][1], 1.00 * 0.25);  //Stress[4];
	rOutput[0][2] = factor * CalculateInternalMoment(Stress[5], Stress[11], Load[1][1], 1.00 * 0.25);  //Stress[5];
        //rOutput[0][2] = factor * CalculateInternalMoment(Stress[5], Stress[1], Load[1][1], rVariables.Length * 0.25); 


        rOutput[1][0] = factor * CalculateInternalMoment(Stress[3], Stress[9], Load[1][1], 1.00 * 0.5);
        rOutput[1][1] = factor * CalculateInternalMoment(Stress[4], Stress[10], Load[1][1], 1.00 * 0.5);
        rOutput[1][2] = factor * CalculateInternalMoment(Stress[5], Stress[11], Load[1][1], 1.00 * 0.5);
        //rOutput[1][2] = factor * CalculateInternalMoment(Stress[5], Stress[1], Load[1][1], rVariables.Length * 0.5);


        rOutput[2][0] = factor * CalculateInternalMoment(Stress[3], Stress[9], Load[1][1], 3.00 * 0.25);
        rOutput[2][1] = factor * CalculateInternalMoment(Stress[4], Stress[10], Load[1][1], 3.00 * 0.25);
        rOutput[2][2] = factor * CalculateInternalMoment(Stress[5], Stress[11], Load[1][1], 3.00 * 0.25);
        //rOutput[2][2] = factor * CalculateInternalMoment(Stress[5], Stress[1], Load[1][1], 3.00 * rVariables.Length * 0.25);

    }

    //only force in x and y axis (global?)
    if(rVariable==FORCE)
    {
        rOutput[0][0] = factor * CalculateInternalAxil(Stress[0], Load[1][0], rVariables.Length * 0.25);
        rOutput[0][1] = factor * CalculateInternalShear(Stress[1], Load[1][1], rVariables.Length * 0.25);
        rOutput[0][2] = 0.00;

        rOutput[1][0] = factor * CalculateInternalAxil(Stress[0], Load[1][0], rVariables.Length * 0.5);
        rOutput[1][1] = factor * CalculateInternalShear(Stress[1], Load[1][1], rVariables.Length * 0.5);
        rOutput[1][2] = 0.00;

        rOutput[2][0] = factor * CalculateInternalAxil(Stress[0], Load[1][0], 3.00 * rVariables.Length * 0.25);
        rOutput[2][1] = factor * CalculateInternalShear(Stress[1], Load[1][1], 3.00 * rVariables.Length * 0.25);
        rOutput[2][2] = 0.00;
    }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************


double LargeDisplacementBeamElement::CalculateInternalMoment(const double& M1, const double& M2, const double& ShearLoad, const double& X)
{
    //return Mo - Vo*X + 0.5 * ShearLoad * X * X;
    return M1*(1-X) - M2*X;
}

//************************************************************************************
//************************************************************************************


double LargeDisplacementBeamElement::CalculateInternalShear(const double& Q, const double& ShearLoad, const double& X)
{
    return  -Q + ShearLoad * X;
}

//************************************************************************************
//************************************************************************************


double LargeDisplacementBeamElement::CalculateInternalAxil(const double& N, const double& AxialLoad, const double& X)
{
    return  -N + AxialLoad * X;
}


//************************************************************************************
//************************************************************************************


void LargeDisplacementBeamElement::CalculateLocalNodalStress(Vector& Stress, Vector& rVolumeForce)
{
    KRATOS_TRY

    array_1d<double, 12 > CurrentDisplacement;
    CurrentDisplacement(0)		=   GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_X);
    CurrentDisplacement(1)		=   GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Y);
    CurrentDisplacement(2)		=   GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Z);
    CurrentDisplacement(3)		=   GetGeometry()[0].GetSolutionStepValue(ROTATION_X);
    CurrentDisplacement(4)		=   GetGeometry()[0].GetSolutionStepValue(ROTATION_Y);
    CurrentDisplacement(5)		=   GetGeometry()[0].GetSolutionStepValue(ROTATION_Z);
    CurrentDisplacement(6)		=   GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_X);
    CurrentDisplacement(7)		=   GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Y);
    CurrentDisplacement(8)		=   GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Z);
    CurrentDisplacement(9)		=   GetGeometry()[1].GetSolutionStepValue(ROTATION_X);
    CurrentDisplacement(10)	        =   GetGeometry()[1].GetSolutionStepValue(ROTATION_Y);
    CurrentDisplacement(11)	        =   GetGeometry()[1].GetSolutionStepValue(ROTATION_Z);

    Matrix Rotation = ZeroMatrix(12);
    this->CalculateTransformationMatrix(Rotation);

    Matrix LocalMatrix = ZeroMatrix(12);
    this->CalculateLocalStiffnessMatrix(LocalMatrix);

    array_1d<double, 12 > LocalDisplacement;
    noalias(LocalDisplacement) = prod(Matrix(trans(Rotation)), CurrentDisplacement);

    Vector LocalForceVector  = ZeroVector(12);
    this->CalculateLocalBodyForce( LocalForceVector, rVolumeForce );

    if( Stress.size() != 12 )
      Stress.resize(12, false);

    // K·u - fext = fint;  where fint is named: Local Stress
    noalias(Stress) = prod(LocalMatrix, LocalDisplacement); 
    Stress -= LocalForceVector; 

    //std::cout<<" Stress "<<Stress<<std::endl;

    return;

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
      KRATOS_ERROR( std::invalid_argument, "This element works only in 3D and with 2 noded linear elements", "")
    }

    //verify that the variables are correctly initialized
    if(VELOCITY.Key() == 0)
        KRATOS_ERROR( std::invalid_argument,"VELOCITY has Key zero! (check if the application is correctly registered", "" )
    if(DISPLACEMENT.Key() == 0)
        KRATOS_ERROR( std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered", "" )
    if(ACCELERATION.Key() == 0)
        KRATOS_ERROR( std::invalid_argument,"ACCELERATION has Key zero! (check if the application is correctly registered", "" )
    if(DENSITY.Key() == 0)
        KRATOS_ERROR( std::invalid_argument,"DENSITY has Key zero! (check if the application is correctly registered", "" )
     if(VOLUME_ACCELERATION.Key() == 0)
        // KRATOS_ERROR( std::invalid_argument,"VOLUME_ACCELERATION has Key zero! (check if the application is correctly registered", "" )
    if(CROSS_AREA.Key() == 0)
        KRATOS_ERROR( std::invalid_argument,"CROSS_AREA has Key zero! (check if the application is correctly registered", "" )
    if(LOCAL_INERTIA.Key() == 0)
        KRATOS_ERROR( std::invalid_argument,"LOCAL_INERTIA has Key zero! (check if the application is correctly registered", "" )
    if(ROTATION.Key() == 0)
        KRATOS_ERROR( std::invalid_argument,"ROTATION has Key zero! (check if the application is correctly registered", "" )

    //verify that the dofs exist
    for(unsigned int i=0; i<this->GetGeometry().size(); i++)
    {
        if(this->GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false)
            KRATOS_ERROR( std::invalid_argument,"missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id() )
        if(this->GetGeometry()[i].HasDofFor(DISPLACEMENT_X) == false || this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false || this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Z) == false)
            KRATOS_ERROR( std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id() )
    }

    //verify that the area is given by properties
    if (this->GetProperties().Has(CROSS_AREA)==false)
    {
        if( GetValue(AREA) == 0.0 )
            KRATOS_ERROR( std::logic_error,"CROSS_AREA not provided for this element", this->Id() )
    }

    //verify that the inertia is given by properties
    if (this->GetProperties().Has(LOCAL_INERTIA)==false)
    {
        if( GetValue(INERTIA)(0,0) == 0.0 )
	  KRATOS_ERROR( std::logic_error,"INERTIA not provided for this element ", this->Id() )
    }


    return 0;

    KRATOS_CATCH( "" )
}




} // Namespace Kratos


