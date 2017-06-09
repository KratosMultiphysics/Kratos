//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_conditions/line_load_3D_condition.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//***********************************************************************************
//***********************************************************************************
LineLoad3DCondition::LineLoad3DCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : ForceLoadCondition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//***********************************************************************************
//***********************************************************************************
LineLoad3DCondition::LineLoad3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : ForceLoadCondition(NewId, pGeometry, pProperties)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
LineLoad3DCondition::LineLoad3DCondition( LineLoad3DCondition const& rOther )
    : ForceLoadCondition(rOther)     
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer LineLoad3DCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new LineLoad3DCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}


//************************************CLONE*******************************************
//************************************************************************************
Condition::Pointer LineLoad3DCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
  LineLoad3DCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

  NewCondition.SetData(this->GetData());
  NewCondition.SetFlags(this->GetFlags());

  //-----------//      
  return Condition::Pointer( new LineLoad3DCondition(NewCondition) );
}


//***********************************************************************************
//***********************************************************************************
LineLoad3DCondition::~LineLoad3DCondition()
{
}

//************* GETTING METHODS

//************************************************************************************
//************************************************************************************

void LineLoad3DCondition::InitializeGeneralVariables(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{

  ForceLoadCondition::InitializeGeneralVariables(rVariables, rCurrentProcessInfo);

  //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
  rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );
  
  //Calculate Delta Position
  //rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);

  //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
  //rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

  //Calculate Total Delta Position
  rVariables.DeltaPosition = CalculateTotalDeltaPosition(rVariables.DeltaPosition);

  //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_0/d£]
  rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );
  
}

//*************************COMPUTE DELTA POSITION*************************************
//************************************************************************************


Matrix& LineLoad3DCondition::CalculateDeltaPosition(Matrix & rDeltaPosition)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    
    rDeltaPosition.resize(number_of_nodes , dimension, false);
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

//*************************COMPUTE TOTAL DELTA POSITION*******************************
//************************************************************************************

Matrix& LineLoad3DCondition::CalculateTotalDeltaPosition(Matrix & rDeltaPosition)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rDeltaPosition.resize(number_of_nodes , dimension, false);
    rDeltaPosition = zero_matrix<double>( number_of_nodes , dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
 
        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rDeltaPosition(i,j) = CurrentDisplacement[j];
        }
    }

    return rDeltaPosition;

    KRATOS_CATCH( "" )
}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************

void LineLoad3DCondition::CalculateKinematics(GeneralVariables& rVariables,
					     const double& rPointNumber)
{
    KRATOS_TRY

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //get first vector of the plane
    rVariables.Tangent1[0] = rVariables.J[rPointNumber](0, 0); // x_1,e
    rVariables.Tangent1[1] = rVariables.J[rPointNumber](1, 0); // x_2,e
    rVariables.Tangent1[2] = rVariables.J[rPointNumber](2, 0); // x_3,e

    //normal in the x-y plane (must be generalized)
    rVariables.Normal[0] = -rVariables.J[rPointNumber](1, 0); //-x_2,e
    rVariables.Normal[1] =  rVariables.J[rPointNumber](0, 0); // x_1,e
    rVariables.Normal[2] =  rVariables.J[rPointNumber](2, 0); // x_3,e

    //Jacobian to the last known configuration
    rVariables.Jacobian = norm_2(rVariables.Tangent1);

    //Set Shape Functions Values for this integration point
    rVariables.N =row( Ncontainer, rPointNumber);

    //Get domain size
    rVariables.DomainSize = GetGeometry().Length();

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

Vector& LineLoad3DCondition::CalculateVectorForce(Vector& rVectorForce, GeneralVariables& rVariables)
{

    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    
    if( rVectorForce.size() != dimension )
      rVectorForce.resize(dimension,false);

    noalias(rVectorForce) = ZeroVector(dimension);
    
    //PRESSURE CONDITION:
    rVectorForce = rVariables.Normal;
    rVariables.Pressure = 0.0;

    //defined on condition
    if( this->Has( NEGATIVE_FACE_PRESSURE ) ){
      double& NegativeFacePressure = this->GetValue( NEGATIVE_FACE_PRESSURE );
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
	rVariables.Pressure += rVariables.N[i] * NegativeFacePressure;
    }

    if( this->Has( POSITIVE_FACE_PRESSURE ) ){
      double& PositiveFacePressure = this->GetValue( POSITIVE_FACE_PRESSURE );
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
	rVariables.Pressure -= rVariables.N[i] * PositiveFacePressure;
    }

    //defined on condition nodes
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      if( GetGeometry()[i].SolutionStepsDataHas( NEGATIVE_FACE_PRESSURE) ) 
	rVariables.Pressure += rVariables.N[i] * ( GetGeometry()[i].FastGetSolutionStepValue( NEGATIVE_FACE_PRESSURE ) );
      if( GetGeometry()[i].SolutionStepsDataHas( POSITIVE_FACE_PRESSURE) ) 
	rVariables.Pressure -= rVariables.N[i] * ( GetGeometry()[i].FastGetSolutionStepValue( POSITIVE_FACE_PRESSURE ) );     
    }
    
    rVectorForce *= rVariables.Pressure;
   
    //FORCE CONDITION:
    
    //defined on condition
    if( this->Has( LINE_LOAD ) ){
      array_1d<double, 3 > & LineLoad = this->GetValue( LINE_LOAD );
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
	{
	  for( unsigned int k = 0; k < dimension; k++ )
	    rVectorForce[k] += rVariables.N[i] * LineLoad[k];
	}
    }

    //defined on condition nodes
    if( this->Has( LINE_LOADS_VECTOR ) ){
      Vector& LineLoads = this->GetValue( LINE_LOADS_VECTOR );
      unsigned int counter = 0;
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
	{
	  counter = i*3;
	  for( unsigned int k = 0; k < dimension; k++ )
	    {
	      rVectorForce[k] += rVariables.N[i] * LineLoads[counter+k];
	    }
	  
	}
    }
    
    //defined on geometry nodes
    for (unsigned int i = 0; i < number_of_nodes; i++)
      {
	if( GetGeometry()[i].SolutionStepsDataHas( LINE_LOAD ) ){
	  array_1d<double, 3 > & LineLoad = GetGeometry()[i].FastGetSolutionStepValue( LINE_LOAD );
	  for( unsigned int k = 0; k < dimension; k++ )
	    rVectorForce[k] += rVariables.N[i] * LineLoad[k];
	}
      }

	//follower forces
	double testFollower = 0.00;
	for (int i = 0; i < number_of_nodes; ++i)
	{
		if (GetGeometry()[i].SolutionStepsDataHas(FOLLOWER_LINE_LOAD)) {
			array_1d<double, 3 > & LineLoad = GetGeometry()[i].FastGetSolutionStepValue(FOLLOWER_LINE_LOAD);
			const double checkFollower = MathUtils<double>::Norm(LineLoad);
			testFollower += checkFollower;
			if (checkFollower != 0.00)
			{
				for (unsigned int k = 0; k < dimension; k++)
				{
					rVectorForce[k] += rVariables.N[i] * LineLoad[k];
				}
			}
		}
	}
	if (testFollower != 0 && number_of_nodes == 2 && dimension == 3)
	{
		this->CalculateFollowerForceDirection(rVectorForce);
	}


    return rVectorForce;

    KRATOS_CATCH( "" )
}


//calculate load direction for Follower Force

void LineLoad3DCondition::CalculateFollowerForceDirection(Vector& rVectorForce)
{
	KRATOS_TRY;
	const int dimension = GetGeometry().WorkingSpaceDimension();
	const int number_of_nodes = 2;

	Vector LineLoadDir = ZeroVector(dimension);
	for (int i = 0; i < dimension; ++i) {
		LineLoadDir[i] = rVectorForce[i];
	}

	const double VectorNormLoad = MathUtils<double>::Norm(LineLoadDir);
	if (VectorNormLoad != 0.00) LineLoadDir /= VectorNormLoad;

	//////////////////////////////////////////////////////
	//initial CS
	//////////////////////////////////////////////////////
	const int size = number_of_nodes * dimension;
	const int local_size = size * 2;

	if (this->mIsFirstStep == true)
	{
		array_1d<double, 3> DirectionVectorX = ZeroVector(dimension);
		array_1d<double, 3> DirectionVectorY = ZeroVector(dimension);
		array_1d<double, 3> DirectionVectorZ = ZeroVector(dimension);
		Vector ReferenceCoordinates = ZeroVector(size);

		ReferenceCoordinates[0] = this->GetGeometry()[0].X0();
		ReferenceCoordinates[1] = this->GetGeometry()[0].Y0();
		ReferenceCoordinates[2] = this->GetGeometry()[0].Z0();
		ReferenceCoordinates[3] = this->GetGeometry()[1].X0();
		ReferenceCoordinates[4] = this->GetGeometry()[1].Y0();
		ReferenceCoordinates[5] = this->GetGeometry()[1].Z0();

		for (int i = 0; i < dimension; ++i)
		{
			DirectionVectorX[i] = (ReferenceCoordinates[i + dimension]
				- ReferenceCoordinates[i]);
		}

		//use orientation class 1st constructor
		OrientationGeometry element_axis(DirectionVectorX, 0.00);
		element_axis.CalculateBasisVectors(DirectionVectorX, DirectionVectorY,
			DirectionVectorZ);

		this->mQuaternionVEC_A = ZeroVector(dimension);
		this->mQuaternionVEC_B = ZeroVector(dimension);
		this->mQuaternionSCA_A = 1.00;
		this->mQuaternionSCA_B = 1.00;


		this->mDirectionVectorXOriginal = ZeroVector(dimension);
		this->mDirectionVectorYOriginal = ZeroVector(dimension);
		this->mDirectionVectorZOriginal = ZeroVector(dimension);

		this->mDirectionVectorXOriginal = DirectionVectorX;
		this->mDirectionVectorYOriginal = DirectionVectorY;
		this->mDirectionVectorZOriginal = DirectionVectorZ;


		this->mIsFirstStep = false;
	}


	//////////////////////////////////////////////////////
	//update CS
	//////////////////////////////////////////////////////


	Vector ActualDeformation = ZeroVector(local_size);
	Vector LastStepDeformation = ZeroVector(local_size);
	Vector IncrementDeformation = ZeroVector(local_size);
	this->GetValuesVector(ActualDeformation, 0);
	this->GetValuesVector(LastStepDeformation, 1);
	IncrementDeformation = ActualDeformation - LastStepDeformation;

	//test
	//IncrementDeformation = ZeroVector(local)

	Vector dPhiA = ZeroVector(dimension);
	Vector dPhiB = ZeroVector(dimension);

	for (int i = 0; i < dimension; ++i) {
		dPhiA[i] = IncrementDeformation[i + 3];
		dPhiB[i] = IncrementDeformation[i + 9];
	}

	//calculating quaternions
	Vector drA_vec = ZeroVector(dimension);
	Vector drB_vec = ZeroVector(dimension);
	double drA_sca, drB_sca;

	drA_vec = 0.50 * dPhiA;
	drB_vec = 0.50 * dPhiB;

	drA_sca = 0.00;
	drB_sca = 0.00;
	for (int i = 0; i < dimension; ++i) {
		drA_sca += drA_vec[i] * drA_vec[i];
		drB_sca += drB_vec[i] * drB_vec[i];
	}
	drA_sca = 1.00 - drA_sca;
	drB_sca = 1.00 - drB_sca;

	drA_sca = sqrt(drA_sca);
	drB_sca = sqrt(drB_sca);


	Vector tempVec = ZeroVector(dimension);
	double tempSca = 0.00;

	//Node A
	tempVec = this->mQuaternionVEC_A;
	tempSca = this->mQuaternionSCA_A;

	this->mQuaternionSCA_A = drA_sca *tempSca;
	for (int i = 0; i < dimension; ++i) {
		this->mQuaternionSCA_A -= drA_vec[i] * tempVec[i];
	}
	this->mQuaternionVEC_A = drA_sca*tempVec;
	this->mQuaternionVEC_A += tempSca * drA_vec;
	this->mQuaternionVEC_A += MathUtils<double>::CrossProduct(drA_vec, tempVec);

	//Node B
	tempVec = this->mQuaternionVEC_B;
	tempSca = this->mQuaternionSCA_B;

	this->mQuaternionSCA_B = drB_sca *tempSca;
	for (int i = 0; i < dimension; ++i) {
		this->mQuaternionSCA_B -= drB_vec[i] * tempVec[i];
	}

	this->mQuaternionVEC_B = drB_sca*tempVec;
	this->mQuaternionVEC_B += tempSca * drB_vec;
	this->mQuaternionVEC_B += MathUtils<double>::CrossProduct(drB_vec, tempVec);


	//scalar part of difference quaternion
	double scalar_diff;
	scalar_diff = (this->mQuaternionSCA_A + this->mQuaternionSCA_B) *
		(this->mQuaternionSCA_A + this->mQuaternionSCA_B);

	tempVec = this->mQuaternionVEC_A + this->mQuaternionVEC_B;
	scalar_diff += MathUtils<double>::Norm(tempVec) *
		MathUtils<double>::Norm(tempVec);

	scalar_diff = 0.50 * sqrt(scalar_diff);

	//mean rotation quaternion
	double meanRotationScalar;
	meanRotationScalar = (this->mQuaternionSCA_A + this->mQuaternionSCA_B) * 0.50;
	meanRotationScalar = meanRotationScalar / scalar_diff;

	Vector meanRotationVector = ZeroVector(dimension);
	meanRotationVector = (this->mQuaternionVEC_A + this->mQuaternionVEC_B) * 0.50;
	meanRotationVector = meanRotationVector / scalar_diff;

	//vector part of difference quaternion
	Vector vector_diff = ZeroVector(dimension);
	vector_diff = this->mQuaternionSCA_A * this->mQuaternionVEC_B;
	vector_diff -= this->mQuaternionSCA_B * this->mQuaternionVEC_A;
	vector_diff += MathUtils<double>::CrossProduct(this->mQuaternionVEC_A,
		this->mQuaternionVEC_B);

	vector_diff = 0.50 * vector_diff / scalar_diff;

	//rotate inital element basis
	const double r0 = meanRotationScalar;
	const double r1 = meanRotationVector[0];
	const double r2 = meanRotationVector[1];
	const double r3 = meanRotationVector[2];

	Quaternion<double> q(r0, r1, r2, r3);
	Vector rotatedNX0 = this->mDirectionVectorXOriginal;
	Vector rotatedNY0 = this->mDirectionVectorYOriginal;
	Vector rotatedNZ0 = this->mDirectionVectorZOriginal;
	q.RotateVector3(rotatedNX0);
	q.RotateVector3(rotatedNY0);
	q.RotateVector3(rotatedNZ0);


	Matrix RotatedCS = ZeroMatrix(dimension, dimension);
	for (int i = 0; i < dimension; ++i) {
		RotatedCS(i, 0) = rotatedNX0[i];
		RotatedCS(i, 1) = rotatedNY0[i];
		RotatedCS(i, 2) = rotatedNZ0[i];
	}

	//rotate basis to element axis + redefine R
	Vector n_bisectrix = ZeroVector(dimension);
	Vector deltaX = ZeroVector(dimension);
	double VectorNorm;

	deltaX[0] = this->GetGeometry()[1].X0() + ActualDeformation[6] -
		(this->GetGeometry()[0].X0() + ActualDeformation[0]);
	deltaX[1] = this->GetGeometry()[1].Y0() + ActualDeformation[7] -
		(this->GetGeometry()[0].Y0() + ActualDeformation[1]);
	deltaX[2] = this->GetGeometry()[1].Z0() + ActualDeformation[8] -
		(this->GetGeometry()[0].Z0() + ActualDeformation[2]);

	VectorNorm = MathUtils<double>::Norm(deltaX);
	deltaX /= VectorNorm;

	n_bisectrix = rotatedNX0 + deltaX;
	VectorNorm = MathUtils<double>::Norm(n_bisectrix);
	n_bisectrix /= VectorNorm;

	Matrix n_xyz = ZeroMatrix(dimension);
	for (int i = 0; i < dimension; ++i) {
		n_xyz(i, 0) = -1.0 * RotatedCS(i, 0);
		n_xyz(i, 1) = 1.0 * RotatedCS(i, 1);
		n_xyz(i, 2) = 1.0 * RotatedCS(i, 2);
	}

	Matrix Identity = ZeroMatrix(dimension);
	for (int i = 0; i < dimension; ++i) Identity(i, i) = 1.0;
	Identity -= 2.0 * outer_prod(n_bisectrix, n_bisectrix);
	n_xyz = prod(Identity, n_xyz);

	rVectorForce = prod(n_xyz, rVectorForce);
	KRATOS_CATCH("");
}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

int LineLoad3DCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}

//***********************************************************************************
//***********************************************************************************

void LineLoad3DCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ForceLoadCondition )
}

void LineLoad3DCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ForceLoadCondition )
}


OrientationGeometry::OrientationGeometry(array_1d<double, 3>& v1, const double theta) {

	KRATOS_TRY
		//!!!!!!!!!! if crossproduct with array_1d type switch input order !!!!!!!
		//If only direction of v1 is given -> Default case
		const int number_of_nodes = 2;
	const int dimension = 3;
	const int size = number_of_nodes * dimension;
	const int MatSize = 2 * size;

	array_1d<double, 3> GlobalZ = ZeroVector(dimension);
	GlobalZ[2] = 1.0;

	array_1d<double, 3> v2 = ZeroVector(dimension);
	array_1d<double, 3> v3 = ZeroVector(dimension);

	double VectorNorm;
	VectorNorm = MathUtils<double>::Norm(v1);
	if (VectorNorm != 0) v1 /= VectorNorm;

	if (v1[2] == 1.00) {
		v2[1] = 1.0;
		v3[0] = -1.0;
	}

	if (v1[2] == -1.00) {
		v2[1] = 1.0;
		v3[0] = 1.0;
	}

	if (fabs(v1[2]) != 1.00) {

		v2 = MathUtils<double>::CrossProduct(v1, GlobalZ);
		VectorNorm = MathUtils<double>::Norm(v2);
		if (VectorNorm != 0) v2 /= VectorNorm;

		v3 = MathUtils<double>::CrossProduct(v2, v1);
		VectorNorm = MathUtils<double>::Norm(v3);
		if (VectorNorm != 0) v3 /= VectorNorm;
	}

	//manual rotation around the beam axis
	if (theta != 0) {
		const Vector nz_temp = v3;
		const Vector ny_temp = v2;
		const double CosTheta = cos(theta);
		const double SinTheta = sin(theta);

		v2 = ny_temp * CosTheta + nz_temp * SinTheta;
		VectorNorm = MathUtils<double>::Norm(v2);
		if (VectorNorm != 0) v2 /= VectorNorm;

		v3 = nz_temp * CosTheta - ny_temp * SinTheta;
		VectorNorm = MathUtils<double>::Norm(v3);
		if (VectorNorm != 0) v3 /= VectorNorm;
	}

	Matrix RotationMatrix = ZeroMatrix(dimension);
	for (int i = 0; i < dimension; ++i) {
		RotationMatrix(i, 0) = v1[i];
		RotationMatrix(i, 1) = v2[i];
		RotationMatrix(i, 2) = v3[i];
	}

	this->GetQuaternion() = Quaternion<double>::FromRotationMatrix(RotationMatrix);

	KRATOS_CATCH("")
}

OrientationGeometry::OrientationGeometry(array_1d<double, 3>& v1, array_1d<double, 3>& v2) {

	KRATOS_TRY
		//If the user defines an aditional direction v2
		const int number_of_nodes = 2;
	const int dimension = 3;
	const int size = number_of_nodes * dimension;
	const int MatSize = 2 * size;

	array_1d<double, 3> v3 = ZeroVector(dimension);

	double VectorNorm;
	VectorNorm = MathUtils<double>::Norm(v1);
	if (VectorNorm != 0) v1 /= VectorNorm;

	VectorNorm = MathUtils<double>::Norm(v2);
	if (VectorNorm != 0) v2 /= VectorNorm;

	v3 = MathUtils<double>::CrossProduct(v2, v1);
	VectorNorm = MathUtils<double>::Norm(v3);
	if (VectorNorm != 0) v3 /= VectorNorm;


	Matrix RotationMatrix = ZeroMatrix(dimension);
	for (int i = 0; i < dimension; ++i) {
		RotationMatrix(i, 0) = v1[i];
		RotationMatrix(i, 1) = v2[i];
		RotationMatrix(i, 2) = v3[i];
	}

	this->GetQuaternion() = Quaternion<double>::FromRotationMatrix(RotationMatrix);

	KRATOS_CATCH("")
}

void OrientationGeometry::CalculateRotationMatrix(bounded_matrix<double, 3, 3>& R) {

	KRATOS_TRY
		if (R.size1() != 3 || R.size2() != 3) R.resize(3, 3, false);
	const Quaternion<double> q = this->GetQuaternion();
	q.ToRotationMatrix(R);
	KRATOS_CATCH("")
}

void OrientationGeometry::CalculateBasisVectors(array_1d<double, 3>& v1,
	array_1d<double, 3>& v2,
	array_1d<double, 3>& v3) {

	KRATOS_TRY
		const Quaternion<double> q = this->GetQuaternion();
	Matrix R = ZeroMatrix(3);
	q.ToRotationMatrix(R);
	if (v1.size() != 3) v1.resize(3, false);
	if (v2.size() != 3) v2.resize(3, false);
	if (v3.size() != 3) v3.resize(3, false);

	for (int i = 0; i < 3; ++i) {
		v1[i] = R(i, 0);
		v2[i] = R(i, 1);
		v3[i] = R(i, 2);
	}
	KRATOS_CATCH("")
}


} // Namespace Kratos.
