// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Massimo Petracca
//

#include "custom_elements/rocco_net_element_3D4N.hpp"
#include "custom_utilities/shellq4_corotational_coordinate_transformation.hpp"
#include "structural_mechanics_application_variables.h"
#include "includes/define.h"

#include "geometries/quadrilateral_3d_4.h"

#include <string>
#include <iomanip>

namespace Kratos
{

//SYSTEM CALCULATION
void RoccoNetElement3D4N::Initialize()
{
	KRATOS_TRY
	this->mKb = this->GetProperties()[KB_ROCCO];
	this->mNumberWindings = this->GetProperties()[WIRE_WINDINGS_ROCCO];
	this->mThicknessWire = this->GetProperties()[WIRE_THICKNESS_ROCCO];
	this->mDiameter = this->GetProperties()[DIAMETER_ROCCO];
	this->mNodalMass_custom = this->GetProperties()[NODAL_MASS_ROCCO];
	const int Kt_multiplyer = this->GetProperties()[KT_MULT_ROCCO];


	//this is calculated from input variables
	this->mKt = this->mKb * Kt_multiplyer;
	this->mCircumference = this->mDiameter * PI;
	this->mLengthDiagonalReference = this->mDiameter -
		(this->mThicknessWire*sqrt(this->mNumberWindings));
    KRATOS_CATCH("")
}
void RoccoNetElement3D4N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY
	const int number_of_nodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int element_size = number_of_nodes * dimension;


	this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
	this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
	KRATOS_CATCH("")
}
void RoccoNetElement3D4N::CalculateRightHandSide(VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY
	const int number_of_nodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int element_size = number_of_nodes * dimension;

	rRightHandSideVector = ZeroVector(element_size);
	//RHS
	Vector TempInternalForces = ZeroVector(element_size);
	this->AddDiagonalRHSForces(TempInternalForces);
	rRightHandSideVector -= TempInternalForces;

	TempInternalForces = ZeroVector(element_size);
	this->AddCableRHSForces(TempInternalForces);
	rRightHandSideVector -= TempInternalForces;
	KRATOS_CATCH("")
}
void RoccoNetElement3D4N::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
	ProcessInfo& rCurrentProcessInfo) 
{
	KRATOS_TRY
	const int NumNodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int LocalSize = NumNodes * dimension;

	rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);

	//LHS
	Matrix TempStiffness = ZeroMatrix(LocalSize, LocalSize);
	this->AddDiagonalStiffnessContribution(TempStiffness);
	rLeftHandSideMatrix += TempStiffness;

	TempStiffness = ZeroMatrix(LocalSize, LocalSize);
	this->AddCableStiffnessContribution(TempStiffness);
	rLeftHandSideMatrix += TempStiffness;

	//save LHS to pass to CalculateLHS()
	this->mLHS = ZeroMatrix(LocalSize, LocalSize);
	this->mLHS = rLeftHandSideMatrix;
	KRATOS_CATCH("")
}



//RHS CALCUALTION
void RoccoNetElement3D4N::AddDiagonalRHSForces(Vector& rInternalForces)
{
	KRATOS_TRY;
	const int number_of_nodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int element_size = number_of_nodes * dimension;
	const int spring_size = dimension * 2;


	rInternalForces = ZeroVector(element_size);
	//use distance of input nodes as ref. diagonal to get u
	//assume same refLength for both diagonals
	const double ReferenceDiagonalLength = this->CalculateDiagonalLengthRefNodes(0, 2);
	const double LengthDiagonal1 = this->CalculateActualDiagonalLength(0, 2);
	const double LengthDiagonal2 = this->CalculateActualDiagonalLength(1, 3);

	const double u1 = LengthDiagonal1 - ReferenceDiagonalLength;
	const double u2 = LengthDiagonal2 - ReferenceDiagonalLength;
	const double StiffnessKb = this->mKb;

	//check if NormalResistance is reached
	const double ReferenceCableLength = this->mCircumference;
	const double ActualCableLength = this->CalculateActualCableLength();

	if (ActualCableLength <= ReferenceCableLength &&
		this->mNormalResistanceReached == false)
	{
		this->mdN1 = u1;
		this->mdN2 = u2;
	}
	else
	{
		this->mNormalResistanceReached = true; // true
	}

	//create diagonal RHS local   
	double N1, N2;
	// Diagonal 1
	if (u1 <= 0.00) N1 = 0.00;
	else if (u1 > 0.00 && mNormalResistanceReached == false) N1 = StiffnessKb * u1;
	else N1 = StiffnessKb * this->mdN1;

	// Diagonal 1
	if (u2 <= 0.00) N2 = 0.00;
	else if (u2 > 0.00 && mNormalResistanceReached == false) N2 = StiffnessKb * u2;
	else N2 = StiffnessKb * this->mdN2;


	Vector InnerForces1 = ZeroVector(spring_size);
	Vector InnerForces2 = ZeroVector(spring_size);
	InnerForces1[0] = -1.00 * N1;
	InnerForces1[3] = 1.00 * N1;
	InnerForces2[0] = -1.00 * N2;
	InnerForces2[3] = 1.00 * N2;


	//create diagonal RHS global
	Matrix TempRotationmatrix = ZeroMatrix(spring_size, spring_size);
	this->CreateTransformationMatrix(TempRotationmatrix, 0, 2);
	InnerForces1 = prod(TempRotationmatrix, InnerForces1);

	TempRotationmatrix = ZeroMatrix(spring_size, spring_size);
	this->CreateTransformationMatrix(TempRotationmatrix, 1, 3);
	InnerForces2 = prod(TempRotationmatrix, InnerForces2);


	for (int i = 0; i < dimension; ++i)
	{
		rInternalForces[i] += InnerForces1[i];
		rInternalForces[i + spring_size] += InnerForces1[i + dimension];

		rInternalForces[i + dimension] = InnerForces2[i];
		rInternalForces[i + dimension + spring_size] = InnerForces2[i + dimension];
	}

	KRATOS_CATCH("");
}
void RoccoNetElement3D4N::AddCableRHSForces(Vector& rInternalForces)
{
	KRATOS_TRY;
	const int number_of_nodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int element_size = number_of_nodes * dimension;
	const int spring_size = dimension * 2;


	rInternalForces = ZeroVector(element_size);
	//With respect to AXEL Volkwein assign the same N to every cable part

	const double ReferenceCableLength = this->mCircumference;
	const double ActualCableLength = this->CalculateActualCableLength();
	const double U_cable = ActualCableLength - ReferenceCableLength;
	const double StiffnessKt = this->mKt;



	double N_cable = 0.00;
	if (U_cable > 0.00) N_cable = StiffnessKt * U_cable;

	Vector Cable_Inner_Forces = ZeroVector(spring_size);
	Vector Temp_Inner_Forces = ZeroVector(spring_size);
	Cable_Inner_Forces[0] = -1.00 * N_cable;
	Cable_Inner_Forces[3] = 1.00 * N_cable;

	//cable element 1
	Matrix TempRotationMatrix = ZeroMatrix(spring_size, spring_size);
	this->CreateTransformationMatrix(TempRotationMatrix, 0, 1);
	Temp_Inner_Forces = prod(TempRotationMatrix, Cable_Inner_Forces);
	for (int i = 0; i < spring_size; ++i)
	{
		rInternalForces[i] += Temp_Inner_Forces[i];
	}

	//cable element 2
	TempRotationMatrix = ZeroMatrix(spring_size, spring_size);
	this->CreateTransformationMatrix(TempRotationMatrix, 1, 2);
	Temp_Inner_Forces = prod(TempRotationMatrix, Cable_Inner_Forces);
	for (int i = 0; i < spring_size; ++i)
	{
		rInternalForces[i + dimension] += Temp_Inner_Forces[i];
	}

	//cable element 3
	TempRotationMatrix = ZeroMatrix(spring_size, spring_size);
	this->CreateTransformationMatrix(TempRotationMatrix, 2, 3);
	Temp_Inner_Forces = prod(TempRotationMatrix, Cable_Inner_Forces);
	for (int i = 0; i < spring_size; ++i)
	{
		rInternalForces[i + spring_size] += Temp_Inner_Forces[i];
	}

	//cable element 4
	TempRotationMatrix = ZeroMatrix(spring_size, spring_size);
	this->CreateTransformationMatrix(TempRotationMatrix, 3, 0);
	Temp_Inner_Forces = prod(TempRotationMatrix, Cable_Inner_Forces);
	for (int i = 0; i < dimension; ++i)
	{
		rInternalForces[i + spring_size + dimension] += Temp_Inner_Forces[i];
		rInternalForces[i] += Temp_Inner_Forces[i+dimension];
	}

	KRATOS_CATCH("");
}









//LENGTH CALCULATION
double RoccoNetElement3D4N::CalculateActualCableLength()
{
	const int number_of_nodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int element_size = number_of_nodes * dimension;

	Vector VectorCoordinates = ZeroVector(element_size);
	for (int i = 0; i < number_of_nodes; ++i)
	{
		int index = i * dimension;
		VectorCoordinates[index] = this->GetGeometry()[i].X0();
		VectorCoordinates[index + 1] = this->GetGeometry()[i].Y0();
		VectorCoordinates[index + 2] = this->GetGeometry()[i].Z0();
	}
	Vector ActualDeformation = ZeroVector(element_size);
	this->GetValuesVector(ActualDeformation, 0);
	VectorCoordinates += ActualDeformation;

	Vector N1 = ZeroVector(dimension);
	Vector N2 = ZeroVector(dimension);
	Vector N3 = ZeroVector(dimension);
	Vector N4 = ZeroVector(dimension);

	for (int i = 0; i < dimension; ++i) N1[i] = VectorCoordinates[i];
	for (int i = 0; i < dimension; ++i) N2[i] = VectorCoordinates[i+dimension];
	for (int i = 0; i < dimension; ++i) N3[i] = VectorCoordinates[i+(2*dimension)];
	for (int i = 0; i < dimension; ++i) N4[i] = VectorCoordinates[i+(3*dimension)];

	Vector L1 = ZeroVector(dimension);
	Vector L2 = ZeroVector(dimension);
	Vector L3 = ZeroVector(dimension);
	Vector L4 = ZeroVector(dimension);

	L1 = N2 - N1;
	L2 = N3 - N2;
	L3 = N4 - N3;
	L4 = N1 - N4;

	const double norm1 = MathUtils<double>::Norm(L1);
	const double norm2 = MathUtils<double>::Norm(L2);
	const double norm3 = MathUtils<double>::Norm(L3);
	const double norm4 = MathUtils<double>::Norm(L4);

	const double NORM = norm1 + norm2 + norm3 + norm4;
	return NORM;
}
double RoccoNetElement3D4N::CalculateActualDiagonalLength(const int node1,
	const int node2)
{
	const int dimension = this->GetGeometry().WorkingSpaceDimension();


	Vector Coordinates = ZeroVector(dimension);
	const double du = this->GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT_X, 0) -
		this->GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT_X, 0);
	const double dv = this->GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT_Y, 0) -
		this->GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT_Y, 0);
	const double dw = this->GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT_Z, 0) -
		this->GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT_Z, 0);
	const double dx = this->GetGeometry()[node2].X0() - this->GetGeometry()[node1].X0();
	const double dy = this->GetGeometry()[node2].Y0() - this->GetGeometry()[node1].Y0();
	const double dz = this->GetGeometry()[node2].Z0() - this->GetGeometry()[node1].Z0();

	const double length = sqrt((du + dx)*(du + dx) + (dv + dy)*(dv + dy) + (dw + dz)*(dw + dz));
	return length;

}
double RoccoNetElement3D4N::CalculateDiagonalLengthRefNodes(const int node1, const int node2)
{
	const int dimension = this->GetGeometry().WorkingSpaceDimension();


	Vector Coordinates = ZeroVector(dimension);
	const double dx = this->GetGeometry()[node2].X0() - this->GetGeometry()[node1].X0();
	const double dy = this->GetGeometry()[node2].Y0() - this->GetGeometry()[node1].Y0();
	const double dz = this->GetGeometry()[node2].Z0() - this->GetGeometry()[node1].Z0();

	const double length = sqrt((dx * dx) + (dy * dy) + (dz * dz));
	return length;

}
double RoccoNetElement3D4N::CalculateCableLengthRefNodes()
{
	const int number_of_nodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int element_size = number_of_nodes * dimension;

	Vector VectorCoordinates = ZeroVector(element_size);
	for (int i = 0; i < number_of_nodes; ++i)
	{
		int index = i * dimension;
		VectorCoordinates[index] = this->GetGeometry()[i].X0();
		VectorCoordinates[index + 1] = this->GetGeometry()[i].Y0();
		VectorCoordinates[index + 2] = this->GetGeometry()[i].Z0();
	}

	Vector N1 = ZeroVector(dimension);
	Vector N2 = ZeroVector(dimension);
	Vector N3 = ZeroVector(dimension);
	Vector N4 = ZeroVector(dimension);

	for (int i = 0; i < dimension; ++i) N1[i] = VectorCoordinates[i];
	for (int i = 0; i < dimension; ++i) N2[i] = VectorCoordinates[i + dimension];
	for (int i = 0; i < dimension; ++i) N3[i] = VectorCoordinates[i + (2 * dimension)];
	for (int i = 0; i < dimension; ++i) N4[i] = VectorCoordinates[i + (3 * dimension)];

	Vector L1 = ZeroVector(dimension);
	Vector L2 = ZeroVector(dimension);
	Vector L3 = ZeroVector(dimension);
	Vector L4 = ZeroVector(dimension);

	L1 = N2 - N1;
	L2 = N3 - N2;
	L3 = N4 - N3;
	L4 = N1 - N4;

	const double norm1 = MathUtils<double>::Norm(L1);
	const double norm2 = MathUtils<double>::Norm(L2);
	const double norm3 = MathUtils<double>::Norm(L3);
	const double norm4 = MathUtils<double>::Norm(L4);

	const double NORM = norm1 + norm2 + norm3 + norm4;
	return NORM;
}


//STIFFNESS CALCULATION
void RoccoNetElement3D4N::AddDiagonalStiffnessContribution(Matrix& rLeftHandSideMatrix)
{
	KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int element_size = number_of_nodes * dimension;
	const int spring_size = dimension * 2;

	rLeftHandSideMatrix = ZeroMatrix(element_size, element_size);

	Matrix DraftDiagonalK = ZeroMatrix(spring_size, spring_size);
	Matrix TempDiagonalK = ZeroMatrix(spring_size, spring_size);
	Matrix TempTransformationMatrix = ZeroMatrix(spring_size, spring_size);
	Matrix AuxRotationMatrix = ZeroMatrix(spring_size, spring_size);
	Matrix TempLeftHandSide = ZeroMatrix(element_size, element_size);
	double StiffnessKb; //enhance function for this
	DraftDiagonalK(0, 0) = 1.00;
	DraftDiagonalK(0, 3) = -1.00;
	DraftDiagonalK(3, 0) = -1.00;
	DraftDiagonalK(3, 3) = 1.00;


	//Diagonal Node 1 -> 3
	StiffnessKb = this->CalculateDiagonalStiffnessKb(0, 2); //enhance function for this
	TempDiagonalK = DraftDiagonalK * StiffnessKb;
	this->CreateTransformationMatrix(TempTransformationMatrix, 0, 2);
	AuxRotationMatrix = prod(TempTransformationMatrix, TempDiagonalK);
	TempDiagonalK = prod(AuxRotationMatrix,
		Matrix(trans(TempTransformationMatrix)));

	for (int i = 0; i < dimension; ++i)
	{
		for (int j = 0; j < dimension; ++j)
		{
			TempLeftHandSide(i, j) += TempDiagonalK(i, j);
			TempLeftHandSide(i, j + spring_size) +=
				TempDiagonalK(i, j + dimension);
			TempLeftHandSide(i + spring_size, j) +=
				TempDiagonalK(i + dimension, j);
			TempLeftHandSide(i + spring_size, j + spring_size) +=
				TempDiagonalK(i + dimension, j + dimension);
		}
	}
	rLeftHandSideMatrix += TempLeftHandSide;

	//Diagonal Node 2 -> 4
	StiffnessKb = this->CalculateDiagonalStiffnessKb(1, 3); //enhance function for this
	TempDiagonalK = ZeroMatrix(spring_size, spring_size);
	TempTransformationMatrix = ZeroMatrix(spring_size, spring_size);
	AuxRotationMatrix = ZeroMatrix(spring_size, spring_size);
	TempLeftHandSide = ZeroMatrix(element_size, element_size);

	TempDiagonalK = DraftDiagonalK * StiffnessKb;
	this->CreateTransformationMatrix(TempTransformationMatrix, 1, 3);
	AuxRotationMatrix = prod(TempTransformationMatrix, TempDiagonalK);
	TempDiagonalK = prod(AuxRotationMatrix,
		Matrix(trans(TempTransformationMatrix)));

	for (int i = 0; i < dimension; ++i)
	{
		int ii = i + dimension;
		for (int j = 0; j < dimension; ++j)
		{
			int jj = j + dimension;
			TempLeftHandSide(ii, jj) += TempDiagonalK(i, j);
			TempLeftHandSide(ii, jj + spring_size) +=
				TempDiagonalK(i, j + dimension);
			TempLeftHandSide(ii + spring_size, jj) +=
				TempDiagonalK(i + dimension, j);
			TempLeftHandSide(ii + spring_size, jj + spring_size) +=
				TempDiagonalK(i + dimension, j + dimension);
		}
	}
	rLeftHandSideMatrix += TempLeftHandSide;
	KRATOS_CATCH("")
}
void RoccoNetElement3D4N::AddCableStiffnessContribution(Matrix& rLeftHandSideMatrix)
{
	KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int element_size = number_of_nodes * dimension;
	const int spring_size = dimension * 2;

	const double ReferenceCableLength = this->mCircumference;
	const double ActualCableLength = this->CalculateActualCableLength();

	rLeftHandSideMatrix = ZeroMatrix(element_size, element_size);

	double StiffnessKt = this->mKt*4.00;

	if (ReferenceCableLength > ActualCableLength && this->mCableResistanceReached == false)
	{
		this->mCableResistanceReached = false;
	}
	else this->mCableResistanceReached = true;

	if (this->mCableResistanceReached == true)
	{
		Matrix TempDiagonalK = ZeroMatrix(spring_size, spring_size);
		Matrix TempTransformationMatrix = ZeroMatrix(spring_size, spring_size);
		Matrix AuxRotationMatrix = ZeroMatrix(spring_size, spring_size);
		Matrix DraftDiagonalK = ZeroMatrix(spring_size, spring_size);
		Matrix TempLeftHandSide = ZeroMatrix(element_size, element_size);

		DraftDiagonalK(0, 0) = 1.00 * StiffnessKt;
		DraftDiagonalK(0, 3) = -1.00 * StiffnessKt;
		DraftDiagonalK(3, 0) = -1.00 * StiffnessKt;
		DraftDiagonalK(3, 3) = 1.00 * StiffnessKt;

		// Kt relates now to the respective displacement not the total delta L....

		//node 1 -> 2
		this->CreateTransformationMatrix(TempTransformationMatrix, 0, 1);
		AuxRotationMatrix = prod(TempTransformationMatrix, DraftDiagonalK);
		TempDiagonalK = prod(AuxRotationMatrix,
			Matrix(trans(TempTransformationMatrix)));

		for (int i = 0; i < spring_size; ++i)
		{
			for (int j = 0; j < spring_size; ++j)
			{
				TempLeftHandSide(i, j) += TempDiagonalK(i, j);
			}
		}
		rLeftHandSideMatrix += TempLeftHandSide;


		//node 2 -> 3
		TempDiagonalK = ZeroMatrix(spring_size, spring_size);
		TempTransformationMatrix = ZeroMatrix(spring_size, spring_size);
		AuxRotationMatrix = ZeroMatrix(spring_size, spring_size);
		TempLeftHandSide = ZeroMatrix(element_size, element_size);

		this->CreateTransformationMatrix(TempTransformationMatrix, 1, 2);
		AuxRotationMatrix = prod(TempTransformationMatrix, DraftDiagonalK);
		TempDiagonalK = prod(AuxRotationMatrix,
			Matrix(trans(TempTransformationMatrix)));

		for (int i = 0; i < spring_size; ++i)
		{
			for (int j = 0; j < spring_size; ++j)
			{
				TempLeftHandSide(i + dimension, j + dimension) += TempDiagonalK(i, j);
			}
		}
		rLeftHandSideMatrix += TempLeftHandSide;

		//node 3 -> 4
		TempDiagonalK = ZeroMatrix(spring_size, spring_size);
		TempTransformationMatrix = ZeroMatrix(spring_size, spring_size);
		AuxRotationMatrix = ZeroMatrix(spring_size, spring_size);
		TempLeftHandSide = ZeroMatrix(element_size, element_size);

		this->CreateTransformationMatrix(TempTransformationMatrix, 2, 3);
		AuxRotationMatrix = prod(TempTransformationMatrix, DraftDiagonalK);
		TempDiagonalK = prod(AuxRotationMatrix,
			Matrix(trans(TempTransformationMatrix)));

		for (int i = 0; i < spring_size; ++i)
		{
			for (int j = 0; j < spring_size; ++j)
			{
				TempLeftHandSide(i + (2 * dimension), j + (2 * dimension))
					+= TempDiagonalK(i, j);
			}
		}
		rLeftHandSideMatrix += TempLeftHandSide;

		//node 4 -> 1
		TempDiagonalK = ZeroMatrix(spring_size, spring_size);
		TempTransformationMatrix = ZeroMatrix(spring_size, spring_size);
		AuxRotationMatrix = ZeroMatrix(spring_size, spring_size);
		TempLeftHandSide = ZeroMatrix(element_size, element_size);

		this->CreateTransformationMatrix(TempTransformationMatrix, 3, 0);
		AuxRotationMatrix = prod(TempTransformationMatrix, DraftDiagonalK);
		TempDiagonalK = prod(AuxRotationMatrix,
			Matrix(trans(TempTransformationMatrix)));

		for (int i = 0; i < dimension; ++i)
		{
			for (int j = 0; j < dimension; ++j)
			{
				TempLeftHandSide(i, j) += TempDiagonalK(i, j);
				TempLeftHandSide(i + (3 * dimension), j + (3 * dimension)) +=
					TempDiagonalK(i + dimension, j + dimension);

				TempLeftHandSide(i, j + (3 * dimension)) += TempDiagonalK(i, j + dimension);
				TempLeftHandSide(i + (3 * dimension), j) += TempDiagonalK(i + dimension, j);
			}
		}
		rLeftHandSideMatrix += TempLeftHandSide;
	}
	KRATOS_CATCH("")
}
double RoccoNetElement3D4N::CalculateDiagonalStiffnessKb(const int node1, const int node2)
{
	const int number_of_nodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int element_size = number_of_nodes * dimension;

	double StiffnessKb = this->mKb;
	const double length = this->CalculateActualDiagonalLength(node1, node2);
	const double RefLength = this->mLengthDiagonalReference;

	// use dring - tring as ref length
	if (length <= RefLength) StiffnessKb = 0.00;

	return StiffnessKb;
}





void RoccoNetElement3D4N::CreateTransformationMatrix(Matrix& rRotationMatrix,
	const int node1, const int node2) 
	{
	KRATOS_TRY
	const int number_of_nodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int local_size = 2 * dimension;

	//1st calculate transformation matrix
	Vector DirectionVectorX = ZeroVector(dimension);
	Vector DirectionVectorY = ZeroVector(dimension);
	Vector DirectionVectorZ = ZeroVector(dimension);
	Vector ReferenceCoordinates = ZeroVector(local_size);
	Vector GlobalZ = ZeroVector(dimension);
	GlobalZ[2] = 1.0;

	ReferenceCoordinates[0] = this->GetGeometry()[node1].X();
	ReferenceCoordinates[1] = this->GetGeometry()[node1].Y();
	ReferenceCoordinates[2] = this->GetGeometry()[node1].Z();
	ReferenceCoordinates[3] = this->GetGeometry()[node2].X();
	ReferenceCoordinates[4] = this->GetGeometry()[node2].Y();
	ReferenceCoordinates[5] = this->GetGeometry()[node2].Z();

	for (int i = 0; i < dimension; ++i)
	{
		DirectionVectorX[i] = (ReferenceCoordinates[i + dimension] -
			ReferenceCoordinates[i]);
	}
	// local x-axis (e1_local) is the beam axis  (in GID is e3_local)
	double VectorNorm;
	VectorNorm = MathUtils<double>::Norm(DirectionVectorX);
	if (VectorNorm != 0) DirectionVectorX /= VectorNorm;

	if (DirectionVectorX[2] == 1.00) {
		DirectionVectorY[1] = 1.0;
		DirectionVectorZ[0] = -1.0;
	}

	if (DirectionVectorX[2] == -1.00) {
		DirectionVectorY[1] = 1.0;
		DirectionVectorZ[0] = 1.0;
	}

	if (fabs(DirectionVectorX[2]) != 1.00) {

		DirectionVectorY = MathUtils<double>::CrossProduct(GlobalZ,
			DirectionVectorX);
		VectorNorm = MathUtils<double>::Norm(DirectionVectorY);
		if (VectorNorm != 0) DirectionVectorY /= VectorNorm;

		DirectionVectorZ = MathUtils<double>::CrossProduct(DirectionVectorX,
			DirectionVectorY);
		VectorNorm = MathUtils<double>::Norm(DirectionVectorZ);
		if (VectorNorm != 0) DirectionVectorZ /= VectorNorm;
	}

	//2nd fill big rotation matrix
	MatrixType CurrentCS = ZeroMatrix(dimension, dimension);
	for (int i = 0; i < dimension; ++i) {
		CurrentCS(i, 0) = DirectionVectorX[i];
		CurrentCS(i, 1) = DirectionVectorY[i];
		CurrentCS(i, 2) = DirectionVectorZ[i];
	}

	rRotationMatrix = ZeroMatrix(local_size, local_size);
	if (rRotationMatrix.size1() != local_size) {
		rRotationMatrix.resize(local_size, local_size, false);
	}
	//Building the rotation matrix for the local element matrix
	for (int kk = 0; kk < local_size; kk += dimension)
	{
		for (int i = 0; i<dimension; ++i)
		{
			for (int j = 0; j<dimension; ++j)
			{
				rRotationMatrix(i + kk, j + kk) = CurrentCS(i, j);
			}
		}
	}
	KRATOS_CATCH("")
}


RoccoNetElement3D4N::RoccoNetElement3D4N(IndexType NewId,
	GeometryType::Pointer pGeometry)
	: Element(NewId, pGeometry) {}

RoccoNetElement3D4N::RoccoNetElement3D4N(IndexType NewId,
	GeometryType::Pointer pGeometry,
	PropertiesType::Pointer pProperties)
	: Element(NewId, pGeometry, pProperties) {}

RoccoNetElement3D4N::~RoccoNetElement3D4N() {}

Element::Pointer RoccoNetElement3D4N::Create(IndexType NewId,
	NodesArrayType const& rThisNodes,
	PropertiesType::Pointer pProperties) const
{
	const GeometryType& rGeom = this->GetGeometry();
	return BaseType::Pointer(new RoccoNetElement3D4N(
		NewId, rGeom.Create(rThisNodes), pProperties));
}




int RoccoNetElement3D4N::Check(const ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY

		GeometryType& geom = GetGeometry();

	// verify that the variables are correctly initialized
	if (DISPLACEMENT.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "");
	if (ROTATION.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument, "ROTATION has Key zero! (check if the application is correctly registered", "");
	if (VELOCITY.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "");
	if (ACCELERATION.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument, "ACCELERATION has Key zero! (check if the application is correctly registered", "");
	if (DENSITY.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument, "DENSITY has Key zero! (check if the application is correctly registered", "");
	if (SHELL_CROSS_SECTION.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument, "SHELL_CROSS_SECTION has Key zero! (check if the application is correctly registered", "");
	if (THICKNESS.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument, "THICKNESS has Key zero! (check if the application is correctly registered", "");
	if (CONSTITUTIVE_LAW.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument, "CONSTITUTIVE_LAW has Key zero! (check if the application is correctly registered", "");

	// verify that the dofs exist
	for (unsigned int i = 0; i<geom.size(); i++)
	{
		if (geom[i].SolutionStepsDataHas(DISPLACEMENT) == false)
			KRATOS_THROW_ERROR(std::invalid_argument, "missing variable DISPLACEMENT on node ", geom[i].Id());
		if (geom[i].HasDofFor(DISPLACEMENT_X) == false || geom[i].HasDofFor(DISPLACEMENT_Y) == false || geom[i].HasDofFor(DISPLACEMENT_Z) == false)
			KRATOS_THROW_ERROR(std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id());
		if (geom[i].GetBufferSize() < 2)
			KRATOS_THROW_ERROR(std::logic_error, "This Element needs at least a buffer size = 2", "");
	}

	if (this->GetProperties().Has(KB_ROCCO) == false ||
		this->GetProperties()[KB_ROCCO] == 0)
	{
		KRATOS_ERROR << ("KB_ROCCO not provided for this element", this->Id()) << std::endl;
	}

	if (this->GetProperties().Has(WIRE_WINDINGS_ROCCO) == false ||
		this->GetProperties()[WIRE_WINDINGS_ROCCO] == 0)
	{
		KRATOS_ERROR << ("WIRE_WINDINGS_ROCCO not provided for this element", this->Id()) << std::endl;
	}

	if (this->GetProperties().Has(WIRE_THICKNESS_ROCCO) == false ||
		this->GetProperties()[WIRE_THICKNESS_ROCCO] == 0)
	{
		KRATOS_ERROR << ("WIRE_THICKNESS_ROCCO not provided for this element", this->Id()) << std::endl;
	}

	if (this->GetProperties().Has(DIAMETER_ROCCO) == false ||
		this->GetProperties()[DIAMETER_ROCCO] == 0)
	{
		KRATOS_ERROR << ("DIAMETER_ROCCO not provided for this element", this->Id()) << std::endl;
	}

	if (this->GetProperties().Has(NODAL_MASS_ROCCO) == false ||
		this->GetProperties()[NODAL_MASS_ROCCO] == 0)
	{
		KRATOS_ERROR << ("NODAL_MASS_ROCCO not provided for this element", this->Id()) << std::endl;
	}

	if (this->GetProperties().Has(KT_MULT_ROCCO) == false ||
		this->GetProperties()[KT_MULT_ROCCO] == 0)
	{
		KRATOS_ERROR << ("KT_MULT_ROCCO not provided for this element", this->Id()) << std::endl;
	}

	

	return 0;

	KRATOS_CATCH("")
}

void RoccoNetElement3D4N::EquationIdVector(EquationIdVectorType& rResult,
	ProcessInfo& rCurrentProcessInfo)
{
	const int number_of_nodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int local_size = number_of_nodes * dimension;

	if (rResult.size() != local_size) rResult.resize(local_size, false);

	for (int i = 0; i < number_of_nodes; ++i)
	{
		int index = i * dimension;
		rResult[index] = this->GetGeometry()[i]
			.GetDof(DISPLACEMENT_X).EquationId();
		rResult[index + 1] = this->GetGeometry()[i]
			.GetDof(DISPLACEMENT_Y).EquationId();
		rResult[index + 2] = this->GetGeometry()[i]
			.GetDof(DISPLACEMENT_Z).EquationId();
	}
}

void RoccoNetElement3D4N::GetDofList(DofsVectorType& rElementalDofList,
	ProcessInfo& CurrentProcessInfo)
{
	const int number_of_nodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int local_size = number_of_nodes * dimension;

	if (rElementalDofList.size() != local_size) {
		rElementalDofList.resize(local_size);
	}

	for (int i = 0; i < number_of_nodes; ++i)
	{
		int index = i * dimension;
		rElementalDofList[index] = this->GetGeometry()[i]
			.pGetDof(DISPLACEMENT_X);
		rElementalDofList[index + 1] = this->GetGeometry()[i]
			.pGetDof(DISPLACEMENT_Y);
		rElementalDofList[index + 2] = this->GetGeometry()[i]
			.pGetDof(DISPLACEMENT_Z);
	}
}


void RoccoNetElement3D4N::GetValuesVector(Vector& rValues, int Step)
{
	KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int element_size = number_of_nodes * dimension;

	if (rValues.size() != element_size) rValues.resize(element_size, false);

	for (int i = 0; i < number_of_nodes; ++i)
	{
		int index = i * dimension;
		rValues[index] = this->GetGeometry()[i]
			.FastGetSolutionStepValue(DISPLACEMENT_X, Step);
		rValues[index + 1] = this->GetGeometry()[i]
			.FastGetSolutionStepValue(DISPLACEMENT_Y, Step);
		rValues[index + 2] = this->GetGeometry()[i]
			.FastGetSolutionStepValue(DISPLACEMENT_Z, Step);
	}
	KRATOS_CATCH("")
}

void RoccoNetElement3D4N::GetFirstDerivativesVector(Vector& rValues, int Step)
{
	KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int element_size = number_of_nodes * dimension;

	if (rValues.size() != element_size) rValues.resize(element_size, false);

	for (int i = 0; i < number_of_nodes; ++i)
	{
		int index = i * dimension;
		rValues[index] = this->GetGeometry()[i]
			.FastGetSolutionStepValue(VELOCITY_X, Step);
		rValues[index + 1] = this->GetGeometry()[i]
			.FastGetSolutionStepValue(VELOCITY_Y, Step);
		rValues[index + 2] = this->GetGeometry()[i]
			.FastGetSolutionStepValue(VELOCITY_Z, Step);
	}
	KRATOS_CATCH("")
}

void RoccoNetElement3D4N::GetSecondDerivativesVector(Vector& rValues, int Step)
{
	KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int element_size = number_of_nodes * dimension;

	if (rValues.size() != element_size) rValues.resize(element_size, false);

	for (int i = 0; i < number_of_nodes; ++i)
	{
		int index = i * dimension;
		rValues[index] = this->GetGeometry()[i]
			.FastGetSolutionStepValue(ACCELERATION_X, Step);
		rValues[index + 1] = this->GetGeometry()[i]
			.FastGetSolutionStepValue(ACCELERATION_Y, Step);
		rValues[index + 2] = this->GetGeometry()[i]
			.FastGetSolutionStepValue(ACCELERATION_Z, Step);
	}
	KRATOS_CATCH("")
}

void RoccoNetElement3D4N::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int element_size = number_of_nodes * dimension;

	if (rMassMatrix.size1() != element_size)
	{
		rMassMatrix.resize(element_size, element_size, false);
	}

	rMassMatrix = ZeroMatrix(element_size, element_size);

	const double NodalMassRocco = this->mNodalMass_custom;

	for (int i = 0; i < element_size; ++i)
	{
		rMassMatrix(i, i) = NodalMassRocco;
	}
	KRATOS_CATCH("")
}

void RoccoNetElement3D4N::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int element_size = number_of_nodes * dimension;

	if (rDampingMatrix.size1() != element_size)
	{
		rDampingMatrix.resize(element_size, element_size, false);
	}
	///TODO!
	rDampingMatrix = ZeroMatrix(element_size, element_size);
	KRATOS_CATCH("")
}


void RoccoNetElement3D4N::AddExplicitContribution(const VectorType& rRHSVector,
	const Variable<VectorType>& rRHSVariable,
	Variable<array_1d<double, 3> >& rDestinationVariable,
	const ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY;
	const int number_of_nodes = this->GetGeometry().PointsNumber();
	const int dimension = this->GetGeometry().WorkingSpaceDimension();
	const int element_size = number_of_nodes * dimension;

	if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL)
	{

		for (int i = 0; i< number_of_nodes; ++i)
		{
			int index = dimension * i;

			GetGeometry()[i].SetLock();

			array_1d<double, 3 > &ForceResidual =
				GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

			for (int j = 0; j<dimension; ++j)
			{
				ForceResidual[j] += rRHSVector[index + j];
			}

			GetGeometry()[i].UnSetLock();
		}
	}
	KRATOS_CATCH("")
}




} // namespace Kratos.
