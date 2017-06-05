// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//                   
//                   
//
#include "custom_elements/cr_beam_element_3D2N.hpp"
#include "structural_mechanics_application_variables.h"
#include "includes/define.h"



namespace Kratos
{

	CrBeamElement3D2N::CrBeamElement3D2N(IndexType NewId,
		GeometryType::Pointer pGeometry, bool rLinear)
		: Element(NewId, pGeometry)
	{
		this->mIsLinearElement = rLinear;
	}

	CrBeamElement3D2N::CrBeamElement3D2N(IndexType NewId,
		GeometryType::Pointer pGeometry,
		PropertiesType::Pointer pProperties, bool rLinear)
		: Element(NewId, pGeometry, pProperties)
	{
		this->mIsLinearElement = rLinear;
	}

	Element::Pointer CrBeamElement3D2N::Create(IndexType NewId,
		NodesArrayType const& rThisNodes,
		PropertiesType::Pointer pProperties) const
	{
		const GeometryType& rGeom = this->GetGeometry();
		return BaseType::Pointer(new CrBeamElement3D2N(
			NewId, rGeom.Create(rThisNodes), pProperties, this->mIsLinearElement));
	}

	CrBeamElement3D2N::~CrBeamElement3D2N() {}

	void CrBeamElement3D2N::EquationIdVector(EquationIdVectorType& rResult,
		ProcessInfo& rCurrentProcessInfo) {

		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int local_size = number_of_nodes * dimension * 2;

		if (rResult.size() != local_size) rResult.resize(local_size);

		for (int i = 0; i < number_of_nodes; ++i)
		{
			int index = i * number_of_nodes * dimension;
			rResult[index] = this->GetGeometry()[i].GetDof(DISPLACEMENT_X)
				.EquationId();
			rResult[index + 1] = this->GetGeometry()[i].GetDof(DISPLACEMENT_Y)
				.EquationId();
			rResult[index + 2] = this->GetGeometry()[i].GetDof(DISPLACEMENT_Z)
				.EquationId();

			rResult[index + 3] = this->GetGeometry()[i].GetDof(ROTATION_X)
				.EquationId();
			rResult[index + 4] = this->GetGeometry()[i].GetDof(ROTATION_Y)
				.EquationId();
			rResult[index + 5] = this->GetGeometry()[i].GetDof(ROTATION_Z)
				.EquationId();
		}

	}

	void CrBeamElement3D2N::GetDofList(DofsVectorType& rElementalDofList,
		ProcessInfo& rCurrentProcessInfo) {

		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int local_size = number_of_nodes * dimension * 2;

		if (rElementalDofList.size() != local_size) {
			rElementalDofList.resize(local_size);
		}

		for (int i = 0; i < number_of_nodes; ++i)
		{
			int index = i * number_of_nodes * dimension;
			rElementalDofList[index] = this->GetGeometry()[i]
				.pGetDof(DISPLACEMENT_X);
			rElementalDofList[index + 1] = this->GetGeometry()[i]
				.pGetDof(DISPLACEMENT_Y);
			rElementalDofList[index + 2] = this->GetGeometry()[i]
				.pGetDof(DISPLACEMENT_Z);

			rElementalDofList[index + 3] = this->GetGeometry()[i]
				.pGetDof(ROTATION_X);
			rElementalDofList[index + 4] = this->GetGeometry()[i]
				.pGetDof(ROTATION_Y);
			rElementalDofList[index + 5] = this->GetGeometry()[i]
				.pGetDof(ROTATION_Z);
		}
	}

	void CrBeamElement3D2N::Initialize() {

		KRATOS_TRY;
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int local_size = number_of_nodes * dimension;


		this->mPoisson = this->GetProperties()[POISSON_RATIO];
		this->mArea = this->GetProperties()[CROSS_AREA];
		this->mYoungsModulus = this->GetProperties()[YOUNG_MODULUS];
		this->mShearModulus = this->mYoungsModulus / (2.0 * (1.0 + this->mPoisson));
		this->mLength = this->CalculateReferenceLength();
		this->mCurrentLength = this->CalculateCurrentLength();
		this->mDensity = this->GetProperties()[DENSITY];
		this->mInertiaX = this->GetProperties()[IT];
		this->mInertiaY = this->GetProperties()[IY];
		this->mInertiaZ = this->GetProperties()[IZ];
		//effecive shear Area
		if (this->GetProperties().Has(AREA_EFFECTIVE_Y) == true) {
			this->mEffAreaY = GetProperties()[AREA_EFFECTIVE_Y];
		}
		else this->mEffAreaY = 0.00;
		if (this->GetProperties().Has(AREA_EFFECTIVE_Z) == true) {
			this->mEffAreaZ = GetProperties()[AREA_EFFECTIVE_Z];
		}
		else this->mEffAreaZ = 0.00;

		//rotational inertia
		if (this->GetProperties().Has(INERTIA_ROT_Y) == true) {
			this->mRotInertiaY = GetProperties()[INERTIA_ROT_Y];
		}
		else this->mRotInertiaY = this->mInertiaY;
		if (this->GetProperties().Has(INERTIA_ROT_Z) == true) {
			this->mRotInertiaZ = GetProperties()[INERTIA_ROT_Z];
		}
		else this->mRotInertiaZ = this->mInertiaZ;

		//caluclate Psi,y,z
		this->mPsiY = this->CalculatePsi(this->mInertiaY, this->mEffAreaZ);
		this->mPsiZ = this->CalculatePsi(this->mInertiaZ, this->mEffAreaY);
		//manual beam rotation
		if (this->GetProperties().Has(ANG_ROT) == true) {
			this->mtheta = this->GetProperties()[ANG_ROT];
		}
		else this->mtheta = 0.00;
		if (this->mLength == 0.00) {
			KRATOS_ERROR << ("Zero length found in element #", this->Id()) <<
				std::endl;
		}

		if (this->mIterationCount == 0)
		{
			this->mNodalForces = ZeroVector(local_size * 2);
		}
		KRATOS_CATCH("")
	}

	Matrix CrBeamElement3D2N::CreateElementStiffnessMatrix_Material() {

		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int local_size = number_of_nodes * dimension * 2;

		const double E = this->mYoungsModulus;
		const double G = this->mShearModulus;
		const double A = this->mArea;
		const double L = this->mLength;

		const double J = this->mInertiaX;
		const double Iy = this->mInertiaY;
		const double Iz = this->mInertiaZ;

		const double Psi_y = this->mPsiY;
		const double Psi_z = this->mPsiZ;

		Matrix LocalStiffnessMatrix = ZeroMatrix(local_size, local_size);
		const double L3 = L*L*L;
		const double L2 = L*L;


		LocalStiffnessMatrix(0, 0) = E*A / L;
		LocalStiffnessMatrix(6, 0) = -1.0 * LocalStiffnessMatrix(0, 0);
		LocalStiffnessMatrix(0, 6) = LocalStiffnessMatrix(6, 0);
		LocalStiffnessMatrix(6, 6) = LocalStiffnessMatrix(0, 0);

		LocalStiffnessMatrix(1, 1) = 12.0 * E * Iz * Psi_z / L3;
		LocalStiffnessMatrix(1, 7) = -1.0 * LocalStiffnessMatrix(1, 1);
		LocalStiffnessMatrix(1, 5) = 6.0 * E * Iz * Psi_z / L2;
		LocalStiffnessMatrix(1, 11) = LocalStiffnessMatrix(1, 5);

		LocalStiffnessMatrix(2, 2) = 12.0 * E *Iy * Psi_y / L3;
		LocalStiffnessMatrix(2, 8) = -1.0 * LocalStiffnessMatrix(2, 2);
		LocalStiffnessMatrix(2, 4) = -6.0 * E *Iy *Psi_y / L2;
		LocalStiffnessMatrix(2, 10) = LocalStiffnessMatrix(2, 4);

		LocalStiffnessMatrix(4, 2) = LocalStiffnessMatrix(2, 4);
		LocalStiffnessMatrix(5, 1) = LocalStiffnessMatrix(1, 5);
		LocalStiffnessMatrix(3, 3) = G*J / L;
		LocalStiffnessMatrix(4, 4) = E*Iy*(3.0 * Psi_y + 1.0) / L;
		LocalStiffnessMatrix(5, 5) = E*Iz*(3.0 * Psi_z + 1.0) / L;
		LocalStiffnessMatrix(4, 8) = -1.0 * LocalStiffnessMatrix(4, 2);
		LocalStiffnessMatrix(5, 7) = -1.0 * LocalStiffnessMatrix(5, 1);
		LocalStiffnessMatrix(3, 9) = -1.0 * LocalStiffnessMatrix(3, 3);
		LocalStiffnessMatrix(4, 10) = E*Iy*(3.0 * Psi_y - 1) / L;
		LocalStiffnessMatrix(5, 11) = E*Iz*(3.0 * Psi_z - 1) / L;

		LocalStiffnessMatrix(7, 1) = LocalStiffnessMatrix(1, 7);
		LocalStiffnessMatrix(7, 5) = LocalStiffnessMatrix(5, 7);
		LocalStiffnessMatrix(7, 7) = LocalStiffnessMatrix(1, 1);
		LocalStiffnessMatrix(7, 11) = LocalStiffnessMatrix(7, 5);

		LocalStiffnessMatrix(8, 2) = LocalStiffnessMatrix(2, 8);
		LocalStiffnessMatrix(8, 4) = LocalStiffnessMatrix(4, 8);
		LocalStiffnessMatrix(8, 8) = LocalStiffnessMatrix(2, 2);
		LocalStiffnessMatrix(8, 10) = LocalStiffnessMatrix(8, 4);

		LocalStiffnessMatrix(9, 3) = LocalStiffnessMatrix(3, 9);
		LocalStiffnessMatrix(9, 9) = LocalStiffnessMatrix(3, 3);

		LocalStiffnessMatrix(10, 2) = LocalStiffnessMatrix(2, 10);
		LocalStiffnessMatrix(10, 4) = LocalStiffnessMatrix(4, 10);
		LocalStiffnessMatrix(10, 8) = LocalStiffnessMatrix(8, 10);
		LocalStiffnessMatrix(10, 10) = LocalStiffnessMatrix(4, 4);

		LocalStiffnessMatrix(11, 1) = LocalStiffnessMatrix(1, 11);
		LocalStiffnessMatrix(11, 5) = LocalStiffnessMatrix(5, 11);
		LocalStiffnessMatrix(11, 7) = LocalStiffnessMatrix(7, 11);
		LocalStiffnessMatrix(11, 11) = LocalStiffnessMatrix(5, 5);

		return LocalStiffnessMatrix;
		KRATOS_CATCH("")
	}

	Matrix CrBeamElement3D2N::CreateElementStiffnessMatrix_Geometry(
		const Vector qe) {

		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int local_size = number_of_nodes * dimension * 2;

		const double N = qe[6];
		const double Mt = qe[9];
		const double my_A = qe[4];
		const double mz_A = qe[5];
		const double my_B = qe[10];
		const double mz_B = qe[11];

		const double L = this->mCurrentLength;
		const double Qy = -1.00 * (mz_A + mz_B) / L;
		const double Qz = (my_A + my_B) / L;

		Matrix LocalStiffnessMatrix = ZeroMatrix(local_size, local_size);

		LocalStiffnessMatrix(0, 1) = -Qy / L;
		LocalStiffnessMatrix(0, 2) = -Qz / L;
		LocalStiffnessMatrix(0, 7) = -1.0 * LocalStiffnessMatrix(0, 1);
		LocalStiffnessMatrix(0, 8) = -1.0 * LocalStiffnessMatrix(0, 2);

		LocalStiffnessMatrix(1, 0) = LocalStiffnessMatrix(0, 1);

		LocalStiffnessMatrix(1, 1) = 1.2 * N / L;

		LocalStiffnessMatrix(1, 3) = my_A / L;
		LocalStiffnessMatrix(1, 4) = Mt / L;

		LocalStiffnessMatrix(1, 5) = N / 10.0;

		LocalStiffnessMatrix(1, 6) = LocalStiffnessMatrix(0, 7);
		LocalStiffnessMatrix(1, 7) = -1.00 * LocalStiffnessMatrix(1, 1);
		LocalStiffnessMatrix(1, 9) = my_B / L;
		LocalStiffnessMatrix(1, 10) = -1.00 * LocalStiffnessMatrix(1, 4);
		LocalStiffnessMatrix(1, 11) = LocalStiffnessMatrix(1, 5);

		LocalStiffnessMatrix(2, 0) = LocalStiffnessMatrix(0, 2);
		LocalStiffnessMatrix(2, 2) = LocalStiffnessMatrix(1, 1);
		LocalStiffnessMatrix(2, 3) = mz_A / L;
		LocalStiffnessMatrix(2, 4) = -1.00 * LocalStiffnessMatrix(1, 5);
		LocalStiffnessMatrix(2, 5) = LocalStiffnessMatrix(1, 4);
		LocalStiffnessMatrix(2, 6) = LocalStiffnessMatrix(0, 8);
		LocalStiffnessMatrix(2, 8) = LocalStiffnessMatrix(1, 7);
		LocalStiffnessMatrix(2, 9) = mz_B / L;
		LocalStiffnessMatrix(2, 10) = LocalStiffnessMatrix(2, 4);
		LocalStiffnessMatrix(2, 11) = LocalStiffnessMatrix(1, 10);

		for (int i = 0; i < 3; ++i) {
			LocalStiffnessMatrix(3, i) = LocalStiffnessMatrix(i, 3);
		}
		LocalStiffnessMatrix(3, 4) = (-mz_A / 3.00) + (mz_B / 6.00);
		LocalStiffnessMatrix(3, 5) = (my_A / 3.00) - (my_B / 6.00);
		LocalStiffnessMatrix(3, 7) = -my_A / L;
		LocalStiffnessMatrix(3, 8) = -mz_A / L;
		LocalStiffnessMatrix(3, 10) = L*Qy / 6.00;
		LocalStiffnessMatrix(3, 11) = L*Qz / 6.00;

		for (int i = 0; i < 4; ++i) {
			LocalStiffnessMatrix(4, i) = LocalStiffnessMatrix(i, 4);
		}
		LocalStiffnessMatrix(4, 4) = 2.00 * L*N / 15.00;
		LocalStiffnessMatrix(4, 7) = -Mt / L;
		LocalStiffnessMatrix(4, 8) = N / 10.00;
		LocalStiffnessMatrix(4, 9) = LocalStiffnessMatrix(3, 10);
		LocalStiffnessMatrix(4, 10) = -L*N / 30.00;
		LocalStiffnessMatrix(4, 11) = Mt / 2.00;


		for (int i = 0; i < 5; ++i) {
			LocalStiffnessMatrix(5, i) = LocalStiffnessMatrix(i, 5);
		}
		LocalStiffnessMatrix(5, 5) = LocalStiffnessMatrix(4, 4);
		LocalStiffnessMatrix(5, 7) = -N / 10.0;
		LocalStiffnessMatrix(5, 8) = -Mt / L;
		LocalStiffnessMatrix(5, 9) = LocalStiffnessMatrix(3, 11);
		LocalStiffnessMatrix(5, 10) = -1.00 * LocalStiffnessMatrix(4, 11);
		LocalStiffnessMatrix(5, 11) = LocalStiffnessMatrix(4, 10);

		for (int i = 0; i < 6; ++i) {
			LocalStiffnessMatrix(6, i) = LocalStiffnessMatrix(i, 6);
		}
		LocalStiffnessMatrix(6, 7) = LocalStiffnessMatrix(0, 1);
		LocalStiffnessMatrix(6, 8) = LocalStiffnessMatrix(0, 2);

		for (int i = 0; i < 7; ++i) {
			LocalStiffnessMatrix(7, i) = LocalStiffnessMatrix(i, 7);
		}
		LocalStiffnessMatrix(7, 7) = LocalStiffnessMatrix(1, 1);
		LocalStiffnessMatrix(7, 9) = -1.00 * LocalStiffnessMatrix(1, 9);
		LocalStiffnessMatrix(7, 10) = LocalStiffnessMatrix(4, 1);
		LocalStiffnessMatrix(7, 11) = LocalStiffnessMatrix(2, 4);

		for (int i = 0; i < 8; ++i) {
			LocalStiffnessMatrix(8, i) = LocalStiffnessMatrix(i, 8);
		}
		LocalStiffnessMatrix(8, 8) = LocalStiffnessMatrix(1, 1);
		LocalStiffnessMatrix(8, 9) = -1.00 * LocalStiffnessMatrix(2, 9);
		LocalStiffnessMatrix(8, 10) = LocalStiffnessMatrix(1, 5);
		LocalStiffnessMatrix(8, 11) = LocalStiffnessMatrix(1, 4);

		for (int i = 0; i < 9; ++i) {
			LocalStiffnessMatrix(9, i) = LocalStiffnessMatrix(i, 9);
		}
		LocalStiffnessMatrix(9, 10) = (mz_A / 6.00) - (mz_B / 3.00);
		LocalStiffnessMatrix(9, 11) = (-my_A / 6.00) + (my_B / 3.00);

		for (int i = 0; i < 10; ++i) {
			LocalStiffnessMatrix(10, i) = LocalStiffnessMatrix(i, 10);
		}
		LocalStiffnessMatrix(10, 10) = LocalStiffnessMatrix(4, 4);

		for (int i = 0; i < 11; ++i) {
			LocalStiffnessMatrix(11, i) = LocalStiffnessMatrix(i, 11);
		}
		LocalStiffnessMatrix(11, 11) = LocalStiffnessMatrix(4, 4);

		return LocalStiffnessMatrix;
		KRATOS_CATCH("")
	}

	Matrix CrBeamElement3D2N::CalculateDeformationStiffness() {

		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int local_size = number_of_nodes * dimension;

		Matrix Kd = ZeroMatrix(local_size, local_size);
		const double E = this->mYoungsModulus;
		const double G = this->mShearModulus;
		const double A = this->mArea;
		const double L = this->mLength;
		const double J = this->mInertiaX;
		const double Iy = this->mInertiaY;
		const double Iz = this->mInertiaZ;

		const double Psi_y = this->mPsiY;
		const double Psi_z = this->mPsiZ;

		Kd(0, 0) = G * J / L;
		Kd(1, 1) = E * Iy / L;
		Kd(2, 2) = E * Iz / L;
		Kd(3, 3) = E * A / L;
		Kd(4, 4) = 3.0 * E * Iy * Psi_y / L;
		Kd(5, 5) = 3.0 * E * Iz * Psi_z / L;


		//add geometric stiffness part
		if (this->mIsLinearElement == false)
		{
			const double l = this->mCurrentLength;
			const double N = this->mNodalForces[6];

			const double Qy = -1.00 * (this->mNodalForces[5] +
				this->mNodalForces[11]) / l;

			const double Qz = 1.00 * (this->mNodalForces[4] +
				this->mNodalForces[10]) / l;

			const double N1 = l*N / 12.00;
			const double N2 = l*N / 20.00;
			const double Qy1 = -l*Qy / 6.00;
			const double Qz1 = -l*Qz / 6.00;

			Kd(1, 1) += N1;
			Kd(2, 2) += N1;
			Kd(4, 4) += N2;
			Kd(5, 5) += N2;

			Kd(0, 1) += Qy1;
			Kd(0, 2) += Qz1;
			Kd(1, 0) += Qy1;
			Kd(2, 0) += Qz1;

		}
		return Kd;
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateInitialLocalCS() {

		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = number_of_nodes * dimension;
		const int local_size = size * 2;

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
		Orientation element_axis(DirectionVectorX, this->mtheta);
		element_axis.CalculateBasisVectors(DirectionVectorX, DirectionVectorY,
			DirectionVectorZ);
		//save them to update the local axis in every following iter. step
		this->mNX0 = DirectionVectorX;
		this->mNY0 = DirectionVectorY;
		this->mNZ0 = DirectionVectorZ;


		Matrix Temp = ZeroMatrix(dimension);
		this->mRotationMatrix0 = ZeroMatrix(local_size);
		element_axis.CalculateRotationMatrix(Temp);
		this->AssembleSmallInBigMatrix(Temp, this->mRotationMatrix0);
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateTransformationMatrix(Matrix& rRotationMatrix) {

		KRATOS_TRY
			//12x12
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = number_of_nodes * dimension;
		const int MatSize = 2 * size;

		//initialize local CS
		if (this->mIterationCount == 0) this->CalculateInitialLocalCS();

		//update local CS
		Matrix AuxRotationMatrix = ZeroMatrix(dimension);
		AuxRotationMatrix = this->UpdateRotationMatrixLocal();

		if (rRotationMatrix.size1() != MatSize) {
			rRotationMatrix.resize(MatSize, MatSize, false);
		}

		rRotationMatrix = ZeroMatrix(MatSize);
		//Building the rotation matrix for the local element matrix
		this->AssembleSmallInBigMatrix(AuxRotationMatrix, rRotationMatrix);
		KRATOS_CATCH("")
	}

	Matrix CrBeamElement3D2N::CalculateTransformationS() {

		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = number_of_nodes * dimension;
		const int MatSize = 2 * size;

		const double L = this->mCurrentLength;
		Matrix S = ZeroMatrix(MatSize, size);
		S(0, 3) = -1.00;
		S(1, 5) = 2.00 / L;
		S(2, 4) = -2.00 / L;
		S(3, 0) = -1.00;
		S(4, 1) = -1.00;
		S(4, 4) = 1.00;
		S(5, 2) = -1.00;
		S(5, 5) = 1.00;
		S(6, 3) = 1.00;
		S(7, 5) = -2.00 / L;
		S(8, 4) = 2.00 / L;
		S(9, 0) = 1.00;
		S(10, 1) = 1.00;
		S(10, 4) = 1.00;
		S(11, 2) = 1.00;
		S(11, 5) = 1.00;

		return S;
		KRATOS_CATCH("")
	}

	Matrix CrBeamElement3D2N::UpdateRotationMatrixLocal() {

		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = number_of_nodes * dimension;
		const int MatSize = 2 * size;

		Vector dPhiA = ZeroVector(dimension);
		Vector dPhiB = ZeroVector(dimension);
		Vector IncrementDeformation = ZeroVector(MatSize);
		IncrementDeformation = this->mIncrementDeformation;

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


		//1st solution step
		if (mIterationCount == 0) {
			this->mQuaternionVEC_A = ZeroVector(dimension);
			this->mQuaternionVEC_B = ZeroVector(dimension);
			this->mQuaternionSCA_A = 1.00;
			this->mQuaternionSCA_B = 1.00;
		}

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
		Vector rotatedNX0 = this->mNX0;
		Vector rotatedNY0 = this->mNY0;
		Vector rotatedNZ0 = this->mNZ0;
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

		deltaX[0] = this->mTotalNodalPosistion[3] - this->mTotalNodalPosistion[0];
		deltaX[1] = this->mTotalNodalPosistion[4] - this->mTotalNodalPosistion[1];
		deltaX[2] = this->mTotalNodalPosistion[5] - this->mTotalNodalPosistion[2];


		VectorNorm = MathUtils<double>::Norm(deltaX);
		if (VectorNorm != 0.00) deltaX /= VectorNorm;


		n_bisectrix = rotatedNX0 + deltaX;
		VectorNorm = MathUtils<double>::Norm(n_bisectrix);
		if (VectorNorm != 0.00) n_bisectrix /= VectorNorm;

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


		//save current CS for GID OUTPUT
		this->mNX = ZeroVector(dimension);
		this->mNY = ZeroVector(dimension);
		this->mNZ = ZeroVector(dimension);
		for (int i = 0; i < dimension; ++i)
		{
			this->mNX[i] = n_xyz(i, 0);
			this->mNY[i] = n_xyz(i, 1);
			this->mNZ[i] = n_xyz(i, 2);
		}

		//calculating deformation modes
		this->mPhiS = ZeroVector(dimension);
		this->mPhiA = ZeroVector(dimension);
		this->mPhiS = prod(Matrix(trans(n_xyz)), vector_diff);
		this->mPhiS *= 4.00;

		rotatedNX0 = ZeroVector(dimension);
		tempVec = ZeroVector(dimension);
		for (int i = 0; i < dimension; ++i) rotatedNX0[i] = n_xyz(i, 0);
		tempVec = MathUtils<double>::CrossProduct(rotatedNX0, n_bisectrix);
		this->mPhiA = prod(Matrix(trans(n_xyz)), tempVec);
		this->mPhiA *= 4.00;

		if (this->mIterationCount == 0)
		{
			this->mPhiS = ZeroVector(dimension);
			this->mPhiA = ZeroVector(dimension);
		}
		return n_xyz;
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::GetValuesVector(Vector& rValues, int Step) {

		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int element_size = number_of_nodes * dimension * 2;

		if (rValues.size() != element_size) rValues.resize(element_size, false);

		for (int i = 0; i < number_of_nodes; ++i)
		{
			int index = i * dimension * 2;
			rValues[index] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(DISPLACEMENT_X, Step);
			rValues[index + 1] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(DISPLACEMENT_Y, Step);
			rValues[index + 2] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(DISPLACEMENT_Z, Step);

			rValues[index + 3] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ROTATION_X, Step);
			rValues[index + 4] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ROTATION_Y, Step);
			rValues[index + 5] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ROTATION_Z, Step);
		}
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::GetFirstDerivativesVector(Vector& rValues, int Step)
	{

		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int element_size = number_of_nodes * dimension * 2;

		if (rValues.size() != element_size) rValues.resize(element_size, false);

		for (int i = 0; i < number_of_nodes; ++i)
		{
			int index = i * dimension * 2;
			rValues[index] = this->GetGeometry()[i].
				FastGetSolutionStepValue(VELOCITY_X, Step);
			rValues[index + 1] = this->GetGeometry()[i].
				FastGetSolutionStepValue(VELOCITY_Y, Step);
			rValues[index + 2] = this->GetGeometry()[i].
				FastGetSolutionStepValue(VELOCITY_Z, Step);

			rValues[index + 3] = this->GetGeometry()[i].
				FastGetSolutionStepValue(ANGULAR_VELOCITY_X, Step);
			rValues[index + 4] = this->GetGeometry()[i].
				FastGetSolutionStepValue(ANGULAR_VELOCITY_Y, Step);
			rValues[index + 5] = this->GetGeometry()[i].
				FastGetSolutionStepValue(ANGULAR_VELOCITY_Z, Step);
		}

		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::GetSecondDerivativesVector(Vector& rValues, int Step)
	{

		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int element_size = number_of_nodes * dimension * 2;

		if (rValues.size() != element_size) rValues.resize(element_size, false);

		for (int i = 0; i < number_of_nodes; ++i)
		{
			int index = i * dimension * 2;

			rValues[index] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ACCELERATION_X, Step);
			rValues[index + 1] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ACCELERATION_Y, Step);
			rValues[index + 2] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ACCELERATION_Z, Step);

			rValues[index + 3] = this->GetGeometry()[i].
				FastGetSolutionStepValue(ANGULAR_ACCELERATION_X, Step);
			rValues[index + 4] = this->GetGeometry()[i].
				FastGetSolutionStepValue(ANGULAR_ACCELERATION_Y, Step);
			rValues[index + 5] = this->GetGeometry()[i].
				FastGetSolutionStepValue(ANGULAR_ACCELERATION_Z, Step);
		}
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateMassMatrix(MatrixType& rMassMatrix,
		ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		const int number_of_nodes = GetGeometry().PointsNumber();
		const int dimension = GetGeometry().WorkingSpaceDimension();
		const int smallMatSize = number_of_nodes * 2;
		const int MatSize = number_of_nodes * dimension * 2;

		if (rMassMatrix.size1() != MatSize) {
			rMassMatrix.resize(MatSize, MatSize, false);
		}
		rMassMatrix = ZeroMatrix(MatSize, MatSize);

		////////////////////////////
		//TESTING
		this->mIsLumpedMassMatrix = false;
		////////////////////////////


		if (this->mIsLumpedMassMatrix == true)
		{
			this->CalculateLumpedMassMatrix(rMassMatrix, rCurrentProcessInfo);
		}
		else
		{
			this->CalculateConsistentMassMatrix(rMassMatrix, rCurrentProcessInfo);

			Matrix RotationMatrix = ZeroMatrix(MatSize);
			Matrix aux_matrix = ZeroMatrix(MatSize);

			RotationMatrix = this->mRotationMatrix;
			aux_matrix = prod(RotationMatrix, rMassMatrix);
			rMassMatrix = prod(aux_matrix,
				Matrix(trans(RotationMatrix)));
		}
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateDampingMatrix(MatrixType& rDampingMatrix,
		ProcessInfo& rCurrentProcessInfo) {

		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int MatSize = number_of_nodes * dimension * 2;

		if (rDampingMatrix.size1() != MatSize)
		{
			rDampingMatrix.resize(MatSize, MatSize, false);
		}

		rDampingMatrix = ZeroMatrix(MatSize, MatSize);

		Matrix StiffnessMatrix = ZeroMatrix(MatSize, MatSize);

		this->CalculateLeftHandSide(StiffnessMatrix, rCurrentProcessInfo);

		Matrix MassMatrix = ZeroMatrix(MatSize, MatSize);

		this->CalculateMassMatrix(MassMatrix, rCurrentProcessInfo);

		double alpha = 0.0;
		if (this->GetProperties().Has(RAYLEIGH_ALPHA))
		{
			alpha = this->GetProperties()[RAYLEIGH_ALPHA];
		}
		else if (rCurrentProcessInfo.Has(RAYLEIGH_ALPHA))
		{
			alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
		}

		double beta = 0.0;
		if (this->GetProperties().Has(RAYLEIGH_BETA))
		{
			beta = this->GetProperties()[RAYLEIGH_BETA];
		}
		else if (rCurrentProcessInfo.Has(RAYLEIGH_BETA))
		{
			beta = rCurrentProcessInfo[RAYLEIGH_BETA];
		}

		rDampingMatrix += alpha * MassMatrix;
		rDampingMatrix += beta  * StiffnessMatrix;

		KRATOS_CATCH("")
	}


	Vector CrBeamElement3D2N::CalculateBodyForces()
	{
		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int localSize = number_of_nodes * dimension;
		const int MatSize = number_of_nodes * dimension * 2;

		//getting shapefunctionvalues for linear SF
		const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues(
			GeometryData::GI_GAUSS_1);

		Vector EquivalentLineLoad = ZeroVector(dimension);
		Vector BodyForcesGlobal = ZeroVector(MatSize);

		//calculating equivalent line load
		for (int i = 0; i < number_of_nodes; ++i)
		{
			EquivalentLineLoad += this->mArea*this->mDensity*
				this->GetGeometry()[i].
				FastGetSolutionStepValue(VOLUME_ACCELERATION)*Ncontainer(0, i);
		}


		// adding the nodal forces
		for (int i = 0; i < number_of_nodes; ++i)
		{
			int index = i*localSize;
			for (int j = 0; j < dimension; ++j)
			{
				BodyForcesGlobal[j + index] =
					EquivalentLineLoad[j] * Ncontainer(0, i) * this->mCurrentLength;
			}
		}

		// adding the nodal moments
		this->CalculateAndAddWorkEquivalentNodalForcesLineLoad
			(EquivalentLineLoad, BodyForcesGlobal, this->mCurrentLength);


		// return the total ForceVector
		return BodyForcesGlobal;
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateAndAddWorkEquivalentNodalForcesLineLoad(
		const Vector ForceInput, VectorType& rRightHandSideVector,
		const double GeometryLength)
	{
		KRATOS_TRY;
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		//calculate orthogonal load vector
		Vector GeometricOrientation = ZeroVector(dimension);
		GeometricOrientation[0] = this->GetGeometry()[1].X()
			- this->GetGeometry()[0].X();
		GeometricOrientation[1] = this->GetGeometry()[1].Y()
			- this->GetGeometry()[0].Y();
		if (dimension == 3)
		{
			GeometricOrientation[2] = this->GetGeometry()[1].Z()
				- this->GetGeometry()[0].Z();
		}

		const double VectorNormA = MathUtils<double>::Norm(GeometricOrientation);
		if (VectorNormA != 0.00) GeometricOrientation /= VectorNormA;

		Vector LineLoadDir = ZeroVector(dimension);
		for (int i = 0; i < dimension; ++i)
		{
			LineLoadDir[i] = ForceInput[i];
		}

		const double VectorNormB = MathUtils<double>::Norm(LineLoadDir);
		if (VectorNormB != 0.00) LineLoadDir /= VectorNormB;

		double cosAngle = 0.00;
		for (int i = 0; i < dimension; ++i)
		{
			cosAngle += LineLoadDir[i] * GeometricOrientation[i];
		}

		const double sinAngle = sqrt(1.00 - (cosAngle*cosAngle));
		const double NormForceVectorOrth = sinAngle * VectorNormB;


		Vector NodeA = ZeroVector(dimension);
		NodeA[0] = this->GetGeometry()[0].X();
		NodeA[1] = this->GetGeometry()[0].Y();
		if (dimension == 3)	NodeA[2] = this->GetGeometry()[0].Z();

		Vector NodeB = ZeroVector(dimension);
		NodeB = NodeA + LineLoadDir;

		Vector NodeC = ZeroVector(dimension);
		NodeC = NodeA + (GeometricOrientation*cosAngle);

		Vector LoadOrthogonalDir = ZeroVector(dimension);
		LoadOrthogonalDir = NodeB - NodeC;
		const double VectorNormC = MathUtils<double>::Norm(LoadOrthogonalDir);
		if (VectorNormC != 0.00) LoadOrthogonalDir /= VectorNormC;



		// now caluclate respective work equivilent nodal moments

		const double CustomMoment = NormForceVectorOrth *
			GeometryLength*GeometryLength / 12.00;

		Vector MomentNodeA = ZeroVector(dimension);
		MomentNodeA = MathUtils<double>::CrossProduct(GeometricOrientation,
			LoadOrthogonalDir);
		MomentNodeA *= CustomMoment;

		for (int i = 0; i < dimension; ++i)
		{
			rRightHandSideVector[(1 * dimension) + i] += MomentNodeA[i];
			rRightHandSideVector[(3 * dimension) + i] -= MomentNodeA[i];
		}

		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo) {

		KRATOS_TRY
			const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = NumNodes * dimension;
		const int LocalSize = NumNodes * dimension * 2;

		this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

		//Nodal element forces global
		Vector nodalForcesGlobal_q = ZeroVector(LocalSize);
		nodalForcesGlobal_q = prod(this->mRotationMatrix,
			this->mNodalForces);

		//create+compute RHS
		//update Residual
		rRightHandSideVector = ZeroVector(LocalSize);
		rRightHandSideVector -= nodalForcesGlobal_q;


		//LINEAR BEAM ELEMENT
		if (this->mIsLinearElement == true)
		{
			Vector NodalDeformation = ZeroVector(LocalSize);
			this->GetValuesVector(NodalDeformation);
			rRightHandSideVector = ZeroVector(LocalSize);
			rRightHandSideVector -= prod(rLeftHandSideMatrix, NodalDeformation);
		}
		//assign global element variables
		this->mLHS = rLeftHandSideMatrix;
		//add bodyforces 
		rRightHandSideVector += this->CalculateBodyForces();
		this->mIterationCount++;
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateRightHandSide(
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo) {

		KRATOS_TRY;
		const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = NumNodes * dimension;
		const int LocalSize = NumNodes * dimension * 2;
		rRightHandSideVector = ZeroVector(LocalSize);
		this->UpdateIncrementDeformation();
		Matrix TransformationMatrix = ZeroMatrix(LocalSize);
		this->CalculateTransformationMatrix(TransformationMatrix);
		Vector elementForces_t = ZeroVector(size);
		elementForces_t = this->CalculateElementForces();
		Vector nodalForcesLocal_qe = ZeroVector(LocalSize);
		Matrix TransformationMatrixS = ZeroMatrix(LocalSize, size);
		TransformationMatrixS = this->CalculateTransformationS();
		nodalForcesLocal_qe = prod(TransformationMatrixS,
			elementForces_t);
		//save local nodal forces
		this->mNodalForces = ZeroVector(LocalSize);
		this->mNodalForces = nodalForcesLocal_qe;

		Vector nodalForcesGlobal_q = ZeroVector(LocalSize);
		nodalForcesGlobal_q = prod(TransformationMatrix, nodalForcesLocal_qe);
		rRightHandSideVector -= nodalForcesGlobal_q;


		//LINEAR BEAM ELEMENT
		if (this->mIsLinearElement == true)
		{
			Matrix LeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
			this->CalculateLeftHandSide(LeftHandSideMatrix,rCurrentProcessInfo);
			Vector NodalDeformation = ZeroVector(LocalSize);
			this->GetValuesVector(NodalDeformation);
			rRightHandSideVector = ZeroVector(LocalSize);
			rRightHandSideVector -= prod(LeftHandSideMatrix, NodalDeformation);
		}

		//add bodyforces 
		rRightHandSideVector += this->CalculateBodyForces();
		KRATOS_CATCH("")

	}

	void CrBeamElement3D2N::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
		ProcessInfo& rCurrentProcessInfo) {

		KRATOS_TRY
			const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = NumNodes * dimension;
		const int LocalSize = NumNodes * dimension * 2;

		//update displacement_delta
		this->UpdateIncrementDeformation();

		//calculate Transformation Matrix
		Matrix TransformationMatrix = ZeroMatrix(LocalSize);
		this->CalculateTransformationMatrix(TransformationMatrix);
		this->mRotationMatrix = ZeroMatrix(LocalSize);
		this->mRotationMatrix = TransformationMatrix;

		//deformation modes
		Vector elementForces_t = ZeroVector(size);
		elementForces_t = this->CalculateElementForces();

		//Nodal element forces local
		Vector nodalForcesLocal_qe = ZeroVector(LocalSize);
		Matrix TransformationMatrixS = ZeroMatrix(LocalSize, size);
		TransformationMatrixS = this->CalculateTransformationS();
		nodalForcesLocal_qe = prod(TransformationMatrixS, elementForces_t);

		//save local nodal forces
		this->mNodalForces = ZeroVector(LocalSize);
		this->mNodalForces = nodalForcesLocal_qe;

		//resizing the matrices + create memory for LHS
		rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		//creating LHS
		rLeftHandSideMatrix +=
			this->CreateElementStiffnessMatrix_Material();
		rLeftHandSideMatrix +=
			this->CreateElementStiffnessMatrix_Geometry(nodalForcesLocal_qe);


		Matrix aux_matrix = ZeroMatrix(LocalSize);
		aux_matrix = prod(TransformationMatrix, rLeftHandSideMatrix);
		rLeftHandSideMatrix = prod(aux_matrix,
			Matrix(trans(TransformationMatrix)));

		//LINEAR BEAM ELEMENT
		if (this->mIsLinearElement == true)
		{
			TransformationMatrix = this->mRotationMatrix0;
			rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
			rLeftHandSideMatrix +=
				this->CreateElementStiffnessMatrix_Material();
			aux_matrix = ZeroMatrix(LocalSize);
			aux_matrix = prod(TransformationMatrix, rLeftHandSideMatrix);
			rLeftHandSideMatrix = prod(aux_matrix,
				Matrix(trans(TransformationMatrix)));
		}
		KRATOS_CATCH("")
	}

	Vector CrBeamElement3D2N::CalculateElementForces() {

		KRATOS_TRY
			const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int LocalSize = NumNodes * dimension;

		Vector deformation_modes_total_V = ZeroVector(LocalSize);
		deformation_modes_total_V[3] = this->mCurrentLength - this->mLength;
		for (int i = 0; i < 3; ++i) deformation_modes_total_V[i] = this->mPhiS[i];
		for (int i = 0; i < 2; ++i) deformation_modes_total_V[i + 4] = this->mPhiA[i + 1];
		//calculate element forces
		Vector element_forces_t = ZeroVector(LocalSize);
		Matrix deformation_stiffness_Kd = ZeroMatrix(LocalSize);

		deformation_stiffness_Kd = this->CalculateDeformationStiffness();
		element_forces_t = prod(deformation_stiffness_Kd,
			deformation_modes_total_V);

		return element_forces_t;
		KRATOS_CATCH("")
	}

	double CrBeamElement3D2N::CalculateCurrentLength() {

		KRATOS_TRY
			const double du = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
		const double dv = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
		const double dw = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
		const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
		const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
		const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
		const double l = sqrt((du + dx)*(du + dx) + (dv + dy)*(dv + dy) +
			(dw + dz)*(dw + dz));
		return l;
		KRATOS_CATCH("")

	}
	double CrBeamElement3D2N::CalculatePsi(const double I, const double A_eff) {

		KRATOS_TRY
			const double E = this->mYoungsModulus;
		const double L = this->mCurrentLength;
		const double G = this->mShearModulus;

		const double phi = (12.0 * E * I) / (L*L * G*A_eff);
		double psi;
		//interpret input A_eff == 0 as shearstiff -> psi = 1.0
		if (A_eff == 0.00) psi = 1.00;
		else psi = 1.0 / (1.0 + phi);

		return psi;
		KRATOS_CATCH("")
	}

	double CrBeamElement3D2N::CalculateReferenceLength() {

		KRATOS_TRY
			const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
		const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
		const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
		const double L = sqrt(dx*dx + dy*dy + dz*dz);
		return L;
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::UpdateIncrementDeformation() {

		KRATOS_TRY
			const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = NumNodes * dimension;
		const int LocalSize = NumNodes * dimension * 2;

		Vector actualDeformation = ZeroVector(LocalSize);
		Vector total_nodal_def = ZeroVector(LocalSize);
		Vector total_nodal_pos = ZeroVector(size);
		this->mIncrementDeformation = ZeroVector(LocalSize);

		if (mIterationCount == 0) this->mTotalNodalDeformation = ZeroVector(LocalSize);
		this->GetValuesVector(actualDeformation, 0);
		this->mIncrementDeformation = actualDeformation
			- this->mTotalNodalDeformation;

		this->mTotalNodalDeformation = ZeroVector(LocalSize);
		this->mTotalNodalDeformation = actualDeformation;

		this->mTotalNodalPosistion = ZeroVector(size);
		this->mTotalNodalPosistion[0] = this->GetGeometry()[0].X0()
			+ actualDeformation[0];
		this->mTotalNodalPosistion[1] = this->GetGeometry()[0].Y0()
			+ actualDeformation[1];
		this->mTotalNodalPosistion[2] = this->GetGeometry()[0].Z0()
			+ actualDeformation[2];

		this->mTotalNodalPosistion[3] = this->GetGeometry()[1].X0()
			+ actualDeformation[6];
		this->mTotalNodalPosistion[4] = this->GetGeometry()[1].Y0()
			+ actualDeformation[7];
		this->mTotalNodalPosistion[5] = this->GetGeometry()[1].Z0()
			+ actualDeformation[8];

		this->mCurrentLength = this->CalculateCurrentLength();
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateOnIntegrationPoints(
		const Variable<array_1d<double, 3 > >& rVariable,
		std::vector< array_1d<double, 3 > >& rOutput,
		const ProcessInfo& rCurrentProcessInfo) {

		KRATOS_TRY
			const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = NumNodes * dimension;
		const int LocalSize = NumNodes * dimension * 2;


		//element with two nodes can only represent results at one node 
		const int&  write_points_number = GetGeometry()
			.IntegrationPointsNumber(Kratos::GeometryData::GI_GAUSS_3);
		if (rOutput.size() != write_points_number) {
			rOutput.resize(write_points_number);
		}


		this->UpdateIncrementDeformation();
		//calculate Transformation Matrix
		Matrix TransformationMatrix = ZeroMatrix(LocalSize);
		this->CalculateTransformationMatrix(TransformationMatrix);
		//deformation modes
		Vector elementForces_t = ZeroVector(size);
		elementForces_t = this->CalculateElementForces();
		Vector Stress = ZeroVector(LocalSize);
		Matrix TransformationMatrixS = ZeroMatrix(LocalSize, size);
		TransformationMatrixS = this->CalculateTransformationS();
		Stress = prod(TransformationMatrixS, elementForces_t);

		//LINEAR BEAM ELEMENT
		if (this->mIsLinearElement == true)
		{
			Matrix LeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
			LeftHandSideMatrix = this->mLHS;

			Vector NodalDeformation = ZeroVector(LocalSize);
			this->GetValuesVector(NodalDeformation);
			Stress = ZeroVector(LocalSize);
			Stress = prod(LeftHandSideMatrix, NodalDeformation);
			Matrix TransformationMatrix = ZeroMatrix(LocalSize);
			TransformationMatrix = this->mRotationMatrix;
			Stress = prod(Matrix(trans(TransformationMatrix)), Stress);
		}


		//rOutput[GP 1,2,3][x,y,z]

		if (rVariable == MOMENT)
		{
			rOutput[0][0] = -1.0 *Stress[3] * 0.75 + Stress[9] * 0.25;
			rOutput[1][0] = -1.0 *Stress[3] * 0.50 + Stress[9] * 0.50;
			rOutput[2][0] = -1.0 *Stress[3] * 0.25 + Stress[9] * 0.75;

			rOutput[0][1] = -1.0 *Stress[4] * 0.75 + Stress[10] * 0.25;
			rOutput[1][1] = -1.0 *Stress[4] * 0.50 + Stress[10] * 0.50;
			rOutput[2][1] = -1.0 *Stress[4] * 0.25 + Stress[10] * 0.75;

			rOutput[0][2] = -1.0 *Stress[5] * 0.75 + Stress[11] * 0.25;
			rOutput[1][2] = -1.0 *Stress[5] * 0.50 + Stress[11] * 0.50;
			rOutput[2][2] = -1.0 *Stress[5] * 0.25 + Stress[11] * 0.75;

		}
		if (rVariable == FORCE)
		{
			rOutput[0][0] = -1.0 * Stress[0] * 0.75 + Stress[6] * 0.25;
			rOutput[1][0] = -1.0 * Stress[0] * 0.50 + Stress[6] * 0.50;
			rOutput[2][0] = -1.0 * Stress[0] * 0.25 + Stress[6] * 0.75;

			rOutput[0][1] = -1.0 * Stress[1] * 0.75 + Stress[7] * 0.25;
			rOutput[1][1] = -1.0 *Stress[1] * 0.50 + Stress[7] * 0.50;
			rOutput[2][1] = -1.0 *Stress[1] * 0.25 + Stress[7] * 0.75;

			rOutput[0][2] = -1.0 *Stress[2] * 0.75 + Stress[8] * 0.25;
			rOutput[1][2] = -1.0 *Stress[2] * 0.50 + Stress[8] * 0.50;
			rOutput[2][2] = -1.0 *Stress[2] * 0.25 + Stress[8] * 0.75;

		}

		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::GetValueOnIntegrationPoints(
		const Variable<array_1d<double, 3 > >& rVariable,
		std::vector< array_1d<double, 3 > >& rOutput,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		this->CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
		KRATOS_CATCH("")
	}



	void CrBeamElement3D2N::CalculateOnIntegrationPoints(const Variable<Vector >& rVariable,
		std::vector< Vector >& rOutput,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		if (rVariable == LOCAL_AXES_VECTOR)
		{
			rOutput.resize(3);
			for (int i = 0; i < 3; ++i) rOutput[i] = ZeroVector(3);

			if (this->mIsLinearElement == true)
			{
				rOutput[0] = this->mNX0;
				rOutput[1] = this->mNY0;
				rOutput[2] = this->mNZ0;
			}
			else
			{
				rOutput[0] = this->mNX;
				rOutput[1] = this->mNY;
				rOutput[2] = this->mNZ;
			}
		}

		KRATOS_CATCH("");
	}

	void CrBeamElement3D2N::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
		std::vector<Vector>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		this->CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::AssembleSmallInBigMatrix(Matrix SmallMatrix,
		Matrix& BigMatrix) {

		KRATOS_TRY
			const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = number_of_nodes * dimension;
		const int MatSize = 2 * size;

		for (int kk = 0; kk < MatSize; kk += dimension)
		{
			for (int i = 0; i<dimension; ++i)
			{
				for (int j = 0; j<dimension; ++j)
				{
					BigMatrix(i + kk, j + kk) = SmallMatrix(i, j);
				}
			}
		}
		KRATOS_CATCH("")
	}


	void CrBeamElement3D2N::BuildSingleMassMatrix(MatrixType& rMassMatrix,
		double Phi, double CT, double CR, double L)
	{
		KRATOS_TRY;
		const int number_of_nodes = GetGeometry().PointsNumber();
		const int MatSize = number_of_nodes * 2;

		if (rMassMatrix.size1() != MatSize) {
			rMassMatrix.resize(MatSize, MatSize, false);
		}
		rMassMatrix = ZeroMatrix(MatSize, MatSize);
		Matrix TempMassMatrix = ZeroMatrix(MatSize, MatSize);
		const double Phi2 = Phi * Phi;
		const double L2 = L*L;


		TempMassMatrix(0, 0) = (13.00 / 35.00) + (7.00 / 10.00)*Phi
			+ (1.00 / 3.00)*Phi2;
		TempMassMatrix(0, 1) = ((11.00 / 210.00) + (11.00 / 210.00)*Phi
			+ (1.00 / 24.00)*Phi2)*L;
		TempMassMatrix(0, 2) = (9.00 / 70.00) + (3.00 / 10.00)*Phi
			+ (1.00 / 6.00)*Phi2;
		TempMassMatrix(0, 3) = -((13.00 / 420.00) + (3.00 / 40.00)*Phi
			+ (1.00 / 24.00)*Phi2)*L;
		TempMassMatrix(1, 0) = TempMassMatrix(0, 1);
		TempMassMatrix(1, 1) = ((1.00 / 105.00) + (1.00 / 60.00)*Phi
			+ (1.00 / 120.00)*Phi2)*L2;
		TempMassMatrix(1, 2) = ((13.00 / 420.00) + (3.00 / 40.00)*Phi
			+ (1.00 / 24.00)*Phi2)*L;
		TempMassMatrix(1, 3) = -((1.00 / 140.00) + (1.00 / 60.00)*Phi
			+ (1.00 / 120.00)*Phi2)*L2;
		TempMassMatrix(2, 0) = TempMassMatrix(0, 2);
		TempMassMatrix(2, 1) = TempMassMatrix(1, 2);
		TempMassMatrix(2, 2) = (13.00 / 35.00) + (7.00 / 10.00)*Phi
			+ (1.00 / 3.00)*Phi2;
		TempMassMatrix(2, 3) = -((11.00 / 210.00) + (11.00 / 210.00)*Phi
			+ (1.00 / 24.00)*Phi2)*L;
		TempMassMatrix(3, 0) = TempMassMatrix(0, 3);
		TempMassMatrix(3, 1) = TempMassMatrix(1, 3);
		TempMassMatrix(3, 2) = TempMassMatrix(2, 3);
		TempMassMatrix(3, 3) = ((1.00 / 105.00) + (1.00 / 60.00)*Phi
			+ (1.00 / 120.00)*Phi2)*L2;

		TempMassMatrix *= CT;
		rMassMatrix += TempMassMatrix;


		TempMassMatrix = ZeroMatrix(MatSize, MatSize);

		TempMassMatrix(0, 0) = 6.00 / 5.00;
		TempMassMatrix(0, 1) = ((1.00 / 10.00) - (1.00 / 2.00)*Phi)*L;
		TempMassMatrix(0, 2) = -6.00 / 5.00;
		TempMassMatrix(0, 3) = ((1.00 / 10.00) - (1.00 / 2.00)*Phi)*L;
		TempMassMatrix(1, 0) = TempMassMatrix(0, 1);
		TempMassMatrix(1, 1) = ((2.00 / 15.00) + (1.00 / 6.00)*Phi
			+ (1.00 / 3.00)*Phi2)*L2;
		TempMassMatrix(1, 2) = ((-1.00 / 10.00) + (1.00 / 2.00)*Phi)*L;
		TempMassMatrix(1, 3) = -((1.00 / 30.00) + (1.00 / 6.00)*Phi
			- (1.00 / 6.00)*Phi2)*L2;
		TempMassMatrix(2, 0) = TempMassMatrix(0, 2);
		TempMassMatrix(2, 1) = TempMassMatrix(1, 2);
		TempMassMatrix(2, 2) = 6.00 / 5.00;
		TempMassMatrix(2, 3) = ((-1.00 / 10.00) + (1.00 / 2.00)*Phi)*L;
		TempMassMatrix(3, 0) = TempMassMatrix(0, 3);
		TempMassMatrix(3, 1) = TempMassMatrix(1, 3);
		TempMassMatrix(3, 2) = TempMassMatrix(2, 3);
		TempMassMatrix(3, 3) = ((2.00 / 15.00) + (1.00 / 6.00)*Phi
			+ (1.00 / 3.00)*Phi2)*L2;

		TempMassMatrix *= CR;
		rMassMatrix += TempMassMatrix;
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateConsistentMassMatrix(MatrixType& rMassMatrix,
		ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		const int number_of_nodes = GetGeometry().PointsNumber();
		const int dimension = GetGeometry().WorkingSpaceDimension();
		const int smallMatSize = number_of_nodes * 2;
		const int MatSize = number_of_nodes * dimension * 2;

		if (rMassMatrix.size1() != MatSize) {
			rMassMatrix.resize(MatSize, MatSize, false);
		}
		rMassMatrix = ZeroMatrix(MatSize, MatSize);

		const double L = this->mLength;
		const double L2 = L * L;
		const double rho = this->mDensity;
		const double A = this->mArea;
		const double E = this->mYoungsModulus;
		const double Iy = this->mInertiaY;
		const double Iz = this->mInertiaZ;
		const double J = this->mInertiaX;
		const double G = this->mShearModulus;
		const double Ay = this->mEffAreaY;
		const double Az = this->mEffAreaZ;
		const double IRy = this->mRotInertiaY;
		const double IRz = this->mRotInertiaZ;

		double Phiy = 0.00;
		double Phiz = 0.00;

		if (Ay != 0.00) Phiz = (12.00 * E * Iz) / (L2*G*Ay);
		if (Az != 0.00) Phiy = (12.00 * E * Iy) / (L2*G*Az);

		const double CTy = (rho * A * L) / ((1 + Phiy)*(1 + Phiy));
		const double CTz = (rho * A * L) / ((1 + Phiz)*(1 + Phiz));

		const double CRy = (rho*IRy) / ((1 + Phiy)*(1 + Phiy)*L);
		const double CRz = (rho*IRz) / ((1 + Phiz)*(1 + Phiz)*L);

		//longitudinal forces + torsional moment
		const double M00 = (1.00 / 3.00)*A*rho*L;
		const double M06 = M00 / 2.00;
		rMassMatrix(0, 0) = M00;
		rMassMatrix(0, 6) = M06;
		rMassMatrix(6, 6) = M00;
		rMassMatrix(3, 3) = M00;
		rMassMatrix(3, 9) = M06;
		rMassMatrix(9, 9) = M00;

		Matrix TempBendingMassMatrix = ZeroMatrix(smallMatSize, smallMatSize);
		this->BuildSingleMassMatrix(TempBendingMassMatrix, Phiz, CTz, CRz, L);

		rMassMatrix(1, 1) = TempBendingMassMatrix(0, 0);
		rMassMatrix(1, 5) = TempBendingMassMatrix(0, 1);
		rMassMatrix(1, 7) = TempBendingMassMatrix(0, 2);
		rMassMatrix(1, 11) = TempBendingMassMatrix(0, 3);
		rMassMatrix(5, 5) = TempBendingMassMatrix(1, 1);
		rMassMatrix(5, 7) = TempBendingMassMatrix(1, 2);
		rMassMatrix(5, 11) = TempBendingMassMatrix(1, 3);
		rMassMatrix(7, 7) = TempBendingMassMatrix(2, 2);
		rMassMatrix(7, 11) = TempBendingMassMatrix(2, 3);
		rMassMatrix(11, 11) = TempBendingMassMatrix(3, 3);

		TempBendingMassMatrix = ZeroMatrix(smallMatSize, smallMatSize);
		this->BuildSingleMassMatrix(TempBendingMassMatrix, Phiy, CTy, CRy, L);

		rMassMatrix(2, 2) = TempBendingMassMatrix(0, 0);
		rMassMatrix(2, 4) = TempBendingMassMatrix(0, 1);
		rMassMatrix(2, 8) = TempBendingMassMatrix(0, 2);
		rMassMatrix(2, 10) = TempBendingMassMatrix(0, 3);
		rMassMatrix(4, 4) = TempBendingMassMatrix(1, 1);
		rMassMatrix(4, 8) = TempBendingMassMatrix(1, 2);
		rMassMatrix(4, 10) = TempBendingMassMatrix(1, 3);
		rMassMatrix(8, 8) = TempBendingMassMatrix(2, 2);
		rMassMatrix(8, 10) = TempBendingMassMatrix(2, 3);
		rMassMatrix(10, 10) = TempBendingMassMatrix(3, 3);


		for (int j = 1; j < 12; ++j)
		{
			for (int i = 0; i < j; ++i)
			{
				rMassMatrix(j, i) = rMassMatrix(i, j);
			}
		}

		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateLumpedMassMatrix(MatrixType& rMassMatrix,
		ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		const int number_of_nodes = GetGeometry().PointsNumber();
		const int dimension = GetGeometry().WorkingSpaceDimension();
		const int MatSize = number_of_nodes * dimension * 2;

		if (rMassMatrix.size1() != MatSize) {
			rMassMatrix.resize(MatSize, MatSize, false);
		}
		rMassMatrix = ZeroMatrix(MatSize, MatSize);

		const double TotalMass = this->mArea * this->mLength * this->mDensity;
		const double temp = 0.50 * TotalMass;

		//translatonal mass	
		for (int i = 0; i < number_of_nodes; ++i)
		{
			for (int j = 0; j < dimension; ++j)
			{
				int index = i * (dimension * 2) + j;
				rMassMatrix(index, index) = temp;
			}
		}
		//rotaional mass neglected alpha = 0
		KRATOS_CATCH("")
	}

	CrBeamElement3D2N::IntegrationMethod
		CrBeamElement3D2N::GetIntegrationMethod() const
	{
		//do this to have 3GP as an output in GID
		return Kratos::GeometryData::GI_GAUSS_3;
	}

	void CrBeamElement3D2N::AddExplicitContribution(const VectorType& rRHSVector,
		const Variable<VectorType>& rRHSVariable,
		Variable<array_1d<double, 3> >& rDestinationVariable,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		/////////////////////////////////////////////////
		//NOT FINISHED YET !!!!
		/////////////////////////////////////////////////
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int element_size = number_of_nodes * dimension;
		const int mat_size = element_size * 2;

		if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL)
		{

			for (int i = 0; i< number_of_nodes; ++i)
			{
				int index = element_size * i;

				GetGeometry()[i].SetLock();

				array_1d<double, 3> &ForceResidual =
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


	int  CrBeamElement3D2N::Check(const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

			if (GetGeometry().WorkingSpaceDimension() != 3 || GetGeometry().size() != 2)
			{
				KRATOS_ERROR <<
					("The truss element works only in 3D and with 2 noded elements", "")
					<< std::endl;
			}
		//verify that the variables are correctly initialized
		if (VELOCITY.Key() == 0) {
			KRATOS_ERROR <<
				("VELOCITY has Key zero! (check if the application is correctly registered", "")
				<< std::endl;
		}
		if (DISPLACEMENT.Key() == 0) {
			KRATOS_ERROR <<
				("DISPLACEMENT has Key zero! (check if the application is correctly registered", "")
				<< std::endl;
		}
		if (ACCELERATION.Key() == 0) {
			KRATOS_ERROR <<
				("ACCELERATION has Key zero! (check if the application is correctly registered", "")
				<< std::endl;
		}
		if (DENSITY.Key() == 0) {
			KRATOS_ERROR <<
				("DENSITY has Key zero! (check if the application is correctly registered", "")
				<< std::endl;
		}
		if (CROSS_AREA.Key() == 0) {
			KRATOS_ERROR <<
				("CROSS_AREA has Key zero! (check if the application is correctly registered", "")
				<< std::endl;
		}
		//verify that the dofs exist
		for (int i = 0; i<this->GetGeometry().size(); ++i)
		{
			if (this->GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false) {
				KRATOS_ERROR <<
					("missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id())
					<< std::endl;
			}
			if (this->GetGeometry()[i].HasDofFor(DISPLACEMENT_X) == false ||
				this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false ||
				this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Z) == false) {
				KRATOS_ERROR <<
					("missing one of the dofs for the variable DISPLACEMENT on node ",
						GetGeometry()[i].Id()) << std::endl;
			}
		}



		if (this->GetProperties().Has(CROSS_AREA) == false ||
			this->GetProperties()[CROSS_AREA] == 0)
		{
			KRATOS_ERROR << ("CROSS_AREA not provided for this element", this->Id())
				<< std::endl;
		}

		if (this->GetProperties().Has(YOUNG_MODULUS) == false ||
			this->GetProperties()[YOUNG_MODULUS] == 0)
		{
			KRATOS_ERROR << ("YOUNG_MODULUS not provided for this element", this->Id())
				<< std::endl;
		}
		if (this->GetProperties().Has(DENSITY) == false)
		{
			KRATOS_ERROR << ("DENSITY not provided for this element", this->Id())
				<< std::endl;
		}


		return 0;

		KRATOS_CATCH("")
	}

	Orientation::Orientation(array_1d<double, 3>& v1, const double theta) {

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

	Orientation::Orientation(array_1d<double, 3>& v1, array_1d<double, 3>& v2) {

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

	void Orientation::CalculateRotationMatrix(Matrix& R) {

		KRATOS_TRY
			if (R.size1() != 3 || R.size2() != 3) R.resize(3, 3, false);
		const Quaternion<double> q = this->GetQuaternion();
		q.ToRotationMatrix(R);
		KRATOS_CATCH("")
	}

	void Orientation::CalculateBasisVectors(array_1d<double, 3>& v1,
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


	void CrBeamElement3D2N::save(Serializer& rSerializer) const
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
		rSerializer.save("DeformationModes", this->mDeformationModes);
		rSerializer.save("NodalPosition", this->mTotalNodalPosistion);
		rSerializer.save("NodalDeformation", this->mTotalNodalDeformation);
		rSerializer.save("IterationCounter", this->mIterationCount);
		rSerializer.save("NodalForces", this->mNodalForces);


		rSerializer.save("LocalInitalAxisX", this->mNX0);
		rSerializer.save("LocalInitalAxisY", this->mNY0);
		rSerializer.save("LocalInitalAxisZ", this->mNZ0);


		rSerializer.save("QuaternionVecA", this->mQuaternionVEC_A);
		rSerializer.save("QuaternionVecB", this->mQuaternionVEC_B);
		rSerializer.save("QuaternionScaA", this->mQuaternionSCA_A);
		rSerializer.save("QuaternionScaB", this->mQuaternionSCA_B);

		rSerializer.save("Poisson", this->mPoisson);
		rSerializer.save("Area", this->mArea);
		rSerializer.save("YoungsModulus", this->mYoungsModulus);
		rSerializer.save("ShearModulus", this->mShearModulus);
		rSerializer.save("LengthRef", this->mLength);
		rSerializer.save("LengthCurr", this->mCurrentLength);
		rSerializer.save("Density", this->mDensity);
		rSerializer.save("InertiaIT", this->mInertiaX);
		rSerializer.save("InertiaIY", this->mInertiaY);
		rSerializer.save("InertiaIZ", this->mInertiaZ);
		rSerializer.save("EffAreaY", this->mEffAreaY);
		rSerializer.save("EffAreaZ", this->mEffAreaZ);
		rSerializer.save("RotInertiaY", this->mRotInertiaY);
		rSerializer.save("RotInertiaZ", this->mRotInertiaZ);
		rSerializer.save("PsiY", this->mPsiY);
		rSerializer.save("PsiZ", this->mPsiZ);
		rSerializer.save("ThetaInitialRot", this->mtheta);
		rSerializer.save("mIsLumpedMassMatrix", this->mIsLumpedMassMatrix);
	}

	void CrBeamElement3D2N::load(Serializer& rSerializer)
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
		rSerializer.load("DeformationModes", this->mDeformationModes);
		rSerializer.load("NodalPosition", this->mTotalNodalPosistion);
		rSerializer.load("NodalDeformation", this->mTotalNodalDeformation);
		rSerializer.load("IterationCounter", this->mIterationCount);
		rSerializer.load("NodalForces", this->mNodalForces);

		rSerializer.load("LocalInitalAxisX", this->mNX0);
		rSerializer.load("LocalInitalAxisY", this->mNY0);
		rSerializer.load("LocalInitalAxisZ", this->mNZ0);


		rSerializer.load("QuaternionVecA", this->mQuaternionVEC_A);
		rSerializer.load("QuaternionVecB", this->mQuaternionVEC_B);
		rSerializer.load("QuaternionScaA", this->mQuaternionSCA_A);
		rSerializer.load("QuaternionScaB", this->mQuaternionSCA_B);

		rSerializer.load("Poisson", this->mPoisson);
		rSerializer.load("Area", this->mArea);
		rSerializer.load("YoungsModulus", this->mYoungsModulus);
		rSerializer.load("ShearModulus", this->mShearModulus);
		rSerializer.load("LengthRef", this->mLength);
		rSerializer.load("LengthCurr", this->mCurrentLength);
		rSerializer.load("Density", this->mDensity);
		rSerializer.load("InertiaIT", this->mInertiaX);
		rSerializer.load("InertiaIY", this->mInertiaY);
		rSerializer.load("InertiaIZ", this->mInertiaZ);
		rSerializer.load("EffAreaY", this->mEffAreaY);
		rSerializer.load("EffAreaZ", this->mEffAreaZ);
		rSerializer.load("RotInertiaY", this->mRotInertiaY);
		rSerializer.load("RotInertiaZ", this->mRotInertiaZ);
		rSerializer.load("PsiY", this->mPsiY);
		rSerializer.load("PsiZ", this->mPsiZ);
		rSerializer.load("ThetaInitialRot", this->mtheta);
		rSerializer.load("mIsLumpedMassMatrix", this->mIsLumpedMassMatrix);
	}

} // namespace Kratos.


