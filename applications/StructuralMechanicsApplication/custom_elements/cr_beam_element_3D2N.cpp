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
								GeometryType::Pointer pGeometry)
								: Element(NewId, pGeometry) {}

	CrBeamElement3D2N::CrBeamElement3D2N(IndexType NewId,
								GeometryType::Pointer pGeometry,
								PropertiesType::Pointer pProperties) 
								: Element(NewId, pGeometry, pProperties) {}

	Element::Pointer CrBeamElement3D2N::Create(IndexType NewId,
								NodesArrayType const& rThisNodes,
								PropertiesType::Pointer pProperties) const{
		const GeometryType& rGeom = this->GetGeometry();
		return BaseType::Pointer(new CrBeamElement3D2N(
				NewId, rGeom.Create(rThisNodes), pProperties));
	}

	CrBeamElement3D2N::~CrBeamElement3D2N() {}

	void CrBeamElement3D2N::EquationIdVector(EquationIdVectorType& rResult,
								ProcessInfo& rCurrentProcessInfo){

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
								ProcessInfo& rCurrentProcessInfo){

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

	void CrBeamElement3D2N::Initialize(){

		KRATOS_TRY
		this->mPoisson = this->GetProperties()[POISSON_RATIO];
		this->mArea = this->GetProperties()[CROSS_AREA];
		this->mYoungsModulus = this->GetProperties()[YOUNG_MODULUS];
		this->mShearModulus = this->mYoungsModulus / (2.0 * (1.0 + this->mPoisson));
		this->mLength = this->CalculateReferenceLength();
		this->mCurrentLength = this->CalculateCurrentLength();
		this->mDensity = this->GetProperties()[DENSITY];
		this->mInertiaX = this->GetProperties()[IX];
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

		KRATOS_CATCH("")
	}

	CrBeamElement3D2N::MatrixType 
					CrBeamElement3D2N::CreateElementStiffnessMatrix_Material(){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int local_size = number_of_nodes * dimension * 2;

		const double E = this->mYoungsModulus;
		const double G = this->mShearModulus;
		const double A = this->mArea;
		const double L = this->mCurrentLength;
		const double J = this->mInertiaX;
		const double Iy = this->mInertiaY;
		const double Iz = this->mInertiaZ;

		const double Psi_y = this->mPsiY;
		const double Psi_z = this->mPsiZ;

		MatrixType LocalStiffnessMatrix = ZeroMatrix(local_size, local_size);
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

	CrBeamElement3D2N::MatrixType 
		CrBeamElement3D2N::CreateElementStiffnessMatrix_Geometry(VectorType qe){

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

		MatrixType LocalStiffnessMatrix = ZeroMatrix(local_size, local_size);

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

		for (int i = 0; i < 3; i++) {
			LocalStiffnessMatrix(3, i) = LocalStiffnessMatrix(i, 3);
		}
		LocalStiffnessMatrix(3, 4) = (-mz_A / 3.00) + (mz_B / 6.00);
		LocalStiffnessMatrix(3, 5) = (my_A / 3.00) - (my_B / 6.00);
		LocalStiffnessMatrix(3, 10) = L*Qy / 6.00;
		LocalStiffnessMatrix(3, 11) = L*Qz / 6.00;

		for (int i = 0; i < 4; i++) {
			LocalStiffnessMatrix(4, i) = LocalStiffnessMatrix(i, 4);
		}
		LocalStiffnessMatrix(4, 4) = 2.00 * L*N / 15.00;
		LocalStiffnessMatrix(4, 6) = -1.00 * LocalStiffnessMatrix(3, 1);
		LocalStiffnessMatrix(4, 7) = LocalStiffnessMatrix(2, 11);
		LocalStiffnessMatrix(4, 8) = LocalStiffnessMatrix(2, 10);
		LocalStiffnessMatrix(4, 9) = LocalStiffnessMatrix(3, 10);
		LocalStiffnessMatrix(4, 10) = -L*N / 30.00;
		LocalStiffnessMatrix(4, 11) = Mt / 2.00;


		for (int i = 0; i < 5; i++) {
			LocalStiffnessMatrix(5, i) = LocalStiffnessMatrix(i, 5);
		}
		LocalStiffnessMatrix(5, 5) = LocalStiffnessMatrix(4, 4);
		LocalStiffnessMatrix(5, 6) = -1.00 * LocalStiffnessMatrix(2, 3);
		LocalStiffnessMatrix(5, 7) = LocalStiffnessMatrix(1, 5);
		LocalStiffnessMatrix(5, 8) = LocalStiffnessMatrix(4, 7);
		LocalStiffnessMatrix(5, 9) = LocalStiffnessMatrix(3, 11);
		LocalStiffnessMatrix(5, 10) = -1.00 * LocalStiffnessMatrix(4, 11);
		LocalStiffnessMatrix(5, 11) = LocalStiffnessMatrix(4, 10);

		for (int i = 0; i < 6; i++) {
			LocalStiffnessMatrix(6, i) = LocalStiffnessMatrix(i, 6);
		}
		LocalStiffnessMatrix(6, 7) = LocalStiffnessMatrix(0, 1);
		LocalStiffnessMatrix(6, 8) = LocalStiffnessMatrix(0, 2);

		for (int i = 0; i < 7; i++) {
			LocalStiffnessMatrix(7, i) = LocalStiffnessMatrix(i, 7);
		}
		LocalStiffnessMatrix(7, 7) = LocalStiffnessMatrix(1, 1);
		LocalStiffnessMatrix(7, 9) = -1.00 * LocalStiffnessMatrix(1, 9);
		LocalStiffnessMatrix(7, 10) = LocalStiffnessMatrix(4, 1);
		LocalStiffnessMatrix(7, 11) = LocalStiffnessMatrix(2, 4);

		for (int i = 0; i < 8; i++) {
			LocalStiffnessMatrix(8, i) = LocalStiffnessMatrix(i, 8);
		}
		LocalStiffnessMatrix(8, 8) = LocalStiffnessMatrix(1, 1);
		LocalStiffnessMatrix(8, 9) = -1.00 * LocalStiffnessMatrix(2, 9);
		LocalStiffnessMatrix(8, 10) = LocalStiffnessMatrix(1, 5);
		LocalStiffnessMatrix(8, 11) = LocalStiffnessMatrix(1, 4);

		for (int i = 0; i < 9; i++) {
			LocalStiffnessMatrix(9, i) = LocalStiffnessMatrix(i, 9);
		}
		LocalStiffnessMatrix(9, 10) = (mz_A / 6.00) - (mz_B / 3.00);
		LocalStiffnessMatrix(9, 11) = (-my_A / 6.00) + (my_B / 3.00);

		for (int i = 0; i < 10; i++) {
			LocalStiffnessMatrix(10, i) = LocalStiffnessMatrix(i, 10);
		}
		LocalStiffnessMatrix(10, 10) = LocalStiffnessMatrix(4, 4);

		for (int i = 0; i < 11; i++) {
			LocalStiffnessMatrix(11, i) = LocalStiffnessMatrix(i, 11);
		}
		LocalStiffnessMatrix(11, 11) = LocalStiffnessMatrix(4, 4);

		return LocalStiffnessMatrix;
		KRATOS_CATCH("")
	}

	CrBeamElement3D2N::MatrixType CrBeamElement3D2N::CalculateMaterialStiffness(){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int local_size = number_of_nodes * dimension;

		MatrixType Kd = ZeroMatrix(local_size, local_size);

		const double E = this->mYoungsModulus;
		const double G = this->mShearModulus;
		const double A = this->mArea;
		const double L = this->mCurrentLength;
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

		return Kd;
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateInitialLocalCS(){

		KRATOS_TRY
		//Transformation matrix T = [e1_local, e2_local, e3_local]
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = number_of_nodes * dimension;
		const int MatSize = 2 * size;

		Vector DirectionVectorX = ZeroVector(3);
		Vector DirectionVectorY = ZeroVector(3);
		Vector DirectionVectorZ = ZeroVector(3);
		Vector ReferenceCoordinates = ZeroVector(size);
		Vector GlobalZ = ZeroVector(3);
		GlobalZ[2] = 1.0;

		ReferenceCoordinates[0] = this->GetGeometry()[0].X0();
		ReferenceCoordinates[1] = this->GetGeometry()[0].Y0();
		ReferenceCoordinates[2] = this->GetGeometry()[0].Z0();
		ReferenceCoordinates[3] = this->GetGeometry()[1].X0();
		ReferenceCoordinates[4] = this->GetGeometry()[1].Y0();
		ReferenceCoordinates[5] = this->GetGeometry()[1].Z0();

		for (unsigned int i = 0; i < dimension; i++)
		{
			DirectionVectorX[i] = (ReferenceCoordinates[i + 3] 
				- ReferenceCoordinates[i]);
		}
		// local x-axis (e1_local) is the beam axis  (in GID is e3_local)
		double VectorNorm;
		VectorNorm = MathUtils<double>::Norm(DirectionVectorX);
		if (VectorNorm != 0) DirectionVectorX /= VectorNorm;

		const double tolerance = 1.0 / 1000.0;
		if ((fabs(DirectionVectorX[0]) < tolerance) &&
			(fabs(DirectionVectorX[1]) < tolerance)) {
			DirectionVectorX = ZeroVector(3);
			DirectionVectorX[2] = 1.0;
			DirectionVectorY[1] = 1.0;
			DirectionVectorZ[0] = -1.0;
			if (ReferenceCoordinates[2] > ReferenceCoordinates[5]) {
				DirectionVectorX *= -1.0;
				DirectionVectorZ *= -1.0;
			}
		}
		else {
			DirectionVectorY = MathUtils<double>::CrossProduct(GlobalZ,
				DirectionVectorX);
			VectorNorm = MathUtils<double>::Norm(DirectionVectorY);
			DirectionVectorY /= VectorNorm;

			DirectionVectorZ = MathUtils<double>::CrossProduct(DirectionVectorX,
				DirectionVectorY);
			VectorNorm = MathUtils<double>::Norm(DirectionVectorZ);
			DirectionVectorZ /= VectorNorm;
		}

		//manual rotation around the beam axis
		double theta = this->mtheta;
		Vector nz = DirectionVectorZ;
		Vector ny = DirectionVectorY;
		if (theta != 0) {
			DirectionVectorY = ny * cos(theta) + nz * sin(theta);
			DirectionVectorY /= MathUtils<double>::Norm(DirectionVectorY);

			DirectionVectorZ = nz * cos(theta) - ny * sin(theta);
			DirectionVectorZ /= MathUtils<double>::Norm(DirectionVectorZ);
		}

		this->mNX0 = ZeroVector(3);
		this->mNY0 = ZeroVector(3);
		this->mNZ0 = ZeroVector(3);
		this->mNX0 = DirectionVectorX;
		this->mNY0 = DirectionVectorY;
		this->mNZ0 = DirectionVectorZ;

		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateTransformationMatrix(Matrix& rRotationMatrix){

		KRATOS_TRY
		//12x12
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = number_of_nodes * dimension;
		const int MatSize = 2 * size;

		//initialize local CS
		if(this->mIterationCount == 0) this->CalculateInitialLocalCS();

		//update local CS
		Matrix AuxRotationMatrix = ZeroMatrix(dimension);
		AuxRotationMatrix = this->UpdateRotationMatrixLocal();

		if (rRotationMatrix.size1() != MatSize) {
			rRotationMatrix.resize(MatSize, MatSize, false);
		}

		rRotationMatrix = ZeroMatrix(MatSize);

		//Building the rotation matrix for the local element matrix
		for (int kk = 0; kk < MatSize; kk += dimension)
		{
			for (int i = 0; i<dimension; i++)
			{
				for (int j = 0; j<dimension; j++)
				{
					rRotationMatrix(i + kk, j + kk) = AuxRotationMatrix(i, j);
				}
			}
		}
		mRotationMatrix = rRotationMatrix;
		KRATOS_CATCH("")
	}

	CrBeamElement3D2N::MatrixType CrBeamElement3D2N::CalculateTransformationS(){

		KRATOS_TRY
		VectorType n_x = ZeroVector(3);
		VectorType n_y = ZeroVector(3);
		VectorType n_z = ZeroVector(3);
		n_x[0] = 1.0;
		n_y[1] = 1.0;
		n_z[2] = 1.0;

		const double L = this->mCurrentLength;
		MatrixType S = ZeroMatrix(12, 6);
		for (int i = 0; i < 3; i++) S(3 + i, 0) = -1.0* n_x[i];
		for (int i = 0; i < 3; i++) S(9 + i, 0) = n_x[i];
		for (int i = 0; i < 3; i++) S(3 + i, 1) = -1.0* n_y[i];
		for (int i = 0; i < 3; i++) S(9 + i, 1) = n_y[i];
		for (int i = 0; i < 3; i++) S(3 + i, 2) = -1.0* n_z[i];
		for (int i = 0; i < 3; i++) S(9 + i, 2) = n_z[i];
		for (int i = 0; i < 3; i++) S(0 + i, 3) = -1.0* n_x[i];
		for (int i = 0; i < 3; i++) S(6 + i, 3) = n_x[i];
		for (int i = 0; i < 3; i++) S(3 + i, 4) = n_y[i];
		for (int i = 0; i < 3; i++) S(9 + i, 4) = n_y[i];
		for (int i = 0; i < 3; i++) S(3 + i, 5) = n_z[i];
		for (int i = 0; i < 3; i++) S(9 + i, 5) = n_z[i];

		for (int i = 0; i < 3; i++) S(0 + i, 4) = -2.0 * n_z[i] / L;
		for (int i = 0; i < 3; i++) S(6 + i, 4) = 2.0 * n_z[i] / L;
		for (int i = 0; i < 3; i++) S(0 + i, 5) = 2.0 * n_y[i] / L;
		for (int i = 0; i < 3; i++) S(6 + i, 5) = -2.0 * n_y[i] / L;

		return S;
		KRATOS_CATCH("")
	}

	CrBeamElement3D2N::MatrixType CrBeamElement3D2N::UpdateRotationMatrixLocal(){

		KRATOS_TRY
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		VectorType dPhiA = ZeroVector(dimension);
		VectorType dPhiB = ZeroVector(dimension);
		Vector IncrementDeformation = ZeroVector(12);
		IncrementDeformation = this->mIncrementDeformation;

		for (int i = 0; i < dimension; i++) {
			dPhiA[i] = IncrementDeformation[i + 3];
			dPhiB[i] = IncrementDeformation[i + 9];
		}

		//calculating quaternions
		VectorType drA_vec = ZeroVector(dimension);
		VectorType drB_vec = ZeroVector(dimension);
		double drA_sca, drB_sca;

		drA_vec = 0.50 * dPhiA;
		drB_vec = 0.50 * dPhiB;

		drA_sca = 0.00;
		drB_sca = 0.00;
		for (int i = 0; i < dimension; i++) {
			drA_sca += drA_vec[i] * drA_vec[i];
			drB_sca += drB_vec[i] * drB_vec[i];
		}
		drA_sca = 1.00 - drA_sca;
		drB_sca = 1.00 - drB_sca;

		drA_sca = sqrt(drA_sca);
		drB_sca = sqrt(drB_sca);

		//1st solution step
		if (mIterationCount == 0) {
			mQuaternionVEC_A = ZeroVector(dimension);
			mQuaternionVEC_B = ZeroVector(dimension);
			mQuaternionSCA_A = 1.00;
			mQuaternionSCA_B = 1.00;
		}

		VectorType tempVec = ZeroVector(dimension);
		double tempSca = 0.00;

		//Node A
		tempVec = mQuaternionVEC_A;
		tempSca = mQuaternionSCA_A;

		mQuaternionSCA_A = drA_sca *tempSca;
		for (int i = 0; i < dimension; i++) {
			mQuaternionSCA_A -= drA_vec[i] * tempVec[i];
		}
		mQuaternionVEC_A = drA_sca*tempVec;
		mQuaternionVEC_A += tempSca * drA_vec;
		mQuaternionVEC_A += MathUtils<double>::CrossProduct(drA_vec, tempVec);

		//Node B
		tempVec = mQuaternionVEC_B;
		tempSca = mQuaternionSCA_B;

		mQuaternionSCA_B = drB_sca *tempSca;
		for (int i = 0; i < dimension; i++) {
			mQuaternionSCA_B -= drB_vec[i] * tempVec[i];
		}

		mQuaternionVEC_B = drB_sca*tempVec;
		mQuaternionVEC_B += tempSca * drB_vec;
		mQuaternionVEC_B += MathUtils<double>::CrossProduct(drB_vec, tempVec);


		//scalar part of difference quaternion
		double scalar_diff;
		scalar_diff = (mQuaternionSCA_A + mQuaternionSCA_B) *
			(mQuaternionSCA_A + mQuaternionSCA_B);

		tempVec = mQuaternionVEC_A + mQuaternionVEC_B;
		scalar_diff += MathUtils<double>::Norm(tempVec) *
			MathUtils<double>::Norm(tempVec);

		scalar_diff = 0.50 * sqrt(scalar_diff);

		//mean rotation quaternion
		double meanRotationScalar;
		meanRotationScalar = (mQuaternionSCA_A + mQuaternionSCA_B) * 0.50;
		meanRotationScalar = meanRotationScalar / scalar_diff;

		VectorType meanRotationVector = ZeroVector(dimension);
		meanRotationVector = (mQuaternionVEC_A + mQuaternionVEC_B) * 0.50;
		meanRotationVector = meanRotationVector / scalar_diff;

		//vector part of difference quaternion
		VectorType vector_diff = ZeroVector(dimension);
		vector_diff = mQuaternionSCA_A * mQuaternionVEC_B;
		vector_diff -= mQuaternionSCA_B * mQuaternionVEC_A;
		vector_diff += MathUtils<double>::CrossProduct(mQuaternionVEC_A,
			mQuaternionVEC_B);

		vector_diff = 0.50 * vector_diff / scalar_diff;


		//rotate inital element basis
		const double r0 = meanRotationScalar;
		const double r1 = meanRotationVector[0];
		const double r2 = meanRotationVector[1];
		const double r3 = meanRotationVector[2];

		Quaternion<double> q(r0, r1, r2, r3);
		VectorType rotatedNX0 = this->mNX0;
		VectorType rotatedNY0 = this->mNY0;
		VectorType rotatedNZ0 = this->mNZ0;
		q.RotateVector3(rotatedNX0);
		q.RotateVector3(rotatedNY0);
		q.RotateVector3(rotatedNZ0);

		MatrixType RotatedCS = ZeroMatrix(dimension, dimension);
		for (int i = 0; i < dimension; i++) {
			RotatedCS(i, 0) = rotatedNX0[i];
			RotatedCS(i, 1) = rotatedNY0[i];
			RotatedCS(i, 2) = rotatedNZ0[i];
		}
		

		//rotate basis to element axis + redefine R
		VectorType n_bisectrix = ZeroVector(dimension);
		VectorType nx = ZeroVector(dimension);
		for (int i = 0; i < dimension; i++) nx[i] = RotatedCS(i, 0);

		VectorType deltaX = ZeroVector(dimension);
		deltaX[0] = this->mTotalNodalPosistion[3] - this->mTotalNodalPosistion[0];
		deltaX[1] = this->mTotalNodalPosistion[4] - this->mTotalNodalPosistion[1];
		deltaX[2] = this->mTotalNodalPosistion[5] - this->mTotalNodalPosistion[2];

		double VectorNorm;
		n_bisectrix = nx + (deltaX / this->mCurrentLength);
		VectorNorm = MathUtils<double>::Norm(n_bisectrix);
		n_bisectrix /= VectorNorm;

		MatrixType n_xyz = ZeroMatrix(dimension);
		for (int i = 0; i < dimension; i++) {
			n_xyz(i, 0) = -1.0 * RotatedCS(i, 0);
			n_xyz(i, 1) = 1.0 * RotatedCS(i, 1);
			n_xyz(i, 2) = 1.0 * RotatedCS(i, 2);
		}

		MatrixType Identity = ZeroMatrix(dimension);
		for (int i = 0; i < dimension; i++) Identity(i, i) = 1.0;
		Identity -= 2.0 * outer_prod(n_bisectrix, n_bisectrix);
		n_xyz = prod(Identity, n_xyz);

		//calculating deformation modes
		mPhiS = ZeroVector(dimension);
		mPhiA = ZeroVector(dimension);
		mPhiS = 4.00 * prod(Matrix(trans(n_xyz)), vector_diff);

		nx = ZeroVector(dimension);
		tempVec = ZeroVector(dimension);
		for (int i = 0; i < dimension; i++) nx[i] = n_xyz(i, 0);
		tempVec = MathUtils<double>::CrossProduct(nx, n_bisectrix);
		mPhiA = 4.00 * prod(Matrix(trans(n_xyz)), tempVec);

		return n_xyz;
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::GetValuesVector(Vector& rValues, int Step){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int element_size = number_of_nodes * dimension * 2;

		if (rValues.size() != element_size) rValues.resize(element_size, false);

		for (int i = 0; i < number_of_nodes; i++)
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

	void CrBeamElement3D2N::GetFirstDerivativesVector(Vector& rValues, int Step){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int element_size = number_of_nodes * dimension * 2;

		if (rValues.size() != element_size) rValues.resize(element_size, false);

		for (int i = 0; i < number_of_nodes; i++)
		{
			int index = i * dimension * 2;
			rValues[index] = this->GetGeometry()[i].
				FastGetSolutionStepValue(VELOCITY_X, Step);
			rValues[index + 1] = this->GetGeometry()[i].
				FastGetSolutionStepValue(VELOCITY_Y, Step);
			rValues[index + 2] = this->GetGeometry()[i].
				FastGetSolutionStepValue(VELOCITY_Z, Step);

			//rotational dofs negelcted
			rValues[index + 3] = 0.00;
			rValues[index + 4] = 0.00;
			rValues[index + 5] = 0.00;
		}

		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::GetSecondDerivativesVector(Vector& rValues, int Step){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int element_size = number_of_nodes * dimension * 2;

		if (rValues.size() != element_size) rValues.resize(element_size, false);

		for (int i = 0; i < number_of_nodes; i++)
		{
			int index = i * dimension * 2;

			rValues[index] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ACCELERATION_X, Step);
			rValues[index + 1] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ACCELERATION_Y, Step);
			rValues[index + 2] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ACCELERATION_Z, Step);

			//rotational dofs negelcted
			rValues[index + 3] = 0.00;
			rValues[index + 4] = 0.00;
			rValues[index + 5] = 0.00;
		}
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateMassMatrix(MatrixType& rMassMatrix,
								ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int number_of_nodes = GetGeometry().PointsNumber();
		const int dimension = GetGeometry().WorkingSpaceDimension();
		const int MatSize = number_of_nodes * dimension * 2;

		if (rMassMatrix.size1() != MatSize) {
			rMassMatrix.resize(MatSize, MatSize, false);
		}
		rMassMatrix = ZeroMatrix(MatSize, MatSize);

		const double TotalMass = this->mArea * this->mLength * this->mDensity;

		Vector LumpFact = ZeroVector(number_of_nodes);
		LumpFact = GetGeometry().LumpingFactors(LumpFact);

		//translatonal mass	
		for (int i = 0; i < number_of_nodes; i++)
		{
			double temp = LumpFact[i] * TotalMass;

			for (int j = 0; j < dimension; j++)
			{
				int index = i * (dimension * 2) + j;
				rMassMatrix(index, index) = temp;
			}
		}
		//rotaional mass neglected
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateDampingMatrix(MatrixType& rDampingMatrix,
								ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int MatSize = number_of_nodes * dimension * 2;

		if (rDampingMatrix.size1() != MatSize)
		{
			rDampingMatrix.resize(MatSize, MatSize, false);
		}

		noalias(rDampingMatrix) = ZeroMatrix(MatSize, MatSize);

		MatrixType StiffnessMatrix = ZeroMatrix(MatSize, MatSize);

		this->CalculateLeftHandSide(StiffnessMatrix, rCurrentProcessInfo);

		MatrixType MassMatrix = ZeroMatrix(MatSize, MatSize);

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

		noalias(rDampingMatrix) += alpha * MassMatrix;
		noalias(rDampingMatrix) += beta  * StiffnessMatrix;

		KRATOS_CATCH("")
	}


	CrBeamElement3D2N::VectorType CrBeamElement3D2N::CalculateBodyForces(){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int MatSize = number_of_nodes * dimension * 2;

		//getting shapefunctionvalues 
		const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues();

		//creating necessary values 
		const double TotalMass = this->mArea * this->mLength * this->mDensity;
		VectorType BodyForcesNode = ZeroVector(dimension);
		VectorType BodyForcesGlobal = ZeroVector(MatSize);

		//assemble global Vector
		for (int i = 0; i < number_of_nodes; i++) {
			BodyForcesNode = TotalMass*this->GetGeometry()[i].
				FastGetSolutionStepValue(VOLUME_ACCELERATION)*Ncontainer(0, i);

			for (int j = 0; j < dimension; j++) {
				BodyForcesGlobal[(i*dimension * 2) + j] = BodyForcesNode[j];
			}
		}
		return BodyForcesGlobal;
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
								VectorType& rRightHandSideVector,
								ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = NumNodes * dimension;
		const int LocalSize = NumNodes * dimension * 2;

		//update displacement_delta
		this->Initialize();
		this->UpdateIncrementDeformation();

		//calculate Transformation Matrix
		Matrix TransformationMatrix = ZeroMatrix(LocalSize);
		this->CalculateTransformationMatrix(TransformationMatrix);

		//deformation modes
		VectorType elementForces_t = ZeroVector(size);
		elementForces_t = this->CalculateElementForces();

		//Nodal element forces local
		VectorType nodalForcesLocal_qe = ZeroVector(LocalSize);
		Matrix TransformationMatrixS = ZeroMatrix(LocalSize,size);
		TransformationMatrixS = this->CalculateTransformationS();
		nodalForcesLocal_qe = prod(TransformationMatrixS, elementForces_t);

		//resizing the matrices + create memory for LHS
		if (rLeftHandSideMatrix.size1() != LocalSize) {
			rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		}
		rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		//creating LHS
		noalias(rLeftHandSideMatrix) += 
			this->CreateElementStiffnessMatrix_Material();
		noalias(rLeftHandSideMatrix) += 
			this->CreateElementStiffnessMatrix_Geometry(nodalForcesLocal_qe);

		Matrix aux_matrix = ZeroMatrix(LocalSize);
		noalias(aux_matrix) = prod(TransformationMatrix, rLeftHandSideMatrix);
		noalias(rLeftHandSideMatrix) = prod(aux_matrix,
			Matrix(trans(TransformationMatrix)));

		//Nodal element forces global
		VectorType nodalForcesGlobal_q = ZeroVector(LocalSize);
		noalias(nodalForcesGlobal_q) = prod(TransformationMatrix,
			nodalForcesLocal_qe);

		//create+compute RHS
		//update Residual
		if (rRightHandSideVector.size() != LocalSize) {
			rRightHandSideVector = ZeroVector(LocalSize);
		}
		rRightHandSideVector = ZeroVector(LocalSize);
		noalias(rRightHandSideVector) -= nodalForcesGlobal_q;


		//assign global element variables
		this->mLHS = rLeftHandSideMatrix;

		//add bodyforces 
		this->mBodyForces = ZeroVector(LocalSize);
		this->mBodyForces = this->CalculateBodyForces();
		noalias(rRightHandSideVector) += this->mBodyForces;

		this->mIterationCount++;
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateRightHandSide(
								VectorType& rRightHandSideVector,
								ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = NumNodes * dimension;
		const int LocalSize = NumNodes * dimension * 2;

		if (rRightHandSideVector.size() != LocalSize) {
			rRightHandSideVector = ZeroVector(LocalSize);
		}
		rRightHandSideVector = ZeroVector(LocalSize);
		
		this->UpdateIncrementDeformation();
		Matrix TransformationMatrix = ZeroMatrix(LocalSize);
		this->CalculateTransformationMatrix(TransformationMatrix);
		VectorType elementForces_t = ZeroVector(size);
		elementForces_t = this->CalculateElementForces();
		VectorType nodalForcesLocal_qe = ZeroVector(LocalSize);
		Matrix TransformationMatrixS = ZeroMatrix(LocalSize, size);
		TransformationMatrixS = this->CalculateTransformationS();
		noalias(nodalForcesLocal_qe) = prod(TransformationMatrixS,
			elementForces_t);
		VectorType nodalForcesGlobal_q = ZeroVector(LocalSize);
		nodalForcesGlobal_q = prod(TransformationMatrix, nodalForcesLocal_qe);
		noalias(rRightHandSideVector) -= nodalForcesGlobal_q;

		//add bodyforces 
		noalias(rRightHandSideVector) += mBodyForces;
		KRATOS_CATCH("")

	}

	void CrBeamElement3D2N::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
								ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int LocalSize = NumNodes * dimension * 2;

		if (rLeftHandSideMatrix.size1() != LocalSize) {
			rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		}
		rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		noalias(rLeftHandSideMatrix) = this->mLHS;
		KRATOS_CATCH("")
	}

	CrBeamElement3D2N::VectorType CrBeamElement3D2N::CalculateElementForces(){

		KRATOS_TRY
		Vector deformation_modes_total_V = ZeroVector(6);
		deformation_modes_total_V[3] = this->mCurrentLength - this->mLength;
		for (int i = 0; i < 3; i++) deformation_modes_total_V[i] = mPhiS[i];
		for (int i = 0; i < 2; i++) deformation_modes_total_V[i + 4] = mPhiA[i + 1];
		//calculate element forces
		VectorType element_forces_t = ZeroVector(6);
		MatrixType material_stiffness_Kd = ZeroMatrix(6);

		material_stiffness_Kd = this->CalculateMaterialStiffness();
		noalias(element_forces_t) = prod(material_stiffness_Kd,
			deformation_modes_total_V);

		return element_forces_t;
		KRATOS_CATCH("")
	}

	double CrBeamElement3D2N::CalculateGreenLagrangeStrain(){

		KRATOS_TRY
		const double l = this->mCurrentLength;
		const double L = this->mLength;
		const double e = ((l * l - L * L) / (2 * L * L));
		return e;
		KRATOS_CATCH("")
	}

	double CrBeamElement3D2N::CalculateCurrentLength(){

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
	double CrBeamElement3D2N::CalculatePsi(const double I, const double A_eff){

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

	double CrBeamElement3D2N::CalculateReferenceLength(){

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

		VectorType actualDeformation = ZeroVector(LocalSize);
		VectorType total_nodal_def = ZeroVector(LocalSize);
		VectorType total_nodal_pos = ZeroVector(size);

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
								const ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int size = NumNodes * dimension;
		const int LocalSize = NumNodes * dimension * 2;

		const int&  write_points_number = GetGeometry()
			.IntegrationPointsNumber(GeometryData::GI_GAUSS_3);
		if (rOutput.size() != write_points_number) {
			rOutput.resize(write_points_number);
		}

		this->UpdateIncrementDeformation();
		//calculate Transformation Matrix
		Matrix TransformationMatrix = ZeroMatrix(LocalSize);
		this->CalculateTransformationMatrix(TransformationMatrix);
		//deformation modes
		VectorType elementForces_t = ZeroVector(size);
		elementForces_t = this->CalculateElementForces();
		VectorType Stress = ZeroVector(LocalSize);
		Matrix TransformationMatrixS = ZeroMatrix(LocalSize, size);
		TransformationMatrixS = this->CalculateTransformationS();
		Stress = prod(TransformationMatrixS, elementForces_t);

		////calculate Body Forces here ngelected atm
		//Vector BodyForces = ZeroVector(12);
		//this->CalculateLocalBodyForce(LocalForceVector, rVolumeForce);
		//Stress -= BodyForces;

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
		KRATOS_TRY
		this->CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
		KRATOS_CATCH("")
	}

	int  CrBeamElement3D2N::Check(const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

			if (GetGeometry().WorkingSpaceDimension() != 3 || GetGeometry().size() != 2)
			{
				KRATOS_THROW_ERROR(std::invalid_argument, "The truss element works only in 3D and with 2 noded elements", "")
			}
		//verify that the variables are correctly initialized
		if (VELOCITY.Key() == 0)
			KRATOS_THROW_ERROR(std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "")
			if (DISPLACEMENT.Key() == 0)
				KRATOS_THROW_ERROR(std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "")
				if (ACCELERATION.Key() == 0)
					KRATOS_THROW_ERROR(std::invalid_argument, "ACCELERATION has Key zero! (check if the application is correctly registered", "")
					if (DENSITY.Key() == 0)
						KRATOS_THROW_ERROR(std::invalid_argument, "DENSITY has Key zero! (check if the application is correctly registered", "")
						if (CROSS_AREA.Key() == 0)
							KRATOS_THROW_ERROR(std::invalid_argument, "CROSS_AREA has Key zero! (check if the application is correctly registered", "")
							//verify that the dofs exist
							for (unsigned int i = 0; i<this->GetGeometry().size(); i++)
							{
								if (this->GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false)
									KRATOS_THROW_ERROR(std::invalid_argument, "missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id())
									if (this->GetGeometry()[i].HasDofFor(DISPLACEMENT_X) == false || this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false || this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Z) == false)
										KRATOS_THROW_ERROR(std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id())
							}



		if (this->GetProperties().Has(CROSS_AREA) == false || this->GetProperties()[CROSS_AREA] == 0)
		{
			KRATOS_THROW_ERROR(std::logic_error, "CROSS_AREA not provided for this element", this->Id())
		}

		if (this->GetProperties().Has(YOUNG_MODULUS) == false || this->GetProperties()[YOUNG_MODULUS] == 0)
		{
			KRATOS_THROW_ERROR(std::logic_error, "YOUNG_MODULUS not provided for this element", this->Id())
		}
		if (this->GetProperties().Has(DENSITY) == false)
		{
			KRATOS_THROW_ERROR(std::logic_error, "DENSITY not provided for this element", this->Id())
		}


		return 0;

		KRATOS_CATCH("")
	}

} // namespace Kratos.


