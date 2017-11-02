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
// System includes

// External includes

// Project includes
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
		if (rResult.size() != msElementSize) rResult.resize(msElementSize);

		for (int i = 0; i < msNumberOfNodes; ++i)
		{
			int index = i * msNumberOfNodes * msDimension;
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

		if (rElementalDofList.size() != msElementSize) {
			rElementalDofList.resize(msElementSize);
		}

		for (int i = 0; i < msNumberOfNodes; ++i)
		{
			int index = i * msNumberOfNodes * msDimension;
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
		if (this->mIterationCount == 0)
		{
			this->mNodalForces = ZeroVector(msElementSize);
			this->CalculateInitialLocalCS();
		}
		KRATOS_CATCH("")
	}

	bounded_matrix<double,CrBeamElement3D2N::msElementSize,CrBeamElement3D2N::msElementSize> 
	CrBeamElement3D2N::CreateElementStiffnessMatrix_Material() {

		KRATOS_TRY;
		const double E = this->GetProperties()[YOUNG_MODULUS];
		const double G = this->CalculateShearModulus();
		const double A = this->GetProperties()[CROSS_AREA];
		const double L = this->CalculateReferenceLength();

		bounded_vector<double,msDimension> inertia = this->GetProperties()[LOCAL_INERTIA_VECTOR];
		const double J = inertia[0];
		const double Iy = inertia[1];
		const double Iz = inertia[2];

		double Ay = 0.00;
		if (this->GetProperties().Has(AREA_EFFECTIVE_Y)) {
			Ay = GetProperties()[AREA_EFFECTIVE_Y];
		}

		double Az = 0.00;
		if (this->GetProperties().Has(AREA_EFFECTIVE_Z)) {
			Az = GetProperties()[AREA_EFFECTIVE_Z];
		}
		const double Psi_y = this->CalculatePsi(Iy, Az);
		const double Psi_z = this->CalculatePsi(Iz, Ay);



		bounded_matrix<double,msElementSize,msElementSize> 
		LocalStiffnessMatrix = ZeroMatrix(msElementSize, msElementSize);
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

	bounded_matrix<double,CrBeamElement3D2N::msElementSize,CrBeamElement3D2N::msElementSize> 
 	CrBeamElement3D2N::CreateElementStiffnessMatrix_Geometry() {

		KRATOS_TRY;
		//deformation modes
		Vector elementForces_t = ZeroVector(msLocalSize);
		elementForces_t = this->CalculateElementForces();

		//Nodal element forces local
		Vector nodalForcesLocal_qe = ZeroVector(msElementSize);
		Matrix TransformationMatrixS = ZeroMatrix(msElementSize, msLocalSize);
		TransformationMatrixS = this->CalculateTransformationS();
		nodalForcesLocal_qe = prod(TransformationMatrixS, elementForces_t);

		//save local nodal forces
		this->mNodalForces = ZeroVector(msElementSize);
		this->mNodalForces = nodalForcesLocal_qe;



		const double N = nodalForcesLocal_qe[6];
		const double Mt = nodalForcesLocal_qe[9];
		const double my_A = nodalForcesLocal_qe[4];
		const double mz_A = nodalForcesLocal_qe[5];
		const double my_B = nodalForcesLocal_qe[10];
		const double mz_B = nodalForcesLocal_qe[11];

		const double L = this->CalculateCurrentLength();
		const double Qy = -1.00 * (mz_A + mz_B) / L;
		const double Qz = (my_A + my_B) / L;

		bounded_matrix<double,msElementSize,msElementSize> 
		 LocalStiffnessMatrix = ZeroMatrix(msElementSize, msElementSize);

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

	bounded_matrix<double,CrBeamElement3D2N::msLocalSize,CrBeamElement3D2N::msLocalSize>  
	CrBeamElement3D2N::CalculateDeformationStiffness() {

		KRATOS_TRY
		bounded_matrix<double,msLocalSize,msLocalSize>
		 Kd = ZeroMatrix(msLocalSize, msLocalSize);
		const double E = this->GetProperties()[YOUNG_MODULUS];
		const double G = this->CalculateShearModulus();
		const double A = this->GetProperties()[CROSS_AREA];
		const double L = this->CalculateReferenceLength();
		
		bounded_vector<double,msDimension> inertia = this->GetProperties()[LOCAL_INERTIA_VECTOR];
		const double J = inertia[0];
		const double Iy = inertia[1];
		const double Iz = inertia[2];

		double Ay = 0.00;
		if (this->GetProperties().Has(AREA_EFFECTIVE_Y)) {
			Ay = GetProperties()[AREA_EFFECTIVE_Y];
		}

		double Az = 0.00;
		if (this->GetProperties().Has(AREA_EFFECTIVE_Z)) {
			Az = GetProperties()[AREA_EFFECTIVE_Z];
		}
		const double Psi_y = this->CalculatePsi(Iy, Az);
		const double Psi_z = this->CalculatePsi(Iz, Ay);

		Kd(0, 0) = G * J / L;
		Kd(1, 1) = E * Iy / L;
		Kd(2, 2) = E * Iz / L;
		Kd(3, 3) = E * A / L;
		Kd(4, 4) = 3.0 * E * Iy * Psi_y / L;
		Kd(5, 5) = 3.0 * E * Iz * Psi_z / L;


		//add geometric stiffness part
		if (this->mIsLinearElement == false)
		{
			const double l = this->CalculateCurrentLength();
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
		const double numerical_limit = std::numeric_limits<double>::epsilon();
		array_1d<double, msDimension> DirectionVectorX = ZeroVector(msDimension);
		array_1d<double, msDimension> DirectionVectorY = ZeroVector(msDimension);
		array_1d<double, msDimension> DirectionVectorZ = ZeroVector(msDimension);
		array_1d<double, msLocalSize> ReferenceCoordinates = ZeroVector(msLocalSize);

		ReferenceCoordinates[0] = this->GetGeometry()[0].X0();
		ReferenceCoordinates[1] = this->GetGeometry()[0].Y0();
		ReferenceCoordinates[2] = this->GetGeometry()[0].Z0();
		ReferenceCoordinates[3] = this->GetGeometry()[1].X0();
		ReferenceCoordinates[4] = this->GetGeometry()[1].Y0();
		ReferenceCoordinates[5] = this->GetGeometry()[1].Z0();

		for (unsigned int i = 0; i < msDimension; ++i)
		{
			DirectionVectorX[i] = (ReferenceCoordinates[i + msDimension]
				- ReferenceCoordinates[i]);
		}

		Matrix Temp = ZeroMatrix(msDimension);
		this->mRotationMatrix0 = ZeroMatrix(msElementSize);

		// take user defined local axis 2 from GID input
		if (this->Has(LOCAL_AXIS_2)) 
		{
			double VectorNorm = MathUtils<double>::Norm(DirectionVectorX);
			if (VectorNorm > numerical_limit) DirectionVectorX /= VectorNorm;

			DirectionVectorY = this->GetValue(LOCAL_AXIS_2);
			DirectionVectorZ[0] = DirectionVectorX[1]*DirectionVectorY[2]-DirectionVectorX[2]*DirectionVectorY[1];
			DirectionVectorZ[1] = DirectionVectorX[2]*DirectionVectorY[0]-DirectionVectorX[0]*DirectionVectorY[2];
			DirectionVectorZ[2] = DirectionVectorX[0]*DirectionVectorY[1]-DirectionVectorX[1]*DirectionVectorY[0];

			VectorNorm = MathUtils<double>::Norm(DirectionVectorZ);
			if (VectorNorm > numerical_limit) DirectionVectorZ /= VectorNorm;
			else KRATOS_ERROR << "LOCAL_AXIS_3 has length 0 for element " << this->Id() << std::endl;

			for (int i = 0; i < msDimension; ++i)
			{
				Temp(i, 0) = DirectionVectorX[i];
				Temp(i, 1) = DirectionVectorY[i];
				Temp(i, 2) = DirectionVectorZ[i];
			}
		}

		// if no user defined local axis 2 input available use this
		else
		{
			//use orientation class 1st constructor
			double theta_costum = 0.00;
			if (this->GetProperties().Has(ANG_ROT)) theta_costum = this->GetProperties()[ANG_ROT];
			
			Orientation element_axis(DirectionVectorX, theta_costum);
			element_axis.CalculateBasisVectors(DirectionVectorX, DirectionVectorY,
				DirectionVectorZ);
			element_axis.CalculateRotationMatrix(Temp);
		}

		//save them to update the local axis in every following iter. step
		this->mNX0 = DirectionVectorX;
		this->mNY0 = DirectionVectorY;
		this->mNZ0 = DirectionVectorZ;

		//Create big rotation Matrix
		this->AssembleSmallInBigMatrix(Temp, this->mRotationMatrix0);

		//provide Initial Rotation Matrix for strategies that dont call 'CalculateLocalSystem'
		this->mRotationMatrix = ZeroMatrix(msElementSize);
		this->mRotationMatrix = this->mRotationMatrix0;

		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateTransformationMatrix(bounded_matrix<double,
		CrBeamElement3D2N::msElementSize,CrBeamElement3D2N::msElementSize>& rRotationMatrix) {

		KRATOS_TRY
		//initialize local CS
		if (this->mIterationCount == 0) this->CalculateInitialLocalCS();

		//update local CS
		Matrix AuxRotationMatrix = ZeroMatrix(msDimension);
		AuxRotationMatrix = this->UpdateRotationMatrixLocal();

		rRotationMatrix = ZeroMatrix(msElementSize);
		//Building the rotation matrix for the local element matrix
		this->AssembleSmallInBigMatrix(AuxRotationMatrix, rRotationMatrix);
		KRATOS_CATCH("")
	}

	bounded_matrix<double,CrBeamElement3D2N::msElementSize,CrBeamElement3D2N::msLocalSize>
	CrBeamElement3D2N::CalculateTransformationS() {

		KRATOS_TRY
		const double L = this->CalculateCurrentLength();
		bounded_matrix<double,msElementSize,msLocalSize> S = ZeroMatrix(msElementSize, msLocalSize);
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

	bounded_matrix<double,CrBeamElement3D2N::msDimension,CrBeamElement3D2N::msDimension>
	CrBeamElement3D2N::UpdateRotationMatrixLocal() {

		KRATOS_TRY
		const double numerical_limit = std::numeric_limits<double>::epsilon();
		bounded_vector<double,msDimension> dPhiA = ZeroVector(msDimension);
		bounded_vector<double,msDimension> dPhiB = ZeroVector(msDimension);
		bounded_vector<double,msElementSize>  IncrementDeformation = ZeroVector(msElementSize);
		IncrementDeformation = this->mIncrementDeformation;

		for (unsigned int i = 0; i < msDimension; ++i) {
			dPhiA[i] = IncrementDeformation[i + 3];
			dPhiB[i] = IncrementDeformation[i + 9];
		}

		//calculating quaternions
		Vector drA_vec = ZeroVector(msDimension);
		Vector drB_vec = ZeroVector(msDimension);
		double drA_sca, drB_sca;

		drA_vec = 0.50 * dPhiA;
		drB_vec = 0.50 * dPhiB;

		drA_sca = 0.00;
		drB_sca = 0.00;
		for (unsigned int i = 0; i < msDimension; ++i) {
			drA_sca += drA_vec[i] * drA_vec[i];
			drB_sca += drB_vec[i] * drB_vec[i];
		}
		drA_sca = 1.00 - drA_sca;
		drB_sca = 1.00 - drB_sca;

		drA_sca = sqrt(drA_sca);
		drB_sca = sqrt(drB_sca);


		//1st solution step
		if (mIterationCount == 0) {
			this->mQuaternionVEC_A = ZeroVector(msDimension);
			this->mQuaternionVEC_B = ZeroVector(msDimension);
			this->mQuaternionSCA_A = 1.00;
			this->mQuaternionSCA_B = 1.00;
		}

		Vector tempVec = ZeroVector(msDimension);
		double tempSca = 0.00;

		//Node A
		tempVec = this->mQuaternionVEC_A;
		tempSca = this->mQuaternionSCA_A;

		this->mQuaternionSCA_A = drA_sca *tempSca;
		for (unsigned int i = 0; i < msDimension; ++i) {
			this->mQuaternionSCA_A -= drA_vec[i] * tempVec[i];
		}
		this->mQuaternionVEC_A = drA_sca*tempVec;
		this->mQuaternionVEC_A += tempSca * drA_vec;
		this->mQuaternionVEC_A += MathUtils<double>::CrossProduct(drA_vec, tempVec);

		//Node B
		tempVec = this->mQuaternionVEC_B;
		tempSca = this->mQuaternionSCA_B;

		this->mQuaternionSCA_B = drB_sca *tempSca;
		for (unsigned int i = 0; i < msDimension; ++i) {
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

		bounded_vector<double,msDimension>  meanRotationVector = ZeroVector(msDimension);
		meanRotationVector = (this->mQuaternionVEC_A + this->mQuaternionVEC_B) * 0.50;
		meanRotationVector = meanRotationVector / scalar_diff;

		//vector part of difference quaternion
		bounded_vector<double,msDimension>  vector_diff = ZeroVector(msDimension);
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

		bounded_matrix<double,msDimension,msDimension> RotatedCS = ZeroMatrix(msDimension, msDimension);
		for (unsigned int i = 0; i < msDimension; ++i) {
			RotatedCS(i, 0) = rotatedNX0[i];
			RotatedCS(i, 1) = rotatedNY0[i];
			RotatedCS(i, 2) = rotatedNZ0[i];
		}

		//rotate basis to element axis + redefine R
		Vector n_bisectrix = ZeroVector(msDimension);
		Vector deltaX = ZeroVector(msDimension);
		double VectorNorm;

		deltaX[0] = this->mTotalNodalPosistion[3] - this->mTotalNodalPosistion[0];
		deltaX[1] = this->mTotalNodalPosistion[4] - this->mTotalNodalPosistion[1];
		deltaX[2] = this->mTotalNodalPosistion[5] - this->mTotalNodalPosistion[2];


		VectorNorm = MathUtils<double>::Norm(deltaX);
		if (VectorNorm > numerical_limit) deltaX /= VectorNorm;


		n_bisectrix = rotatedNX0 + deltaX;
		VectorNorm = MathUtils<double>::Norm(n_bisectrix);
		if (VectorNorm > numerical_limit) n_bisectrix /= VectorNorm;

		bounded_matrix<double,msDimension,msDimension> n_xyz = ZeroMatrix(msDimension);
		for (unsigned int i = 0; i < msDimension; ++i) {
			n_xyz(i, 0) = -1.0 * RotatedCS(i, 0);
			n_xyz(i, 1) = 1.0 * RotatedCS(i, 1);
			n_xyz(i, 2) = 1.0 * RotatedCS(i, 2);
		}

		bounded_matrix<double,msDimension,msDimension> Identity = ZeroMatrix(msDimension);
		for (unsigned int i = 0; i < msDimension; ++i) Identity(i, i) = 1.0;
		Identity -= 2.0 * outer_prod(n_bisectrix, n_bisectrix);
		n_xyz = prod(Identity, n_xyz);


		//save current CS for GID OUTPUT
		this->mNX = ZeroVector(msDimension);
		this->mNY = ZeroVector(msDimension);
		this->mNZ = ZeroVector(msDimension);
		for (unsigned int i = 0; i < msDimension; ++i)
		{
			this->mNX[i] = n_xyz(i, 0);
			this->mNY[i] = n_xyz(i, 1);
			this->mNZ[i] = n_xyz(i, 2);
		}

		//calculating deformation modes
		this->mPhiS = ZeroVector(msDimension);
		this->mPhiA = ZeroVector(msDimension);
		this->mPhiS = prod(Matrix(trans(n_xyz)), vector_diff);
		this->mPhiS *= 4.00;

		rotatedNX0 = ZeroVector(msDimension);
		tempVec = ZeroVector(msDimension);
		for (unsigned int i = 0; i < msDimension; ++i) rotatedNX0[i] = n_xyz(i, 0);
		tempVec = MathUtils<double>::CrossProduct(rotatedNX0, n_bisectrix);
		this->mPhiA = prod(Matrix(trans(n_xyz)), tempVec);
		this->mPhiA *= 4.00;

		if (this->mIterationCount == 0)
		{
			this->mPhiS = ZeroVector(msDimension);
			this->mPhiA = ZeroVector(msDimension);
		}
		return n_xyz;
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::GetValuesVector(Vector& rValues, int Step) {

		KRATOS_TRY
		if (rValues.size() != msElementSize) rValues.resize(msElementSize, false);

		for (int i = 0; i < msNumberOfNodes; ++i)
		{
			int index = i * msDimension * 2;
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
		if (rValues.size() != msElementSize) rValues.resize(msElementSize, false);

		for (int i = 0; i < msNumberOfNodes; ++i)
		{
			int index = i * msDimension * 2;
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
		if (rValues.size() != msElementSize) rValues.resize(msElementSize, false);

		for (int i = 0; i < msNumberOfNodes; ++i)
		{
			int index = i * msDimension * 2;

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
		if (rMassMatrix.size1() != msElementSize) {
			rMassMatrix.resize(msElementSize, msElementSize, false);
		}
		rMassMatrix = ZeroMatrix(msElementSize, msElementSize);



		if (this->GetProperties().Has(LUMPED_MASS_MATRIX)) {
			this->mIsLumpedMassMatrix = GetProperties()[LUMPED_MASS_MATRIX];
		}
		else this->mIsLumpedMassMatrix = false;



		if (this->mIsLumpedMassMatrix)
		{
			this->CalculateLumpedMassMatrix(rMassMatrix, rCurrentProcessInfo);
		}
		else
		{
			this->CalculateConsistentMassMatrix(rMassMatrix, rCurrentProcessInfo);

			bounded_matrix<double,msElementSize,msElementSize> RotationMatrix = ZeroMatrix(msElementSize,msElementSize);
			bounded_matrix<double,msElementSize,msElementSize> aux_matrix = ZeroMatrix(msElementSize,msElementSize);

			if (this->mIsLinearElement)
			{
				RotationMatrix = this->mRotationMatrix0;
			}
			else
			{
				RotationMatrix = this->mRotationMatrix;
			}
			aux_matrix = prod(RotationMatrix, rMassMatrix);
			rMassMatrix = prod(aux_matrix,
				Matrix(trans(RotationMatrix)));
		}
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateDampingMatrix(MatrixType& rDampingMatrix,
		ProcessInfo& rCurrentProcessInfo) {

		KRATOS_TRY
		if (rDampingMatrix.size1() != msElementSize)
		{
			rDampingMatrix.resize(msElementSize, msElementSize, false);
		}

		rDampingMatrix = ZeroMatrix(msElementSize, msElementSize);

		Matrix StiffnessMatrix = ZeroMatrix(msElementSize, msElementSize);

		this->CalculateLeftHandSide(StiffnessMatrix, rCurrentProcessInfo);

		Matrix MassMatrix = ZeroMatrix(msElementSize, msElementSize);

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


	bounded_vector<double,CrBeamElement3D2N::msElementSize> CrBeamElement3D2N::CalculateBodyForces()
	{
		KRATOS_TRY
		//getting shapefunctionvalues for linear SF
		const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues(
			GeometryData::GI_GAUSS_1);

		bounded_vector<double,msDimension> EquivalentLineLoad = ZeroVector(msDimension);
		bounded_vector<double,msElementSize> BodyForcesGlobal = ZeroVector(msElementSize);

		const double A = this->GetProperties()[CROSS_AREA];
		const double l = this->CalculateCurrentLength();
		const double rho = this->GetProperties()[DENSITY];

		//calculating equivalent line load
		for (int i = 0; i < msNumberOfNodes; ++i)
		{
			EquivalentLineLoad += A * rho*
				this->GetGeometry()[i].
				FastGetSolutionStepValue(VOLUME_ACCELERATION)*Ncontainer(0, i);
		}


		// adding the nodal forces
		for (int i = 0; i < msNumberOfNodes; ++i)
		{
			int index = i*msLocalSize;
			for (int j = 0; j < msDimension; ++j)
			{
				BodyForcesGlobal[j + index] =
					EquivalentLineLoad[j] * Ncontainer(0, i) * l;
			}
		}

		// adding the nodal moments
		this->CalculateAndAddWorkEquivalentNodalForcesLineLoad
			(EquivalentLineLoad, BodyForcesGlobal, l);


		// return the total ForceVector
		return BodyForcesGlobal;
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateAndAddWorkEquivalentNodalForcesLineLoad(
		const bounded_vector<double,CrBeamElement3D2N::msDimension> ForceInput,
		bounded_vector<double,CrBeamElement3D2N::msElementSize>& rRightHandSideVector,
		const double GeometryLength)
	{
		KRATOS_TRY;
		//calculate orthogonal load vector
		const double numerical_limit = std::numeric_limits<double>::epsilon();
		Vector GeometricOrientation = ZeroVector(msDimension);
		GeometricOrientation[0] = this->GetGeometry()[1].X()
			- this->GetGeometry()[0].X();
		GeometricOrientation[1] = this->GetGeometry()[1].Y()
			- this->GetGeometry()[0].Y();
		if (msDimension == 3)
		{
			GeometricOrientation[2] = this->GetGeometry()[1].Z()
				- this->GetGeometry()[0].Z();
		}

		const double VectorNormA = MathUtils<double>::Norm(GeometricOrientation);
		if (VectorNormA > numerical_limit) GeometricOrientation /= VectorNormA;

		Vector LineLoadDir = ZeroVector(msDimension);
		for (int i = 0; i < msDimension; ++i)
		{
			LineLoadDir[i] = ForceInput[i];
		}

		const double VectorNormB = MathUtils<double>::Norm(LineLoadDir);
		if (VectorNormB > numerical_limit) LineLoadDir /= VectorNormB;

		double cosAngle = 0.00;
		for (int i = 0; i < msDimension; ++i)
		{
			cosAngle += LineLoadDir[i] * GeometricOrientation[i];
		}

		const double sinAngle = sqrt(1.00 - (cosAngle*cosAngle));
		const double NormForceVectorOrth = sinAngle * VectorNormB;


		Vector NodeA = ZeroVector(msDimension);
		NodeA[0] = this->GetGeometry()[0].X();
		NodeA[1] = this->GetGeometry()[0].Y();
		if (msDimension == 3)	NodeA[2] = this->GetGeometry()[0].Z();

		Vector NodeB = ZeroVector(msDimension);
		NodeB = NodeA + LineLoadDir;

		Vector NodeC = ZeroVector(msDimension);
		NodeC = NodeA + (GeometricOrientation*cosAngle);

		Vector LoadOrthogonalDir = ZeroVector(msDimension);
		LoadOrthogonalDir = NodeB - NodeC;
		const double VectorNormC = MathUtils<double>::Norm(LoadOrthogonalDir);
		if (VectorNormC > numerical_limit) LoadOrthogonalDir /= VectorNormC;



		// now caluclate respective work equivilent nodal moments

		const double CustomMoment = NormForceVectorOrth *
			GeometryLength*GeometryLength / 12.00;

		Vector MomentNodeA = ZeroVector(msDimension);
		MomentNodeA = MathUtils<double>::CrossProduct(GeometricOrientation,
			LoadOrthogonalDir);
		MomentNodeA *= CustomMoment;

		for (int i = 0; i < msDimension; ++i)
		{
			rRightHandSideVector[(1 * msDimension) + i] += MomentNodeA[i];
			rRightHandSideVector[(3 * msDimension) + i] -= MomentNodeA[i];
		}

		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo) {

		KRATOS_TRY
		this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

		//Nodal element forces global
		bounded_vector<double,msElementSize> nodalForcesGlobal_q = ZeroVector(msElementSize);
		nodalForcesGlobal_q = prod(this->mRotationMatrix,
			this->mNodalForces);

		//create+compute RHS
		//update Residual
		rRightHandSideVector = ZeroVector(msElementSize);
		rRightHandSideVector -= nodalForcesGlobal_q;


		//LINEAR BEAM ELEMENT
		if (this->mIsLinearElement)
		{
			Vector NodalDeformation = ZeroVector(msElementSize);
			this->GetValuesVector(NodalDeformation);
			rRightHandSideVector = ZeroVector(msElementSize);
			rRightHandSideVector -= prod(rLeftHandSideMatrix, NodalDeformation);
		}
		//add bodyforces 
		rRightHandSideVector += this->CalculateBodyForces();
		this->mIterationCount++;
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateRightHandSide(
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		rRightHandSideVector = ZeroVector(msElementSize);

		if (this->mIsLinearElement == false)
		{
			this->UpdateIncrementDeformation();
			bounded_matrix<double,msElementSize,msElementSize> TransformationMatrix = ZeroMatrix(msElementSize);
			this->CalculateTransformationMatrix(TransformationMatrix);
			bounded_vector<double,msLocalSize> elementForces_t = ZeroVector(msLocalSize);
			elementForces_t = this->CalculateElementForces();
			bounded_vector<double,msElementSize> nodalForcesLocal_qe = ZeroVector(msElementSize);
			bounded_matrix<double,msElementSize,msLocalSize>  TransformationMatrixS = ZeroMatrix(msElementSize, msLocalSize);
			TransformationMatrixS = this->CalculateTransformationS();
			nodalForcesLocal_qe = prod(TransformationMatrixS,
				elementForces_t);
			//save local nodal forces
			this->mNodalForces = ZeroVector(msElementSize);
			this->mNodalForces = nodalForcesLocal_qe;

			bounded_vector<double,msElementSize> nodalForcesGlobal_q = ZeroVector(msElementSize);
			nodalForcesGlobal_q = prod(TransformationMatrix, nodalForcesLocal_qe);
			rRightHandSideVector -= nodalForcesGlobal_q;
		}

		//LINEAR BEAM ELEMENT
		if (this->mIsLinearElement)
		{
			Matrix LeftHandSideMatrix = ZeroMatrix(msElementSize, msElementSize);
			this->CalculateLeftHandSide(LeftHandSideMatrix, rCurrentProcessInfo);
			Vector NodalDeformation = ZeroVector(msElementSize);
			this->GetValuesVector(NodalDeformation);
			rRightHandSideVector = ZeroVector(msElementSize);
			rRightHandSideVector -= prod(LeftHandSideMatrix, NodalDeformation);
		}

		//add bodyforces 
		rRightHandSideVector += this->CalculateBodyForces();
		KRATOS_CATCH("")

	}

	void CrBeamElement3D2N::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
		ProcessInfo& rCurrentProcessInfo) {

		KRATOS_TRY;
		//update displacement_delta
		this->UpdateIncrementDeformation();

		//calculate Transformation Matrix
		bounded_matrix<double,msElementSize,msElementSize> TransformationMatrix = ZeroMatrix(msElementSize);
		this->CalculateTransformationMatrix(TransformationMatrix);
		this->mRotationMatrix = ZeroMatrix(msElementSize);
		this->mRotationMatrix = TransformationMatrix;



		//resizing the matrices + create memory for LHS
		rLeftHandSideMatrix = ZeroMatrix(msElementSize, msElementSize);
		//creating LHS
		rLeftHandSideMatrix +=
			this->CreateElementStiffnessMatrix_Material();
		rLeftHandSideMatrix +=
			this->CreateElementStiffnessMatrix_Geometry();


		bounded_matrix<double,msElementSize,msElementSize> aux_matrix = ZeroMatrix(msElementSize);
		aux_matrix = prod(TransformationMatrix, rLeftHandSideMatrix);
		rLeftHandSideMatrix = prod(aux_matrix,
			Matrix(trans(TransformationMatrix)));

		//LINEAR BEAM ELEMENT
		if (this->mIsLinearElement)
		{
			TransformationMatrix = this->mRotationMatrix0;
			rLeftHandSideMatrix = ZeroMatrix(msElementSize, msElementSize);
			rLeftHandSideMatrix +=
				this->CreateElementStiffnessMatrix_Material();
			bounded_matrix<double,msElementSize,msElementSize> aux_matrix = ZeroMatrix(msElementSize);
			aux_matrix = prod(TransformationMatrix, rLeftHandSideMatrix);
			rLeftHandSideMatrix = prod(aux_matrix,
				Matrix(trans(TransformationMatrix)));
		}
		//assign global element variables
		this->mLHS = rLeftHandSideMatrix;
		KRATOS_CATCH("")
	}

	bounded_vector<double,CrBeamElement3D2N::msLocalSize> CrBeamElement3D2N::CalculateElementForces() {

		KRATOS_TRY;
		bounded_vector<double,msLocalSize> deformation_modes_total_V = ZeroVector(msLocalSize);
		const double L = this->CalculateReferenceLength();
		const double l = this->CalculateCurrentLength();

		deformation_modes_total_V[3] = l - L;
		for (int i = 0; i < 3; ++i) deformation_modes_total_V[i] = this->mPhiS[i];
		for (int i = 0; i < 2; ++i) deformation_modes_total_V[i + 4] = this->mPhiA[i + 1];
		//calculate element forces
		bounded_vector<double,msLocalSize> element_forces_t = ZeroVector(msLocalSize);
		bounded_matrix<double,msLocalSize,msLocalSize> deformation_stiffness_Kd = ZeroMatrix(msLocalSize);

		deformation_stiffness_Kd = this->CalculateDeformationStiffness();
		element_forces_t = prod(deformation_stiffness_Kd,
			deformation_modes_total_V);

		return element_forces_t;
		KRATOS_CATCH("")
	}

	double CrBeamElement3D2N::CalculateCurrentLength() {

		KRATOS_TRY;
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

		KRATOS_TRY;
		const double E = this->GetProperties()[YOUNG_MODULUS];
		const double L = this->CalculateCurrentLength();
		const double G = this->CalculateShearModulus();

		const double phi = (12.0 * E * I) / (L*L * G*A_eff);
		double psi;
		//interpret input A_eff == 0 as shearstiff -> psi = 1.0
		if (A_eff == 0.00) psi = 1.00;
		else psi = 1.0 / (1.0 + phi);

		return psi;
		KRATOS_CATCH("")
	}

	double CrBeamElement3D2N::CalculateReferenceLength() {

		KRATOS_TRY;
		const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
		const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
		const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
		const double L = sqrt(dx*dx + dy*dy + dz*dz);
		return L;
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::UpdateIncrementDeformation() {

		KRATOS_TRY
		Vector actualDeformation = ZeroVector(msElementSize);
		this->mIncrementDeformation = ZeroVector(msElementSize);

		if (mIterationCount == 0) this->mTotalNodalDeformation = ZeroVector(msElementSize);
		this->GetValuesVector(actualDeformation, 0);

		this->mIncrementDeformation = actualDeformation
			- this->mTotalNodalDeformation;

		this->mTotalNodalDeformation = ZeroVector(msElementSize);
		this->mTotalNodalDeformation = actualDeformation;

		this->mTotalNodalPosistion = ZeroVector(msLocalSize);
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
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateOnIntegrationPoints(
		const Variable<array_1d<double, 3 > >& rVariable,
		std::vector< array_1d<double, 3 > >& rOutput,
		const ProcessInfo& rCurrentProcessInfo) {

		KRATOS_TRY
		//element with two nodes can only represent results at one node 
		const unsigned int&  write_points_number = GetGeometry()
			.IntegrationPointsNumber(Kratos::GeometryData::GI_GAUSS_3);
		if (rOutput.size() != write_points_number) {
			rOutput.resize(write_points_number);
		}


		this->UpdateIncrementDeformation();
		//calculate Transformation Matrix
		bounded_matrix<double,msElementSize,msElementSize> TransformationMatrix = ZeroMatrix(msElementSize);
		this->CalculateTransformationMatrix(TransformationMatrix);
		//deformation modes
		bounded_vector<double,msLocalSize> elementForces_t = ZeroVector(msLocalSize);
		elementForces_t = this->CalculateElementForces();
		Vector Stress = ZeroVector(msElementSize);
		bounded_matrix<double,msElementSize,msLocalSize>  TransformationMatrixS = ZeroMatrix(msElementSize, msLocalSize);
		TransformationMatrixS = this->CalculateTransformationS();
		Stress = prod(TransformationMatrixS, elementForces_t);

		//LINEAR BEAM ELEMENT
		if (this->mIsLinearElement)
		{
			Matrix LeftHandSideMatrix = ZeroMatrix(msElementSize, msElementSize);
			LeftHandSideMatrix = this->mLHS;

			Vector NodalDeformation = ZeroVector(msElementSize);
			this->GetValuesVector(NodalDeformation);
			Stress = ZeroVector(msElementSize);
			Stress = prod(LeftHandSideMatrix, NodalDeformation);
			bounded_matrix<double,msElementSize,msElementSize> TransformationMatrix = ZeroMatrix(msElementSize);
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

			if (this->mIsLinearElement)
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
		bounded_matrix<double,
		CrBeamElement3D2N::msElementSize,CrBeamElement3D2N::msElementSize>& BigMatrix) {

		KRATOS_TRY
		for (unsigned int kk = 0; kk < msElementSize; kk += msDimension)
		{
			for (int i = 0; i<msDimension; ++i)
			{
				for (int j = 0; j<msDimension; ++j)
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
		const SizeType MatSize = msNumberOfNodes * 2;

		if (rMassMatrix.size1() != MatSize) {
			rMassMatrix.resize(MatSize, MatSize, false);
		}
		rMassMatrix = ZeroMatrix(MatSize, MatSize);
		bounded_matrix<double,MatSize,MatSize> TempMassMatrix = ZeroMatrix(MatSize, MatSize);
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
		const int smallMatSize = msNumberOfNodes * 2;

		if (rMassMatrix.size1() != msElementSize) {
			rMassMatrix.resize(msElementSize, msElementSize, false);
		}
		rMassMatrix = ZeroMatrix(msElementSize, msElementSize);

		const double L = this->CalculateReferenceLength();
		const double L2 = L * L;
		const double rho = this->GetProperties()[DENSITY];
		const double A = this->GetProperties()[CROSS_AREA];
		const double E = this->GetProperties()[YOUNG_MODULUS];
		Vector inertia = this->GetProperties()[LOCAL_INERTIA_VECTOR];
		const double J = inertia[0];
		const double Iy = inertia[1];
		const double Iz = inertia[2];
		const double G = this->CalculateShearModulus();

		double Ay = 0.00;
		if (this->GetProperties().Has(AREA_EFFECTIVE_Y)) {
			Ay = GetProperties()[AREA_EFFECTIVE_Y];
		}

		double Az = 0.00;
		if (this->GetProperties().Has(AREA_EFFECTIVE_Z)) {
			Az = GetProperties()[AREA_EFFECTIVE_Z];
		}

		double IRy = Iy;
		if (this->GetProperties().Has(INERTIA_ROT_Y)) {
			IRy = GetProperties()[INERTIA_ROT_Y];
		}

		double IRz = Iz;
		if (this->GetProperties().Has(INERTIA_ROT_Y)) {
			IRz = GetProperties()[INERTIA_ROT_Z];
		}

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
		const double M33 = (J*L*rho) / 3.00;
		const double M39 = M33 / 2.00;

		rMassMatrix(0, 0) = M00;
		rMassMatrix(0, 6) = M06;
		rMassMatrix(6, 6) = M00;
		rMassMatrix(3, 3) = M33;
		rMassMatrix(3, 9) = M39;
		rMassMatrix(9, 9) = M33;

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
		if (rMassMatrix.size1() != msElementSize) {
			rMassMatrix.resize(msElementSize, msElementSize, false);
		}
		rMassMatrix = ZeroMatrix(msElementSize, msElementSize);
		const double A = this->GetProperties()[CROSS_AREA];
		const double L = this->CalculateReferenceLength();
		const double rho = this->GetProperties()[DENSITY];

		const double TotalMass = A * L * rho;
		const double temp = 0.50 * TotalMass;

		//translatonal mass	
		for (int i = 0; i < msNumberOfNodes; ++i)
		{
			for (int j = 0; j < msDimension; ++j)
			{
				int index = i * (msDimension * 2) + j;
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
		if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL)
		{

			for (int i = 0; i< msNumberOfNodes; ++i)
			{
				int index = msLocalSize * i;

				GetGeometry()[i].SetLock();

				array_1d<double, 3> &ForceResidual =
					GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

				for (int j = 0; j<msDimension; ++j)
				{
					ForceResidual[j] += rRHSVector[index + j];
				}

				GetGeometry()[i].UnSetLock();
			}
		}


		if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == MOMENT_RESIDUAL)
		{

			for (int i = 0; i< msNumberOfNodes; ++i)
			{
				int index = (msLocalSize * i) + msDimension;

				GetGeometry()[i].SetLock();

				array_1d<double, 3> &MomentResidual =
					GetGeometry()[i].FastGetSolutionStepValue(MOMENT_RESIDUAL);

				for (int j = 0; j<msDimension; ++j)
				{
					MomentResidual[j] += rRHSVector[index + j];
				}

				GetGeometry()[i].UnSetLock();
			}
		}

		KRATOS_CATCH("")
	}

	double CrBeamElement3D2N::CalculateShearModulus() {
		KRATOS_TRY;
		const double nu = this->GetProperties()[POISSON_RATIO];
		const double E = this->GetProperties()[YOUNG_MODULUS];
		const double G = E / (2.0 * (1.0 + nu));
		return G;
		KRATOS_CATCH("")
	}

	int CrBeamElement3D2N::Check(const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

			if (GetGeometry().WorkingSpaceDimension() != 3 || GetGeometry().size() != 2)
			{
				KRATOS_ERROR <<
					"The beam element works only in 3D and with 2 noded elements" << ""
					<< std::endl;
			}
		//verify that the variables are correctly initialized
		if (VELOCITY.Key() == 0) {
			KRATOS_ERROR <<
				"VELOCITY has Key zero! (check if the application is correctly registered" << ""
				<< std::endl;
		}
		if (DISPLACEMENT.Key() == 0) {
			KRATOS_ERROR <<
				"DISPLACEMENT has Key zero! (check if the application is correctly registered"<< ""
				<< std::endl;
		}
		if (ACCELERATION.Key() == 0) {
			KRATOS_ERROR <<
				"ACCELERATION has Key zero! (check if the application is correctly registered" << ""
				<< std::endl;
		}
		if (DENSITY.Key() == 0) {
			KRATOS_ERROR <<
				"DENSITY has Key zero! (check if the application is correctly registered" << ""
				<< std::endl;
		}
		if (CROSS_AREA.Key() == 0) {
			KRATOS_ERROR <<
				"CROSS_AREA has Key zero! (check if the application is correctly registered" << ""
				<< std::endl;
		}
		//verify that the dofs exist
		for (unsigned int i = 0; i<this->GetGeometry().size(); ++i)
		{
			if (this->GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false) {
				KRATOS_ERROR <<
					"missing variable DISPLACEMENT on node " << this->GetGeometry()[i].Id()
					<< std::endl;
			}
			if (this->GetGeometry()[i].HasDofFor(DISPLACEMENT_X) == false ||
				this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false ||
				this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Z) == false) {
				KRATOS_ERROR <<
					"missing one of the dofs for the variable DISPLACEMENT on node " <<
						GetGeometry()[i].Id() << std::endl;
			}
		}



		if (this->GetProperties().Has(CROSS_AREA) == false ||
			this->GetProperties()[CROSS_AREA] == 0)
		{
			KRATOS_ERROR << "CROSS_AREA not provided for this element" << this->Id()
				<< std::endl;
		}

		if (this->GetProperties().Has(YOUNG_MODULUS) == false ||
			this->GetProperties()[YOUNG_MODULUS] == 0)
		{
			KRATOS_ERROR << "YOUNG_MODULUS not provided for this element" << this->Id()
				<< std::endl;
		}
		if (this->GetProperties().Has(DENSITY) == false)
		{
			KRATOS_ERROR << "DENSITY not provided for this element" << this->Id()
				<< std::endl;
		}

		if (this->GetProperties().Has(POISSON_RATIO) == false)
		{
			KRATOS_ERROR << "POISSON_RATIO not provided for this element" << this->Id()
				<< std::endl;
		}

		if (this->GetProperties().Has(LOCAL_INERTIA_VECTOR) == false)
		{
			KRATOS_ERROR << "LOCAL_INERTIA_VECTOR not provided for this element" << this->Id()
				<< std::endl;
		}
		else
		{
			Vector inertia = this->GetProperties()[LOCAL_INERTIA_VECTOR];
			if (inertia.size() != 3)
			{
				KRATOS_ERROR << "LOCAL_INERTIA_VECTOR must be of size 3 in element:" << this->Id()
				<< std::endl;
			}
		}

		return 0;

		KRATOS_CATCH("")
	}

	Orientation::Orientation(array_1d<double, Orientation::msDimension>& v1, const double theta) {

		KRATOS_TRY
		//!!!!!!!!!! if crossproduct with array_1d type switch input order !!!!!!!
		//If only direction of v1 is given -> Default case
		const double numerical_limit = std::numeric_limits<double>::epsilon();
		array_1d<double, msDimension> GlobalZ = ZeroVector(msDimension);
		GlobalZ[2] = 1.0;

		array_1d<double, msDimension> v2 = ZeroVector(msDimension);
		array_1d<double, msDimension> v3 = ZeroVector(msDimension);

		double VectorNorm;
		VectorNorm = MathUtils<double>::Norm(v1);
		if (VectorNorm > numerical_limit) v1 /= VectorNorm;

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
			if (VectorNorm > numerical_limit) v2 /= VectorNorm;

			v3 = MathUtils<double>::CrossProduct(v2, v1);
			VectorNorm = MathUtils<double>::Norm(v3);
			if (VectorNorm > numerical_limit) v3 /= VectorNorm;
		}

		//manual rotation around the beam axis
		if (theta != 0) {
			const Vector nz_temp = v3;
			const Vector ny_temp = v2;
			const double CosTheta = cos(theta);
			const double SinTheta = sin(theta);

			v2 = ny_temp * CosTheta + nz_temp * SinTheta;
			VectorNorm = MathUtils<double>::Norm(v2);
			if (VectorNorm > numerical_limit) v2 /= VectorNorm;

			v3 = nz_temp * CosTheta - ny_temp * SinTheta;
			VectorNorm = MathUtils<double>::Norm(v3);
			if (VectorNorm > numerical_limit) v3 /= VectorNorm;
		}

		Matrix RotationMatrix = ZeroMatrix(msDimension);
		for (int i = 0; i < msDimension; ++i) {
			RotationMatrix(i, 0) = v1[i];
			RotationMatrix(i, 1) = v2[i];
			RotationMatrix(i, 2) = v3[i];
		}

		this->GetQuaternion() = Quaternion<double>::FromRotationMatrix(RotationMatrix);

		KRATOS_CATCH("")
	}

	Orientation::Orientation(array_1d<double, Orientation::msDimension>& v1, array_1d<double, Orientation::msDimension>& v2) {

		KRATOS_TRY
		//If the user defines an aditional direction v2
		const double numerical_limit = std::numeric_limits<double>::epsilon();
		array_1d<double, msDimension> v3 = ZeroVector(msDimension);

		double VectorNorm;
		VectorNorm = MathUtils<double>::Norm(v1);
		if (VectorNorm > numerical_limit) v1 /= VectorNorm;

		VectorNorm = MathUtils<double>::Norm(v2);
		if (VectorNorm > numerical_limit) v2 /= VectorNorm;

		v3 = MathUtils<double>::CrossProduct(v2, v1);
		VectorNorm = MathUtils<double>::Norm(v3);
		if (VectorNorm > numerical_limit) v3 /= VectorNorm;


		Matrix RotationMatrix = ZeroMatrix(msDimension);
		for (int i = 0; i < msDimension; ++i) {
			RotationMatrix(i, 0) = v1[i];
			RotationMatrix(i, 1) = v2[i];
			RotationMatrix(i, 2) = v3[i];
		}

		this->GetQuaternion() = Quaternion<double>::FromRotationMatrix(RotationMatrix);

		KRATOS_CATCH("")
	}

	void Orientation::CalculateRotationMatrix(Matrix& R) {

		KRATOS_TRY
			if (R.size1() != msDimension || R.size2() != msDimension) R.resize(msDimension, msDimension, false);
		const Quaternion<double> q = this->GetQuaternion();
		q.ToRotationMatrix(R);
		KRATOS_CATCH("")
	}

	void Orientation::CalculateBasisVectors(array_1d<double, Orientation::msDimension>& v1,
		array_1d<double, Orientation::msDimension>& v2,
		array_1d<double, Orientation::msDimension>& v3) {

		KRATOS_TRY
			const Quaternion<double> q = this->GetQuaternion();
		Matrix R = ZeroMatrix(msDimension);
		q.ToRotationMatrix(R);
		if (v1.size() != msDimension) v1.resize(msDimension, false);
		if (v2.size() != msDimension) v2.resize(msDimension, false);
		if (v3.size() != msDimension) v3.resize(msDimension, false);

		for (int i = 0; i < msDimension; ++i) {
			v1[i] = R(i, 0);
			v2[i] = R(i, 1);
			v3[i] = R(i, 2);
		}
		KRATOS_CATCH("")
	}


	void CrBeamElement3D2N::CalculateGeometricStiffnessMatrix(MatrixType& rGeometricStiffnessMatrix,
		ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		rGeometricStiffnessMatrix = ZeroMatrix(msElementSize, msElementSize);
		rGeometricStiffnessMatrix = this->CreateElementStiffnessMatrix_Geometry();
		KRATOS_CATCH("")
	}

	void CrBeamElement3D2N::CalculateElasticStiffnessMatrix(MatrixType& rElasticStiffnessMatrix,
		ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		rElasticStiffnessMatrix = ZeroMatrix(msElementSize, msElementSize);
		rElasticStiffnessMatrix = this->CreateElementStiffnessMatrix_Material();
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
		rSerializer.load("mIsLumpedMassMatrix", this->mIsLumpedMassMatrix);
	}

} // namespace Kratos.


