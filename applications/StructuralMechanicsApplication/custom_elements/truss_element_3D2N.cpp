// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors:
//                   
//                   
//
#include "custom_elements/truss_element_3D2N.hpp"
#include "structural_mechanics_application_variables.h"
#include "includes/define.h"


namespace Kratos
{
	TrussElement3D2N::TrussElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry) {}
	TrussElement3D2N::TrussElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Element(NewId, pGeometry, pProperties) {}
	Element::Pointer TrussElement3D2N::Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const
	{
		const GeometryType& rGeom = this->GetGeometry();
		return BaseType::Pointer(new TrussElement3D2N(NewId, rGeom.Create(rThisNodes),pProperties));
	}
	TrussElement3D2N::~TrussElement3D2N(){}
	void TrussElement3D2N::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		uint number_of_nodes = this->GetGeometry().PointsNumber();
		uint dimension = this->GetGeometry().WorkingSpaceDimension();
		uint dim = number_of_nodes * dimension;

		if (rResult.size() != dim)
			rResult.resize(dim);

		for (uint i = 0; i < number_of_nodes; i++)
		{
			int index = i * 3;
			rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
			rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
			rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
		}

	}
	void TrussElement3D2N::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo)
	{
		uint number_of_nodes = this->GetGeometry().PointsNumber();
		uint dimension = this->GetGeometry().WorkingSpaceDimension();
		uint dim = number_of_nodes * dimension;

		if (rElementalDofList.size() != dim)
			rElementalDofList.resize(dim);

		for (uint i = 0; i < GetGeometry().size(); i++)
		{
			int Index = i * 3;
			rElementalDofList[Index] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
			rElementalDofList[Index + 1] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
			rElementalDofList[Index + 2] = GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
		}
	}
	void TrussElement3D2N::Initialize() {
		KRATOS_TRY

		this->mArea = GetProperties()[CROSS_AREA];
		this->mYoungsModulus = GetProperties()[YOUNG_MODULUS];
		this->mLength = this->CalculateReferenceLength();
		this->mCurrentLength = this->CalculateCurrentLength();
		this->mDensity = GetProperties()[DENSITY];

		if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2) == false) this->mPreStress = 0.00;
		else this->mPreStress = GetProperties()[TRUSS_PRESTRESS_PK2];

		if (this->mLength == 0.00) KRATOS_THROW_ERROR(std::invalid_argument, "Zero length found in element #", this->Id());

		KRATOS_CATCH("")
	}
	TrussElement3D2N::MatrixType TrussElement3D2N::CreateElementStiffnessMatrix()
	{
		KRATOS_TRY
		const double E = this->mYoungsModulus;
		const double A = this->mArea;
		const double S_pre = this->mPreStress;
		MatrixType LocalStiffnessMatrix = ZeroMatrix(6, 6);

		double du, dv, dw, dx, dy, dz, l, L, e_GL, L2, L4, K_1;
		du = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X) - GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
		dv = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y) - GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
		dw = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z) - GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
		dx = GetGeometry()[1].X0() - GetGeometry()[0].X0();
		dy = GetGeometry()[1].Y0() - GetGeometry()[0].Y0();
		dz = GetGeometry()[1].Z0() - GetGeometry()[0].Z0();
		L = mLength;
		l = mCurrentLength;
		e_GL = (l*l - L*L) / (2.00 * L*L);
		L2 = L*L;
		L4 = L2*L2;
		//using PK2 -> GreenLagrange
		K_1 = e_GL*E + S_pre;

		LocalStiffnessMatrix(0, 0) = A*L*(K_1 / L2 + E*(dx + du)*(dx + du) / L4); //k11
		LocalStiffnessMatrix(3, 3) = LocalStiffnessMatrix(0, 0); //k44

		LocalStiffnessMatrix(1, 1) = A*L*(K_1 / L2 + E*(dy + dv)*(dy + dv) / L4); //k22
		LocalStiffnessMatrix(4, 4) = LocalStiffnessMatrix(1, 1); //k55

		LocalStiffnessMatrix(2, 2) = A*L*(K_1 / L2 + E*(dz + dw)*(dz + dw) / L4); //k33
		LocalStiffnessMatrix(5, 5) = LocalStiffnessMatrix(2, 2); //k66

		LocalStiffnessMatrix(0, 1) = A*L*((dx + du)*(dy + dv)*E / L4);		//k12
		LocalStiffnessMatrix(1, 0) = LocalStiffnessMatrix(0, 1);			 //k21

		LocalStiffnessMatrix(0, 2) = A*L*((dx + du)*(dz + dw)*E / L4); //k13
		LocalStiffnessMatrix(2, 0) = LocalStiffnessMatrix(0, 2); //k31

		LocalStiffnessMatrix(0, 3) = A*L*(-K_1 / L2 - E*(dx + du)*(dx + du) / L4); //k14
		LocalStiffnessMatrix(3, 0) = LocalStiffnessMatrix(0, 3); //k41

		LocalStiffnessMatrix(0, 4) = A*L*((-1.00)*(dx + du)*(dy + dv)*E / L4); //k15
		LocalStiffnessMatrix(4, 0) = LocalStiffnessMatrix(0, 4); //k51

		LocalStiffnessMatrix(0, 5) = A*L*((-1.00)*(dx + du)*(dz + dw)*E / L4); //k16
		LocalStiffnessMatrix(5, 0) = LocalStiffnessMatrix(0, 5); //k61

		LocalStiffnessMatrix(1, 2) = A*L*((dy + dv)*(dz + dw)*E / L4); //k23
		LocalStiffnessMatrix(2, 1) = LocalStiffnessMatrix(1, 2); //k32

		LocalStiffnessMatrix(1, 3) = A*L*((-1.00)*(dy + dv)*(dx + du)*E / L4); //k24
		LocalStiffnessMatrix(3, 1) = LocalStiffnessMatrix(1, 3); //k42

		LocalStiffnessMatrix(1, 4) = A*L*(-K_1 / L2 - E*(dy + dv)*(dy + dv) / L4);  //k25
		LocalStiffnessMatrix(4, 1) = LocalStiffnessMatrix(1, 4); //k52

		LocalStiffnessMatrix(1, 5) = A*L*((-1.00)*(dy + dv)*(dz + dw)*E / L4); //k26
		LocalStiffnessMatrix(5, 1) = LocalStiffnessMatrix(1, 5); //k62

		LocalStiffnessMatrix(2, 3) = A*L*((-1.00)*(dw + dz)*(dx + du)*E / L4); //k34
		LocalStiffnessMatrix(3, 2) = LocalStiffnessMatrix(2, 3); //k43

		LocalStiffnessMatrix(2, 4) = A*L*((-1.00)*(dw + dz)*(dy + dv)*E / L4);  //k35
		LocalStiffnessMatrix(4, 2) = LocalStiffnessMatrix(2, 4); //k53

		LocalStiffnessMatrix(2, 5) = A*L*(-K_1 / L2 - E*(dz + dw)*(dz + dw) / L4); //k36
		LocalStiffnessMatrix(5, 2) = LocalStiffnessMatrix(2, 5); //k63

		LocalStiffnessMatrix(3, 4) = A*L*((dx + du)*(dy + dv)*E / L4); //k45
		LocalStiffnessMatrix(4, 3) = LocalStiffnessMatrix(3, 4); //k54

		LocalStiffnessMatrix(3, 5) = A*L*((dx + du)*(dz + dw)*E / L4); //k46
		LocalStiffnessMatrix(5, 3) = LocalStiffnessMatrix(3, 5); //k64

		LocalStiffnessMatrix(4, 5) = A*L*((dy + dv)*(dz + dw)*E / L4); //k56
		LocalStiffnessMatrix(5, 4) = LocalStiffnessMatrix(4, 5); //k65

		return LocalStiffnessMatrix;
		KRATOS_CATCH("")
	}
	void TrussElement3D2N::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		// 0.-Initialize the DampingMatrix:
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

		// Resizing as needed the LHS
		const unsigned int MatSize = number_of_nodes * dimension;

		if (rDampingMatrix.size1() != MatSize)
		{
			rDampingMatrix.resize(MatSize, MatSize, false);
		}

		noalias(rDampingMatrix) = ZeroMatrix(MatSize, MatSize);

		// 1.-Calculate StiffnessMatrix:

		MatrixType StiffnessMatrix = ZeroMatrix(MatSize, MatSize);

		this->CalculateLeftHandSide(StiffnessMatrix, rCurrentProcessInfo);

		// 2.-Calculate MassMatrix:

		MatrixType MassMatrix = ZeroMatrix(MatSize, MatSize);

		this->CalculateMassMatrix(MassMatrix, rCurrentProcessInfo);

		// 3.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
		double alpha = 0.0;
		if (GetProperties().Has(RAYLEIGH_ALPHA))
		{
			alpha = GetProperties()[RAYLEIGH_ALPHA];
		}
		else if (rCurrentProcessInfo.Has(RAYLEIGH_ALPHA))
		{
			alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
		}

		double beta = 0.0;
		if (GetProperties().Has(RAYLEIGH_BETA))
		{
			beta = GetProperties()[RAYLEIGH_BETA];
		}
		else if (rCurrentProcessInfo.Has(RAYLEIGH_BETA))
		{
			beta = rCurrentProcessInfo[RAYLEIGH_BETA];
		}

		// 4.-Compose the Damping Matrix:

		// Rayleigh Damping Matrix: alpha*M + beta*K
		noalias(rDampingMatrix) += alpha * MassMatrix;
		noalias(rDampingMatrix) += beta  * StiffnessMatrix;

		KRATOS_CATCH("")
	}
	void TrussElement3D2N::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		const uint number_of_nodes = GetGeometry().size();
		const uint dimension = GetGeometry().WorkingSpaceDimension();
		uint MatSize = number_of_nodes * dimension;

		if (rMassMatrix.size1() != MatSize)
			rMassMatrix.resize(MatSize, MatSize, false);

		rMassMatrix = ZeroMatrix(MatSize, MatSize);

		double TotalMass = 0;
		TotalMass = this->mArea * this->mLength * this->mDensity;

		Vector LumpFact = ZeroVector(number_of_nodes);

		LumpFact = GetGeometry().LumpingFactors(LumpFact);

		

		for (uint i = 0; i < number_of_nodes; i++)
		{
			double temp = LumpFact[i] * TotalMass;

			for (uint j = 0; j < dimension; j++)
			{
				uint index = i *dimension + j;

				rMassMatrix(index, index) = temp;
			}
		}
		KRATOS_CATCH("")
	}
	TrussElement3D2N::VectorType TrussElement3D2N::CalculateBodyForces()
	{
		KRATOS_TRY
		const uint number_of_nodes = GetGeometry().size();
		const uint dimension = GetGeometry().WorkingSpaceDimension();
		uint MatSize = number_of_nodes * dimension;

		//getting shapefunctionvalues 
		const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

		//creating necessary values 
		double TotalMass = this->mArea * this->mLength * this->mDensity;
		VectorType BodyForcesNode = ZeroVector(dimension);
		VectorType BodyForcesGlobal = ZeroVector(MatSize);

		//assemble global Vector
		for (uint i = 0; i < number_of_nodes; i++) {
			BodyForcesNode = TotalMass*this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION)*Ncontainer(0,i);

			for (uint j = 0; j < dimension; j++) {
				BodyForcesGlobal[(i*dimension) + j] = BodyForcesNode[j];
			}
		}
		
		return BodyForcesGlobal;
		KRATOS_CATCH("")
	}
	void TrussElement3D2N::GetValuesVector(Vector& rValues, int Step)
	{
		KRATOS_TRY
		GeometryType& rGeometry = this->GetGeometry();
		const uint number_of_nodes = GetGeometry().PointsNumber();
		const uint dimension = GetGeometry().WorkingSpaceDimension();
		uint      element_size = number_of_nodes * dimension;

		if (rValues.size() != element_size) rValues.resize(element_size, false);

		for (uint i = 0; i < number_of_nodes; i++)
		{
			int index = i * dimension;
			rValues[index] = rGeometry[i].FastGetSolutionStepValue(DISPLACEMENT_X, Step);
			rValues[index + 1] = rGeometry[i].FastGetSolutionStepValue(DISPLACEMENT_Y, Step);
			rValues[index + 2] = rGeometry[i].FastGetSolutionStepValue(DISPLACEMENT_Z, Step);
		}
		KRATOS_CATCH("")
	}
	void TrussElement3D2N::GetFirstDerivativesVector(Vector& rValues, int Step)
	{
		KRATOS_TRY
		GeometryType& rGeometry = this->GetGeometry();
		const uint number_of_nodes = GetGeometry().PointsNumber();
		const uint dimension = GetGeometry().WorkingSpaceDimension();
		uint      element_size = number_of_nodes * dimension;

		if (rValues.size() != element_size) rValues.resize(element_size, false);

		for (uint i = 0; i < number_of_nodes; i++)
		{
			int index = i * dimension;
			rValues[index] = rGeometry[i].FastGetSolutionStepValue(VELOCITY_X, Step);
			rValues[index + 1] = rGeometry[i].FastGetSolutionStepValue(VELOCITY_Y, Step);
			rValues[index + 2] = rGeometry[i].FastGetSolutionStepValue(VELOCITY_Z, Step);
		}
		KRATOS_CATCH("")
	}
	void TrussElement3D2N::GetSecondDerivativesVector(Vector& rValues,int Step)
	{
		KRATOS_TRY
		GeometryType& rGeometry = this->GetGeometry();
		const uint number_of_nodes = GetGeometry().PointsNumber();
		const uint dimension = GetGeometry().WorkingSpaceDimension();
		uint       element_size = number_of_nodes * dimension;

		if (rValues.size() != element_size) rValues.resize(element_size, false);

		for (uint i = 0; i < number_of_nodes; i++)
		{
			uint index = i * dimension;
			rValues[index] = rGeometry[i].FastGetSolutionStepValue(ACCELERATION_X, Step);
			rValues[index + 1] = rGeometry[i].FastGetSolutionStepValue(ACCELERATION_Y, Step);
			rValues[index + 2] = rGeometry[i].FastGetSolutionStepValue(ACCELERATION_Z, Step);
		}

		KRATOS_CATCH("")
	}
	void TrussElement3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		const SizeType NumNodes = this->GetGeometry().PointsNumber();
		const uint dimension = this->GetGeometry().WorkingSpaceDimension();
		const SizeType LocalSize = NumNodes * dimension;

		this->Initialize();

		//resizing the matrices + create memory for LHS
		if (rLeftHandSideMatrix.size1() != LocalSize) rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		//creating LHS
		noalias(rLeftHandSideMatrix) = this->CreateElementStiffnessMatrix();
		mLHS = rLeftHandSideMatrix;

		//create+compute RHS
		VectorType currentDisp = ZeroVector(LocalSize);
		this->GetValuesVector(currentDisp);
		if (rRightHandSideVector.size() != LocalSize) rRightHandSideVector = ZeroVector(LocalSize);
		rRightHandSideVector = ZeroVector(LocalSize);
		//update Residual
		VectorType internalForces = ZeroVector(6);
		this->UpdateInternalForces(internalForces);
		noalias(rRightHandSideVector) -= internalForces;
		//add bodyforces 
		noalias(rRightHandSideVector) += this->CalculateBodyForces();
		this->mIterCount++;
		KRATOS_CATCH("")
	}
	void TrussElement3D2N::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		const SizeType NumNodes = this->GetGeometry().PointsNumber();
		const uint dimension = this->GetGeometry().WorkingSpaceDimension();
		const SizeType LocalSize = NumNodes * dimension;

		if (rRightHandSideVector.size() != LocalSize) rRightHandSideVector = ZeroVector(LocalSize);
		rRightHandSideVector = ZeroVector(LocalSize);

		VectorType currentDisp = ZeroVector(LocalSize);
		this->GetValuesVector(currentDisp, 0);
		rRightHandSideVector -= prod(this->mLHS, currentDisp);

		//add bodyforces 
		rRightHandSideVector += this->CalculateBodyForces();
		KRATOS_CATCH("")
	}
	void TrussElement3D2N::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		const SizeType NumNodes = this->GetGeometry().PointsNumber();
		const uint dimension = this->GetGeometry().WorkingSpaceDimension();
		const SizeType LocalSize = NumNodes * dimension;

		////resizing the matrices + create memory for LHS
		if (rLeftHandSideMatrix.size1() != LocalSize) rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		////creating LHS
		rLeftHandSideMatrix = this->mLHS;
		KRATOS_CATCH("")
	}

	void TrussElement3D2N::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
		if (rOutput.size() != integration_points.size()) rOutput.resize(integration_points.size());
		if (rVariable == TRUSS_PRESTRESS_PK2) rOutput[0] = this->mPreStress;
		KRATOS_CATCH("")
	}
	void TrussElement3D2N::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
		if (rOutput.size() != integration_points.size()) rOutput.resize(integration_points.size());
		if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR)
		{
			Vector Strain = ZeroVector(3);
			Strain[0] = this->CalculateGreenLagrangeStrain();
			Strain[1] = 0;
			Strain[2] = 0;
			rOutput[0] = Strain;
		}
		KRATOS_CATCH("")
	}
	void TrussElement3D2N::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		this->CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		KRATOS_CATCH("")
	}
	void TrussElement3D2N::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		this->CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		KRATOS_CATCH("")
	}




	///////////////////////////// additional dynamic functions /////////////////////////////

	void TrussElement3D2N::CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		const SizeType NumNodes = this->GetGeometry().PointsNumber();
		const uint dimension = this->GetGeometry().WorkingSpaceDimension();
		const SizeType LocalSize = NumNodes * dimension;

		if (rLeftHandSideMatrix.size1() != LocalSize) rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
		rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		this->CalculateMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);

		if (rRightHandSideVector.size() != LocalSize) rRightHandSideVector.resize(LocalSize);
		rRightHandSideVector = ZeroVector(LocalSize);

		double AlphaM = 0.0;
		Vector CurrentAccelerationVector = ZeroVector(LocalSize);
		this->GetSecondDerivativesVector(CurrentAccelerationVector, 0);
		if (rCurrentProcessInfo.Has(BOSSAK_ALPHA)) {
			AlphaM = rCurrentProcessInfo[BOSSAK_ALPHA];
			Vector PreviousAccelerationVector = ZeroVector(LocalSize);
			this->GetSecondDerivativesVector(PreviousAccelerationVector, 1);
			CurrentAccelerationVector *= (1.0 - AlphaM);
			CurrentAccelerationVector += AlphaM * (PreviousAccelerationVector);
		}

		rRightHandSideVector = prod(rLeftHandSideMatrix, CurrentAccelerationVector);	
		KRATOS_CATCH("")
	}

	void TrussElement3D2N::CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		const SizeType NumNodes = this->GetGeometry().PointsNumber();
		const uint dimension = this->GetGeometry().WorkingSpaceDimension();
		const SizeType LocalSize = NumNodes * dimension;

		MatrixType LeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		this->CalculateMassMatrix(LeftHandSideMatrix, rCurrentProcessInfo);

		if (rRightHandSideVector.size() != LocalSize) rRightHandSideVector.resize(LocalSize);
		rRightHandSideVector = ZeroVector(LocalSize);

		Vector CurrentAccelerationVector = ZeroVector(LocalSize);
		this->GetSecondDerivativesVector(CurrentAccelerationVector, 0);

		double AlphaM = 0.0;
		if (rCurrentProcessInfo.Has(BOSSAK_ALPHA)) {
			AlphaM = rCurrentProcessInfo[BOSSAK_ALPHA];
			Vector PreviousAccelerationVector = ZeroVector(LocalSize);
			this->GetSecondDerivativesVector(PreviousAccelerationVector, 1);
			CurrentAccelerationVector *= (1.0 - AlphaM);
			CurrentAccelerationVector += AlphaM * (PreviousAccelerationVector);
		}

		rRightHandSideVector = prod(LeftHandSideMatrix, CurrentAccelerationVector);

		KRATOS_CATCH("")
	}

	void TrussElement3D2N::CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		const SizeType NumNodes = this->GetGeometry().PointsNumber();
		const uint dimension = this->GetGeometry().WorkingSpaceDimension();
		const SizeType LocalSize = NumNodes * dimension;

		if (rLeftHandSideMatrix.size1() != LocalSize) rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
		rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		this->CalculateMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);		
		KRATOS_CATCH("")
	}






	int  TrussElement3D2N::Check(const ProcessInfo& rCurrentProcessInfo)
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
	double TrussElement3D2N::CalculateGreenLagrangeStrain()
	{
		KRATOS_TRY
		double l = this->mCurrentLength;
		double L = this->mLength;
		//longitudinal green lagrange strain
		return ((l * l - L * L) / (2.00 * L * L));
		KRATOS_CATCH("")
	}

	double TrussElement3D2N::CalculateReferenceLength()
	{
		KRATOS_TRY

		double dx, dy, dz, L;
		dx = GetGeometry()[1].X0() - GetGeometry()[0].X0();
		dy = GetGeometry()[1].Y0() - GetGeometry()[0].Y0();
		dz = GetGeometry()[1].Z0() - GetGeometry()[0].Z0();
		L = sqrt(dx*dx + dy*dy + dz*dz);
		return L;
		KRATOS_CATCH("")
	}
	double TrussElement3D2N::CalculateCurrentLength()
	{
		KRATOS_TRY
		double du, dv, dw, dx, dy, dz, l;
		du = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X) - GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
		dv = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y) - GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
		dw = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z) - GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
		dx = GetGeometry()[1].X0() - GetGeometry()[0].X0();
		dy = GetGeometry()[1].Y0() - GetGeometry()[0].Y0();
		dz = GetGeometry()[1].Z0() - GetGeometry()[0].Z0();
		l = sqrt((du + dx)*(du + dx) + (dv + dy)*(dv + dy) + (dw + dz)*(dw + dz));
		return l;
		KRATOS_CATCH("")
	}
	void TrussElement3D2N::UpdateInternalForces(VectorType& rinternalForces)
	{
		KRATOS_TRY
		MatrixType TransformationMatrix = ZeroMatrix(6, 6);
		this->CreateTransformationMatrix(TransformationMatrix);
		if (mIterCount == 0) mInternalStrainGL = 0.00;
		else mInternalStrainGL = this->CalculateGreenLagrangeStrain();
		double l = mCurrentLength;
		double L0 = mLength;
		double E = mYoungsModulus;
		double A = mArea;
		double S_pre = mPreStress;
		double N = ((E*mInternalStrainGL + S_pre) * l * A) / L0;
		//internal force vectors
		VectorType f_local = ZeroVector(6);
		f_local[0] = -1.00 * N;
		f_local[3] = 1.00 * N;

		rinternalForces = ZeroVector(6);
		rinternalForces = prod(TransformationMatrix, f_local);
		KRATOS_CATCH("");
	}

	void TrussElement3D2N::CreateTransformationMatrix(Matrix& rRotationMatrix)
	{
		KRATOS_TRY
		//1st calculate transformation matrix
		Vector DirectionVectorX = ZeroVector(3);
		Vector DirectionVectorY = ZeroVector(3);
		Vector DirectionVectorZ = ZeroVector(3);
		Vector ReferenceCoordinates = ZeroVector(6);
		Vector GlobalZ = ZeroVector(3);
		GlobalZ[2] = 1.0;

		ReferenceCoordinates[0] = GetGeometry()[0].X();
		ReferenceCoordinates[1] = GetGeometry()[0].Y();
		ReferenceCoordinates[2] = GetGeometry()[0].Z();
		ReferenceCoordinates[3] = GetGeometry()[1].X();
		ReferenceCoordinates[4] = GetGeometry()[1].Y();
		ReferenceCoordinates[5] = GetGeometry()[1].Z();

		for (unsigned int i = 0; i < 3; i++)
		{
			DirectionVectorX[i] = (ReferenceCoordinates[i + 3] - ReferenceCoordinates[i]);
		}
		// local x-axis (e1_local) is the beam axis  (in GID is e3_local)
		double VectorNorm;
		VectorNorm = MathUtils<double>::Norm(DirectionVectorX);
		if (VectorNorm != 0) DirectionVectorX /= VectorNorm;

		double tolerance = 1.0 / 1000.0;
		if ((fabs(DirectionVectorX[0]) < tolerance) && (fabs(DirectionVectorX[1]) < tolerance)) {
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
			DirectionVectorY = MathUtils<double>::CrossProduct(GlobalZ, DirectionVectorX);
			VectorNorm = MathUtils<double>::Norm(DirectionVectorY);
			DirectionVectorY /= VectorNorm;

			DirectionVectorZ = MathUtils<double>::CrossProduct(DirectionVectorX, DirectionVectorY);
			VectorNorm = MathUtils<double>::Norm(DirectionVectorZ);
			DirectionVectorZ /= VectorNorm;
		}

		//2nd fill big rotation matrix
		MatrixType IninitialCS = ZeroMatrix(3, 3);
		for (uint i = 0; i < 3; i++) {
			IninitialCS(i, 0) = DirectionVectorX[i];
			IninitialCS(i, 1) = DirectionVectorY[i];
			IninitialCS(i, 2) = DirectionVectorZ[i];
		}

		rRotationMatrix = ZeroMatrix(6, 6);
		if (rRotationMatrix.size1() != 6) rRotationMatrix.resize(6, 6, false);
		//Building the rotation matrix for the local element matrix
		for (unsigned int kk = 0; kk < 6; kk += 3)
		{
			for (unsigned int i = 0; i<3; i++)
			{
				for (unsigned int j = 0; j<3; j++)
				{
					rRotationMatrix(i + kk, j + kk) = IninitialCS(i, j);
				}
			}
		}
		KRATOS_CATCH("")

	}	
} 


