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
#include "custom_elements/truss_element_3D2N.hpp"
#include "structural_mechanics_application_variables.h"
#include "includes/define.h"


namespace Kratos
{
	TrussElement3D2N::TrussElement3D2N(IndexType NewId, 
									GeometryType::Pointer pGeometry)
									: Element(NewId, pGeometry) {}

	TrussElement3D2N::TrussElement3D2N(IndexType NewId,
									GeometryType::Pointer pGeometry,
									PropertiesType::Pointer pProperties) 
									: Element(NewId, pGeometry, pProperties) {}

	Element::Pointer TrussElement3D2N::Create(IndexType NewId,
									NodesArrayType const& rThisNodes,
									PropertiesType::Pointer pProperties) const{
		const GeometryType& rGeom = this->GetGeometry();
		return BaseType::Pointer(new TrussElement3D2N(
			NewId, rGeom.Create(rThisNodes),pProperties));
	}

	TrussElement3D2N::~TrussElement3D2N(){}

	void TrussElement3D2N::EquationIdVector(EquationIdVectorType& rResult,
									ProcessInfo& rCurrentProcessInfo){

		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int local_size = number_of_nodes * dimension;

		if (rResult.size() != local_size) rResult.resize(local_size);

		for (int i = 0; i < number_of_nodes; ++i)
		{
			int index = i * 3;
			rResult[index] = this->GetGeometry()[i].GetDof(DISPLACEMENT_X)
				.EquationId();
			rResult[index + 1] = this->GetGeometry()[i].GetDof(DISPLACEMENT_Y)
				.EquationId();
			rResult[index + 2] = this->GetGeometry()[i].GetDof(DISPLACEMENT_Z)
				.EquationId();
		}

	}
	void TrussElement3D2N::GetDofList(DofsVectorType& rElementalDofList,
									ProcessInfo& rCurrentProcessInfo){

		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int local_size = number_of_nodes * dimension;

		if (rElementalDofList.size() != local_size) {
			rElementalDofList.resize(local_size);
		}

		for (int i = 0; i < GetGeometry().size(); ++i)
		{
			int index = i * 3;
			rElementalDofList[index] = this->GetGeometry()[i]
				.pGetDof(DISPLACEMENT_X);
			rElementalDofList[index + 1] = this->GetGeometry()[i]
				.pGetDof(DISPLACEMENT_Y);
			rElementalDofList[index + 2] = this->GetGeometry()[i]
				.pGetDof(DISPLACEMENT_Z);
		}
	}

	void TrussElement3D2N::Initialize() {

		KRATOS_TRY
		this->mArea = this->GetProperties()[CROSS_AREA];
		this->mYoungsModulus = this->GetProperties()[YOUNG_MODULUS];
		this->mLength = this->CalculateReferenceLength();
		this->mCurrentLength = this->CalculateCurrentLength();
		this->mDensity = this->GetProperties()[DENSITY];

		if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2) == false) {
			this->mPreStress = 0.00;
		}
		else this->mPreStress = this->GetProperties()[TRUSS_PRESTRESS_PK2];

		if (this->GetProperties().Has(TRUSS_IS_CABLE) == false) {
			this->mIsCable = false;
		}
		else this->mIsCable = this->GetProperties()[TRUSS_IS_CABLE];

		if (this->mLength == 0.00) {
			KRATOS_ERROR << ("Zero length found in element #", this->Id()) <<
				std::endl;
		}

		KRATOS_CATCH("")
	}

	TrussElement3D2N::MatrixType TrussElement3D2N::CreateElementStiffnessMatrix(){

		KRATOS_TRY
		const double E = this->mYoungsModulus;
		double A = this->mArea;
		const double S_pre = this->mPreStress;
		MatrixType LocalStiffnessMatrix = ZeroMatrix(6, 6);

		// du... delta displacement in x-direction
		// dv... delta displacement in y-direction
		// dw... delta displacement in z-direction
		// L... inital member length
		// l... deformed member length
		// e_gl... green_lagrange strain

		const double du = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
		const double dv = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
		const double dw = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
		const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
		const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
		const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
		const double L = this->mLength;
		const double l = this->mCurrentLength;
		const double e_gL = (l*l - L*L) / (2.00 * L*L);
		const double L2 = L*L;
		const double L4 = L2*L2;

		const double K_1 = e_gL*E + S_pre;

		//if cable + compressed -> no contribution to global K
		if (this->mIsCable == true && this->mIsCompressed == true) A = 0;

		LocalStiffnessMatrix(0, 0) = A*L*(K_1 / L2 + E*(dx + du)*(dx + du) / L4); 
		LocalStiffnessMatrix(3, 3) = LocalStiffnessMatrix(0, 0); 

		LocalStiffnessMatrix(1, 1) = A*L*(K_1 / L2 + E*(dy + dv)*(dy + dv) / L4); 
		LocalStiffnessMatrix(4, 4) = LocalStiffnessMatrix(1, 1); 

		LocalStiffnessMatrix(2, 2) = A*L*(K_1 / L2 + E*(dz + dw)*(dz + dw) / L4); 
		LocalStiffnessMatrix(5, 5) = LocalStiffnessMatrix(2, 2); 

		LocalStiffnessMatrix(0, 1) = A*L*((dx + du)*(dy + dv)*E / L4);		
		LocalStiffnessMatrix(1, 0) = LocalStiffnessMatrix(0, 1);			 

		LocalStiffnessMatrix(0, 2) = A*L*((dx + du)*(dz + dw)*E / L4); 
		LocalStiffnessMatrix(2, 0) = LocalStiffnessMatrix(0, 2); 

		LocalStiffnessMatrix(0, 3) = A*L*(-K_1 / L2 - E*(dx + du)*(dx + du) / L4); 
		LocalStiffnessMatrix(3, 0) = LocalStiffnessMatrix(0, 3); 

		LocalStiffnessMatrix(0, 4) = A*L*((-1.00)*(dx + du)*(dy + dv)*E / L4); 
		LocalStiffnessMatrix(4, 0) = LocalStiffnessMatrix(0, 4); 

		LocalStiffnessMatrix(0, 5) = A*L*((-1.00)*(dx + du)*(dz + dw)*E / L4); 
		LocalStiffnessMatrix(5, 0) = LocalStiffnessMatrix(0, 5); 

		LocalStiffnessMatrix(1, 2) = A*L*((dy + dv)*(dz + dw)*E / L4); 
		LocalStiffnessMatrix(2, 1) = LocalStiffnessMatrix(1, 2); 

		LocalStiffnessMatrix(1, 3) = A*L*((-1.00)*(dy + dv)*(dx + du)*E / L4); 
		LocalStiffnessMatrix(3, 1) = LocalStiffnessMatrix(1, 3); 

		LocalStiffnessMatrix(1, 4) = A*L*(-K_1 / L2 - E*(dy + dv)*(dy + dv) / L4);  
		LocalStiffnessMatrix(4, 1) = LocalStiffnessMatrix(1, 4); 

		LocalStiffnessMatrix(1, 5) = A*L*((-1.00)*(dy + dv)*(dz + dw)*E / L4); 
		LocalStiffnessMatrix(5, 1) = LocalStiffnessMatrix(1, 5); 

		LocalStiffnessMatrix(2, 3) = A*L*((-1.00)*(dw + dz)*(dx + du)*E / L4); 
		LocalStiffnessMatrix(3, 2) = LocalStiffnessMatrix(2, 3); 

		LocalStiffnessMatrix(2, 4) = A*L*((-1.00)*(dw + dz)*(dy + dv)*E / L4);  
		LocalStiffnessMatrix(4, 2) = LocalStiffnessMatrix(2, 4); 

		LocalStiffnessMatrix(2, 5) = A*L*(-K_1 / L2 - E*(dz + dw)*(dz + dw) / L4); 
		LocalStiffnessMatrix(5, 2) = LocalStiffnessMatrix(2, 5); 

		LocalStiffnessMatrix(3, 4) = A*L*((dx + du)*(dy + dv)*E / L4); 
		LocalStiffnessMatrix(4, 3) = LocalStiffnessMatrix(3, 4); 

		LocalStiffnessMatrix(3, 5) = A*L*((dx + du)*(dz + dw)*E / L4); 
		LocalStiffnessMatrix(5, 3) = LocalStiffnessMatrix(3, 5); 

		LocalStiffnessMatrix(4, 5) = A*L*((dy + dv)*(dz + dw)*E / L4); 
		LocalStiffnessMatrix(5, 4) = LocalStiffnessMatrix(4, 5); 

		return LocalStiffnessMatrix;
		KRATOS_CATCH("")
	}

	void TrussElement3D2N::CalculateDampingMatrix(MatrixType& rDampingMatrix,
									ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().size();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int MatSize = number_of_nodes * dimension;

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

	void TrussElement3D2N::CalculateMassMatrix(MatrixType& rMassMatrix,
									ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().size();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int MatSize = number_of_nodes * dimension;

		if (rMassMatrix.size1() != MatSize) {
			rMassMatrix.resize(MatSize, MatSize, false);
		}

		rMassMatrix = ZeroMatrix(MatSize, MatSize);

		double TotalMass = 0.00;
		TotalMass = this->mArea * this->mLength * this->mDensity;

		Vector LumpFact = ZeroVector(number_of_nodes);

		LumpFact = this->GetGeometry().LumpingFactors(LumpFact);

		for (int i = 0; i < number_of_nodes; ++i)
		{
			double temp = LumpFact[i] * TotalMass;

			for (int j = 0; j < dimension; ++j)
			{
				int index = i *dimension + j;

				rMassMatrix(index, index) = temp;
			}
		}
		KRATOS_CATCH("")
	}

	TrussElement3D2N::VectorType TrussElement3D2N::CalculateBodyForces(){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().size();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int MatSize = number_of_nodes * dimension;

		//getting shapefunctionvalues 
		const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues();

		//creating necessary values 
		double TotalMass = this->mArea * this->mLength * this->mDensity;
		VectorType BodyForcesNode = ZeroVector(dimension);
		VectorType BodyForcesGlobal = ZeroVector(MatSize);

		//assemble global Vector
		for (int i = 0; i < number_of_nodes; ++i) {
			BodyForcesNode = TotalMass*this->GetGeometry()[i]
				.FastGetSolutionStepValue(VOLUME_ACCELERATION)*Ncontainer(0,i);

			for (int j = 0; j < dimension; ++j) {
				BodyForcesGlobal[(i*dimension) + j] = BodyForcesNode[j];
			}
		}
		
		return BodyForcesGlobal;
		KRATOS_CATCH("")
	}

	void TrussElement3D2N::GetValuesVector(Vector& rValues, int Step){

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

	void TrussElement3D2N::GetFirstDerivativesVector(Vector& rValues, int Step){

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

	void TrussElement3D2N::GetSecondDerivativesVector(Vector& rValues,int Step){

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

	void TrussElement3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
										VectorType& rRightHandSideVector,
										ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int LocalSize = NumNodes * dimension;

		this->Initialize();
		//calculate internal forces
		VectorType InternalForces = ZeroVector(6);
		this->UpdateInternalForces(InternalForces);
		//resizing the matrices + create memory for LHS
		if (rLeftHandSideMatrix.size1() != LocalSize) {
			rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		}
		rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		//creating LHS
		noalias(rLeftHandSideMatrix) = this->CreateElementStiffnessMatrix();
		this->mLHS = rLeftHandSideMatrix;

		//create+compute RHS
		if (rRightHandSideVector.size() != LocalSize) {
			rRightHandSideVector = ZeroVector(LocalSize);
		}
		rRightHandSideVector = ZeroVector(LocalSize);
		//update Residual
		noalias(rRightHandSideVector) -= InternalForces;
		//add bodyforces 
		noalias(rRightHandSideVector) += this->CalculateBodyForces();

		if (this->mIsCable == true && this->mIsCompressed == true) {
			rRightHandSideVector = ZeroVector(LocalSize);
		}

		this->mIterCount++;
		KRATOS_CATCH("")
	}

	void TrussElement3D2N::CalculateRightHandSide(VectorType& rRightHandSideVector,
										ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int LocalSize = NumNodes * dimension;

		if (rRightHandSideVector.size() != LocalSize) {
			rRightHandSideVector = ZeroVector(LocalSize);
		}
		rRightHandSideVector = ZeroVector(LocalSize);

		VectorType InternalForces = ZeroVector(6);
		this->UpdateInternalForces(InternalForces);
		noalias(rRightHandSideVector) -= InternalForces;

		//add bodyforces 
		noalias(rRightHandSideVector) += this->CalculateBodyForces();
		if (this->mIsCable == true && this->mIsCompressed == true) {
			rRightHandSideVector = ZeroVector(LocalSize);
		}
		KRATOS_CATCH("")
	}

	void TrussElement3D2N::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
										ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int LocalSize = NumNodes * dimension;

		////resizing the matrices + create memory for LHS
		if (rLeftHandSideMatrix.size1() != LocalSize) {
			rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		}
		rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		////creating LHS
		noalias(rLeftHandSideMatrix) = this->mLHS;
		KRATOS_CATCH("")
	}

	void TrussElement3D2N::CalculateOnIntegrationPoints(
										const Variable<double>& rVariable,
										std::vector<double>& rOutput,
										const ProcessInfo& rCurrentProcessInfo){
		KRATOS_TRY
		const GeometryType::IntegrationPointsArrayType& integration_points =
			GetGeometry().IntegrationPoints();

		if (rOutput.size() != integration_points.size()) {
			rOutput.resize(integration_points.size());
		}
		if (rVariable == TRUSS_PRESTRESS_PK2) rOutput[0] = this->mPreStress;
		KRATOS_CATCH("")
	}

	void TrussElement3D2N::CalculateOnIntegrationPoints(
										const Variable<Vector>& rVariable,
										std::vector<Vector>& rOutput,
										const ProcessInfo& rCurrentProcessInfo){
		KRATOS_TRY
		const GeometryType::IntegrationPointsArrayType& integration_points =
			GetGeometry().IntegrationPoints();
		if (rOutput.size() != integration_points.size()) {
			rOutput.resize(integration_points.size());
		}
		if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR)
		{
			Vector Strain = ZeroVector(3);
			Strain[0] = this->CalculateGreenLagrangeStrain();
			Strain[1] = 0.00;
			Strain[2] = 0.00;
			rOutput[0] = Strain;
		}
		KRATOS_CATCH("")
	}

	void TrussElement3D2N::GetValueOnIntegrationPoints(
										const Variable<double>& rVariable,
										std::vector<double>& rValues,
										const ProcessInfo& rCurrentProcessInfo){
		KRATOS_TRY
		this->CalculateOnIntegrationPoints(rVariable, rValues,
										rCurrentProcessInfo);
		KRATOS_CATCH("")
	}
	void TrussElement3D2N::GetValueOnIntegrationPoints(
										const Variable<Vector>& rVariable,
										std::vector<Vector>& rValues,
										const ProcessInfo& rCurrentProcessInfo){
		KRATOS_TRY
		this->CalculateOnIntegrationPoints(rVariable, rValues,
										rCurrentProcessInfo);
		KRATOS_CATCH("")
	}


	int  TrussElement3D2N::Check(const ProcessInfo& rCurrentProcessInfo){
		KRATOS_TRY

		if (GetGeometry().WorkingSpaceDimension() != 3 || GetGeometry().size() != 2)
			{
				KRATOS_THROW_ERROR(std::invalid_argument,
					"The truss element works only in 3D and with 2 noded elements", "")
			}
		//verify that the variables are correctly initialized
		if (VELOCITY.Key() == 0)
			KRATOS_ERROR << ("VELOCITY has Key zero! (check if the application is correctly registered", "") << std::endl;
		if (DISPLACEMENT.Key() == 0)
			KRATOS_ERROR << ("DISPLACEMENT has Key zero! (check if the application is correctly registered", "") << std::endl;
		if (ACCELERATION.Key() == 0)
			KRATOS_ERROR << ("ACCELERATION has Key zero! (check if the application is correctly registered", "") << std::endl;
		if (DENSITY.Key() == 0)
			KRATOS_ERROR << ("DENSITY has Key zero! (check if the application is correctly registered", "") << std::endl;
		if (CROSS_AREA.Key() == 0)
			KRATOS_ERROR << ("CROSS_AREA has Key zero! (check if the application is correctly registered", "") << std::endl;
		//verify that the dofs exist
		for (unsigned int i = 0; i<this->GetGeometry().size(); ++i)
			{
			if (this->GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false)
				KRATOS_ERROR << ("missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id()) << std::endl;
			if (this->GetGeometry()[i].HasDofFor(DISPLACEMENT_X) == false || this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false || this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Z) == false)
				KRATOS_ERROR << ("missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id()) << std::endl;
			}


		
		if (this->GetProperties().Has(CROSS_AREA) == false ||
			this->GetProperties()[CROSS_AREA] == 0)
		{
			KRATOS_ERROR << ( "CROSS_AREA not provided for this element", this->Id()) << std::endl;
		}

		if (this->GetProperties().Has(YOUNG_MODULUS) == false ||
			this->GetProperties()[YOUNG_MODULUS] == 0)
		{
			KRATOS_ERROR << ("YOUNG_MODULUS not provided for this element", this->Id()) << std::endl;
		}
		if (this->GetProperties().Has(DENSITY) == false)
		{
			KRATOS_ERROR << ("DENSITY not provided for this element", this->Id()) << std::endl;
		}


		return 0;

		KRATOS_CATCH("")
	}
	double TrussElement3D2N::CalculateGreenLagrangeStrain(){

		KRATOS_TRY
		const double l = this->CalculateCurrentLength();
		const double L = this->mLength;
		const double e = ((l * l - L * L) / (2.00f * L * L));
		return e;
		KRATOS_CATCH("")
	}

	double TrussElement3D2N::CalculateReferenceLength(){

		KRATOS_TRY
		const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
		const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
		const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
		const double L = sqrt(dx*dx + dy*dy + dz*dz);
		return L;
		KRATOS_CATCH("")
	}
	double TrussElement3D2N::CalculateCurrentLength(){

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
		const double l = sqrt((du + dx)*(du + dx) + (dv + dy)*(dv + dy) 
			+ (dw + dz)*(dw + dz));
		return l;
		KRATOS_CATCH("")
	}
	void TrussElement3D2N::UpdateInternalForces(VectorType& rinternalForces){

		KRATOS_TRY
		MatrixType TransformationMatrix = ZeroMatrix(6, 6);
		this->CreateTransformationMatrix(TransformationMatrix);
		const double InternalStrainGL = this->CalculateGreenLagrangeStrain();
		const double l = this->CalculateCurrentLength();
		const double L0 = this->mLength;
		const double E = this->mYoungsModulus;
		const double A = this->mArea;
		const double S_pre = this->mPreStress;
		const double N = ((E*InternalStrainGL + S_pre) * l * A) / L0;

		if (N < 0.00) this->mIsCompressed = true;
		else this->mIsCompressed = false;

		//internal force vectors
		VectorType f_local = ZeroVector(6);
		f_local[0] = -1.00 * N;
		f_local[3] = 1.00 * N;

		rinternalForces = ZeroVector(6);
		noalias(rinternalForces) = prod(TransformationMatrix, f_local);
		KRATOS_CATCH("");
	}

	void TrussElement3D2N::CreateTransformationMatrix(Matrix& rRotationMatrix){
		KRATOS_TRY

		//1st calculate transformation matrix
		Vector DirectionVectorX = ZeroVector(3);
		Vector DirectionVectorY = ZeroVector(3);
		Vector DirectionVectorZ = ZeroVector(3);
		Vector ReferenceCoordinates = ZeroVector(6);
		Vector GlobalZ = ZeroVector(3);
		GlobalZ[2] = 1.0;

		ReferenceCoordinates[0] = this->GetGeometry()[0].X();
		ReferenceCoordinates[1] = this->GetGeometry()[0].Y();
		ReferenceCoordinates[2] = this->GetGeometry()[0].Z();
		ReferenceCoordinates[3] = this->GetGeometry()[1].X();
		ReferenceCoordinates[4] = this->GetGeometry()[1].Y();
		ReferenceCoordinates[5] = this->GetGeometry()[1].Z();

		for (unsigned int i = 0; i < 3; ++i)
		{
			DirectionVectorX[i] = (ReferenceCoordinates[i + 3] -
				ReferenceCoordinates[i]);
		}
		// local x-axis (e1_local) is the beam axis  (in GID is e3_local)
		double VectorNorm;
		VectorNorm = MathUtils<double>::Norm(DirectionVectorX);
		if (VectorNorm != 0) DirectionVectorX /= VectorNorm;

		double tolerance = 1.0 / 1000.0;
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

		//2nd fill big rotation matrix
		MatrixType CurrentCS = ZeroMatrix(3, 3);
		for (int i = 0; i < 3; ++i) {
			CurrentCS(i, 0) = DirectionVectorX[i];
			CurrentCS(i, 1) = DirectionVectorY[i];
			CurrentCS(i, 2) = DirectionVectorZ[i];
		}

		rRotationMatrix = ZeroMatrix(6, 6);
		if (rRotationMatrix.size1() != 6) rRotationMatrix.resize(6, 6, false);
		//Building the rotation matrix for the local element matrix
		for (unsigned int kk = 0; kk < 6; kk += 3)
		{
			for (unsigned int i = 0; i<3; ++i)
			{
				for (unsigned int j = 0; j<3; ++j)
				{
					rRotationMatrix(i + kk, j + kk) = CurrentCS(i, j);
				}
			}
		}
		KRATOS_CATCH("")

	}	

} // namespace Kratos.


