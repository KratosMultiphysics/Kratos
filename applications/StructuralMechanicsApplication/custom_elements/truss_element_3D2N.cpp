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
#include "custom_elements/truss_element_3D2N.hpp"
#include "structural_mechanics_application_variables.h"
#include "includes/define.h"


namespace Kratos
{
	TrussElement3D2N::TrussElement3D2N(IndexType NewId, 
									GeometryType::Pointer pGeometry,
									bool rLinear)
									: Element(NewId, pGeometry)
	{
		this->mIsLinearElement = rLinear;
	}

	TrussElement3D2N::TrussElement3D2N(IndexType NewId,
									GeometryType::Pointer pGeometry,
									PropertiesType::Pointer pProperties,
									bool rLinear) 
									: Element(NewId, pGeometry, pProperties) 
	{
		this->mIsLinearElement = rLinear;
	}

	Element::Pointer TrussElement3D2N::Create(IndexType NewId,
									NodesArrayType const& rThisNodes,
									PropertiesType::Pointer pProperties) const
	{
		const GeometryType& rGeom = this->GetGeometry();
		return BaseType::Pointer(new TrussElement3D2N(
			NewId, rGeom.Create(rThisNodes),pProperties, this->mIsLinearElement));
	}

	TrussElement3D2N::~TrussElement3D2N(){}

	void TrussElement3D2N::EquationIdVector(EquationIdVectorType& rResult,
									ProcessInfo& rCurrentProcessInfo){

		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int local_size = number_of_nodes * dimension;

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
		const unsigned int local_size = number_of_nodes * dimension;

		if (rElementalDofList.size() != local_size) {
			rElementalDofList.resize(local_size);
		}

		for (int i = 0; i < number_of_nodes; ++i)
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
		KRATOS_CATCH("")
	}

	TrussElement3D2N::MatrixType TrussElement3D2N::CreateElementStiffnessMatrix(
		ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int local_size = number_of_nodes * dimension;

		MatrixType LocalStiffnessMatrix = ZeroMatrix(local_size, local_size);



		if (this->ReturnIfIsCable() == true && this->mIsCompressed == true)
		{
			LocalStiffnessMatrix = ZeroMatrix(local_size, local_size);
		}

		else
		{
			this->CalculateElasticStiffnessMatrix(LocalStiffnessMatrix,
				rCurrentProcessInfo);
			if (this->mIsLinearElement == false) 
			{
				MatrixType K_geo = ZeroMatrix(local_size, local_size);;
				this->CalculateGeometricStiffnessMatrix(K_geo, rCurrentProcessInfo);

				LocalStiffnessMatrix += K_geo;
			}
		}

		return LocalStiffnessMatrix;
		KRATOS_CATCH("")
	}

	void TrussElement3D2N::CalculateDampingMatrix(MatrixType& rDampingMatrix,
									ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int MatSize = number_of_nodes * dimension;

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

		rDampingMatrix += alpha * MassMatrix;
		rDampingMatrix += beta  * StiffnessMatrix;

		KRATOS_CATCH("")
	}

	void TrussElement3D2N::CalculateMassMatrix(MatrixType& rMassMatrix,
									ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int MatSize = number_of_nodes * dimension;

		if (rMassMatrix.size1() != MatSize) {
			rMassMatrix.resize(MatSize, MatSize, false);
		}

		rMassMatrix = ZeroMatrix(MatSize, MatSize);

		const double A = this->GetProperties()[CROSS_AREA];
		const double L = this->CalculateReferenceLength();
		const double rho = this->GetProperties()[DENSITY];

		const double TotalMass = A * L * rho;

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
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int MatSize = number_of_nodes * dimension;

		//getting shapefunctionvalues 
		const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues(
			GeometryData::GI_GAUSS_1);

		//creating necessary values 
		const double A = this->GetProperties()[CROSS_AREA];
		const double L = this->CalculateReferenceLength();
		const double rho = this->GetProperties()[DENSITY];

		double TotalMass = A * L * rho;
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
		const unsigned int element_size = number_of_nodes * dimension;

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
		const unsigned int element_size = number_of_nodes * dimension;

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
		const unsigned int element_size = number_of_nodes * dimension;

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

		//calculate internal forces
		VectorType InternalForces = ZeroVector(LocalSize);
		this->UpdateInternalForces(InternalForces);
		//resizing the matrices + create memory for LHS
		rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		//creating LHS
		noalias(rLeftHandSideMatrix) = this->CreateElementStiffnessMatrix(
			rCurrentProcessInfo);

		//create+compute RHS
		rRightHandSideVector = ZeroVector(LocalSize);
		//update Residual
		rRightHandSideVector -= InternalForces;
		//add bodyforces 
		rRightHandSideVector += this->CalculateBodyForces();


		if (this->mIsLinearElement == true)
		{
			Vector NodalDeformation = ZeroVector(LocalSize);
			this->GetValuesVector(NodalDeformation, 0);
			rRightHandSideVector = ZeroVector(LocalSize);
			rRightHandSideVector -= prod(rLeftHandSideMatrix, NodalDeformation);
			rRightHandSideVector += this->CalculateBodyForces();
		}

		if (this->ReturnIfIsCable() == true && this->mIsCompressed == true) {
			rRightHandSideVector = ZeroVector(LocalSize);
		}
		KRATOS_CATCH("")
	}

	void TrussElement3D2N::CalculateRightHandSide(VectorType& rRightHandSideVector,
										ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int LocalSize = NumNodes * dimension;

		rRightHandSideVector = ZeroVector(LocalSize);

		VectorType InternalForces = ZeroVector(LocalSize);
		this->UpdateInternalForces(InternalForces);
		rRightHandSideVector -= InternalForces;


		if (this->mIsLinearElement == true)
		{
			Matrix LeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
			this->CalculateLeftHandSide(LeftHandSideMatrix, rCurrentProcessInfo);
			Vector NodalDeformation = ZeroVector(LocalSize);
			this->GetValuesVector(NodalDeformation);
			rRightHandSideVector = ZeroVector(LocalSize);
			rRightHandSideVector -= prod(LeftHandSideMatrix, NodalDeformation);
		}

		//add bodyforces 
		rRightHandSideVector += this->CalculateBodyForces();

		if (this->ReturnIfIsCable() == true && this->mIsCompressed == true) {
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

		//resizing the matrices + create memory for LHS
		rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);
		//creating LHS
		noalias(rLeftHandSideMatrix) = this->CreateElementStiffnessMatrix(
			 rCurrentProcessInfo);
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
		if (rVariable == TRUSS_PRESTRESS_PK2) {
			rOutput[0] = 0.00;
			if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2) == true) {
				rOutput[0] = this->GetProperties()[TRUSS_PRESTRESS_PK2];
			}
		}
		KRATOS_CATCH("")
	}

	void TrussElement3D2N::CalculateOnIntegrationPoints(
										const Variable<Vector>& rVariable,
										std::vector<Vector>& rOutput,
										const ProcessInfo& rCurrentProcessInfo){
		KRATOS_TRY
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const GeometryType::IntegrationPointsArrayType& integration_points =
			GetGeometry().IntegrationPoints();
		if (rOutput.size() != integration_points.size()) {
			rOutput.resize(integration_points.size());
		}
		if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR)
		{
			Vector Strain = ZeroVector(dimension);
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


	bool TrussElement3D2N::ReturnIfIsCable()
	{
		KRATOS_TRY;
		bool IsCable = false;
		if (this->GetProperties().Has(TRUSS_IS_CABLE) == true) {
			IsCable = this->GetProperties()[TRUSS_IS_CABLE];
		}
		return IsCable;
		KRATOS_CATCH("")
	}

	int  TrussElement3D2N::Check(const ProcessInfo& rCurrentProcessInfo){
		KRATOS_TRY

		if (this->GetGeometry().WorkingSpaceDimension() != 3 || this->GetGeometry().PointsNumber() != 2)
			{
				KRATOS_THROW_ERROR(std::invalid_argument,
					"The truss element works only in 3D and with 2 noded elements", "")
			}
		//verify that the variables are correctly initialized
		if (VELOCITY.Key() == 0)
			KRATOS_ERROR << "VELOCITY has Key zero! (check if the application is correctly registered" << "" << std::endl;
		if (DISPLACEMENT.Key() == 0)
			KRATOS_ERROR << "DISPLACEMENT has Key zero! (check if the application is correctly registered" << "" << std::endl;
		if (ACCELERATION.Key() == 0)
			KRATOS_ERROR << "ACCELERATION has Key zero! (check if the application is correctly registered" << "" << std::endl;
		if (DENSITY.Key() == 0)
			KRATOS_ERROR << "DENSITY has Key zero! (check if the application is correctly registered" << "" << std::endl;
		if (CROSS_AREA.Key() == 0)
			KRATOS_ERROR << "CROSS_AREA has Key zero! (check if the application is correctly registered" << "" << std::endl;
		//verify that the dofs exist
		for (unsigned int i = 0; i<this->GetGeometry().PointsNumber(); ++i)
			{
			if (this->GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false)
				KRATOS_ERROR << "missing variable DISPLACEMENT on node " << this->GetGeometry()[i].Id() << std::endl;
			if (this->GetGeometry()[i].HasDofFor(DISPLACEMENT_X) == false || this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false || this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Z) == false)
				KRATOS_ERROR << "missing one of the dofs for the variable DISPLACEMENT on node " << GetGeometry()[i].Id() << std::endl;
			}


		
		if (this->GetProperties().Has(CROSS_AREA) == false ||
			this->GetProperties()[CROSS_AREA] == 0)
		{
			KRATOS_ERROR <<  "CROSS_AREA not provided for this element" << this->Id() << std::endl;
		}

		if (this->GetProperties().Has(YOUNG_MODULUS) == false ||
			this->GetProperties()[YOUNG_MODULUS] == 0)
		{
			KRATOS_ERROR << "YOUNG_MODULUS not provided for this element" << this->Id() << std::endl;
		}
		if (this->GetProperties().Has(DENSITY) == false)
		{
			KRATOS_ERROR << "DENSITY not provided for this element" << this->Id() << std::endl;
		}


		return 0;

		KRATOS_CATCH("")
	}
	double TrussElement3D2N::CalculateGreenLagrangeStrain(){

		KRATOS_TRY
		const double l = this->CalculateCurrentLength();
		const double L = this->CalculateReferenceLength();
		const double e = ((l * l - L * L) / (2.00 * L * L));
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
		const int NumNodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const int LocalSize = NumNodes * dimension;

		MatrixType TransformationMatrix = ZeroMatrix(LocalSize, LocalSize);
		this->CreateTransformationMatrix(TransformationMatrix);
		const double InternalStrainGL = this->CalculateGreenLagrangeStrain();
		const double l = this->CalculateCurrentLength();
		const double L0 = this->CalculateReferenceLength();
		const double E = this->GetProperties()[YOUNG_MODULUS];
		const double A = this->GetProperties()[CROSS_AREA];

		double S_pre = 0.00;
		if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2) == true) {
			S_pre = this->GetProperties()[TRUSS_PRESTRESS_PK2];
		}

		const double N = ((E*InternalStrainGL + S_pre) * l * A) / L0;

		if (N < 0.00) this->mIsCompressed = true;
		else this->mIsCompressed = false;

		//internal force vectors
		VectorType f_local = ZeroVector(LocalSize);
		f_local[0] = -1.00 * N;
		f_local[3] = 1.00 * N;
		rinternalForces = ZeroVector(LocalSize);
		noalias(rinternalForces) = prod(TransformationMatrix, f_local);
		KRATOS_CATCH("");
	}

	void TrussElement3D2N::CreateTransformationMatrix(Matrix& rRotationMatrix){

		KRATOS_TRY
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int local_size = number_of_nodes * dimension;

		//1st calculate transformation matrix
		Vector DirectionVectorX = ZeroVector(dimension);
		Vector DirectionVectorY = ZeroVector(dimension);
		Vector DirectionVectorZ = ZeroVector(dimension);
		Vector ReferenceCoordinates = ZeroVector(local_size);
		Vector GlobalZ = ZeroVector(dimension);
		GlobalZ[2] = 1.0;

		ReferenceCoordinates[0] = this->GetGeometry()[0].X();
		ReferenceCoordinates[1] = this->GetGeometry()[0].Y();
		ReferenceCoordinates[2] = this->GetGeometry()[0].Z();
		ReferenceCoordinates[3] = this->GetGeometry()[1].X();
		ReferenceCoordinates[4] = this->GetGeometry()[1].Y();
		ReferenceCoordinates[5] = this->GetGeometry()[1].Z();

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
		for (unsigned int kk = 0; kk < local_size; kk += dimension)
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



	void TrussElement3D2N::AddExplicitContribution(const VectorType& rRHSVector,
		const Variable<VectorType>& rRHSVariable,
		Variable<array_1d<double, 3> >& rDestinationVariable,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		
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


	void TrussElement3D2N::CalculateGeometricStiffnessMatrix(MatrixType& rGeometricStiffnessMatrix,
		ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int local_size = number_of_nodes * dimension;

		const double E = this->GetProperties()[YOUNG_MODULUS];
		double A = this->GetProperties()[CROSS_AREA];

		double S_pre = 0.00;
		if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2) == true) {
			S_pre = this->GetProperties()[TRUSS_PRESTRESS_PK2];
		}

		rGeometricStiffnessMatrix = ZeroMatrix(local_size, local_size);

		// du... delta displacement in x-direction
		// dv... delta displacement in y-direction
		// dw... delta displacement in z-direction
		// L... inital member length
		// l... deformed member length
		// e_gl... green_lagrange strain

		double du = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
		double dv = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
		double dw = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
		const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
		const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
		const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
		const double L = this->CalculateReferenceLength();
		const double l = this->CalculateCurrentLength();
		double e_gL = (l*l - L*L) / (2.00 * L*L);
		const double L3 = L * L * L;


		const double K_sigma = ((E*A*e_gL) / L) + ((S_pre*A) / L);
		const double K_uij = (E*A) / L3;

		rGeometricStiffnessMatrix(0, 0) = K_sigma + K_uij * (2 * du*dx + du*du);
		rGeometricStiffnessMatrix(3, 3) = rGeometricStiffnessMatrix(0, 0);

		rGeometricStiffnessMatrix(1, 1) = K_sigma + K_uij * (2 * dv*dy + dv*dv);
		rGeometricStiffnessMatrix(4, 4) = rGeometricStiffnessMatrix(1, 1);

		rGeometricStiffnessMatrix(2, 2) = K_sigma + K_uij * (2 * dw*dz + dw*dw);
		rGeometricStiffnessMatrix(5, 5) = rGeometricStiffnessMatrix(2, 2);

		rGeometricStiffnessMatrix(0, 1) = K_uij * (dx*dv + dy*du + du*dv);
		rGeometricStiffnessMatrix(1, 0) = rGeometricStiffnessMatrix(0, 1);

		rGeometricStiffnessMatrix(0, 2) = K_uij * (dx*dw + dz*du + du*dw);
		rGeometricStiffnessMatrix(2, 0) = rGeometricStiffnessMatrix(0, 2);

		rGeometricStiffnessMatrix(0, 3) = -rGeometricStiffnessMatrix(0, 0);
		rGeometricStiffnessMatrix(3, 0) = rGeometricStiffnessMatrix(0, 3);

		rGeometricStiffnessMatrix(0, 4) = -rGeometricStiffnessMatrix(0, 1);
		rGeometricStiffnessMatrix(4, 0) = rGeometricStiffnessMatrix(0, 4);

		rGeometricStiffnessMatrix(0, 5) = -rGeometricStiffnessMatrix(0, 2);
		rGeometricStiffnessMatrix(5, 0) = rGeometricStiffnessMatrix(0, 5);

		rGeometricStiffnessMatrix(1, 2) = K_uij * (dy*dw + dz*dv + dv*dw);
		rGeometricStiffnessMatrix(2, 1) = rGeometricStiffnessMatrix(1, 2);

		rGeometricStiffnessMatrix(1, 3) = rGeometricStiffnessMatrix(0, 4);
		rGeometricStiffnessMatrix(3, 1) = rGeometricStiffnessMatrix(1, 3);

		rGeometricStiffnessMatrix(1, 4) = -rGeometricStiffnessMatrix(1, 1);
		rGeometricStiffnessMatrix(4, 1) = rGeometricStiffnessMatrix(1, 4);

		rGeometricStiffnessMatrix(1, 5) = -rGeometricStiffnessMatrix(1, 2);
		rGeometricStiffnessMatrix(5, 1) = rGeometricStiffnessMatrix(1, 5);

		rGeometricStiffnessMatrix(2, 3) = -rGeometricStiffnessMatrix(0, 2);
		rGeometricStiffnessMatrix(3, 2) = rGeometricStiffnessMatrix(2, 3);

		rGeometricStiffnessMatrix(2, 4) = -rGeometricStiffnessMatrix(1, 2);
		rGeometricStiffnessMatrix(4, 2) = rGeometricStiffnessMatrix(2, 4);

		rGeometricStiffnessMatrix(2, 5) = -rGeometricStiffnessMatrix(2, 2);
		rGeometricStiffnessMatrix(5, 2) = rGeometricStiffnessMatrix(2, 5);

		rGeometricStiffnessMatrix(3, 4) = rGeometricStiffnessMatrix(0, 1);
		rGeometricStiffnessMatrix(4, 3) = rGeometricStiffnessMatrix(3, 4);

		rGeometricStiffnessMatrix(3, 5) = rGeometricStiffnessMatrix(0, 2);
		rGeometricStiffnessMatrix(5, 3) = rGeometricStiffnessMatrix(3, 5);

		rGeometricStiffnessMatrix(4, 5) = rGeometricStiffnessMatrix(1, 2);
		rGeometricStiffnessMatrix(5, 4) = rGeometricStiffnessMatrix(4, 5);
		KRATOS_CATCH("")
	}

	void TrussElement3D2N::CalculateElasticStiffnessMatrix(MatrixType& rElasticStiffnessMatrix,
		ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		const int number_of_nodes = this->GetGeometry().PointsNumber();
		const int dimension = this->GetGeometry().WorkingSpaceDimension();
		const unsigned int local_size = number_of_nodes * dimension;

		const double E = this->GetProperties()[YOUNG_MODULUS];
		double A = this->GetProperties()[CROSS_AREA];

		rElasticStiffnessMatrix = ZeroMatrix(local_size, local_size);

		const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
		const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
		const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
		const double L = this->CalculateReferenceLength();
		const double L3 = L*L*L;

		const double EA = E*A;

		rElasticStiffnessMatrix(0, 0) = (EA*dx*dx) / L3;
		rElasticStiffnessMatrix(3, 3) = rElasticStiffnessMatrix(0, 0);

		rElasticStiffnessMatrix(1, 1) = (EA*dy*dy) / L3;
		rElasticStiffnessMatrix(4, 4) = rElasticStiffnessMatrix(1, 1);

		rElasticStiffnessMatrix(2, 2) = (EA*dz*dz) / L3;
		rElasticStiffnessMatrix(5, 5) = rElasticStiffnessMatrix(2, 2);

		rElasticStiffnessMatrix(0, 1) = (EA*dx*dy) / L3;
		rElasticStiffnessMatrix(1, 0) = rElasticStiffnessMatrix(0, 1);

		rElasticStiffnessMatrix(0, 2) = (EA*dx*dz) / L3;
		rElasticStiffnessMatrix(2, 0) = rElasticStiffnessMatrix(0, 2);

		rElasticStiffnessMatrix(0, 3) = -rElasticStiffnessMatrix(0, 0);
		rElasticStiffnessMatrix(3, 0) = rElasticStiffnessMatrix(0, 3);

		rElasticStiffnessMatrix(0, 4) = -rElasticStiffnessMatrix(0, 1);
		rElasticStiffnessMatrix(4, 0) = rElasticStiffnessMatrix(0, 4);

		rElasticStiffnessMatrix(0, 5) = -rElasticStiffnessMatrix(0, 2);
		rElasticStiffnessMatrix(5, 0) = rElasticStiffnessMatrix(0, 5);

		rElasticStiffnessMatrix(1, 2) = (EA*dy*dz) / L3;
		rElasticStiffnessMatrix(2, 1) = rElasticStiffnessMatrix(1, 2);

		rElasticStiffnessMatrix(1, 3) = rElasticStiffnessMatrix(0, 4);
		rElasticStiffnessMatrix(3, 1) = rElasticStiffnessMatrix(1, 3);

		rElasticStiffnessMatrix(1, 4) = -rElasticStiffnessMatrix(1, 1);
		rElasticStiffnessMatrix(4, 1) = rElasticStiffnessMatrix(1, 4);

		rElasticStiffnessMatrix(1, 5) = -rElasticStiffnessMatrix(1, 2);
		rElasticStiffnessMatrix(5, 1) = rElasticStiffnessMatrix(1, 5);

		rElasticStiffnessMatrix(2, 3) = -rElasticStiffnessMatrix(0, 2);
		rElasticStiffnessMatrix(3, 2) = rElasticStiffnessMatrix(2, 3);

		rElasticStiffnessMatrix(2, 4) = -rElasticStiffnessMatrix(1, 2);
		rElasticStiffnessMatrix(4, 2) = rElasticStiffnessMatrix(2, 4);

		rElasticStiffnessMatrix(2, 5) = -rElasticStiffnessMatrix(2, 2);
		rElasticStiffnessMatrix(5, 2) = rElasticStiffnessMatrix(2, 5);

		rElasticStiffnessMatrix(3, 4) = rElasticStiffnessMatrix(0, 1);
		rElasticStiffnessMatrix(4, 3) = rElasticStiffnessMatrix(3, 4);

		rElasticStiffnessMatrix(3, 5) = rElasticStiffnessMatrix(0, 2);
		rElasticStiffnessMatrix(5, 3) = rElasticStiffnessMatrix(3, 5);

		rElasticStiffnessMatrix(4, 5) = rElasticStiffnessMatrix(1, 2);
		rElasticStiffnessMatrix(5, 4) = rElasticStiffnessMatrix(4, 5);
		KRATOS_CATCH("")
	}


	void TrussElement3D2N::save(Serializer& rSerializer) const
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
		rSerializer.save("mIscompressed", this->mIsCompressed);
		rSerializer.save("LinerEle", this->mIsLinearElement);

	}
	void TrussElement3D2N::load(Serializer& rSerializer)
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
		rSerializer.load("mIscompressed", this->mIsCompressed);
		rSerializer.load("LinerEle", this->mIsLinearElement);
	}
} // namespace Kratos.


