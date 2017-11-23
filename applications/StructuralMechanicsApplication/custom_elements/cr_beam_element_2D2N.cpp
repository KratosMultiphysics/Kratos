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
#include "custom_elements/cr_beam_element_2D2N.hpp"
#include "structural_mechanics_application_variables.h"
#include "includes/define.h"



namespace Kratos
{

	CrBeamElement2D2N::CrBeamElement2D2N(IndexType NewId,
		GeometryType::Pointer pGeometry, bool rLinear)
		: Element(NewId, pGeometry)
	{
		this->mIsLinearElement = rLinear;
	}

	CrBeamElement2D2N::CrBeamElement2D2N(IndexType NewId,
		GeometryType::Pointer pGeometry,
		PropertiesType::Pointer pProperties, bool rLinear)
		: Element(NewId, pGeometry, pProperties)
	{
		this->mIsLinearElement = rLinear;
	}

	Element::Pointer CrBeamElement2D2N::Create(IndexType NewId,
		NodesArrayType const& rThisNodes,
		PropertiesType::Pointer pProperties) const
	{
		const GeometryType& rGeom = this->GetGeometry();
		return BaseType::Pointer(new CrBeamElement2D2N(
			NewId, rGeom.Create(rThisNodes), pProperties, this->mIsLinearElement));
	}

	CrBeamElement2D2N::~CrBeamElement2D2N() {}

	void CrBeamElement2D2N::EquationIdVector(EquationIdVectorType& rResult,
		ProcessInfo& rCurrentProcessInfo) {
		if (rResult.size() != msElementSize) rResult.resize(msElementSize);

		for (int i = 0; i < msNumberOfNodes; ++i)
		{
			int index = i * msLocalSize;
			rResult[index] = this->GetGeometry()[i].GetDof(DISPLACEMENT_X)
				.EquationId();
			rResult[index + 1] = this->GetGeometry()[i].GetDof(DISPLACEMENT_Y)
				.EquationId();

			rResult[index + 2] = this->GetGeometry()[i].GetDof(ROTATION_Z)
				.EquationId();
		}

	}

	void CrBeamElement2D2N::GetDofList(DofsVectorType& rElementalDofList,
		ProcessInfo& rCurrentProcessInfo) {

		if (rElementalDofList.size() != msElementSize) {
			rElementalDofList.resize(msElementSize);
		}

		for (int i = 0; i < msNumberOfNodes; ++i)
		{
			int index = i * msLocalSize;
			rElementalDofList[index] = this->GetGeometry()[i]
				.pGetDof(DISPLACEMENT_X);
			rElementalDofList[index + 1] = this->GetGeometry()[i]
				.pGetDof(DISPLACEMENT_Y);

			rElementalDofList[index + 2] = this->GetGeometry()[i]
				.pGetDof(ROTATION_Z);

		}
	}

	void CrBeamElement2D2N::Initialize() {

		KRATOS_TRY;
		KRATOS_CATCH("")
	}

	void CrBeamElement2D2N::GetValuesVector(Vector& rValues, int Step) {

		KRATOS_TRY
		if (rValues.size() != msElementSize) rValues.resize(msElementSize, false);

		for (int i = 0; i < msNumberOfNodes; ++i)
		{
			int index = i * msLocalSize;
			rValues[index] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(DISPLACEMENT_X, Step);
			rValues[index + 1] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(DISPLACEMENT_Y, Step);
			rValues[index + 2] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ROTATION_Z, Step);
		}
		KRATOS_CATCH("")
	}

	void CrBeamElement2D2N::GetFirstDerivativesVector(Vector& rValues, int Step)
	{

		KRATOS_TRY
		if (rValues.size() != msElementSize) rValues.resize(msElementSize, false);

		for (int i = 0; i < msNumberOfNodes; ++i)
		{
			int index = i * msLocalSize;
			rValues[index] = this->GetGeometry()[i].
				FastGetSolutionStepValue(VELOCITY_X, Step);
			rValues[index + 1] = this->GetGeometry()[i].
				FastGetSolutionStepValue(VELOCITY_Y, Step);
			rValues[index + 2] = this->GetGeometry()[i].
				FastGetSolutionStepValue(ANGULAR_VELOCITY_Z, Step);
		}

		KRATOS_CATCH("")
	}

	void CrBeamElement2D2N::GetSecondDerivativesVector(Vector& rValues, int Step)
	{

		KRATOS_TRY
		if (rValues.size() != msElementSize) rValues.resize(msElementSize, false);

		for (int i = 0; i < msNumberOfNodes; ++i)
		{
			int index = i * msLocalSize;

			rValues[index] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ACCELERATION_X, Step);
			rValues[index + 1] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ACCELERATION_Y, Step);
			rValues[index + 2] = this->GetGeometry()[i]
				.FastGetSolutionStepValue(ANGULAR_ACCELERATION_Z, Step);
		}
		KRATOS_CATCH("")
	}

	void CrBeamElement2D2N::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
		{
			std::cout << "LocalSystem" << std::endl;
		}

	void CrBeamElement2D2N::CalculateRightHandSide(VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo)
		{
			std::cout << "RHS" << std::endl;
		}

	void CrBeamElement2D2N::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
		ProcessInfo& rCurrentProcessInfo) 
		{
			std::cout << "LHS" << std::endl;
		}

	/////////////////////////////////////////////////
	///////////// CUSTOM FUNCTIONS --->>
	/////////////////////////////////////////////////
	double CrBeamElement2D2N::CalculateShearModulus() {
		KRATOS_TRY;
		const double nu = this->GetProperties()[POISSON_RATIO];
		const double E = this->GetProperties()[YOUNG_MODULUS];
		const double G = E / (2.0 * (1.0 + nu));
		return G;
		KRATOS_CATCH("")
	}

	double CrBeamElement2D2N::CalculatePsi(const double I, const double A_eff) {

		KRATOS_TRY;
		const double E = this->GetProperties()[YOUNG_MODULUS];
		const double L =this->CalculateCurrentLength();
		const double G = this->CalculateShearModulus();

		const double phi = (12.0 * E * I) / (L*L * G*A_eff);
		double psi;
		//interpret input A_eff == 0 as shearstiff -> psi = 1.0
		if (A_eff == 0.00) psi = 1.00;
		else psi = 1.0 / (1.0 + phi);

		return psi;
		KRATOS_CATCH("")
	}

	double CrBeamElement2D2N::CalculateInitialElementAngle()
	{
		KRATOS_TRY;
		const double numerical_limit = std::numeric_limits<double>::epsilon();
		const double dx_0 = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
		const double dy_0 = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();

		const double norm = std::sqrt((dx_0*dx_0) + (dy_0*dy_0));
		double dx_0_normalized = 0.00;
		if (std::abs(dx_0) > numerical_limit) dx_0_normalized = dx_0 / norm;
		const double phi_0 = std::acos(dx_0_normalized);
		return phi_0;
		KRATOS_CATCH("")
	}

	double CrBeamElement2D2N::CalculateDeformedElementAngle()
	{
		KRATOS_TRY;
		const double numerical_limit = std::numeric_limits<double>::epsilon();

		Vector CurrentDisplacement = ZeroVector(msElementSize);
		this->GetValuesVector(CurrentDisplacement,0);


		const double dx = (this->GetGeometry()[1].X0()+CurrentDisplacement[3])
		 - (this->GetGeometry()[0].X0()+CurrentDisplacement[0]);
		const double dy = (this->GetGeometry()[1].Y0()+CurrentDisplacement[4])
		 - (this->GetGeometry()[0].Y0()+CurrentDisplacement[1]);

		const double norm = std::sqrt((dx*dx) + (dy*dy));
		double dx_normalized = 0.00;
		if (std::abs(dx) > numerical_limit) dx_normalized = dx / norm;
		const double phi = std::acos(dx_normalized);
		return phi;
		KRATOS_CATCH("")
	}



	double CrBeamElement2D2N::CalculateCurrentLength() 
	{
		KRATOS_TRY;
		const double du = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
		const double dv = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);

		const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
		const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();

		const double l = std::sqrt((du + dx)*(du + dx) + (dv + dy)*(dv + dy));
		return l;
		KRATOS_CATCH("")
	}

	double CrBeamElement2D2N::CalculateReferenceLength() {
		KRATOS_TRY;
		const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
		const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
		const double L = std::sqrt((dx*dx) + (dy*dy));
		return L;
		KRATOS_CATCH("")
	}

	bounded_matrix<double,CrBeamElement2D2N::msElementSize,
	CrBeamElement2D2N::msLocalSize> CrBeamElement2D2N::CalculateTransformationS() 
	{
		KRATOS_TRY;
		const double L = this->CalculateCurrentLength();
		bounded_matrix<double,msElementSize,msLocalSize> S = ZeroMatrix(msElementSize, msLocalSize);
		S(0, 0) = -1.00;
		S(1, 2) = 2.00 / L;
		S(2, 1) = -1.00;
		S(2, 2) = 1.00;
		S(3, 0) = 1.00;
		S(4, 2) = -2.00/L;
		S(5, 1) = 1.00;
		S(5, 2) = 1.00;
		return S;
		KRATOS_CATCH("")
	}


	bounded_matrix<double,CrBeamElement2D2N::msLocalSize,
	CrBeamElement2D2N::msLocalSize> CrBeamElement2D2N::CreateElementStiffnessMatrix_Kd_mat()
	{
		KRATOS_TRY	
		// element properties
		const double E = this->GetProperties()[YOUNG_MODULUS];
		const double A = this->GetProperties()[CROSS_AREA];
		const double L = this->CalculateCurrentLength();

		const double Iz = this->GetProperties()[I33];

		double Ay = 0.00;
		if (this->GetProperties().Has(AREA_EFFECTIVE_Y)) {
			Ay = GetProperties()[AREA_EFFECTIVE_Y];
		}

		const double Psi = this->CalculatePsi(Iz, Ay);

		// element material stiffness matrix
		bounded_matrix<double,msLocalSize,msLocalSize> Kd = ZeroMatrix(msLocalSize, msLocalSize);

		Kd(0, 0) = E*A / L;
		Kd(1, 1) = E*Iz / L;
		Kd(2, 2) = 3.00 * Psi * E * Iz;
		return Kd;
		KRATOS_CATCH("")
	}

	bounded_matrix<double,CrBeamElement2D2N::msLocalSize,
	CrBeamElement2D2N::msLocalSize> CrBeamElement2D2N::CreateElementStiffnessMatrix_Kd_geo()
	{
		KRATOS_TRY	
		// element properties
		const double L = this->CalculateCurrentLength();
		const double N = this->DeformationForces[0];

		// element material stiffness matrix
		bounded_matrix<double,msLocalSize,msLocalSize> Kd = ZeroMatrix(msLocalSize, msLocalSize);
		
		Kd(1, 1) = N*L / 12.00;
		Kd(2, 2) = N*L / 20.00;
		return Kd;
		KRATOS_CATCH("")
	}	


	bounded_matrix<double,CrBeamElement2D2N::msElementSize,CrBeamElement2D2N::msElementSize>
	 CrBeamElement2D2N::CreateElementStiffnessMatrix_Kr()
	 {
		KRATOS_TRY	
		// element properties
		const double L = this->CalculateCurrentLength();
		const double N = this->DeformationForces[0];
		const double Q = (-2.00 / L) * this->DeformationForces[2];

		// element material stiffness matrix
		bounded_matrix<double,msElementSize,msElementSize> Kr = ZeroMatrix(msElementSize, msElementSize);

		Kr(0, 1) = -Q;
		Kr(0, 4) = Q;
		Kr(1, 0) = -Q;
		Kr(1, 1) = N;
		Kr(1, 3) = Q;
		Kr(1, 4) = -N;

		Kr(3, 1) = Q;
		Kr(3, 4) = -Q;
		Kr(4, 0) = Q;
		Kr(4, 1) = -N;
		Kr(4, 3) = -Q;
		Kr(4, 4) = N;
		return Kr;
		KRATOS_CATCH("")
	 }



	bounded_matrix<double,CrBeamElement2D2N::msElementSize,CrBeamElement2D2N::msElementSize>
	 CrBeamElement2D2N::CreateElementStiffnessMatrix_Total()
	 {
		KRATOS_TRY	
		// co-rotating K
		bounded_matrix<double,msElementSize,msElementSize> K_r = this->CreateElementStiffnessMatrix_Kr();

		// element K (mat+geo)
		bounded_matrix<double,msLocalSize,msLocalSize> K_d_mat = this->CreateElementStiffnessMatrix_Kd_mat();
		bounded_matrix<double,msLocalSize,msLocalSize> K_d_geo = this->CreateElementStiffnessMatrix_Kd_geo();
		bounded_matrix<double,msLocalSize,msLocalSize> K_d = K_d_mat+K_d_geo;

		bounded_matrix<double,msElementSize,msLocalSize> S = this->CalculateTransformationS();
		bounded_matrix<double,msElementSize,msElementSize> K_d_element = prod(K_d,Matrix(trans(S)));
		K_d_element = prod(S,K_d_element);

		// total K
		bounded_matrix<double,msElementSize,msElementSize> K_total = ZeroMatrix(msElementSize, msElementSize);
		K_total += K_r;
		K_total += K_d_element;
		
		return K_total;
		KRATOS_CATCH("")
	 }


	bounded_vector<double,CrBeamElement2D2N::msLocalSize> CrBeamElement2D2N::CalculateDeformationParameters()
	{
		KRATOS_TRY;

		Vector CurrentDisplacement = ZeroVector(msElementSize);
		this->GetValuesVector(CurrentDisplacement,0);

		bounded_vector<double,msLocalSize> DeformationParameters = ZeroVector(msLocalSize);
		DeformationParameters[0] = this->CalculateCurrentLength() - this->CalculateReferenceLength();
		DeformationParameters[1] = CurrentDisplacement[5] - CurrentDisplacement[2];
		DeformationParameters[2] = CurrentDisplacement[5] + CurrentDisplacement[2];
		DeformationParameters[2] -= 2.00 * (this->CalculateDeformedElementAngle()
		 - this->CalculateInitialElementAngle());

		return DeformationParameters;
		KRATOS_CATCH("")
	}


	bounded_vector<double,CrBeamElement2D2N::msLocalSize>
	 CrBeamElement2D2N::CalculateInternalStresses_DeformationModes()
	{
		KRATOS_TRY;

		bounded_vector<double,msLocalSize> DeformationStresses = ZeroVector(msLocalSize);

		bounded_vector<double,msLocalSize> DeformationModes = this->CalculateDeformationParameters();

		bounded_matrix<double,msLocalSize,msLocalSize> K_d_mat = this->CreateElementStiffnessMatrix_Kd_mat();
		bounded_matrix<double,msLocalSize,msLocalSize> K_d_geo = this->CreateElementStiffnessMatrix_Kd_geo();
		bounded_matrix<double,msLocalSize,msLocalSize> K_d = K_d_mat+K_d_geo;



		DeformationStresses = prod(K_d,DeformationModes);

		return DeformationStresses;
		KRATOS_CATCH("")
	}

	bounded_matrix<double,CrBeamElement2D2N::msElementSize,CrBeamElement2D2N::msElementSize>
	 CrBeamElement2D2N::CreateRotationMatrix()
	{
		KRATOS_TRY;
		const double current_element_angle = this->CalculateDeformedElementAngle();
		const double c = std::cos(current_element_angle);
		const double s = std::sin(current_element_angle);

		bounded_matrix<double,msElementSize,msElementSize> RotationMatrix = ZeroMatrix(msElementSize,msElementSize);

		RotationMatrix(0, 0) = c;
		RotationMatrix(0, 1) = -s;
		RotationMatrix(1, 0) = s;
		RotationMatrix(2, 2) = c;

		RotationMatrix(3, 3) = c;
		RotationMatrix(3, 4) = -s;
		RotationMatrix(4, 3) = s;
		RotationMatrix(5, 5) = c;

		return RotationMatrix;
		KRATOS_CATCH("")
	}


	void CrBeamElement2D2N::GlobalizeMatrix(Matrix &A)
	{
		KRATOS_TRY;
		bounded_matrix<double,msElementSize,msElementSize> R = this->CreateRotationMatrix();
		A = prod(A,Matrix(trans(R)));
		A = prod(R,A);
		KRATOS_CATCH("")
	}

	void CrBeamElement2D2N::GlobalizeVector(Vector &A)
	{
		KRATOS_TRY;
		bounded_matrix<double,msElementSize,msElementSize> R = this->CreateRotationMatrix();
		A = prod(R,A);
		KRATOS_CATCH("")
	}
	
	bounded_vector<double,CrBeamElement2D2N::msElementSize>
	 CrBeamElement2D2N::ReturnElementForces_Local()
	 {
		KRATOS_TRY;
		bounded_matrix<double,msElementSize,msLocalSize> S = this->CalculateTransformationS();
		bounded_vector<double,msLocalSize> t = this->CalculateDeformationParameters();

		bounded_vector<double,msElementSize> qe = prod(S,t);
		return qe;
		KRATOS_CATCH("")
	 }


} // namespace Kratos.


