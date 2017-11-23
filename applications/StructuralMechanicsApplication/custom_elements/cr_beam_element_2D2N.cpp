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
		if (norm > numerical_limit) dx_0_normalized = dx_0 / norm;
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
		if (norm > numerical_limit) dx_normalized = dx / norm;
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
	CrBeamElement2D2N::msLocalSize> CrBeamElement2D2N::CreateElementStiffnessMatrix_Kd()
	{
		KRATOS_TRY	
		bounded_matrix<double,msLocalSize,msLocalSize> Kd = ZeroMatrix(msLocalSize, msLocalSize);

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
		bounded_matrix<double,msLocalSize,msLocalSize> 
		LocalStiffnessMatrix = ZeroMatrix(msLocalSize, msLocalSize);

		LocalStiffnessMatrix(0, 0) = E*A / L;
		LocalStiffnessMatrix(1, 1) = E*Iz / L;
		LocalStiffnessMatrix(2, 2) = 3.00 * Psi * E * Iz;
		return LocalStiffnessMatrix;
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




} // namespace Kratos.


