// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//  license: 	 structural_mechanics_application/license.txt
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

	void CrBeamElement2D2N::CalculateMassMatrix(MatrixType& rMassMatrix,
		ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		if (rMassMatrix.size1() != msElementSize) {
			rMassMatrix.resize(msElementSize, msElementSize, false);
		}
		rMassMatrix = ZeroMatrix(msElementSize, msElementSize);

		const double L = this->CalculateCurrentLength();
		const double A = this->GetProperties()[CROSS_AREA];
		const double rho = this->GetProperties()[DENSITY];

		const double pre_beam = (rho * A * L) / 420.00;
		const double pre_bar = (rho * A * L) / 6.00;

		// bar part
		rMassMatrix(0, 0) = 2.00 * pre_bar;
		rMassMatrix(0, 3) = 1.00 * pre_bar;
		rMassMatrix(3, 0) = 1.00 * pre_bar;
		rMassMatrix(3, 3) = 2.00 * pre_bar;

		// beam part

		rMassMatrix(1, 1) = pre_beam * 156.00;
		rMassMatrix(1, 2) = pre_beam * 22.00 * L;
		rMassMatrix(1, 4) = pre_beam * 54.00;
		rMassMatrix(1, 5) = pre_beam * (-13.00) * L;

		rMassMatrix(2, 1) = pre_beam * 22.00 * L;
		rMassMatrix(2, 2) = pre_beam * 4.00 * L * L;
		rMassMatrix(2, 4) = pre_beam * 13.00 * L;
		rMassMatrix(2, 5) = pre_beam * (-3.00) * L * L;

		rMassMatrix(4, 1) = pre_beam * 54.00;
		rMassMatrix(4, 2) = pre_beam * 13.00 * L;
		rMassMatrix(4, 4) = pre_beam * 156.00;
		rMassMatrix(4, 5) = pre_beam * (-22.00) * L;

		rMassMatrix(5, 1) = pre_beam * (-13.00) * L;
		rMassMatrix(5, 2) = pre_beam * (-3.00) * L * L;
		rMassMatrix(5, 4) = pre_beam * (-22.00) * L;
		rMassMatrix(5, 5) = pre_beam * (4.00) * L * L;		
		
		this->GlobalizeMatrix(rMassMatrix);
		
		KRATOS_CATCH("")
	}

	void CrBeamElement2D2N::CalculateDampingMatrix(MatrixType& rDampingMatrix,
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

	void CrBeamElement2D2N::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
		{
			KRATOS_TRY;
			// t 
			this->DeformationForces = this->CalculateInternalStresses_DeformationModes();
			if(this->mIsLinearElement) this->DeformationForces = ZeroVector(msLocalSize);

			// qe
			Vector NodalForces = ZeroVector(msElementSize);
			NodalForces = this->ReturnElementForces_Local();
			// q
			this->GlobalizeVector(NodalForces);
			this->F_int_global = NodalForces;
			// Kt
			this->CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
				
			// residual >>> r = f_ext - f_int
			if(this->mIsLinearElement) this->CalculateRightHandSideLinear(rRightHandSideVector,rLeftHandSideMatrix);
			else
			{
				rRightHandSideVector = ZeroVector(msElementSize);
				rRightHandSideVector -= NodalForces;
			}

			rRightHandSideVector += this->CalculateBodyForces();

			KRATOS_CATCH("")
		}

	void CrBeamElement2D2N::CalculateRightHandSide(VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo)
		{
			KRATOS_TRY;
			if(this->mIsLinearElement) this->CalculateRightHandSideLinear(rRightHandSideVector,this->K_master);
			else
			{
				rRightHandSideVector = ZeroVector(msElementSize);
				rRightHandSideVector -= this->F_int_global;
			}
			rRightHandSideVector += this->CalculateBodyForces();
			KRATOS_CATCH("")
		}

	void CrBeamElement2D2N::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
		ProcessInfo& rCurrentProcessInfo) 
		{
			KRATOS_TRY;
			rLeftHandSideMatrix = this->CreateElementStiffnessMatrix_Total();
			this->GlobalizeMatrix(rLeftHandSideMatrix);
			this->K_master = rLeftHandSideMatrix;
			KRATOS_CATCH("")
		}

	/////////////////////////////////////////////////
	///////////// CUSTOM FUNCTIONS --->>
	/////////////////////////////////////////////////

	bounded_vector<double,CrBeamElement2D2N::msElementSize> CrBeamElement2D2N::CalculateBodyForces()
	{
		KRATOS_TRY
		//getting shapefunctionvalues for linear SF
		const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues(
			GeometryData::GI_GAUSS_1);

		bounded_vector<double,3> EquivalentLineLoad = ZeroVector(3); 
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

	void CrBeamElement2D2N::CalculateAndAddWorkEquivalentNodalForcesLineLoad(
		const bounded_vector<double,3> ForceInput,
		bounded_vector<double,CrBeamElement2D2N::msElementSize>& rRightHandSideVector,
		const double GeometryLength)
	{
		KRATOS_TRY;
		//calculate orthogonal load vector
		const double numerical_limit = std::numeric_limits<double>::epsilon();
		Vector GeometricOrientation = ZeroVector(3);
		GeometricOrientation[0] = this->GetGeometry()[1].X()
			- this->GetGeometry()[0].X();
		GeometricOrientation[1] = this->GetGeometry()[1].Y()
			- this->GetGeometry()[0].Y();
		GeometricOrientation[2] = 0.000;


		const double VectorNormA = MathUtils<double>::Norm(GeometricOrientation);
		if (VectorNormA > numerical_limit) GeometricOrientation /= VectorNormA;

		Vector LineLoadDir = ZeroVector(3);
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

		const double sinAngle = std::sqrt(1.00 - (cosAngle*cosAngle));
		const double NormForceVectorOrth = sinAngle * VectorNormB;


		Vector NodeA = ZeroVector(3);
		NodeA[0] = this->GetGeometry()[0].X();
		NodeA[1] = this->GetGeometry()[0].Y();
		NodeA[2] = 0.00;


		Vector NodeB = ZeroVector(3);
		NodeB = NodeA + LineLoadDir;

		Vector NodeC = ZeroVector(3);
		NodeC = NodeA + (GeometricOrientation*cosAngle);

		Vector LoadOrthogonalDir = ZeroVector(3);
		LoadOrthogonalDir = NodeB - NodeC;
		const double VectorNormC = MathUtils<double>::Norm(LoadOrthogonalDir);
		if (VectorNormC > numerical_limit) LoadOrthogonalDir /= VectorNormC;



		// now caluclate respective work equivilent nodal moments

		const double CustomMoment = NormForceVectorOrth *
			GeometryLength*GeometryLength / 12.00;

		Vector MomentNodeA = ZeroVector(3);
		MomentNodeA = MathUtils<double>::CrossProduct(GeometricOrientation,
			LoadOrthogonalDir);
		MomentNodeA *= CustomMoment;


		rRightHandSideVector[msDimension] += MomentNodeA[2];
		rRightHandSideVector[(2*msDimension)+1] -= MomentNodeA[2];

		KRATOS_CATCH("")
	}
	void CrBeamElement2D2N::CalculateRightHandSideLinear(
		VectorType& rRightHandSideVector, MatrixType rLeftHandSideMatrix)
		{
			KRATOS_TRY;
			Vector NodalDeformation = ZeroVector(msElementSize);
			this->GetValuesVector(NodalDeformation);
			rRightHandSideVector = ZeroVector(msElementSize);
			rRightHandSideVector -= prod(rLeftHandSideMatrix, NodalDeformation);
			KRATOS_CATCH("")
		}




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

		const double dx = this->GetGeometry()[1].X0()-this->GetGeometry()[0].X0();
		const double dy = this->GetGeometry()[1].Y0()-this->GetGeometry()[0].Y0();

		const double norm = std::sqrt((dx*dx) + (dy*dy));

		double phi;
		if ((dx > numerical_limit) & (std::abs(dy) < numerical_limit)) phi = 0.00; // dy = 0 and dx > 0
		else if ((dx < -numerical_limit) & (std::abs(dy) < numerical_limit)) phi = Globals::Pi; // dy = 0 and dx < 0
		else if (std::abs(dx) < numerical_limit)
		{
			phi = Globals::Pi / 2.00; // dy > 0 and dx = 0
			if (dy < -numerical_limit) phi = 1.500 * Globals::Pi;  // dy < 0 and dx = 0
		}
		else
		{
			phi = (norm - dx) / dy;
			phi = std::atan(phi);
			phi = 2.00 * phi;
		}

		return phi;
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

		double phi;
		if ((dx > numerical_limit) & (std::abs(dy) < numerical_limit)) phi = 0.00; // dy = 0 and dx > 0
		else if ((dx < -numerical_limit) & (std::abs(dy) < numerical_limit)) phi = Globals::Pi; // dy = 0 and dx < 0
		else if (std::abs(dx) < numerical_limit)
		{
			phi = Globals::Pi / 2.00; // dy > 0 and dx = 0
			if (dy < -numerical_limit) phi = 1.500 * Globals::Pi;  // dy < 0 and dx = 0
		}
		else
		{
			phi = (norm - dx) / dy;
			phi = std::atan(phi);
			phi = 2.00 * phi;
		}

		return phi;
		KRATOS_CATCH("")
	}

	double CrBeamElement2D2N::CalculateCurrentLength() 
	{
		KRATOS_TRY;
		const double numerical_limit = std::numeric_limits<double>::epsilon();
		const double du = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
		const double dv = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y)
			- this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);

		const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
		const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();

		const double l = std::sqrt((du + dx)*(du + dx) + (dv + dy)*(dv + dy));

		KRATOS_ERROR_IF(l<numerical_limit) << "length 0 for element " << this->Id() << std::endl;
		return l;
		KRATOS_CATCH("")
	}

	double CrBeamElement2D2N::CalculateReferenceLength() {
		KRATOS_TRY;
		const double numerical_limit = std::numeric_limits<double>::epsilon();
		const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
		const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
		const double L = std::sqrt((dx*dx) + (dy*dy));

		KRATOS_ERROR_IF(L<numerical_limit) << "length 0 for element " << this->Id() << std::endl;
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
		Kd(2, 2) = 3.00 * Psi * E * Iz / L;
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
		//calculate v

		Vector CurrentDisplacement = ZeroVector(msElementSize);
		this->GetValuesVector(CurrentDisplacement,0);

		bounded_vector<double,msLocalSize> DeformationParameters = ZeroVector(msLocalSize);
		DeformationParameters[0] = this->CalculateCurrentLength() - this->CalculateReferenceLength();
		DeformationParameters[1] = CurrentDisplacement[5] - CurrentDisplacement[2];
		DeformationParameters[2] = CurrentDisplacement[5] + CurrentDisplacement[2];
		DeformationParameters[2] -= 2.00 * (this->CalculateDeformedElementAngle()
		 - this->CalculateInitialElementAngle());

		//calculate modulus 2PI for phi_a
		DeformationParameters[2] = this->Modulus2Pi(DeformationParameters[2] + Globals::Pi) - Globals::Pi;

		return DeformationParameters;
		KRATOS_CATCH("")
	}


	bounded_vector<double,CrBeamElement2D2N::msLocalSize>
	 CrBeamElement2D2N::CalculateInternalStresses_DeformationModes()
	{
		KRATOS_TRY;
		//calculate t

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
		RotationMatrix(1, 1) = c;
		RotationMatrix(2, 2) = 1.00;

		RotationMatrix(3, 3) = c;
		RotationMatrix(3, 4) = -s;
		RotationMatrix(4, 3) = s;
		RotationMatrix(4, 4) = c;
		RotationMatrix(5, 5) = 1.00;



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
		// calculate qe

		bounded_matrix<double,msElementSize,msLocalSize> S = this->CalculateTransformationS();
		bounded_vector<double,msLocalSize> t = this->CalculateInternalStresses_DeformationModes();

		bounded_vector<double,msElementSize> qe = prod(S,t);
		return qe;
		KRATOS_CATCH("")
	 }



	 double CrBeamElement2D2N::Modulus2Pi(double A)
	 {
		KRATOS_TRY;	
		const int B = A / (2.00 * Globals::Pi);
		const double C = A - (B*2.00*Globals::Pi);
		return C;
		KRATOS_CATCH("")
	 }

	int CrBeamElement2D2N::Check(const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		const double numerical_limit = std::numeric_limits<double>::epsilon();
			if (GetGeometry().WorkingSpaceDimension() != 2 || GetGeometry().size() != 2)
			{
				KRATOS_ERROR <<
					"The beam element works only in 2D and with 2 noded elements" << ""
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
				this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false) {
				KRATOS_ERROR <<
					"missing one of the dofs for the variable DISPLACEMENT on node " <<
						GetGeometry()[i].Id() << std::endl;
			}
		}



		if (this->GetProperties().Has(CROSS_AREA) == false ||
			this->GetProperties()[CROSS_AREA] <= numerical_limit)
		{
			KRATOS_ERROR << "CROSS_AREA not provided for this element" << this->Id()
				<< std::endl;
		}

		if (this->GetProperties().Has(YOUNG_MODULUS) == false ||
			this->GetProperties()[YOUNG_MODULUS] <= numerical_limit)
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

		if (this->GetProperties().Has(I33) == false)
		{
			KRATOS_ERROR << "I33 not provided for this element" << this->Id()
				<< std::endl;
		}
		return 0;

		KRATOS_CATCH("")
		}


	void CrBeamElement2D2N::save(Serializer& rSerializer) const
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
	}

	void CrBeamElement2D2N::load(Serializer& rSerializer)
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
	}
} // namespace Kratos.


