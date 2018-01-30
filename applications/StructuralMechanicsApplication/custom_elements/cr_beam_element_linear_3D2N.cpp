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
#include "custom_elements/cr_beam_element_linear_3D2N.hpp"
#include "structural_mechanics_application_variables.h"
#include "includes/define.h"



namespace Kratos
{

	CrBeamElementLinear3D2N::CrBeamElementLinear3D2N(IndexType NewId,
		GeometryType::Pointer pGeometry)
		: CrBeamElement3D2N(NewId, pGeometry)
	{}

	CrBeamElementLinear3D2N::CrBeamElementLinear3D2N(IndexType NewId,
		GeometryType::Pointer pGeometry,
		PropertiesType::Pointer pProperties)
		: CrBeamElement3D2N(NewId, pGeometry, pProperties)
	{}

	Element::Pointer CrBeamElementLinear3D2N::Create(IndexType NewId,
		NodesArrayType const& rThisNodes,
		PropertiesType::Pointer pProperties) const
	{
		const GeometryType& rGeom = this->GetGeometry();
		return Kratos::make_shared<CrBeamElementLinear3D2N>(
			NewId, rGeom.Create(rThisNodes), pProperties);
	}

	CrBeamElementLinear3D2N::~CrBeamElementLinear3D2N() {}


	void CrBeamElementLinear3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo) {

		KRATOS_TRY
		this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

		Vector NodalDeformation = ZeroVector(msElementSize);
		this->GetValuesVector(NodalDeformation);
		rRightHandSideVector = ZeroVector(msElementSize);
		rRightHandSideVector -= prod(rLeftHandSideMatrix, NodalDeformation);

		//add bodyforces 
		rRightHandSideVector += this->CalculateBodyForces();
		this->IncrementIterationCounter();
		KRATOS_CATCH("")
	}


	void CrBeamElementLinear3D2N::CalculateRightHandSide(
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
		rRightHandSideVector = ZeroVector(msElementSize);

		Matrix LeftHandSideMatrix = ZeroMatrix(msElementSize, msElementSize);
		this->CalculateLeftHandSide(LeftHandSideMatrix, rCurrentProcessInfo);
		Vector NodalDeformation = ZeroVector(msElementSize);
		this->GetValuesVector(NodalDeformation);
		rRightHandSideVector = ZeroVector(msElementSize);
		rRightHandSideVector -= prod(LeftHandSideMatrix, NodalDeformation);

		//add bodyforces 
		rRightHandSideVector += this->CalculateBodyForces();
		KRATOS_CATCH("")

	}

	void CrBeamElementLinear3D2N::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
		ProcessInfo& rCurrentProcessInfo) {

		KRATOS_TRY;
		bounded_matrix<double,msElementSize,msElementSize>
		 TransformationMatrix = this->GetReferenceRotationMatrix();
		rLeftHandSideMatrix = ZeroMatrix(msElementSize, msElementSize);
		rLeftHandSideMatrix +=
			this->CreateElementStiffnessMatrix_Material();
		bounded_matrix<double,msElementSize,msElementSize> aux_matrix = ZeroMatrix(msElementSize);
		aux_matrix = prod(TransformationMatrix, rLeftHandSideMatrix);
		rLeftHandSideMatrix = prod(aux_matrix,
			Matrix(trans(TransformationMatrix)));
		KRATOS_CATCH("")
	}






	bounded_matrix<double,CrBeamElement3D2N::msLocalSize,CrBeamElement3D2N::msLocalSize>  
		CrBeamElementLinear3D2N::CalculateDeformationStiffness() {

		KRATOS_TRY
		bounded_matrix<double,msLocalSize,msLocalSize>
		 Kd = ZeroMatrix(msLocalSize, msLocalSize);
		const double E = this->GetProperties()[YOUNG_MODULUS];
		const double G = this->CalculateShearModulus();
		const double A = this->GetProperties()[CROSS_AREA];
		const double L = this->CalculateReferenceLength();
		
		const double J = this->GetProperties()[TORSIONAL_INERTIA];
		const double Iy = this->GetProperties()[I22];
		const double Iz = this->GetProperties()[I33];

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

		return Kd;
		KRATOS_CATCH("")
	}


	void CrBeamElementLinear3D2N::CalculateOnIntegrationPoints(
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

			Matrix LeftHandSideMatrix = CreateElementStiffnessMatrix_Material();

			Vector NodalDeformation = ZeroVector(msElementSize);
			this->GetValuesVector(NodalDeformation);

			bounded_matrix<double,msElementSize,msElementSize> TransformationMatrix = this->GetReferenceRotationMatrix();
			NodalDeformation = prod(Matrix(trans(TransformationMatrix)),NodalDeformation);

			Vector Stress = prod(LeftHandSideMatrix, NodalDeformation); 


			//rOutput[GP 1,2,3][x,y,z]

			if (rVariable == MOMENT)
			{
				rOutput[0][0] = -1.0 *Stress[3] * 0.75 + Stress[9] * 0.25;
				rOutput[1][0] = -1.0 *Stress[3] * 0.50 + Stress[9] * 0.50;
				rOutput[2][0] = -1.0 *Stress[3] * 0.25 + Stress[9] * 0.75;

				rOutput[0][1] = -1.0 *Stress[4] * 0.75 + Stress[10] * 0.25;
				rOutput[1][1] = -1.0 *Stress[4] * 0.50 + Stress[10] * 0.50;
				rOutput[2][1] = -1.0 *Stress[4] * 0.25 + Stress[10] * 0.75;

				rOutput[0][2] = 1.0 *Stress[5] * 0.75 - Stress[11] * 0.25;
				rOutput[1][2] = 1.0 *Stress[5] * 0.50 - Stress[11] * 0.50;
				rOutput[2][2] = 1.0 *Stress[5] * 0.25 - Stress[11] * 0.75;

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

	void CrBeamElementLinear3D2N::CalculateOnIntegrationPoints(const Variable<Vector >& rVariable,
		std::vector< Vector >& rOutput,
		const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		if (rVariable == LOCAL_AXES_VECTOR)
		{
			rOutput.resize(3);
			for (int i = 0; i < 3; ++i) rOutput[i] = ZeroVector(3);
			rOutput[0] = this->GetNX0();
			rOutput[1] = this->GetNY0();
			rOutput[2] = this->GetNZ0();
		}

		KRATOS_CATCH("");
	}

	void CrBeamElementLinear3D2N::save(Serializer& rSerializer) const
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
	}

	void CrBeamElementLinear3D2N::load(Serializer& rSerializer)
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
	}

} // namespace Kratos.


