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
#include "custom_elements/cable_element_3D2N.hpp"
#include "structural_mechanics_application_variables.h"
#include "includes/define.h"


namespace Kratos
{
	CableElement3D2N::CableElement3D2N(IndexType NewId, 
									GeometryType::Pointer pGeometry)
									: TrussElement3D2N(NewId, pGeometry)
	{}

	CableElement3D2N::CableElement3D2N(IndexType NewId,
									GeometryType::Pointer pGeometry,
									PropertiesType::Pointer pProperties) 
									: TrussElement3D2N(NewId, pGeometry, pProperties) 
	{}

	Element::Pointer CableElement3D2N::Create(IndexType NewId,
									NodesArrayType const& rThisNodes,
									PropertiesType::Pointer pProperties) const
	{
		const GeometryType& rGeom = this->GetGeometry();
		return Kratos::make_shared<CableElement3D2N>(
			NewId, rGeom.Create(rThisNodes),pProperties);
	}

	CableElement3D2N::~CableElement3D2N(){}



	bounded_matrix<double,TrussElement3D2N::msLocalSize,
	TrussElement3D2N::msLocalSize> CableElement3D2N::CreateElementStiffnessMatrix(
		ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		bounded_matrix<double,msLocalSize,msLocalSize>
		 LocalStiffnessMatrix = ZeroMatrix(msLocalSize, msLocalSize);

		if (this->mIsCompressed)
		{
			LocalStiffnessMatrix = ZeroMatrix(msLocalSize, msLocalSize);
		}

		else
		{
			this->CalculateElasticStiffnessMatrix(LocalStiffnessMatrix,
				rCurrentProcessInfo);
			bounded_matrix<double,msLocalSize,msLocalSize> K_geo = ZeroMatrix(msLocalSize, msLocalSize);
			this->CalculateGeometricStiffnessMatrix(K_geo, rCurrentProcessInfo);

			LocalStiffnessMatrix += K_geo;
		}

		return LocalStiffnessMatrix;
		KRATOS_CATCH("")
	}


	void CableElement3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
										VectorType& rRightHandSideVector,
										ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		//calculate internal forces
		bounded_vector<double,msLocalSize> InternalForces = ZeroVector(msLocalSize);
		this->UpdateInternalForces(InternalForces);
		//resizing the matrices + create memory for LHS
		rLeftHandSideMatrix = ZeroMatrix(msLocalSize, msLocalSize);
		//creating LHS
		noalias(rLeftHandSideMatrix) = this->CreateElementStiffnessMatrix(
			rCurrentProcessInfo);

		if (this->mIsCompressed) {
			rRightHandSideVector = ZeroVector(msLocalSize);
		}
		else {
		//create+compute RHS
		rRightHandSideVector = ZeroVector(msLocalSize);
		//update Residual
		rRightHandSideVector -= InternalForces;
		//add bodyforces 
		rRightHandSideVector += this->CalculateBodyForces();
		}

		KRATOS_CATCH("")
	}

	void CableElement3D2N::CalculateRightHandSide(VectorType& rRightHandSideVector,
										ProcessInfo& rCurrentProcessInfo){

		KRATOS_TRY
		bounded_vector<double,msLocalSize> InternalForces = ZeroVector(msLocalSize);
		this->UpdateInternalForces(InternalForces);

		if (this->mIsCompressed) {
			rRightHandSideVector = ZeroVector(msLocalSize);
		}
		else {
			rRightHandSideVector = ZeroVector(msLocalSize);
			rRightHandSideVector -= InternalForces;
			//add bodyforces 
			rRightHandSideVector += this->CalculateBodyForces();
		}

		KRATOS_CATCH("")
	}

	void CableElement3D2N::UpdateInternalForces(bounded_vector<double,
		TrussElement3D2N::msLocalSize>& rinternalForces){

		KRATOS_TRY
		bounded_matrix<double,msLocalSize,msLocalSize>
		 TransformationMatrix = ZeroMatrix(msLocalSize, msLocalSize);

		this->CreateTransformationMatrix(TransformationMatrix);
		const double InternalStrainGL = this->CalculateGreenLagrangeStrain();
		const double l = this->CalculateCurrentLength();
		const double L0 = this->CalculateReferenceLength();
		const double E = this->GetProperties()[YOUNG_MODULUS];
		const double A = this->GetProperties()[CROSS_AREA];

		double prestress = 0.00;
		if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
			prestress = this->GetProperties()[TRUSS_PRESTRESS_PK2];
		}

		const double N = ((E*InternalStrainGL + prestress) * l * A) / L0;

		if (N < 0.00) this->mIsCompressed = true;
		else this->mIsCompressed = false;

		//internal force vectors
		bounded_vector<double,msLocalSize> f_local = ZeroVector(msLocalSize);
		f_local[0] = -1.00 * N;
		f_local[3] = 1.00 * N;
		rinternalForces = ZeroVector(msLocalSize);
		noalias(rinternalForces) = prod(TransformationMatrix, f_local);
		KRATOS_CATCH("");
	}

	
	void CableElement3D2N::save(Serializer& rSerializer) const
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
		rSerializer.save("mIscompressed", this->mIsCompressed);

	}
	void CableElement3D2N::load(Serializer& rSerializer)
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
		rSerializer.load("mIscompressed", this->mIsCompressed);
	}
} // namespace Kratos.


