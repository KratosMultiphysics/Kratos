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

#if !defined(KRATOS_CR_BEAM_ELEMENT_2D2N_H_INCLUDED )
#define  KRATOS_CR_BEAM_ELEMENT_2D2N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/serializer.h"

namespace Kratos
{

	class CrBeamElement2D2N : public Element
	{
	private:
		//const values
		static constexpr int msNumberOfNodes = 2;
		static constexpr int msDimension = 2;
		static constexpr unsigned int msLocalSize = 3;
		static constexpr unsigned int msElementSize = msLocalSize * 2;

	public:
		KRATOS_CLASS_POINTER_DEFINITION(CrBeamElement2D2N);


		typedef Element BaseType;
		typedef BaseType::GeometryType GeometryType;
		typedef BaseType::NodesArrayType NodesArrayType;
		typedef BaseType::PropertiesType PropertiesType;
		typedef BaseType::IndexType IndexType;
		typedef BaseType::SizeType SizeType;
		typedef BaseType::MatrixType MatrixType;
		typedef BaseType::VectorType VectorType;
		typedef BaseType::EquationIdVectorType EquationIdVectorType;
		typedef BaseType::DofsVectorType DofsVectorType;


		CrBeamElement2D2N(IndexType NewId, GeometryType::Pointer pGeometry,
						bool rLinear = false);
		CrBeamElement2D2N(IndexType NewId, GeometryType::Pointer pGeometry,
						PropertiesType::Pointer pProperties,
						bool rLinear = false);


		~CrBeamElement2D2N() override;


		BaseType::Pointer Create(
			IndexType NewId,
			NodesArrayType const& rThisNodes,
			PropertiesType::Pointer pProperties) const override;

		void EquationIdVector(
			EquationIdVectorType& rResult,
			ProcessInfo& rCurrentProcessInfo) override;

		void GetDofList(
			DofsVectorType& rElementalDofList,
			ProcessInfo& rCurrentProcessInfo) override;

		void Initialize() override;

		void GetValuesVector(
			Vector& rValues,
			int Step = 0) override;

		void GetSecondDerivativesVector(
			Vector& rValues,
			int Step = 0) override;

		void GetFirstDerivativesVector(
			Vector& rValues,
			int Step = 0) override;

		void CalculateMassMatrix(
			MatrixType& rMassMatrix,
			ProcessInfo& rCurrentProcessInfo) override;

		void CalculateDampingMatrix(
			MatrixType& rDampingMatrix,
			ProcessInfo& rCurrentProcessInfo) override;			

		void CalculateLocalSystem(
			MatrixType& rLeftHandSideMatrix,
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo) override;

		void CalculateRightHandSide(
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo) override;

		void CalculateLeftHandSide(
			MatrixType& rLeftHandSideMatrix,
			ProcessInfo& rCurrentProcessInfo) override;

		int Check(const ProcessInfo& rCurrentProcessInfo) override;

	/////////////////////////////////////////////////
	///////////// CUSTOM FUNCTIONS --->>
	/////////////////////////////////////////////////
		double CalculateShearModulus();
		double CalculatePsi(const double I, const double A_eff);
		double CalculateInitialElementAngle();
		double CalculateDeformedElementAngle();
		
		bounded_vector<double,msElementSize> CalculateBodyForces();  
		void CalculateAndAddWorkEquivalentNodalForcesLineLoad(
			const bounded_vector<double,3> ForceInput,
			bounded_vector<double,msElementSize>& rRightHandSideVector,
			const double GeometryLength);


		bounded_matrix<double,msElementSize,msLocalSize> CalculateTransformationS();
		double CalculateCurrentLength();
		double CalculateReferenceLength();

		bounded_matrix<double,msLocalSize,msLocalSize> CreateElementStiffnessMatrix_Kd_mat();
		bounded_matrix<double,msLocalSize,msLocalSize> CreateElementStiffnessMatrix_Kd_geo();
		bounded_matrix<double,msElementSize,msElementSize> CreateElementStiffnessMatrix_Kr();
		bounded_matrix<double,msElementSize,msElementSize> CreateElementStiffnessMatrix_Total();

		void CalculateRightHandSideLinear(VectorType& rRightHandSideVector, MatrixType rLeftHandSideMatrix);

		void GlobalizeMatrix(Matrix &A);
		void GlobalizeVector(Vector &A);
		double Modulus2Pi(double A);

		bounded_matrix<double,msElementSize,msElementSize> CreateRotationMatrix();

		bounded_vector<double,msLocalSize> CalculateDeformationParameters();
		bounded_vector<double,msLocalSize> CalculateInternalStresses_DeformationModes();
		bounded_vector<double,msElementSize> ReturnElementForces_Local();

	private:

		bool mIsLinearElement = false;
		bounded_vector<double,msLocalSize> DeformationForces = ZeroVector(msLocalSize);
		Vector F_int_global = ZeroVector(msElementSize);
		Matrix K_master = ZeroMatrix(msElementSize,msElementSize);
		CrBeamElement2D2N() {};



		friend class Serializer;
		void save(Serializer& rSerializer) const override;
		void load(Serializer& rSerializer) override;
	};

}

#endif
