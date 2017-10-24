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

#if !defined(KRATOS_CR_BEAM_ELEMENT_3D2N_H_INCLUDED )
#define  KRATOS_CR_BEAM_ELEMENT_3D2N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/serializer.h"

namespace Kratos
{

	class CrBeamElement3D2N : public Element
	{
	private:
		//const values
		static constexpr int msNumberOfNodes = 2;
		static constexpr int msDimension = 3;
		static constexpr unsigned int msLocalSize = msNumberOfNodes * msDimension;
		static constexpr unsigned int msElementSize = msLocalSize * 2;

	public:
		KRATOS_CLASS_POINTER_DEFINITION(CrBeamElement3D2N);


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


		CrBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry,
						bool rLinear = false);
		CrBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry,
						PropertiesType::Pointer pProperties,
						bool rLinear = false);


		~CrBeamElement3D2N() override;


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

		bounded_matrix<double,msElementSize,msElementSize> CreateElementStiffnessMatrix_Material();
		bounded_matrix<double,msElementSize,msElementSize>  CreateElementStiffnessMatrix_Geometry();
		bounded_matrix<double,msLocalSize,msLocalSize> CalculateDeformationStiffness();
		bounded_matrix<double,msElementSize,msLocalSize> CalculateTransformationS();

		bounded_vector<double,msLocalSize> CalculateElementForces();

		void CalculateTransformationMatrix(
			bounded_matrix<double,msElementSize,msElementSize>& rRotationMatrix);

		void CalculateInitialLocalCS();

		bounded_matrix<double,msDimension,msDimension> UpdateRotationMatrixLocal();

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

		void CalculateMassMatrix(
			MatrixType& rMassMatrix,
			ProcessInfo& rCurrentProcessInfo) override;

		void CalculateLumpedMassMatrix(
			MatrixType& rMassMatrix,
			ProcessInfo& rCurrentProcessInfo);

		void CalculateConsistentMassMatrix(
			MatrixType& rMassMatrix,
			ProcessInfo& rCurrentProcessInfo);

		void BuildSingleMassMatrix(
			MatrixType& rMassMatrix,
			const double Phi, const double CT, const double CR, const double L);

		void CalculateDampingMatrix(
			MatrixType& rDampingMatrix,
			ProcessInfo& rCurrentProcessInfo) override;

		void AddExplicitContribution(const VectorType& rRHSVector,
			const Variable<VectorType>& rRHSVariable,
			Variable<array_1d<double, 3> >& rDestinationVariable,
			const ProcessInfo& rCurrentProcessInfo) override;

		void GetValuesVector(
			Vector& rValues,
			int Step = 0) override;

		void GetSecondDerivativesVector(
			Vector& rValues,
			int Step = 0) override;

		void GetFirstDerivativesVector(
			Vector& rValues,
			int Step = 0) override;

		void AssembleSmallInBigMatrix(Matrix SmallMatrix, bounded_matrix<double,
			msElementSize,msElementSize>& BigMatrix);

		int Check(const ProcessInfo& rCurrentProcessInfo) override;


		double CalculateCurrentLength();
		double CalculatePsi(const double I, const double A_eff);
		double CalculateShearModulus();
		double CalculateReferenceLength();
		void UpdateIncrementDeformation();

		bounded_vector<double,msElementSize> CalculateBodyForces();  

		void CalculateOnIntegrationPoints(
			const Variable<array_1d<double, 3 > >& rVariable,
			std::vector< array_1d<double, 3 > >& rOutput,
			const ProcessInfo& rCurrentProcessInfo) override;

		void GetValueOnIntegrationPoints(
			const Variable<array_1d<double, 3 > >& rVariable,
			std::vector< array_1d<double, 3 > >& rOutput,
			const ProcessInfo& rCurrentProcessInfo) override;

		void CalculateOnIntegrationPoints(
			const Variable<Vector >& rVariable,
			std::vector< Vector >& rOutput,
			const ProcessInfo& rCurrentProcessInfo) override;

		void GetValueOnIntegrationPoints(
			const Variable<Vector>& rVariable,
			std::vector<Vector>& rValues,
			const ProcessInfo& rCurrentProcessInfo) override;


		IntegrationMethod GetIntegrationMethod() const override;

		void CalculateAndAddWorkEquivalentNodalForcesLineLoad(
			const bounded_vector<double,msDimension> ForceInput,
			bounded_vector<double,msElementSize>& rRightHandSideVector,
			const double GeometryLength);


		void CalculateGeometricStiffnessMatrix(MatrixType& rGeometricStiffnessMatrix,
			ProcessInfo& rCurrentProcessInfo);

		void CalculateElasticStiffnessMatrix(MatrixType& rElasticStiffnessMatrix,
			ProcessInfo& rCurrentProcessInfo);

	private:
		double mdPhi_x_a, mRotInertiaY, mRotInertiaZ;
		Vector mNX, mNY, mNZ, mRHS, mTotalDef, mTotalPos;
		Vector mTotalNodalDeformation, mTotalNodalPosistion, mBodyForces;
		Vector mDeformationModes, mIncrementDeformation;
		Matrix mLHS, mRotationMatrix;
		bounded_matrix<double,msElementSize,msElementSize> mRotationMatrix0;
		Vector mNX0, mNY0, mNZ0;
		Vector mQuaternionVEC_A, mQuaternionVEC_B;
		double mQuaternionSCA_A, mQuaternionSCA_B;
		Vector mPhiS, mPhiA;
		Vector mNodalForces;

		int mIterationCount = 0;
		bool mIsLinearElement = false;
		bool mIsLumpedMassMatrix = false;

		CrBeamElement3D2N() {};



		friend class Serializer;
		void save(Serializer& rSerializer) const override;
		void load(Serializer& rSerializer) override;
	};


	class Orientation 
	{
	private:
		//const values
		static constexpr int msDimension = 3;
		
	public:
		Orientation(array_1d<double, msDimension>& v1, const double theta = 0.00);
		Orientation(array_1d<double, msDimension>& v1, array_1d<double, msDimension>& v2);


		void CalculateRotationMatrix(Matrix& R);
		void CalculateBasisVectors(array_1d<double, msDimension>& v1,
								   array_1d<double, msDimension>& v2,
								   array_1d<double, msDimension>& v3);

		Quaternion<double>& GetQuaternion() { return morientation; }

	private:
		Quaternion<double> morientation;
	};


}

#endif
