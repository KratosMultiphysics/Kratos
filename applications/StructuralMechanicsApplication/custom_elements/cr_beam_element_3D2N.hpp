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
#if !defined(KRATOS_CR_BEAM_ELEMENT_3D2N_H_INCLUDED )
#define  KRATOS_CR_BEAM_ELEMENT_3D2N_H_INCLUDED


#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos
{

	class CrBeamElement3D2N : public Element
	{
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
		typedef unsigned int uint;


		CrBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry);
		CrBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);


		virtual ~CrBeamElement3D2N();


		BaseType::Pointer Create(
			IndexType NewId,
			NodesArrayType const& rThisNodes,
			PropertiesType::Pointer pProperties) const;

		void EquationIdVector(
			EquationIdVectorType& rResult,
			ProcessInfo& rCurrentProcessInfo);

		void GetDofList(
			DofsVectorType& rElementalDofList,
			ProcessInfo& rCurrentProcessInfo);

		void Initialize();

		MatrixType CreateElementStiffnessMatrix_Material();
		MatrixType CreateElementStiffnessMatrix_Geometry(VectorType qe);
		MatrixType CalculateMaterialStiffness();
		MatrixType CalculateTransformationS();

		VectorType CalculateElementForces();

		void CalculateTransformationMatrix(
			Matrix& rRotationMatrix);

		//'Local' -> 3x3
		void CalculateInitialLocalCS();
		MatrixType UpdateRotationMatrixLocal();

		void CalculateLocalSystem(
			MatrixType& rLeftHandSideMatrix,
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo);

		void CalculateRightHandSide(
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo);

		void CalculateLeftHandSide(
			MatrixType& rLeftHandSideMatrix,
			ProcessInfo& rCurrentProcessInfo);

		void CalculateMassMatrix(
			MatrixType& rMassMatrix,
			ProcessInfo& rCurrentProcessInfo);

		void CalculateDampingMatrix(
			MatrixType& rDampingMatrix,
			ProcessInfo& rCurrentProcessInfo);

		void GetValuesVector(
			Vector& rValues,
			int Step = 0);

		void GetSecondDerivativesVector(
			Vector& rValues,
			int Step = 0);

		void GetFirstDerivativesVector(
			Vector& rValues,
			int Step = 0);

		void CalculateSecondDerivativesContributions(
			MatrixType& rLeftHandSideMatrix,
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo);

		void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
		void CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

		int  Check(const ProcessInfo& rCurrentProcessInfo);


		double CalculateGreenLagrangeStrain();
		double CalculateCurrentLength();
		double CalculatePsi(double I, double A_eff);
		double CalculateReferenceLength();
		void UpdateIncrementDeformation();

		VectorType CalculateBodyForces();  

		void CalculateOnIntegrationPoints(
			const Variable<array_1d<double, 3 > >& rVariable,
			std::vector< array_1d<double, 3 > >& rOutput,
			const ProcessInfo& rCurrentProcessInfo);

		void GetValueOnIntegrationPoints(
			const Variable<array_1d<double, 3 > >& rVariable,
			std::vector< array_1d<double, 3 > >& rOutput,
			const ProcessInfo& rCurrentProcessInfo);


	private:
		double mArea, mYoungsModulus, mLength, mDensity, mInertiaX, mInertiaY, mInertiaZ;
		double mEffAreaY, mEffAreaZ, mShearModulus, mPoisson, mtheta, mCurrentLength, mPsiY, mPsiZ;
		double mdPhi_x_a;
		VectorType mNX, mNY, mNZ, mRHS, mTotalDef, mTotalPos, mTotalNodalDeformation, mTotalNodalPosistion;
		VectorType mDeformationModes, mIncrementDeformation, mNodalForces, mBodyForces;
		MatrixType mLHS, mlocalStiffness, mRotationMatrix;
		VectorType mNX0, mNY0, mNZ0;
		VectorType mQuaternionVEC_A, mQuaternionVEC_B;
		double mQuaternionSCA_A, mQuaternionSCA_B;
		VectorType mPhiS, mPhiA;

		int mIterationCount = 0;

		CrBeamElement3D2N() {};

		friend class Serializer;
	};


}


#endif
