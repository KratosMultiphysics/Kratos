// ==============================================================================
/*
CR_BEAM_ELEMENT_3D2N
Main author: Klaus B. Sautter
klaus.sautter@tum.de

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//==============================================================================

/* ****************************************************************************
*  Projectname:         $CR_BEAM_ELEMENT_3D2N
*  Last Modified by:    $Author: klaus.sautter@tum.de $
*  Date:                $Date: April 2017 $
*  Revision:            $Revision: 1.0 $
* ***************************************************************************/

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
		CrBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry,
						PropertiesType::Pointer pProperties);


		virtual ~CrBeamElement3D2N();


		BaseType::Pointer Create(
			IndexType NewId,
			NodesArrayType const& rThisNodes,
			PropertiesType::Pointer pProperties) const;

		void EquationIdVector(
			EquationIdVectorType& rResult,
			ProcessInfo& rCurrentProcessInfo) override;

		void GetDofList(
			DofsVectorType& rElementalDofList,
			ProcessInfo& rCurrentProcessInfo) override;

		void Initialize();

		MatrixType CreateElementStiffnessMatrix_Material();
		MatrixType CreateElementStiffnessMatrix_Geometry(VectorType qe);
		MatrixType CalculateMaterialStiffness();
		MatrixType CalculateTransformationS();

		VectorType CalculateElementForces();

		void CalculateTransformationMatrix(
			Matrix& rRotationMatrix);

		void CalculateInitialLocalCS();
		MatrixType UpdateRotationMatrixLocal();

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

		void CalculateDampingMatrix(
			MatrixType& rDampingMatrix,
			ProcessInfo& rCurrentProcessInfo) override;

		void GetValuesVector(
			Vector& rValues,
			int Step = 0) override;

		void GetSecondDerivativesVector(
			Vector& rValues,
			int Step = 0) override;

		void GetFirstDerivativesVector(
			Vector& rValues,
			int Step = 0) override;

		int Check(const ProcessInfo& rCurrentProcessInfo);


		double CalculateGreenLagrangeStrain();
		double CalculateCurrentLength();
		double CalculatePsi(double I, double A_eff);
		double CalculateReferenceLength();
		void UpdateIncrementDeformation();

		VectorType CalculateBodyForces();  

		void CalculateOnIntegrationPoints(
			const Variable<array_1d<double, 3 > >& rVariable,
			std::vector< array_1d<double, 3 > >& rOutput,
			const ProcessInfo& rCurrentProcessInfo) override;

		void GetValueOnIntegrationPoints(
			const Variable<array_1d<double, 3 > >& rVariable,
			std::vector< array_1d<double, 3 > >& rOutput,
			const ProcessInfo& rCurrentProcessInfo) override;


	private:
		double mArea, mYoungsModulus, mLength, mDensity, mInertiaX, mInertiaY;
		double mInertiaZ, mCurrentLength, mPsiY, mPsiZ;
		double mEffAreaY, mEffAreaZ, mShearModulus, mPoisson, mtheta;
		double mdPhi_x_a;
		VectorType mNX, mNY, mNZ, mRHS, mTotalDef, mTotalPos;
		VectorType mTotalNodalDeformation, mTotalNodalPosistion, mBodyForces;
		VectorType mDeformationModes, mIncrementDeformation, mNodalForces;
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
