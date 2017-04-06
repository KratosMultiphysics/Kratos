// ==============================================================================
/*
TRUSS_ELEMENT_3D2N
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
*  Projectname:         $TRUSS_ELEMENT_3D2N
*  Last Modified by:    $Author: klaus.sautter@tum.de $
*  Date:                $Date: April 2017 $
*  Revision:            $Revision: 1.0 $
* ***************************************************************************/

#if !defined(KRATOS_TRUSS_ELEMENT_3D2N_H_INCLUDED )
#define  KRATOS_TRUSS_ELEMENT_3D2N_H_INCLUDED


#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos
{

	class TrussElement3D2N : public Element
	{
	public:
		KRATOS_CLASS_POINTER_DEFINITION(TrussElement3D2N);


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


		TrussElement3D2N(IndexType NewId, 
						GeometryType::Pointer pGeometry);
		TrussElement3D2N(IndexType NewId,
						GeometryType::Pointer pGeometry,
						PropertiesType::Pointer pProperties);


		virtual ~TrussElement3D2N();


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

		MatrixType CreateElementStiffnessMatrix();

		void CalculateOnIntegrationPoints(
			const Variable<double>& rVariable,
			std::vector<double>& rOutput,
			const ProcessInfo& rCurrentProcessInfo) override;

		void GetValueOnIntegrationPoints(
			const Variable<double>& rVariable,
			std::vector<double>& rValues,
			const ProcessInfo& rCurrentProcessInfo) override;

		void UpdateInternalForces(
			VectorType& rinternalForces);
		void CreateTransformationMatrix(
			Matrix& rRotationMatrix);
		double CalculateCurrentLength();

		void CalculateOnIntegrationPoints(
			const Variable<Vector>& rVariable,
			std::vector<Vector>& rOutput,
			const ProcessInfo& rCurrentProcessInfo) override;

		void GetValueOnIntegrationPoints(
			const Variable<Vector>& rVariable,
			std::vector<Vector>& rValues,
			const ProcessInfo& rCurrentProcessInfo) override;

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

		int  Check(
			const ProcessInfo& rCurrentProcessInfo);


		double CalculateGreenLagrangeStrain();
		double CalculateReferenceLength();

		VectorType CalculateBodyForces();  


	private:
		double mPreStress, mArea, mYoungsModulus, mLength, mDensity;
		double mCurrentLength;
		MatrixType mLHS;
		int mIterCount = 0; 
		bool mIsCable;
		bool mIsCompressed;

		TrussElement3D2N() {};

		friend class Serializer;
	};


}


#endif
