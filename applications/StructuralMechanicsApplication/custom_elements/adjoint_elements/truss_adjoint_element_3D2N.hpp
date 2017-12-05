// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Long Chen
//                   
//                   
//


#if !defined(KRATOS_TRUSS_ADJOINT_ELEMENT_3D2N_H_INCLUDED )
#define  KRATOS_TRUSS_ADJOINT_ELEMENT_3D2N_H_INCLUDED


#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "custom_elements/truss_element_3D2N.hpp"

namespace Kratos
{

	class TrussAdjointElement3D2N : public TrussElement3D2N
	{
	public:
		KRATOS_CLASS_POINTER_DEFINITION(TrussAdjointElement3D2N);


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


        TrussAdjointElement3D2N(IndexType NewId,
						GeometryType::Pointer pGeometry,
						bool rLinear = false);
        TrussAdjointElement3D2N(IndexType NewId,
						GeometryType::Pointer pGeometry,
						PropertiesType::Pointer pProperties,
						bool rLinear = false);


		~TrussAdjointElement3D2N() override;


		BaseType::Pointer Create(
			IndexType NewId,
			NodesArrayType const& rThisNodes,
			PropertiesType::Pointer pProperties) const override;

        BaseType::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
            PropertiesType::Pointer pProperties) const override;


		void EquationIdVector(
			EquationIdVectorType& rResult,
			ProcessInfo& rCurrentProcessInfo) override;

		void GetDofList(
			DofsVectorType& rElementalDofList,
			ProcessInfo& rCurrentProcessInfo) override;

        // adjoint element functions

        void CalculateRightHandSideUnDist(
            VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo);

        void CalculateRightHandSideDist(
            VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo);


        double GetDisturbanceMeasureCorrectionFactor(const Variable<double>& rDesignVariable);

        double GetDisturbanceMeasureCorrectionFactor(const Variable<array_1d<double, 3>>& rDesignVariable);

        void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
            const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rDesignVariable, Matrix& rOutput,
            const ProcessInfo& rCurrentProcessInfo) override;
        void Calculate(const Variable<Vector >& rVariable,
            Vector& rOutput,
            const ProcessInfo& rCurrentProcessInfo) override;

        void Calculate(const Variable<Matrix >& rVariable,
            Matrix& rOutput,
            const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable, Matrix& rOutput,
            const ProcessInfo& rCurrentProcessInfo);

        void CalculateStressDesignVariableDerivative(const Variable<double>& rDesignVariable,
            const Variable<Vector>& rStressVariable, Matrix& rOutput,
            const ProcessInfo& rCurrentProcessInfo);

        void CalculateStressDesignVariableDerivative(const Variable<array_1d<double, 3>>& rDesignVariable,
            const Variable<Vector>& rStressVariable,
            Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo);

        void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rOutput,
            const ProcessInfo& rCurrentProcessInfo) override;

        void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo) override;

        void GetValuesVector(Vector& rValues, int Step = 0);

        int Check(const ProcessInfo& rCurrentProcessInfo) override;



        /*
		void Initialize() override;

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

		int  Check(
			const ProcessInfo& rCurrentProcessInfo) override;


		double CalculateGreenLagrangeStrain();
		double CalculateReferenceLength();

		VectorType CalculateBodyForces();  

		bool ReturnIfIsCable();
        */

    protected:
        TrussAdjointElement3D2N():TrussElement3D2N()
        {}


	private:
		bool mIsCompressed;
		bool mIsLinearElement = false;
        

		friend class Serializer;
		void save(Serializer& rSerializer) const override;
		void load(Serializer& rSerializer) override;
	};


}


#endif
