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

#if !defined(KRATOS_TRUSS_ELEMENT_3D2N_H_INCLUDED )
#define  KRATOS_TRUSS_ELEMENT_3D2N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos
{
	/** 
     * @class TrussElement3D2N
     * 
     * @brief This is a 3D-2node truss element with 3 translational dofs per node
     * 
     * @author Klaus B Sautter
     */

	class TrussElement3D2N : public Element
	{
	protected:
		//const values
		static constexpr int msNumberOfNodes = 2;
		static constexpr int msDimension = 3;
		static constexpr unsigned int msLocalSize = msNumberOfNodes * msDimension;
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


		TrussElement3D2N() {};
		TrussElement3D2N(IndexType NewId, 
						GeometryType::Pointer pGeometry);
		TrussElement3D2N(IndexType NewId,
						GeometryType::Pointer pGeometry,
						PropertiesType::Pointer pProperties);


		~TrussElement3D2N() override;


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

		/**
         * @brief This function calculates the total stiffness matrix for the element
         */
		virtual bounded_matrix<double,msLocalSize,msLocalSize>
		 CreateElementStiffnessMatrix(ProcessInfo& rCurrentProcessInfo);

		void CalculateOnIntegrationPoints(
			const Variable<double>& rVariable,
			std::vector<double>& rOutput,
			const ProcessInfo& rCurrentProcessInfo) override;

		void GetValueOnIntegrationPoints(
			const Variable<double>& rVariable,
			std::vector<double>& rValues,
			const ProcessInfo& rCurrentProcessInfo) override;

		void GetValueOnIntegrationPoints(
			const Variable<array_1d<double, 3 > >& rVariable,
			std::vector< array_1d<double, 3 > >& rOutput,
			const ProcessInfo& rCurrentProcessInfo) override;

        /**
         * @brief This function updates the internal normal force w.r.t. the current deformations
         * @param rinternalForces The current updated internal forces
         */
		virtual void UpdateInternalForces(bounded_vector<double,msLocalSize>& rinternalForces);

		/**
         * @brief This function calculates the transformation matrix to globalize vectors and/or matrices
         * @param rRotationMatrix The transformation matrix
         */
		void CreateTransformationMatrix(bounded_matrix<double,msLocalSize,msLocalSize>& rRotationMatrix);

		void CalculateOnIntegrationPoints(
			const Variable<Vector>& rVariable,
			std::vector<Vector>& rOutput,
			const ProcessInfo& rCurrentProcessInfo) override;

		void GetValueOnIntegrationPoints(
			const Variable<Vector>& rVariable,
			std::vector<Vector>& rValues,
			const ProcessInfo& rCurrentProcessInfo) override;

		void CalculateOnIntegrationPoints(
			const Variable<array_1d<double, 3 > >& rVariable,
			std::vector< array_1d<double, 3 > >& rOutput,
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

		/**
         * @brief This function calculates the current Green-Lagrange strain
         */
		double CalculateGreenLagrangeStrain();

		/**
         * @brief This function calculates the reference length
         */
		double CalculateReferenceLength();

		/**
         * @brief This function calculates the current length
         */
		double CalculateCurrentLength();

		/**
         * @brief This function calculates self-weight forces
         */
		bounded_vector<double,msLocalSize> CalculateBodyForces();  
		
		/**
         * @brief This function assembles the geometric stiffness part of the total stiffness matrix
         * @param rGeometricStiffnessMatrix The geometric stiffness matrix
         * @param rCurrentProcessInfo The current process information
         */
		void CalculateGeometricStiffnessMatrix(bounded_matrix<double,msLocalSize,msLocalSize>& rGeometricStiffnessMatrix,
			ProcessInfo& rCurrentProcessInfo);

		/**
         * @brief This function assembles the elastic stiffness part of the total stiffness matrix
         * @param rElasticStiffnessMatrix The elastic stiffness matrix
         * @param rCurrentProcessInfo The current process information
         */
		void CalculateElasticStiffnessMatrix(bounded_matrix<double,msLocalSize,msLocalSize>& rElasticStiffnessMatrix,
			ProcessInfo& rCurrentProcessInfo);

		/**
         * @brief This function calculates the current nodal postion for the transformation matrix
         * @param rReferenceCoordinates The current coordinates
         */
		virtual void WriteTransformationCoordinates(
			bounded_vector<double,msLocalSize>& rReferenceCoordinates);


	private:
		friend class Serializer;
		void save(Serializer& rSerializer) const override;
		void load(Serializer& rSerializer) override;
	};


}


#endif
