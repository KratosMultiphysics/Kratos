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
	//const values
	constexpr int number_of_nodes = 2;
	constexpr int dimension = 3;
	constexpr unsigned int local_size = number_of_nodes * dimension;

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
						GeometryType::Pointer pGeometry,
						bool rLinear = false);
		TrussElement3D2N(IndexType NewId,
						GeometryType::Pointer pGeometry,
						PropertiesType::Pointer pProperties,
						bool rLinear = false);


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

		bounded_matrix<double,local_size,local_size> CreateElementStiffnessMatrix(ProcessInfo& rCurrentProcessInfo);

		void CalculateOnIntegrationPoints(
			const Variable<double>& rVariable,
			std::vector<double>& rOutput,
			const ProcessInfo& rCurrentProcessInfo) override;

		void GetValueOnIntegrationPoints(
			const Variable<double>& rVariable,
			std::vector<double>& rValues,
			const ProcessInfo& rCurrentProcessInfo) override;

		void UpdateInternalForces(bounded_vector<double,local_size>& rinternalForces);
		void CreateTransformationMatrix(bounded_matrix<double,local_size,local_size>& rRotationMatrix);
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

		bounded_vector<double,local_size> CalculateBodyForces();  

		bool ReturnIfIsCable();

		void CalculateGeometricStiffnessMatrix(bounded_matrix<double,local_size,local_size>& rGeometricStiffnessMatrix,
			ProcessInfo& rCurrentProcessInfo);

		void CalculateElasticStiffnessMatrix(bounded_matrix<double,local_size,local_size>& rElasticStiffnessMatrix,
			ProcessInfo& rCurrentProcessInfo);


	private:
		bool mIsCompressed;
		bool mIsLinearElement = false;

		TrussElement3D2N() {};

		friend class Serializer;
		void save(Serializer& rSerializer) const override;
		void load(Serializer& rSerializer) override;
	};


}


#endif
