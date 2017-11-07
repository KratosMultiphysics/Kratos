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

#if !defined(KRATOS_GNL_BEAM_LUMPE_3D2N_H_INCLUDED )
#define  KRATOS_GNL_BEAM_LUMPE_3D2N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/serializer.h"

namespace Kratos
{

	class GnlBeamLumpe3D2N : public Element
	{
	private:
		//const values
		static constexpr int msNumberOfNodes = 2;
		static constexpr int msDimension = 3;
		static constexpr unsigned int msLocalSize = msNumberOfNodes * msDimension;
		static constexpr unsigned int msElementSize = msLocalSize * 2;

	public:
		KRATOS_CLASS_POINTER_DEFINITION(GnlBeamLumpe3D2N);


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


		GnlBeamLumpe3D2N(IndexType NewId, GeometryType::Pointer pGeometry);
		GnlBeamLumpe3D2N(IndexType NewId, GeometryType::Pointer pGeometry,
						PropertiesType::Pointer pProperties);


		~GnlBeamLumpe3D2N() override;


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

		void GetValuesVector(
			Vector& rValues,
			int Step = 0) override;

		void GetSecondDerivativesVector(
			Vector& rValues,
			int Step = 0) override;

		void GetFirstDerivativesVector(
			Vector& rValues,
			int Step = 0) override;


		//////////////////// custom functions ////////////////////
		bounded_matrix<double,msLocalSize,msElementSize> Mode_TransformationMatrix();
		bounded_matrix<double,msLocalSize,msLocalSize> Delta_StiffnessMatrix();
		bounded_matrix<double,msElementSize,msElementSize> StiffnessMatrix();
		bounded_matrix<double,msDimension,msDimension> RotationMatrix(bounded_vector<double,msDimension> Psi_k);
		bounded_matrix<double,msDimension,msDimension> RotationMatrix0();
		void AssembleTransformationMatrix(Matrix RotationMatrix, bounded_matrix<double,
			msElementSize,msElementSize>& TransformationMatrix);
		double CalculateReferenceLength();

	private:
		
		double Iy = 0.0001;
		double Iz = 0.0002;
		double It = 0.0003;
		double A = 0.01;
		double E = 210000000000;
		double G = 80800000000;
		double alphay = 1.00;
		double alphaz = 1.00;
		double mu = 1.00;
        
		GnlBeamLumpe3D2N() {};



		friend class Serializer;
	};


}

#endif
