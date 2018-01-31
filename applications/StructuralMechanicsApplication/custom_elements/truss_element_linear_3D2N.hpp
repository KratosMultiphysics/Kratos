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

#if !defined(KRATOS_TRUSS_ELEMENT_LINEAR_3D2N_H_INCLUDED )
#define  KRATOS_TRUSS_ELEMENT_LINEAR_3D2N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/truss_element_3D2N.hpp"
#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos
{
	class TrussElementLinear3D2N : public TrussElement3D2N
	{
	public:
		KRATOS_CLASS_POINTER_DEFINITION(TrussElementLinear3D2N);

		TrussElementLinear3D2N() {};
		TrussElementLinear3D2N(IndexType NewId, 
						GeometryType::Pointer pGeometry);
		TrussElementLinear3D2N(IndexType NewId,
						GeometryType::Pointer pGeometry,
						PropertiesType::Pointer pProperties);


		~TrussElementLinear3D2N() override;


		BaseType::Pointer Create(
			IndexType NewId,
			NodesArrayType const& rThisNodes,
			PropertiesType::Pointer pProperties) const override;

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

		void AddPrestressLinear(VectorType& rRightHandSideVector);

		void CalculateOnIntegrationPoints(
			const Variable<array_1d<double, 3 > >& rVariable,
			std::vector< array_1d<double, 3 > >& rOutput,
			const ProcessInfo& rCurrentProcessInfo) override; 

		bounded_matrix<double,msLocalSize,msLocalSize>
		 CreateElementStiffnessMatrix(ProcessInfo& rCurrentProcessInfo) override;


		void WriteTransformationCoordinates(
			bounded_vector<double,msLocalSize>& rReferenceCoordinates) override;
	};


}


#endif
