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

#if !defined(KRATOS_CABLE_ELEMENT_3D2N_H_INCLUDED )
#define  KRATOS_CABLE_ELEMENT_3D2N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/truss_element_3D2N.hpp"
#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos
{
	class CableElement3D2N : public TrussElement3D2N
	{

	public:
		KRATOS_CLASS_POINTER_DEFINITION(CableElement3D2N);



		CableElement3D2N() {};
		CableElement3D2N(IndexType NewId, 
						GeometryType::Pointer pGeometry);
		CableElement3D2N(IndexType NewId,
						GeometryType::Pointer pGeometry,
						PropertiesType::Pointer pProperties);


		~CableElement3D2N() override;


		BaseType::Pointer Create(
			IndexType NewId,
			NodesArrayType const& rThisNodes,
			PropertiesType::Pointer pProperties) const override;

		void CalculateLocalSystem(
			MatrixType& rLeftHandSideMatrix,
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo) override;

		bounded_matrix<double,msLocalSize,msLocalSize>
		 CreateElementStiffnessMatrix(ProcessInfo& rCurrentProcessInfo) override;

		void CalculateRightHandSide(
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo) override;

		void UpdateInternalForces(bounded_vector<double,msLocalSize>& rinternalForces) override;

	private:
		bool mIsCompressed;

		friend class Serializer;
		void save(Serializer& rSerializer) const override;
		void load(Serializer& rSerializer) override;
	};


}


#endif
