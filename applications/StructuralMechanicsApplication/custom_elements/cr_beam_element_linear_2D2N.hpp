// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//  license: 	 structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//                   
//                   
//

#if !defined(KRATOS_CR_BEAM_ELEMENT_LINEAR_2D2N_H_INCLUDED )
#define  KRATOS_CR_BEAM_ELEMENT_LINEAR_2D2N_H_INCLUDED

// System includes

// External includes


// Project includes
#include "custom_elements/cr_beam_element_2D2N.hpp"
#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/serializer.h"

namespace Kratos
{

	class CrBeamElementLinear2D2N : public CrBeamElement2D2N
	{
	public:
		KRATOS_CLASS_POINTER_DEFINITION(CrBeamElementLinear2D2N);


		CrBeamElementLinear2D2N(IndexType NewId, GeometryType::Pointer pGeometry);
		CrBeamElementLinear2D2N(IndexType NewId, GeometryType::Pointer pGeometry,
						PropertiesType::Pointer pProperties);

		BaseType::Pointer Create(
			IndexType NewId,
			NodesArrayType const& rThisNodes,
			PropertiesType::Pointer pProperties) const override;


		~CrBeamElementLinear2D2N() override;		

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


	/////////////////////////////////////////////////
	///////////// CUSTOM FUNCTIONS --->>
	/////////////////////////////////////////////////
		bounded_matrix<double,msElementSize,msElementSize> CreateRotationMatrix() override;

	private:
		CrBeamElementLinear2D2N() {};
		Matrix K_master = ZeroMatrix(msElementSize,msElementSize);

		
		friend class Serializer;
		void save(Serializer& rSerializer) const override;
		void load(Serializer& rSerializer) override;
	};

}

#endif
