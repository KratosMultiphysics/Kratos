//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"
#include "custom_elements/compressible_potential_flow_element.h"

namespace Kratos {
	namespace Testing {

		typedef ModelPart::IndexType									 IndexType;
		typedef ModelPart::NodeIterator					          NodeIteratorType;

	    /** Checks the CompressiblePotentialFlowElement element.
	     * Checks the LHS and RHS computation.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElement, CompressiblePotentialApplicationFastSuite)
		{

			ModelPart modelPart("Main");

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);		

			// Set the element properties
			Properties::Pointer pElemProp = modelPart.pGetProperties(0);

			// Geometry creation
			modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
			modelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
			std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};
			modelPart.CreateNewElement("CompressiblePotentialFlowElement2D3N", 1, elemNodes, pElemProp);

			Element::Pointer pElement = modelPart.pGetElement(1);

			// Define the nodal values
			Vector potential(3);
			potential(0) = 1.0;
			potential(1) = 2.0;
			potential(2) = 3.0;

			for(unsigned int i=0; i<3; i++){
				pElement->GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE) = potential(i);
			}

			// Compute RHS and LHS
			Vector RHS = ZeroVector(3);
			Matrix LHS = ZeroMatrix(3,3);

			pElement->CalculateLocalSystem(LHS, RHS, modelPart.GetProcessInfo());

			// Check the RHS values (the RHS is computed as the LHS x previous_solution, 
			// hence, it is assumed that if the RHS is correct, the LHS is correct as well)
			KRATOS_CHECK_NEAR(RHS(0), 0.5, 1e-7);
			KRATOS_CHECK_NEAR(RHS(1), 0.0, 1e-7);
			KRATOS_CHECK_NEAR(RHS(2), -0.5, 1e-7);
		}
	} // namespace Testing
}  // namespace Kratos.
