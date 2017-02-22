//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/properties.h"
#include "includes/model_part.h"
#include "custom_elements/navier_stokes.h"

namespace Kratos {
	namespace Testing {

		typedef ModelPart::IndexType	IndexType;

	    /** Checks the NavierStokes2D3N element.
	     * Checks the LHS and RHS computation using a small perturbation.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(TestElementNavierStokes2D3N, FluidDynamicsApplicationFastSuite)
		{

			ModelPart modelPart("Main");

			modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
			modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);

			std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};

			Properties::Pointer pElemProp;

			modelPart.CreateNewElement("NavierStokes2D3N", 1, elemNodes, pElemProp);

			KRATOS_WATCH(modelPart);

	    }

	    /** Checks the NavierStokes3D4N element.
	     * Checks the LHS and RHS computation using a small perturbation.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(TestElementNavierStokes3D4N, FluidDynamicsApplicationFastSuite)
		{

			ModelPart modelPart("Main");

			modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
			modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
			modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);

			std::vector<ModelPart::IndexType> elemNodes {1, 2, 3, 4};

			Properties::Pointer pElemProp;

			modelPart.CreateNewElement("NavierStokes3D4N", 1, elemNodes, pElemProp);

			KRATOS_WATCH(modelPart);

	    }

	} // namespace Testing
}  // namespace Kratos.
