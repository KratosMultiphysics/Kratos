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
#include "spaces/ublas_space.h"
#include "includes/properties.h"
#include "includes/model_part.h"
#include "custom_elements/navier_stokes.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"

namespace Kratos {
	namespace Testing {

		typedef UblasSpace<double, Matrix, Vector> 									 SpaceType;

		typedef typename ModelPart::IndexType										 IndexType;
		typedef typename ModelPart::NodeIterator					          NodeIteratorType;

	    /** Checks the NavierStokes2D3N element.
	     * Checks the LHS and RHS computation using a small perturbation.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(TestElementNavierStokes2D3N, FluidDynamicsApplicationFastSuite)
		{

			ModelPart modelPart("Main");

			modelPart.AddNodalSolutionStepVariable(VELOCITY);
			modelPart.AddNodalSolutionStepVariable(PRESSURE);

			modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
			modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);

			std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};

			Properties::Pointer pElemProp = modelPart.pGetProperties(0);
			pElemProp->SetValue(DENSITY, 1.0);
			pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0);
			Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
			pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

			modelPart.CreateNewElement("NavierStokes2D3N", 1, elemNodes, pElemProp);

			Element::Pointer pElement = modelPart.pGetElement(1);

			// Define the nodal values
			array_1d<double, 3> velocity_values;
			velocity_values[0] = 1.0;
			velocity_values[1] = 1.0;
			velocity_values[2] = 1.0;
			double pressure_value = 1.0;

			// Set the nodal values
			for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node)
			{
				std::cout << *it_node << std::endl;
				it_node->FastGetSolutionStepValue(PRESSURE) = pressure_value;
				it_node->FastGetSolutionStepValue(VELOCITY) = velocity_values;
			}

			// Compute RHS and LHS
			Vector RHS = ZeroVector(9);
			Matrix LHS = ZeroMatrix(9,9);

			std::cout << "Before computing LHS(u) and RHS(u)..." << std::endl;

			pElement->CalculateLocalSystem(LHS, RHS, modelPart.GetProcessInfo());

			std::cout << "Computed LHS(u) and RHS(u)" << std::endl;

			// Compute the error of the perturbation
			double perturbation = 1e-3;
			std::vector<double> error_norms;

			for (unsigned int i=1; i<3; ++i)
			{
				perturbation /= i;

				// Set the perturbed nodal values
				Vector perturbation_vector = ZeroVector(9);
				for (unsigned int i=0; i<perturbation_vector.size(); ++i)
				{
					perturbation_vector[i] = perturbation;
				}

				array_1d<double, 3> velocity_values_perturbation;
				velocity_values[0] = 1.0+perturbation;
				velocity_values[1] = 1.0+perturbation;
				velocity_values[2] = 1.0+perturbation;
				double pressure_value_perturbation = 1.0+perturbation;

				for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node)
				{
					it_node->FastGetSolutionStepValue(PRESSURE) = pressure_value_perturbation;
					it_node->FastGetSolutionStepValue(VELOCITY) = velocity_values_perturbation;
				}

				// Compute perturbed RHS and LHS
				Vector error = ZeroVector(9);
				Vector RHS_obtained = ZeroVector(9);
				Vector RHS_perturbed = ZeroVector(9);
				Vector solution_increment = ZeroVector(9);

				pElement->CalculateRightHandSide(RHS_perturbed, modelPart.GetProcessInfo());

				SpaceType::Mult(LHS, perturbation_vector, solution_increment);
				noalias(RHS_obtained) = RHS + solution_increment;
				noalias(error) = RHS_perturbed - RHS_obtained;
				error_norms.push_back(SpaceType::TwoNorm(error));
			}

	    }

	    /** Checks the NavierStokes3D4N element.
	     * Checks the LHS and RHS computation using a small perturbation.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(TestElementNavierStokes3D4N, FluidDynamicsApplicationFastSuite)
		{

			ModelPart modelPart("Main");

			modelPart.AddNodalSolutionStepVariable(VELOCITY);
			modelPart.AddNodalSolutionStepVariable(PRESSURE);

			modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
			modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
			modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);

			std::vector<ModelPart::IndexType> elemNodes {1, 2, 3, 4};

			Properties::Pointer pElemProp;

			modelPart.CreateNewElement("NavierStokes3D4N", 1, elemNodes, pElemProp);

			Element::Pointer pElement = modelPart.pGetElement(1);
			// Define the nodal values
			array_1d<double, 3> velocity_values;
			velocity_values[0] = 1.0;
			velocity_values[1] = 1.0;
			velocity_values[2] = 1.0;
			double pressure_value = 1.0;

			// Set the nodal values
			for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node)
			{
				it_node->FastGetSolutionStepValue(PRESSURE) = pressure_value;
				it_node->FastGetSolutionStepValue(VELOCITY) = velocity_values;
			}

			// Compute RHS and LHS
			Vector RHS = ZeroVector(12);
			Matrix LHS = ZeroMatrix(12,12);

			pElement->CalculateLocalSystem(LHS, RHS, modelPart.GetProcessInfo());

			// Compute the error of the perturbation
			double perturbation = 1e-3;
			std::vector<double> error_norms;

			for (unsigned int i=1; i<3; ++i)
			{
				perturbation /= i;

				// Set the perturbed nodal values
				Vector perturbation_vector = ZeroVector(12);
				for (unsigned int i=0; i<perturbation_vector.size(); ++i)
				{
					perturbation_vector[i] = perturbation;
				}

				array_1d<double, 3> velocity_values_perturbation;
				velocity_values[0] = 1.0+perturbation;
				velocity_values[1] = 1.0+perturbation;
				velocity_values[2] = 1.0+perturbation;
				double pressure_value_perturbation = 1.0+perturbation;

				for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node)
				{
					it_node->FastGetSolutionStepValue(PRESSURE) = pressure_value_perturbation;
					it_node->FastGetSolutionStepValue(VELOCITY) = velocity_values_perturbation;
				}

				// Compute perturbed RHS and LHS
				Vector error = ZeroVector(12);
				Vector RHS_obtained = ZeroVector(12);
				Vector RHS_perturbed = ZeroVector(12);
				Vector solution_increment = ZeroVector(12);

				pElement->CalculateRightHandSide(RHS_perturbed, modelPart.GetProcessInfo());

				SpaceType::Mult(LHS, perturbation_vector, solution_increment);
				noalias(RHS_obtained) = RHS + solution_increment;
				noalias(error) = RHS_perturbed - RHS_obtained;
				error_norms.push_back(SpaceType::TwoNorm(error));
			}

	    }

	} // namespace Testing
}  // namespace Kratos.
