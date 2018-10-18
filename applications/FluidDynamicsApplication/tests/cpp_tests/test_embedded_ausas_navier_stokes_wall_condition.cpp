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
#include "processes/find_nodal_neighbours_process.h"
#include "custom_elements/embedded_ausas_navier_stokes.h"
#include "custom_conditions/embedded_ausas_navier_stokes_wall_condition.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"

namespace Kratos {
	namespace Testing {

		typedef ModelPart::IndexType									 IndexType;
		typedef ModelPart::NodeIterator					          NodeIteratorType;

	    /** Checks the EmbeddedNavierStokes2D3N element.
	     * Checks the LHS and RHS computation.
	     */
		KRATOS_TEST_CASE_IN_SUITE(EmbeddedAusasNavierStokesWallCondition2D3N, FluidDynamicsApplicationFastSuite)
		{

			ModelPart modelPart("Main");
			modelPart.SetBufferSize(3);

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(DISTANCE);
			modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			modelPart.AddNodalSolutionStepVariable(DENSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
			modelPart.AddNodalSolutionStepVariable(SOUND_VELOCITY);
			modelPart.AddNodalSolutionStepVariable(PRESSURE);
			modelPart.AddNodalSolutionStepVariable(VELOCITY);
			modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
			modelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);

			// Process info creation
			double delta_time = 0.1;
			modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
			modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
			modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
			Vector bdf_coefs(3);
			bdf_coefs[0] = 3.0/(2.0*delta_time);
			bdf_coefs[1] = -2.0/delta_time;
			bdf_coefs[2] = 0.5*delta_time;
			modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

			// Set the element properties
			Properties::Pointer pProp = modelPart.pGetProperties(0);
			pProp->SetValue(DENSITY, 1000.0);
			pProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
			Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
			pProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

			// Geometry creation
			modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
			modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
			std::vector<ModelPart::IndexType> condNodes {1, 3};
			std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};
            modelPart.CreateNewElement("EmbeddedAusasNavierStokes2D3N", 1, elemNodes, pProp);
            modelPart.CreateNewCondition("EmbeddedAusasNavierStokesWallCondition2D2N", 1, condNodes, pProp);

			Element::Pointer pElement = modelPart.pGetElement(1);
			Condition::Pointer pCondition = modelPart.pGetCondition(1);

			// Define the nodal values
			Matrix vel_original(3,2);
			vel_original(0,0) = 0.0; vel_original(0,1) = 0.1;
			vel_original(1,0) = 0.1; vel_original(1,1) = 0.2;
			vel_original(2,0) = 0.2; vel_original(2,1) = 0.3;

			// Set the nodal DENSITY and DYNAMIC_VISCOSITY values
			for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
				it_node->FastGetSolutionStepValue(DENSITY) = pProp->GetValue(DENSITY);
				it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pProp->GetValue(DYNAMIC_VISCOSITY);
			}

			Geometry<Node<3>> &r_geometry = pElement->GetGeometry();

			for(unsigned int i=0; i<3; i++){
				r_geometry[i].FastGetSolutionStepValue(PRESSURE)    = 0.0;
				r_geometry[i].FastGetSolutionStepValue(PRESSURE, 1) = 0.0;
				r_geometry[i].FastGetSolutionStepValue(PRESSURE, 2) = 0.0;
				for(unsigned int k=0; k<2; k++){
					r_geometry[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel_original(i,k);
					r_geometry[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i,k);
					r_geometry[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*vel_original(i,k);
					r_geometry[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
					r_geometry[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
					r_geometry[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
				}
			}
			r_geometry[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
			r_geometry[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
			r_geometry[2].FastGetSolutionStepValue(DISTANCE) =  0.5;

			array_1d<double, 3> distances_vector;
			for (unsigned int i = 0; i < r_geometry.size(); ++i) {
				distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
			}
			modelPart.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

			// Find the nodal neighbours
			FindNodalNeighboursProcess find_nodal_neighbours_process(modelPart);
			find_nodal_neighbours_process.Execute();

            // Compute the condition RHS and LHS
            Vector condRHS = ZeroVector(6);
            Matrix condLHS = ZeroMatrix(6,6);

            pCondition->Initialize();
			pCondition->InitializeSolutionStep(modelPart.GetProcessInfo());
			pCondition->CalculateLocalSystem(condLHS, condRHS, modelPart.GetProcessInfo());

			const double tolerance = 1e-10;
			KRATOS_CHECK_NEAR(condRHS(0), 0.0, tolerance);
			KRATOS_CHECK_NEAR(condRHS(1), 0.0, tolerance);
			KRATOS_CHECK_NEAR(condRHS(2), 0.0, tolerance);
			KRATOS_CHECK_NEAR(condRHS(3), 0.0, tolerance);
			KRATOS_CHECK_NEAR(condRHS(4), 0.0, tolerance);
			KRATOS_CHECK_NEAR(condRHS(5), 2.0/30.0, tolerance);
	    }

	    // /** Checks the EmbeddedNavierStokes3D4N element.
	    //  * Checks the LHS and RHS computation.
	    //  */
	    KRATOS_TEST_CASE_IN_SUITE(EmbeddedAusasNavierStokesWallCondition3D4N, FluidDynamicsApplicationFastSuite)
		{

			ModelPart modelPart("Main");
			modelPart.SetBufferSize(3);

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			modelPart.AddNodalSolutionStepVariable(DENSITY);
			modelPart.AddNodalSolutionStepVariable(REACTION);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
			modelPart.AddNodalSolutionStepVariable(SOUND_VELOCITY);
			modelPart.AddNodalSolutionStepVariable(PRESSURE);
			modelPart.AddNodalSolutionStepVariable(VELOCITY);
			modelPart.AddNodalSolutionStepVariable(DISTANCE);
			modelPart.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);
			modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
			modelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);

			// Process info creation
			double delta_time = 0.1;
			modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
			modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
			modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
			Vector bdf_coefs(3);
			bdf_coefs[0] = 3.0/(2.0*delta_time);
			bdf_coefs[1] = -2.0/delta_time;
			bdf_coefs[2] = 0.5*delta_time;
			modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

			// Set the element properties
			Properties::Pointer pProp = modelPart.pGetProperties(0);
			pProp->SetValue(DENSITY, 1000.0);
			pProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
			Newtonian3DLaw::Pointer pConsLaw(new Newtonian3DLaw());
			pProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

			// Geometry creation
			modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			modelPart.CreateNewNode(2, 0.5, 0.25, 0.25);
			modelPart.CreateNewNode(3, 0.25, 0.5, 0.125);
			modelPart.CreateNewNode(4, 0.3, 0.3, 1.0);
			std::vector<ModelPart::IndexType> condNodes {1, 2, 4};
			std::vector<ModelPart::IndexType> elemNodes {1, 2, 3, 4};
			modelPart.CreateNewElement("EmbeddedAusasNavierStokes3D4N", 1, elemNodes, pProp);
			modelPart.CreateNewCondition("EmbeddedAusasNavierStokesWallCondition3D3N", 1, condNodes, pProp);

			Element::Pointer pElement = modelPart.pGetElement(1);
			Condition::Pointer pCondition = modelPart.pGetCondition(1);

			// Define the nodal values
			Matrix vel_original(4,3);
			vel_original(0,0) = 0.0; vel_original(0,1) = 0.1; vel_original(0,2) = 0.2;
			vel_original(1,0) = 0.1; vel_original(1,1) = 0.2; vel_original(1,2) = 0.3;
			vel_original(2,0) = 0.2; vel_original(2,1) = 0.3; vel_original(2,2) = 0.4;
			vel_original(3,0) = 0.3; vel_original(3,1) = 0.4; vel_original(3,2) = 0.5;

			// Set the nodal DENSITY and DYNAMIC_VISCOSITY values
			for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
                it_node->FastGetSolutionStepValue(DENSITY) = pProp->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pProp->GetValue(DYNAMIC_VISCOSITY);
            }

			Geometry<Node<3>> &r_geometry = pElement->GetGeometry();

			for(unsigned int i=0; i<4; i++){
				r_geometry[i].FastGetSolutionStepValue(PRESSURE)    = 0.0;
				r_geometry[i].FastGetSolutionStepValue(PRESSURE, 1) = 0.0;
				r_geometry[i].FastGetSolutionStepValue(PRESSURE, 2) = 0.0;
				for(unsigned int k=0; k<3; k++){
					r_geometry[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel_original(i,k);
					r_geometry[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i,k);
					r_geometry[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*vel_original(i,k);
					r_geometry[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
					r_geometry[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
					r_geometry[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
				}
			}
			r_geometry[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
			r_geometry[1].FastGetSolutionStepValue(DISTANCE) =  1.0;
			r_geometry[2].FastGetSolutionStepValue(DISTANCE) = -1.0;
			r_geometry[3].FastGetSolutionStepValue(DISTANCE) =  1.0;

			array_1d<double, 4> distances_vector;
			for (unsigned int i = 0; i < r_geometry.size(); ++i) {
				distances_vector(i) = r_geometry[i].FastGetSolutionStepValue(DISTANCE);
			}
			modelPart.Elements()[1].SetValue(ELEMENTAL_DISTANCES, distances_vector);

			// Find the nodal neighbours
			FindNodalNeighboursProcess find_nodal_neighbours_process(modelPart);
			find_nodal_neighbours_process.Execute();

			// Compute condition RHS and LHS
            Vector condRHS = ZeroVector(12);
            Matrix condLHS = ZeroMatrix(12,12);

            pCondition->Initialize();
			pCondition->InitializeSolutionStep(modelPart.GetProcessInfo());
			pCondition->CalculateLocalSystem(condLHS, condRHS, modelPart.GetProcessInfo());

			const double tolerance = 1e-7;
			KRATOS_CHECK_NEAR(condRHS(0),  0.0, tolerance);
			KRATOS_CHECK_NEAR(condRHS(1),  0.0, tolerance);
			KRATOS_CHECK_NEAR(condRHS(2),  0.0, tolerance);
			KRATOS_CHECK_NEAR(condRHS(3),  0.0034375, tolerance);
			KRATOS_CHECK_NEAR(condRHS(4),  0.0, tolerance);
			KRATOS_CHECK_NEAR(condRHS(5),  0.0, tolerance);
			KRATOS_CHECK_NEAR(condRHS(6),  0.0, tolerance);
			KRATOS_CHECK_NEAR(condRHS(7),  0.0115625, tolerance);
			KRATOS_CHECK_NEAR(condRHS(8),  0.0, tolerance);
			KRATOS_CHECK_NEAR(condRHS(9),  0.0, tolerance);
			KRATOS_CHECK_NEAR(condRHS(10), 0.0, tolerance);
			KRATOS_CHECK_NEAR(condRHS(11), 0.0111458, tolerance);
		}
	} // namespace Testing
}  // namespace Kratos.
