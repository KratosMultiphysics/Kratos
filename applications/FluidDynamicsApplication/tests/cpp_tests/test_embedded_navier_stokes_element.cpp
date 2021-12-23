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
#include "containers/model.h"
#include "includes/properties.h"
#include "includes/model_part.h"
#include "custom_elements/embedded_navier_stokes.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"

namespace Kratos {
	namespace Testing {

		typedef ModelPart::IndexType									 IndexType;
		typedef ModelPart::NodeIterator					          NodeIteratorType;

	    /** Checks the EmbeddedNavierStokes2D3N element.
	     * Checks the LHS and RHS computation.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(ElementEmbeddedNavierStokes2D3N, FluidDynamicsApplicationFastSuite)
		{
			Model model;
			ModelPart& modelPart = model.CreateModelPart("Main");
			modelPart.SetBufferSize(3);

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(DISTANCE);
			modelPart.AddNodalSolutionStepVariable(REACTION);
			modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
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
			Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
			pElemProp->SetValue(DENSITY, 1000.0);
			pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
			Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
			pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

			// Geometry creation
			modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
			modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
			std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};
			modelPart.CreateNewElement("EmbeddedNavierStokes2D3N", 1, elemNodes, pElemProp);

			Element::Pointer pElement = modelPart.pGetElement(1);

			// Define the nodal values
			Matrix vel_original(3,2);
			vel_original(0,0) = 0.0; vel_original(0,1) = 0.1;
			vel_original(1,0) = 0.1; vel_original(1,1) = 0.2;
			vel_original(2,0) = 0.2; vel_original(2,1) = 0.3;

			array_1d<double, 3> embedded_vel;
			embedded_vel(0) = 1.0;
			embedded_vel(1) = 2.0;
			embedded_vel(2) = 0.0;
			for(unsigned int i=0; i<3; i++){
				pElement->GetGeometry()[i].SetValue(EMBEDDED_VELOCITY, embedded_vel);
				pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE)    = 0.0;
				pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE, 1) = 0.0;
				pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE, 2) = 0.0;
				for(unsigned int k=0; k<2; k++){
					pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel_original(i,k);
					pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i,k);
					pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*vel_original(i,k);
					pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
					pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
					pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
				}
			}
			pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
			pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
			pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) =  0.5;

			// Compute RHS and LHS
			Vector RHS = ZeroVector(9);
			Matrix LHS = ZeroMatrix(9,9);

			const auto& r_process_info = modelPart.GetProcessInfo();
			pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
			pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

			// Check the RHS values (the RHS is computed as the LHS x previous_solution,
			// hence, it is assumed that if the RHS is correct, the LHS is correct as well)
			KRATOS_CHECK_NEAR(RHS(0), 0.0475309, 1e-7);
			KRATOS_CHECK_NEAR(RHS(1), 0.0975309, 1e-7);
			KRATOS_CHECK_NEAR(RHS(2), -0.0546391, 1e-7);
			KRATOS_CHECK_NEAR(RHS(3), 0.0469136, 1e-7);
			KRATOS_CHECK_NEAR(RHS(4), 0.0969136, 1e-7);
			KRATOS_CHECK_NEAR(RHS(5), 0.0176796, 1e-7);
			KRATOS_CHECK_NEAR(RHS(6), 28.2303, 1e-1);
			KRATOS_CHECK_NEAR(RHS(7), 46.459, 1e-1);
			KRATOS_CHECK_NEAR(RHS(8), 0.0202928, 1e-7);
		}

	    // /** Checks the EmbeddedNavierStokes3D4N element.
	    //  * Checks the LHS and RHS computation.
	    //  */
	    KRATOS_TEST_CASE_IN_SUITE(ElementEmbeddedNavierStokes3D4N, FluidDynamicsApplicationFastSuite)
		{
			Model model;
			ModelPart& modelPart = model.CreateModelPart("Main", 3);

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
			modelPart.AddNodalSolutionStepVariable(SOUND_VELOCITY);
			modelPart.AddNodalSolutionStepVariable(PRESSURE);
			modelPart.AddNodalSolutionStepVariable(VELOCITY);
			modelPart.AddNodalSolutionStepVariable(DISTANCE);
			modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);

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
			Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
			pElemProp->SetValue(DENSITY, 1000.0);
			pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
			Newtonian3DLaw::Pointer pConsLaw(new Newtonian3DLaw());
			pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

			// Geometry creation
			modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			modelPart.CreateNewNode(2, 0.5, 0.25, 0.25);
			modelPart.CreateNewNode(3, 0.25, 0.5, 0.125);
			modelPart.CreateNewNode(4, 0.3, 0.3, 1.0);
			std::vector<ModelPart::IndexType> elemNodes {1, 2, 3, 4};
			modelPart.CreateNewElement("EmbeddedNavierStokes3D4N", 1, elemNodes, pElemProp);

			Element::Pointer pElement = modelPart.pGetElement(1);

			// Define the nodal values
			Matrix vel_original(4,3);
			vel_original(0,0) = 0.0; vel_original(0,1) = 0.1; vel_original(0,2) = 0.2;
			vel_original(1,0) = 0.1; vel_original(1,1) = 0.2; vel_original(1,2) = 0.3;
			vel_original(2,0) = 0.2; vel_original(2,1) = 0.3; vel_original(2,2) = 0.4;
			vel_original(3,0) = 0.3; vel_original(3,1) = 0.4; vel_original(3,2) = 0.5;

			array_1d<double, 3> embedded_vel;
			embedded_vel(0) = 1.0;
			embedded_vel(1) = 2.0;
			embedded_vel(2) = 3.0;
			for(unsigned int i=0; i<4; i++){
				pElement->GetGeometry()[i].SetValue(EMBEDDED_VELOCITY, embedded_vel);
				pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE)    = 0.0;
				pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE, 1) = 0.0;
				pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE, 2) = 0.0;
				for(unsigned int k=0; k<3; k++){
					pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel_original(i,k);
					pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i,k);
					pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*vel_original(i,k);
					pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
					pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
					pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
				}
			}
			pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
			pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) =  1.0;
			pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -1.0;
			pElement->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) =  1.0;

			// Compute RHS and LHS
			Vector RHS = ZeroVector(16);
			Matrix LHS = ZeroMatrix(16,16);

			const auto& r_process_info = modelPart.GetProcessInfo();
			pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
			pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

			// Check the RHS values (the RHS is computed as the LHS x previous_solution,
			// hence, it is assumed that if the RHS is correct, the LHS is correct as well)
			KRATOS_CHECK_NEAR(RHS(0), 0.023845, 1e-6);
			KRATOS_CHECK_NEAR(RHS(1), 0.048607, 1e-6);
			KRATOS_CHECK_NEAR(RHS(2), 0.0733691, 1e-7);
			KRATOS_CHECK_NEAR(RHS(3), -0.00618707, 1e-8);
			KRATOS_CHECK_NEAR(RHS(4), 0.179191, 1e-3);
			KRATOS_CHECK_NEAR(RHS(5), 3.49007, 1e-2);
			KRATOS_CHECK_NEAR(RHS(6), 4.62064, 1e-2);
			KRATOS_CHECK_NEAR(RHS(7), -0.00350587, 1e-8);
			KRATOS_CHECK_NEAR(RHS(8), 0.0229279, 1e-7);
			KRATOS_CHECK_NEAR(RHS(9), 0.0476899, 1e-7);
			KRATOS_CHECK_NEAR(RHS(10), 0.072452, 1e-6);
			KRATOS_CHECK_NEAR(RHS(11), 0.00207051, 1e-8);
			KRATOS_CHECK_NEAR(RHS(12), 2.66708, 1e-3);
			KRATOS_CHECK_NEAR(RHS(13), 4.29637, 1e-2);
			KRATOS_CHECK_NEAR(RHS(14), 5.67408, 1e-2);
			KRATOS_CHECK_NEAR(RHS(15), 0.000903677, 1e-9);

		}
	} // namespace Testing
}  // namespace Kratos.
