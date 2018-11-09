//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez, Ruben Zorrilla
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
#include "utilities/math_utils.h"

#include "custom_elements/two_fluid_navier_stokes.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"
#include "custom_constitutive/newtonian_two_fluid_3d_law.h"

#include "processes/find_nodal_neighbours_process.h"
#include "utilities/normal_calculation_utils.h"

namespace Kratos {
	namespace Testing {

		typedef ModelPart::IndexType									 IndexType;
		typedef ModelPart::NodeIterator					          NodeIteratorType;

	    /** Checks the TwoFluidNavierStokes2D3N element.
	     * Checks the LHS and RHS computation
	     */
	    KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokes2D3N, FluidDynamicsApplicationFastSuite)
		{

			ModelPart modelPart("Main");
			modelPart.SetBufferSize(3);

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			modelPart.AddNodalSolutionStepVariable(DENSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
			modelPart.AddNodalSolutionStepVariable(PRESSURE);
			modelPart.AddNodalSolutionStepVariable(VELOCITY);
			modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
			modelPart.AddNodalSolutionStepVariable(DISTANCE);

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
			Properties::Pointer pElemProp = modelPart.pGetProperties(0);
			pElemProp->SetValue(DENSITY, 1000.0);
			pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
			Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
			pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

			// Geometry creation
			modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
			modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
			std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};
			modelPart.CreateNewElement("TwoFluidNavierStokes2D3N", 1, elemNodes, pElemProp);

			Element::Pointer pElement = modelPart.pGetElement(1);

			// Define the nodal values
			Matrix vel_original(3,2);
			vel_original(0,0) = 0.0; vel_original(0,1) = 0.1;
			vel_original(1,0) = 0.1; vel_original(1,1) = 0.2;
			vel_original(2,0) = 0.2; vel_original(2,1) = 0.3;

			// Set the nodal DENSITY and DYNAMIC_VISCOSITY values
			for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
				it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
				it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
				it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = -9.81;
			}

			for(unsigned int i=0; i<3; i++){
				pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE)    = 0.0;
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

			pElement->Initialize(); // Initialize the element to initialize the constitutive law
			pElement->CalculateLocalSystem(LHS, RHS, modelPart.GetProcessInfo());

			// Check the RHS values (the RHS is computed as the LHS x previous_solution, 
			// hence, it is assumed that if the RHS is correct, the LHS is correct as well)
			KRATOS_CHECK_NEAR(RHS(0), -42.7102445, 1e-7);
			KRATOS_CHECK_NEAR(RHS(1), 29.7509341, 1e-7);
			KRATOS_CHECK_NEAR(RHS(2), -0.2010738, 1e-7);
			KRATOS_CHECK_NEAR(RHS(3), 90.0727138, 1e-7);
			KRATOS_CHECK_NEAR(RHS(4), 94.9469995, 1e-7);
			KRATOS_CHECK_NEAR(RHS(5), 0.0671813, 1e-7);
			KRATOS_CHECK_NEAR(RHS(6), 75.7625307, 1e-7);
			KRATOS_CHECK_NEAR(RHS(7), 146.5520664, 1e-7);
			KRATOS_CHECK_NEAR(RHS(8), -0.0161075, 1e-7);

	    }

	    // /** Checks the TwoFluidNavierStokes3D4N element.
	    //  * Checks the LHS and RHS computation using a small perturbation.
	    //  */
	    KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokes3D4N, FluidDynamicsApplicationFastSuite)
		{

			ModelPart modelPart("Main");
			modelPart.SetBufferSize(3);

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			modelPart.AddNodalSolutionStepVariable(DENSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
			modelPart.AddNodalSolutionStepVariable(PRESSURE);
			modelPart.AddNodalSolutionStepVariable(VELOCITY);
			modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
			modelPart.AddNodalSolutionStepVariable(DISTANCE);

			// Process info creation
			double delta_time = 0.1;
			modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
			modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
			Vector bdf_coefs(3);
			bdf_coefs[0] = 3.0/(2.0*delta_time);
			bdf_coefs[1] = -2.0/delta_time;
			bdf_coefs[2] = 0.5*delta_time;
			modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

			// Set the element properties
			Properties::Pointer pElemProp = modelPart.pGetProperties(0);
			pElemProp->SetValue(DENSITY, 1000.0);
			pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
			NewtonianTwoFluid3DLaw::Pointer pConsLaw(new NewtonianTwoFluid3DLaw());
			pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

			// Geometry creation
			modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
			modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
			modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
			std::vector<ModelPart::IndexType> elemNodes {1, 2, 3, 4};
			modelPart.CreateNewElement("TwoFluidNavierStokes3D4N", 1, elemNodes, pElemProp);

			Element::Pointer pElement = modelPart.pGetElement(1);

			// Define the nodal values
			Matrix vel_original(4,3);
			vel_original(0,0) = 0.0; vel_original(0,1) = 0.1; vel_original(0,2) = 0.2;
			vel_original(1,0) = 0.1; vel_original(1,1) = 0.2; vel_original(1,2) = 0.3;
			vel_original(2,0) = 0.2; vel_original(2,1) = 0.3; vel_original(2,2) = 0.4;
			vel_original(3,0) = 0.3; vel_original(3,1) = 0.4; vel_original(3,2) = 0.5;

			// Set the nodal BODY_FORCE, DENSITY and DYNAMIC_VISCOSITY values
			for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
				it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
				it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
				it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = -9.81;
			}

			for(unsigned int i=0; i<4; i++){
				pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE)    = 0.0;
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

			pElement->Initialize(); // Initialize the element to initialize the constitutive law
			pElement->CalculateLocalSystem(LHS, RHS, modelPart.GetProcessInfo());


			// Check the RHS values (the RHS is computed as the LHS x previous_solution, 
			// hence, it is assumed that if the RHS is correct, the LHS is correct as well)
			KRATOS_CHECK_NEAR(RHS(0), -118.03889688, 1e-7);
			KRATOS_CHECK_NEAR(RHS(1), 14.78243843, 1e-7);
			KRATOS_CHECK_NEAR(RHS(2), -213.62425998, 1e-7);
			KRATOS_CHECK_NEAR(RHS(3), -0.05448084, 1e-7);
			KRATOS_CHECK_NEAR(RHS(4), 22.78371861, 1e-7);
			KRATOS_CHECK_NEAR(RHS(5), 27.42496766, 1e-7);
			KRATOS_CHECK_NEAR(RHS(6), -403.33194690, 1e-7);
			KRATOS_CHECK_NEAR(RHS(7), 0.15778148, 1e-7);
			KRATOS_CHECK_NEAR(RHS(8), 60.36872332, 1e-7);
			KRATOS_CHECK_NEAR(RHS(9), 15.21517182, 1e-7);
			KRATOS_CHECK_NEAR(RHS(10), -419.85282231, 1e-7);
			KRATOS_CHECK_NEAR(RHS(11), 0.00659561, 1e-7);
			KRATOS_CHECK_NEAR(RHS(12), 80.61562161, 1e-7);
			KRATOS_CHECK_NEAR(RHS(13), 37.68158875, 1e-7);
			KRATOS_CHECK_NEAR(RHS(14), -453.71180415, 1e-7);
			KRATOS_CHECK_NEAR(RHS(15), -0.20989625, 1e-7);

		}

		/** Checks the TwoFluidNavierStokes2D3N element in a hydrostatic case.
	     *  Checks the computation of the RHS
	     */
	    KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokes2D3NHydrostatic, FluidDynamicsApplicationFastSuite)
		{
			// std::cout << "\nThis is supposed to represent a hydro-static case with a surface at y = 1\n\n" << std::endl;
			ModelPart modelPart("Main");
			modelPart.SetBufferSize(3);

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			modelPart.AddNodalSolutionStepVariable(DENSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
			modelPart.AddNodalSolutionStepVariable(PRESSURE);
			modelPart.AddNodalSolutionStepVariable(VELOCITY);
			modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
			modelPart.AddNodalSolutionStepVariable(DISTANCE);

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
			Properties::Pointer pElemProp = modelPart.pGetProperties(0);
			pElemProp->SetValue(DENSITY, 1000.0);
			pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-03);
			Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
			pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

			// Geometry creation
			modelPart.CreateNewNode(1, 2.0, 0.0, 0.0);  // 0 = node 1
			modelPart.CreateNewNode(2, 2.0, 2.0, 0.0);	// 1 = node 2
			modelPart.CreateNewNode(3, 0.0, 2.0, 0.0);	// 2 = node 3

			std::vector<ModelPart::IndexType> elemNodes1 {1, 2, 3};
			
			modelPart.CreateNewElement("TwoFluidNavierStokes2D3N", 1, elemNodes1, pElemProp);

			Element::Pointer pElement = modelPart.pGetElement(1);

			// Define the nodal values as 0 for hydrostatic case
			Matrix vel_original(3,2);
			vel_original(0,0) = 0.0; vel_original(0,1) = 0.0;
			vel_original(1,0) = 0.0; vel_original(1,1) = 0.0;
			vel_original(2,0) = 0.0; vel_original(2,1) = 0.0;

			// Setting equal nodal values for DENSITY, DYNAMIC_VISCOSITY, BODY_FORCE
			for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
				it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
				it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
				it_node->FastGetSolutionStepValue(BODY_FORCE_X) = 0.0;
				it_node->FastGetSolutionStepValue(BODY_FORCE_Y) = -10.0;
				it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = 0.0;
			}

			// Setting the density (different for nodes since element cut by surface)
			pElement->GetGeometry()[0].FastGetSolutionStepValue(DENSITY) = 2.0;
			pElement->GetGeometry()[1].FastGetSolutionStepValue(DENSITY) = 1.0;
			pElement->GetGeometry()[2].FastGetSolutionStepValue(DENSITY) = 1.0;

			for(unsigned int i=0; i<3; i++){
				for(unsigned int k=0; k<2; k++){
					pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel_original(i,k);
					pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = vel_original(i,k);
					pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = vel_original(i,k);
					// pElement->GetGeometry()[i].Fix(VELOCITY);
					pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
					pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
					pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
				}
			}

			pElement->GetGeometry()[0].Fix(VELOCITY_X);
			pElement->GetGeometry()[0].Fix(VELOCITY_Y);

			// Setting the density (different for nodes to define the position of the surface)
			pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
			pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
			pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = 1.0;

			// Simon : Setting the pressure
			pElement->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE) = 30.0;
			pElement->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE) = 0.0;
			pElement->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE) = 0.0;

			pElement->GetGeometry()[0].Fix(PRESSURE);

			// Compute RHS and LHS
			Vector RHS = ZeroVector(9);
			Matrix LHS = ZeroMatrix(9,9);

			pElement->Initialize(); // Initialize the element to initialize the constitutive law
			pElement->CalculateLocalSystem(LHS, RHS, modelPart.GetProcessInfo());

			double det;
			MathUtils<double>::InvertMatrix(LHS, LHS, det);

			// auto solVec = RHS;

			const Vector solVec = prod(LHS, RHS);

			// for (int i = 0; i < 9; i++)
			// {
			// 	std::cout << "solVec(" << i << ") = " << solVec(i) << std::endl;
			// }

			// Check the RHS values (the RHS is computed as the LHS x previous_solution, 
			// hence, it is assumed that if the RHS is correct, the LHS is correct as well)

			std::cout << "RHS(0) = " << RHS(0) << std::endl;
			std::cout << "RHS(1) = " << RHS(1) << std::endl;
			std::cout << "RHS(2) = " << RHS(2) << std::endl;
			std::cout << "RHS(3) = " << RHS(3) << std::endl;
			std::cout << "RHS(4) = " << RHS(4) << std::endl;
			std::cout << "RHS(5) = " << RHS(5) << std::endl;
			std::cout << "RHS(6) = " << RHS(6) << std::endl;
			std::cout << "RHS(7) = " << RHS(7) << std::endl;
			std::cout << "RHS(8) = " << RHS(8) << std::endl;

			// for (int i = 0; i < 9; i++)
			// {
			// 	for (int j = 0; j < 9; j++)
			// 	{
			// 		std::cout << "LHS(" << i << "," << j << ") = " << LHS(i,j) << std::endl;
			// 	}
			// }
			
			// The remaining residuals in the velocities have the size of the boundary integrals over the enriched pressure.
			// If the "standard" pressure shape functions are used, the results do not hold.
			KRATOS_CHECK_NEAR(RHS(0), 0.0, 1e-7);		// U_x at node 1
			KRATOS_CHECK_NEAR(RHS(1), -17.5, 1e-7); 	// U_y at node 1  
			KRATOS_CHECK_NEAR(RHS(2), 0.0, 1e-7);		// P   at node 1

			KRATOS_CHECK_NEAR(RHS(3), 7.5, 1e-7);		// U_x at node 2		
			KRATOS_CHECK_NEAR(RHS(4), 0.0, 1e-7);		// U_y at node 2
			KRATOS_CHECK_NEAR(RHS(5), 0.0, 1e-7);		// P   at node 2

			KRATOS_CHECK_NEAR(RHS(6), -7.5, 1e-7);		// U_x at node 3
			KRATOS_CHECK_NEAR(RHS(7), -7.5, 1e-7);		// U_y at node 3
			KRATOS_CHECK_NEAR(RHS(8), 0.0, 1e-7);		// P   at node 3
			// std::cout << "\nThe test case finished. End.\n" << std::endl;
	    }


		/** Checks the TwoFluidNavierStokes2D3N element in a hydrostatic case.
	     *  Checks the computation of the RHS
	     */
		KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokes2D3NHydrostaticBehr, FluidDynamicsApplicationFastSuite){

            ModelPart modelPart("TestPart");
			modelPart.SetBufferSize(3);

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			modelPart.AddNodalSolutionStepVariable(DENSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
			modelPart.AddNodalSolutionStepVariable(PRESSURE);
			modelPart.AddNodalSolutionStepVariable(VELOCITY);
			modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
			modelPart.AddNodalSolutionStepVariable(DISTANCE);
            modelPart.AddNodalSolutionStepVariable(NORMAL);

			// Process info creation
            double delta_time = 0.01;
			modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.0);
			modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
			modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
			Vector bdf_coefs(3);
			bdf_coefs[0] = 3.0/(2.0*delta_time);
			bdf_coefs[1] = -2.0/delta_time;
			bdf_coefs[2] = 0.5*delta_time;
			modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

			// Set the element properties
			Properties::Pointer pElemProp = modelPart.pGetProperties(0);
			pElemProp->SetValue(DENSITY, 1000.0);
			pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
			Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
			pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

			// Geometry creation
			modelPart.CreateNewNode(1, 0.0, 2.0, 0.0);
			modelPart.CreateNewNode(2, 0.0, 1.0, 0.0);
            modelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
            modelPart.CreateNewNode(4, 1.0, 0.0, 0.0);

            // Creation of elements
			std::vector<ModelPart::IndexType> elemNodes1 {1, 2, 3};
            std::vector<ModelPart::IndexType> elemNodes2 {2, 4, 3};
			modelPart.CreateNewElement("TwoFluidNavierStokes2D3N", 1, elemNodes1, pElemProp);
            modelPart.CreateNewElement("TwoFluidNavierStokes2D3N", 2, elemNodes2, pElemProp);
            Element::Pointer pElement1 = modelPart.pGetElement(1);
            Element::Pointer pElement2 = modelPart.pGetElement(2);

            std::vector<ModelPart::IndexType> condNodes1 {1, 2};
            std::vector<ModelPart::IndexType> condNodes2 {2, 4};
            modelPart.CreateNewCondition("BehrWallCondition2D2N", 1, condNodes1, pElemProp);
            modelPart.CreateNewCondition("BehrWallCondition2D2N", 2, condNodes2, pElemProp);
            Condition::Pointer pCondition1 = modelPart.pGetCondition(1);
            Condition::Pointer pCondition2 = modelPart.pGetCondition(2);

            Vector elemRHS1 = ZeroVector(9);
            Vector elemRHS2 = ZeroVector(9);
            Matrix elemLHS = ZeroMatrix(9,9);

            Vector condRHS1 = ZeroVector(6);
            Vector condRHS2 = ZeroVector(6);
            Matrix condLHS = ZeroMatrix(6,6);

			for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
				it_node->FastGetSolutionStepValue(DENSITY) = 1000.0;
				it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
				it_node->FastGetSolutionStepValue(BODY_FORCE_Y) = -10.0;
			}

			for(unsigned int i=0; i<3; i++){
				for(unsigned int k=0; k<2; k++){
					pElement1->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = 0.0;
					pElement1->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.0;
					pElement1->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.0;
					pElement1->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
					pElement1->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
					pElement1->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
				}
			}

            pElement1->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    = 10000.0;
            pElement1->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    = 20000.0;
            pElement1->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)    = 20000.0;

            for(unsigned int i=0; i<3; i++){
				for(unsigned int k=0; k<2; k++){
					pElement2->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = 0.0;
					pElement2->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.0;
					pElement2->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.0;
					pElement2->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
					pElement2->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
					pElement2->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
				}
			}
            
            pElement2->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    = 20000.0;
            pElement2->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    = 30000.0;
            pElement2->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)    = 20000.0;

			pElement1->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
			pElement1->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -2.0;
			pElement1->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -2.0;
			pElement2->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -2.0;
			pElement2->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -3.0;
			pElement2->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -2.0;

            // Initialization 
            pElement1->Initialize();
            pElement2->Initialize();

			FindNodalNeighboursProcess find_nodal_neighbours_process(modelPart);
			find_nodal_neighbours_process.Execute();

            NormalCalculationUtils find_nodal_normal_utility;
            find_nodal_normal_utility.CalculateOnSimplex(modelPart, 2);

            pCondition1->Initialize();
            pCondition2->Initialize();

            // Computing locel contributions
            pCondition1->CalculateLocalSystem(condLHS, condRHS1, modelPart.GetProcessInfo());
            pCondition2->CalculateLocalSystem(condLHS, condRHS2, modelPart.GetProcessInfo());

			pElement1->CalculateLocalSystem(elemLHS, elemRHS1, modelPart.GetProcessInfo());
            pElement2->CalculateLocalSystem(elemLHS, elemRHS2, modelPart.GetProcessInfo());

            // Assembly of the residual for node 2 (node between the conditions)
            Vector contriFromElem1 = ZeroVector(3);
            Vector contriFromElem2 = ZeroVector(3);
            Vector contriFromCond1 = ZeroVector(3);
            Vector contriFromCond2 = ZeroVector(3);

            for (unsigned i = 0; i < 3; i++){
                contriFromElem1[i] = elemRHS1[3+i];
                contriFromElem2[i] = elemRHS2[0+i];

                contriFromCond1[i] = condRHS1[3+i];
                contriFromCond2[i] = condRHS2[0+i];
            }

            // KRATOS_WATCH( contriFromElem1 )
            // KRATOS_WATCH( contriFromElem2 )
            // KRATOS_WATCH( contriFromCond1 )
            // KRATOS_WATCH( contriFromCond2 )

            Vector residualAtNodeTwo = contriFromElem1 + contriFromElem2 + contriFromCond1 + contriFromCond2;
            Vector normalAtNodeTwo = pElement1->GetGeometry()[1].FastGetSolutionStepValue(NORMAL);

			// KRATOS_WATCH( normalAtNodeTwo )

			Vector tangentialComponent;
			array_1d<double,3> residualAtNodeTwoVector;
			array_1d<double,3> normalAtNodeTwoVector;

			tangentialComponent = MathUtils<double>::CrossProduct( residualAtNodeTwo, normalAtNodeTwo );

            // KRATOS_WATCH( sum )
            std::cout << "tangentialComponent = " << tangentialComponent << std::endl;

			KRATOS_CHECK_NEAR( tangentialComponent[0], 0.0, 1e-7);
			KRATOS_CHECK_NEAR( tangentialComponent[1], 0.0, 1e-7);
			KRATOS_CHECK_NEAR( tangentialComponent[2], 0.0, 1e-7);
        }

		
		
		KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokes2D3NMotionBehr, FluidDynamicsApplicationFastSuite){

            ModelPart modelPart("TestPart");
			modelPart.SetBufferSize(3);
			// Motion is applied to the only NON-BOUNDARY node
			// This is done to generate viscous stress, which were not considered so far

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			modelPart.AddNodalSolutionStepVariable(DENSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
			modelPart.AddNodalSolutionStepVariable(PRESSURE);
			modelPart.AddNodalSolutionStepVariable(VELOCITY);
			modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
			modelPart.AddNodalSolutionStepVariable(DISTANCE);
            modelPart.AddNodalSolutionStepVariable(NORMAL);

			// Process info creation
            double delta_time = 0.01;
			modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.0);
			modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
			modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
			Vector bdf_coefs(3);
			bdf_coefs[0] = 3.0/(2.0*delta_time);
			bdf_coefs[1] = -2.0/delta_time;
			bdf_coefs[2] = 0.5*delta_time;
			modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

			// Set the element properties
			Properties::Pointer pElemProp = modelPart.pGetProperties(0);
			pElemProp->SetValue(DENSITY, 1000.0);
			pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
			Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
			pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

			// Geometry creation
			modelPart.CreateNewNode(1, 0.0, 2.0, 0.0);
			modelPart.CreateNewNode(2, 0.0, 1.0, 0.0);
            modelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
            modelPart.CreateNewNode(4, 0.0, 0.0, 0.0);

            // Creation of elements
			std::vector<ModelPart::IndexType> elemNodes1 {1, 2, 3};
            std::vector<ModelPart::IndexType> elemNodes2 {2, 4, 3};
			modelPart.CreateNewElement("TwoFluidNavierStokes2D3N", 1, elemNodes1, pElemProp);
            modelPart.CreateNewElement("TwoFluidNavierStokes2D3N", 2, elemNodes2, pElemProp);
            Element::Pointer pElement1 = modelPart.pGetElement(1);
            Element::Pointer pElement2 = modelPart.pGetElement(2);

            std::vector<ModelPart::IndexType> condNodes1 {1, 2};
            std::vector<ModelPart::IndexType> condNodes2 {2, 4};
            modelPart.CreateNewCondition("BehrWallCondition2D2N", 1, condNodes1, pElemProp);
            modelPart.CreateNewCondition("BehrWallCondition2D2N", 2, condNodes2, pElemProp);
            Condition::Pointer pCondition1 = modelPart.pGetCondition(1);
            Condition::Pointer pCondition2 = modelPart.pGetCondition(2);

            Vector elemRHS1 = ZeroVector(9);
            Vector elemRHS2 = ZeroVector(9);
            Matrix elemLHS = ZeroMatrix(9,9);

            Vector condRHS1 = ZeroVector(6);
            Vector condRHS2 = ZeroVector(6);
            Matrix condLHS = ZeroMatrix(6,6);

			for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
				it_node->FastGetSolutionStepValue(DENSITY) = 1000.0;
				it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
				it_node->FastGetSolutionStepValue(BODY_FORCE_X) = 0.0;
				it_node->FastGetSolutionStepValue(BODY_FORCE_Y) = 0.0;
				it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = 0.0;
			}

			for(unsigned int i=0; i<3; i++){
				for (unsigned int k = 0; k < 2; k++)
				{
					// fixing the mesh
					pElement1->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
					pElement1->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
					pElement1->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;

					pElement2->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
					pElement2->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
					pElement2->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
				}
			}

			for(unsigned int i=0; i<3; i++){
				// setting Vx to 0.0
				pElement1->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[0]    = 0.0;
				pElement1->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[0] = 0.0;
				pElement1->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[0] = 0.0;
				pElement2->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[0]    = 0.0;
				pElement2->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[0] = 0.0;
				pElement2->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[0] = 0.0;
				// setting Vy to 0.0
				pElement1->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[1]    = 1.0;
				pElement1->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[1] = 1.0;
				pElement1->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[1] = 1.0;
				pElement2->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[1]    = 1.0;
				pElement2->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[1] = 1.0;
				pElement2->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[1] = 1.0;
			}


            pElement1->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    = 10000.0;
            pElement1->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    = 10000.0;
            pElement1->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)    = 10000.0;
			pElement2->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    = 10000.0;
            pElement2->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    = 10000.0;
            pElement2->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)    = 10000.0;


			pElement1->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
			pElement1->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -2.0;
			pElement1->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -2.0;
			pElement2->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -2.0;
			pElement2->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -3.0;
			pElement2->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -2.0;

            // Initialization 
            pElement1->Initialize();
            pElement2->Initialize();

			FindNodalNeighboursProcess find_nodal_neighbours_process(modelPart);
			find_nodal_neighbours_process.Execute();

            NormalCalculationUtils find_nodal_normal_utility;
            find_nodal_normal_utility.CalculateOnSimplex(modelPart, 2);

            pCondition1->Initialize();
            pCondition2->Initialize();

            // Computing locel contributions
            pCondition1->CalculateLocalSystem(condLHS, condRHS1, modelPart.GetProcessInfo());
            pCondition2->CalculateLocalSystem(condLHS, condRHS2, modelPart.GetProcessInfo());

			pElement1->CalculateLocalSystem(elemLHS, elemRHS1, modelPart.GetProcessInfo());
            pElement2->CalculateLocalSystem(elemLHS, elemRHS2, modelPart.GetProcessInfo());

            // Assembly of the residual for node 2 (node between the conditions)
            Vector contriFromElem1 = ZeroVector(3);
            Vector contriFromElem2 = ZeroVector(3);
            Vector contriFromCond1 = ZeroVector(3);
            Vector contriFromCond2 = ZeroVector(3);

            for (unsigned i = 0; i < 3; i++){
                contriFromElem1[i] = elemRHS1[3+i];
                contriFromElem2[i] = elemRHS2[0+i];

                contriFromCond1[i] = condRHS1[3+i];
                contriFromCond2[i] = condRHS2[0+i];
            }

            Vector residualAtNodeTwo = contriFromElem1 + contriFromElem2 + contriFromCond1 + contriFromCond2;
            Vector normalAtNodeTwo = pElement1->GetGeometry()[1].FastGetSolutionStepValue(NORMAL);

			// KRATOS_WATCH( normalAtNodeTwo )

			Vector tangentialComponent;
			Vector normalComponent;
			array_1d<double,3> residualAtNodeTwoVector;
			array_1d<double,3> normalAtNodeTwoVector;

			tangentialComponent = MathUtils<double>::CrossProduct( residualAtNodeTwo, normalAtNodeTwo );
			normalComponent = MathUtils<double>::Dot3( residualAtNodeTwo, normalAtNodeTwo ) * normalAtNodeTwo;

            // KRATOS_WATCH( sum )
            std::cout << "tangentialComponent = " << tangentialComponent << std::endl;
			std::cout << "normalComponent = " << normalComponent << std::endl;

			KRATOS_CHECK_NEAR( normalComponent[0], -10000.0, 1e-4);
			KRATOS_CHECK_NEAR( normalComponent[1], 0.0, 1e-4);
			KRATOS_CHECK_NEAR( normalComponent[2], 0.0, 1e-4);
        }

	} // namespace Testing
}  // namespace Kratos.
