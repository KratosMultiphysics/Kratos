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
#include "custom_elements/navier_stokes.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"

namespace Kratos {
	namespace Testing {

		typedef ModelPart::IndexType									 IndexType;
		typedef ModelPart::NodeIterator					          NodeIteratorType;

	    /** Checks the NavierStokes2D3N element.
	     * Checks the LHS and RHS computation using a small perturbation.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(ElementNavierStokes2D3N, FluidDynamicsApplicationFastSuite)
		{
			Model model;
			ModelPart& modelPart = model.CreateModelPart("Main", 3);

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			modelPart.AddNodalSolutionStepVariable(DENSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
			modelPart.AddNodalSolutionStepVariable(SOUND_VELOCITY);
			modelPart.AddNodalSolutionStepVariable(PRESSURE);
			modelPart.AddNodalSolutionStepVariable(VELOCITY);
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
			modelPart.CreateNewElement("NavierStokes2D3N", 1, elemNodes, pElemProp);

			Element::Pointer pElement = modelPart.pGetElement(1);

			// Define the nodal values
			Matrix vel_original(3,2);
			vel_original(0,0) = 0.0; vel_original(0,1) = 0.1;
			vel_original(1,0) = 0.1; vel_original(1,1) = 0.2;
			vel_original(2,0) = 0.2; vel_original(2,1) = 0.3;

			Matrix perturbation_matrix(3,2);
			perturbation_matrix(0,0) = 1.0; perturbation_matrix(0,1) = 0.1;
			perturbation_matrix(1,0) = 0.1; perturbation_matrix(1,1) = 0.2;
			perturbation_matrix(2,0) = 0.2; perturbation_matrix(2,1) = 0.3;

			// Set the nodal DENSITY and DYNAMIC_VISCOSITY values
			for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
				it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
				it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
			}

			for(unsigned int i=0; i<3; i++){
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

			// Compute RHS and LHS
			Vector RHS = ZeroVector(9);
			Matrix LHS = ZeroMatrix(9,9);

			pElement->Initialize(); // Initialize the element to initialize the constitutive law
			pElement->CalculateLocalSystem(LHS, RHS, modelPart.GetProcessInfo());

			// Compute the error of the perturbation
			double perturbation = 2e-2;
			std::vector<double> error_norms;
			std::vector<double> error_vel_norms;
			std::vector<double> error_pres_norms;

			for (unsigned int j=1; j<5; ++j)
			{
				perturbation /= 2;
				Vector perturbation_vector = ZeroVector(9);

				for(unsigned int i=0; i<3; i++){
					for(unsigned int k=0; k<2; k++){
						pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k] = vel_original(i,k) + perturbation*perturbation_matrix(i,k);
						perturbation_vector(i*3+k) = perturbation*perturbation_matrix(i,k);
					}
				}

				// Compute perturbed RHS and LHS
				Vector error = ZeroVector(9);
				Vector RHS_obtained = ZeroVector(9);
				Vector RHS_perturbed = ZeroVector(9);
				Vector solution_increment = ZeroVector(9);

				pElement->CalculateRightHandSide(RHS_perturbed, modelPart.GetProcessInfo());

				solution_increment = prod(LHS, perturbation_vector);
				noalias(RHS_obtained) = RHS - solution_increment;
				noalias(error) = RHS_perturbed - RHS_obtained;
				error_norms.push_back(norm_2(error));
			}

			// // Check quadratic convergence (if FullNR has been selected when generating the element)
			// for(unsigned int i=1; i<error_norms.size(); ++i)
			// 	KRATOS_CHECK_NEAR(error_norms[i-1]/error_norms[i], 4.0, 1e-1);

			// Check quadratic convergence (if Picard has been selected when generating the element)
			for(unsigned int i=1; i<error_norms.size(); ++i)
				KRATOS_CHECK_NEAR(error_norms[i-1]/error_norms[i], 2.0, 2.5e-1);

			// std::cout<<std::endl;
			// for(unsigned int i=0;i<error_norms.size();++i){
			// 	std::cout << "Error norm "<< i << " " << error_norms[i] << std::endl;
			// }

	    }

		/** Checks the NavierStokes2D3N element.
		 * Checks the LHS and RHS stationary solid rigid movements.
		 */
	    KRATOS_TEST_CASE_IN_SUITE(ElementNavierStokes2D3NStationary, FluidDynamicsApplicationFastSuite)
		{
			Model model;
			ModelPart& modelPart = model.CreateModelPart("Main", 3);

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			modelPart.AddNodalSolutionStepVariable(DENSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
			modelPart.AddNodalSolutionStepVariable(SOUND_VELOCITY);
			modelPart.AddNodalSolutionStepVariable(PRESSURE);
			modelPart.AddNodalSolutionStepVariable(VELOCITY);
			modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);

			// Process info creation
			double delta_time = 0.1;
			modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.01);
			modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+03);
			modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
			Vector bdf_coefs(3);
			bdf_coefs[0] = 0.0;
			bdf_coefs[1] = 0.0;
			bdf_coefs[2] = 0.0;
			modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

			// Set the element properties
			Properties::Pointer pElemProp = modelPart.pGetProperties(0);
			pElemProp->SetValue(DENSITY, 1.0);
			pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0);
			Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
			pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

			// Geometry creation
			modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
			modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
			std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};
			modelPart.CreateNewElement("NavierStokes2D3N", 1, elemNodes, pElemProp);

			Element::Pointer pElement = modelPart.pGetElement(1);

			// Define the nodal values
			array_1d<double, 3> velocity_values;
			velocity_values[0] = 0.0;
			velocity_values[1] = 0.0;
			velocity_values[2] = 0.0;
			double pressure_value = 0.0;

			// Set the nodal values
			for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
				it_node->FastGetSolutionStepValue(PRESSURE) = pressure_value;
				it_node->FastGetSolutionStepValue(VELOCITY) = velocity_values;
				it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
				it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
			}

			// Compute RHS and LHS
			Vector RHS = ZeroVector(9);
			Matrix LHS = ZeroMatrix(9,9);

			pElement->Initialize(); // Initialize the element to initialize the constitutive law
			pElement->CalculateLocalSystem(LHS, RHS, modelPart.GetProcessInfo());

			// Check obtained RHS
			double sum_RHS = 0.0;
			for (unsigned int i=0; i<RHS.size(); ++i)
				sum_RHS += RHS[i];
			KRATOS_CHECK_NEAR(sum_RHS, 0.0, 1e-12);

			// Check rigid movement modes
			Vector a(9);
			Vector rhs(9);
			double sum_rhs;
			// Mode 1 check
			a[0] = 1.0; a[1] = 0.0; a[2] = 0.0;	a[3] = 1.0;	a[4] = 0.0;	a[5] = 0.0;	a[6] = 1.0;	a[7] = 0.0;	a[8] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i)
				sum_rhs += rhs[i];
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-12);
			// Mode 2 check
			a[0] = 0.0; a[1] = 1.0; a[2] = 0.0;	a[3] = 0.0;	a[4] = 1.0;	a[5] = 0.0;	a[6] = 0.0;	a[7] = 1.0;	a[8] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i)
				sum_rhs += rhs[i];
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-12);
			// Mode 3 check
			a[0] = 1.0; a[1] = 1.0; a[2] = 0.0;	a[3] = 1.0;	a[4] = 1.0;	a[5] = 0.0;	a[6] = 1.0;	a[7] = 1.0;	a[8] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i)
				sum_rhs += rhs[i];
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-12);

	    }

	    // /** Checks the NavierStokes3D4N element.
	    //  * Checks the LHS and RHS computation using a small perturbation.
	    //  */
	    KRATOS_TEST_CASE_IN_SUITE(ElementNavierStokes3D4N, FluidDynamicsApplicationFastSuite)
		{
			Model model;
			ModelPart& modelPart = model.CreateModelPart("Main", 3);

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			modelPart.AddNodalSolutionStepVariable(DENSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
			modelPart.AddNodalSolutionStepVariable(SOUND_VELOCITY);
			modelPart.AddNodalSolutionStepVariable(PRESSURE);
			modelPart.AddNodalSolutionStepVariable(VELOCITY);
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
			Properties::Pointer pElemProp = modelPart.pGetProperties(0);
			pElemProp->SetValue(DENSITY, 1000.0);
			pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
			Newtonian3DLaw::Pointer pConsLaw(new Newtonian3DLaw());
			pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

			// Geometry creation
			modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
			modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
			modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
			std::vector<ModelPart::IndexType> elemNodes {1, 2, 3, 4};
			modelPart.CreateNewElement("NavierStokes3D4N", 1, elemNodes, pElemProp);

			Element::Pointer pElement = modelPart.pGetElement(1);

			// Define the nodal values
			Matrix vel_original(4,3);
			vel_original(0,0) = 0.0; vel_original(0,1) = 0.1; vel_original(0,2) = 0.2;
			vel_original(1,0) = 0.1; vel_original(1,1) = 0.2; vel_original(1,2) = 0.3;
			vel_original(2,0) = 0.2; vel_original(2,1) = 0.3; vel_original(2,2) = 0.4;
			vel_original(3,0) = 0.3; vel_original(3,1) = 0.4; vel_original(3,2) = 0.5;

			Matrix perturbation_matrix(4,3);
			perturbation_matrix(0,0) = 0.0; perturbation_matrix(0,1) = 0.1; perturbation_matrix(0,2) = 0.5;
			perturbation_matrix(1,0) = 0.1; perturbation_matrix(1,1) = 0.2; perturbation_matrix(1,2) = 0.6;
			perturbation_matrix(2,0) = 0.2; perturbation_matrix(2,1) = 0.3; perturbation_matrix(2,2) = 0.7;
			perturbation_matrix(3,0) = 0.3; perturbation_matrix(3,1) = 0.4; perturbation_matrix(3,2) = 0.8;

			// Set the nodal DENSITY and DYNAMIC_VISCOSITY values
			for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
				it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
				it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
			}

			for(unsigned int i=0; i<4; i++){
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

			// Compute RHS and LHS
			Vector RHS = ZeroVector(16);
			Matrix LHS = ZeroMatrix(16,16);

			pElement->Initialize(); // Initialize the element to initialize the constitutive law
			pElement->CalculateLocalSystem(LHS, RHS, modelPart.GetProcessInfo());

			// Compute the error of the perturbation
			double perturbation = 2e-2;
			std::vector<double> error_norms;
			std::vector<double> error_vel_norms;
			std::vector<double> error_pres_norms;

			for (unsigned int j=1; j<5; ++j)
			{
				perturbation /= 2;
				Vector perturbation_vector = ZeroVector(16);

				for(unsigned int i=0; i<4; i++){
					for(unsigned int k=0; k<3; k++){
						pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k] = vel_original(i,k) + perturbation*perturbation_matrix(i,k);
						perturbation_vector(i*4+k) = perturbation*perturbation_matrix(i,k);
					}
				}

				// Compute perturbed RHS and LHS
				Vector error = ZeroVector(16);
				Vector RHS_obtained = ZeroVector(16);
				Vector RHS_perturbed = ZeroVector(16);
				Vector solution_increment = ZeroVector(16);

				pElement->CalculateRightHandSide(RHS_perturbed, modelPart.GetProcessInfo());

				solution_increment = prod(LHS, perturbation_vector);
				noalias(RHS_obtained) = RHS - solution_increment;
				noalias(error) = RHS_perturbed - RHS_obtained;
				error_norms.push_back(norm_2(error));

			}

			// // Check quadratic convergence (if FullNR has been selected when generating the element)
			// for(unsigned int i=1; i<error_norms.size(); ++i)
			// 	KRATOS_CHECK_NEAR(error_norms[i-1]/error_norms[i], 4.0, 1e-1);

			// Check quadratic convergence (if Picard has been selected when generating the element)
			for(unsigned int i=1; i<error_norms.size(); ++i)
				KRATOS_CHECK_NEAR(error_norms[i-1]/error_norms[i], 2.0, 2.5e-1);

			// std::cout<<std::endl;
			// for(unsigned int i=0; i<error_norms.size(); ++i)
			// 	std::cout << "Error norm "<< i << " " << error_norms[i] << std::endl;

		}

		KRATOS_TEST_CASE_IN_SUITE(ElementNavierStokes3D4NStationary, FluidDynamicsApplicationFastSuite)
		{
			Model model;
			ModelPart& modelPart = model.CreateModelPart("Main", 3);

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			modelPart.AddNodalSolutionStepVariable(DENSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
			modelPart.AddNodalSolutionStepVariable(SOUND_VELOCITY);
			modelPart.AddNodalSolutionStepVariable(PRESSURE);
			modelPart.AddNodalSolutionStepVariable(VELOCITY);
			modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);

			// Process info creation
			double delta_time = 0.1;
			modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.01);
			modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+03);
			modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
			Vector bdf_coefs(3);
			bdf_coefs[0] = 0.0;
			bdf_coefs[1] = 0.0;
			bdf_coefs[2] = 0.0;
			modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

			// Set the element properties
			Properties::Pointer pElemProp = modelPart.pGetProperties(0);
			pElemProp->SetValue(DENSITY, 1.0);
			pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0);
			Newtonian3DLaw::Pointer pConsLaw(new Newtonian3DLaw());
			pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

			// Geometry creation
			modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
			modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
			modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
			modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
			std::vector<ModelPart::IndexType> elemNodes {1, 2, 3, 4};
			modelPart.CreateNewElement("NavierStokes3D4N", 1, elemNodes, pElemProp);

			Element::Pointer pElement = modelPart.pGetElement(1);

			// Define the nodal values
			array_1d<double, 3> velocity_values;
			velocity_values[0] = 0.0;
			velocity_values[1] = 0.0;
			velocity_values[2] = 0.0;
			double pressure_value = 0.0;

			// Set the nodal values
			for (NodeIteratorType it_node=modelPart.NodesBegin(); it_node<modelPart.NodesEnd(); ++it_node){
				it_node->FastGetSolutionStepValue(PRESSURE) = pressure_value;
				it_node->FastGetSolutionStepValue(VELOCITY) = velocity_values;
				it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
				it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
			}

			// Compute RHS and LHS
			Vector RHS = ZeroVector(16);
			Matrix LHS = ZeroMatrix(16,16);

			pElement->Initialize(); // Initialize the element to initialize the constitutive law
			pElement->CalculateLocalSystem(LHS, RHS, modelPart.GetProcessInfo());

			// Check obtained RHS
			double sum_RHS = 0.0;
			for (unsigned int i=0; i<RHS.size(); ++i)
				sum_RHS += RHS[i];
			KRATOS_CHECK_NEAR(sum_RHS, 0.0, 1e-12);

			// Check modes
			Vector a(16);
			Vector rhs(16);
			double sum_rhs;
			// Mode 1 check
			a[0] = 1.0; a[1] = 0.0; a[2] = 0.0;	a[3] = 0.0;	a[4] = 1.0;	a[5] = 0.0;	a[6] = 0.0;	a[7] = 0.0;	a[8] = 1.0; a[9] = 0.0; a[10] = 0.0; a[11] = 0.0; a[12] = 1.0; a[13] = 0.0; a[14] = 0.0; a[15] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i)
				sum_rhs += rhs[i];
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-12);
			// Mode 2 check
			a[0] = 0.0; a[1] = 1.0; a[2] = 0.0;	a[3] = 0.0;	a[4] = 0.0;	a[5] = 1.0;	a[6] = 0.0;	a[7] = 0.0;	a[8] = 0.0; a[9] = 1.0; a[10] = 0.0; a[11] = 0.0; a[12] = 0.0; a[13] = 1.0; a[14] = 0.0; a[15] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i)
				sum_rhs += rhs[i];
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-12);
			// Mode 3 check
			a[0] = 0.0; a[1] = 0.0; a[2] = 1.0;	a[3] = 0.0;	a[4] = 0.0;	a[5] = 0.0;	a[6] = 1.0;	a[7] = 0.0;	a[8] = 0.0; a[9] = 0.0; a[10] = 1.0; a[11] = 0.0; a[12] = 0.0; a[13] = 0.0; a[14] = 1.0; a[15] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i)
				sum_rhs += rhs[i];
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-12);
			// Mode 4 check
			a[0] = 1.0; a[1] = 0.0; a[2] = 1.0;	a[3] = 0.0;	a[4] = 1.0;	a[5] = 0.0;	a[6] = 1.0;	a[7] = 0.0;	a[8] = 1.0; a[9] = 0.0; a[10] = 1.0; a[11] = 0.0; a[12] = 1.0; a[13] = 0.0; a[14] = 1.0; a[15] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i)
				sum_rhs += rhs[i];
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-12);
			// Mode 5 check
			a[0] = 0.0; a[1] = 1.0; a[2] = 1.0;	a[3] = 0.0;	a[4] = 0.0;	a[5] = 1.0;	a[6] = 1.0;	a[7] = 0.0;	a[8] = 0.0; a[9] = 1.0; a[10] = 1.0; a[11] = 0.0; a[12] = 0.0; a[13] = 1.0; a[14] = 1.0; a[15] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i)
				sum_rhs += rhs[i];
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-12);
			// Mode 6 check
			a[0] = 1.0; a[1] = 1.0; a[2] = 0.0;	a[3] = 0.0;	a[4] = 1.0;	a[5] = 1.0;	a[6] = 0.0;	a[7] = 0.0;	a[8] = 1.0; a[9] = 1.0; a[10] = 0.0; a[11] = 0.0; a[12] = 1.0; a[13] = 1.0; a[14] = 0.0; a[15] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i)
				sum_rhs += rhs[i];
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, 1e-12);

		}

	} // namespace Testing
}  // namespace Kratos.
