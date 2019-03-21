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
#include "containers/model.h"
#include "includes/model_part.h"
#include "custom_elements/time_averaged_navier_stokes.h"
#include "custom_conditions/time_averaged_navier_stokes_wall_condition.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"

namespace Kratos {
	namespace Testing {

		typedef ModelPart::IndexType									 IndexType;
		typedef ModelPart::NodeIterator					          NodeIteratorType;

	    /** Checks the TimeAveragedNavierStokes2D3N element.
	     * Checks the LHS and RHS computation using a small perturbation.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(ElementTimeAveragedNavierStokes2D3N, FluidDynamicsApplicationFastSuite)
		{			
			Model model;
			unsigned int buffer_size = 5;
			double exp_factor = 1.0;
			ModelPart& modelPart = model.CreateModelPart("Main",buffer_size);

			// Variables addition
			modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
			modelPart.AddNodalSolutionStepVariable(DENSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
			modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
			modelPart.AddNodalSolutionStepVariable(SOUND_VELOCITY);
			modelPart.AddNodalSolutionStepVariable(TIME_AVERAGED_PRESSURE);
			modelPart.AddNodalSolutionStepVariable(TIME_AVERAGED_VELOCITY);
			modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);

			// Process info creation
			double delta_time = 0.1;
			modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.0);
			modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
			// Initializing time and delta time
			modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
			modelPart.GetProcessInfo().SetValue(TIME,0.0);
			// Initializing bdf coefficients
			ComputeBDFCoefficientsProcess bdf_process(modelPart,2);
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
			modelPart.CreateNewElement("TimeAveragedNavierStokes2D3N", 1, elemNodes, pElemProp);

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

			// Initializing Solution Step Values
			for(unsigned int i=0; i<3; i++){
				pElement->GetGeometry()[i].FastGetSolutionStepValue(TIME_AVERAGED_PRESSURE)    = 0.0;
				for(unsigned int k=0; k<2; k++){
					pElement->GetGeometry()[i].FastGetSolutionStepValue(TIME_AVERAGED_VELOCITY)[k]    = vel_original(i,k);
					pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
				}
			}

			// Fake time advance
			double current_time = 0.0;
			int time_step = 0;
			for (unsigned int i = 0; i < buffer_size; ++i) {
				modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
				modelPart.GetProcessInfo().SetValue(TIME,0.0);
				modelPart.CloneSolutionStep();
			}
			KRATOS_WATCH("BUFFER IS FILLED")

			// Compute RHS and LHS
			Vector RHS = ZeroVector(9);
			Matrix LHS = ZeroMatrix(9,9);

			pElement->Initialize(); // Initialize the element to initialize the constitutive law
			KRATOS_WATCH("ELEMENT INITIALIZED")
			pElement->CalculateLocalSystem(LHS, RHS, modelPart.GetProcessInfo());
			KRATOS_WATCH("LOCAL SYSTEM CALCULATED")

			// Compute the error of the perturbation
			double perturbation = 2e-2;
			std::vector<double> error_norms;
			std::vector<double> error_vel_norms;
			std::vector<double> error_pres_norms;

			for (unsigned int j=0; j<5; ++j)
			{
				modelPart.CloneTimeStep(current_time);
				time_step += 1;
				KRATOS_WATCH(time_step)
				double previous_time_step = modelPart.GetProcessInfo().GetPreviousTimeStepInfo()[DELTA_TIME];
				KRATOS_WATCH(previous_time_step)
				double next_delta_time = previous_time_step*exp_factor;
				KRATOS_WATCH(next_delta_time)
				current_time += previous_time_step;
				KRATOS_WATCH(current_time)
				modelPart.GetProcessInfo().SetValue(TIME,current_time);
				modelPart.GetProcessInfo().SetValue(DELTA_TIME,next_delta_time);
				
				bdf_process.Execute();
				const Vector& BDFVector = modelPart.GetProcessInfo().GetValue(BDF_COEFFICIENTS);

				perturbation /= 2;
				Vector perturbation_vector = ZeroVector(9);

				for(unsigned int i=0; i<3; i++){
					for(unsigned int k=0; k<2; k++){
						pElement->GetGeometry()[i].FastGetSolutionStepValue(TIME_AVERAGED_VELOCITY)[k] = vel_original(i,k) + perturbation*perturbation_matrix(i,k);
						perturbation_vector(i*3+k) = perturbation*perturbation_matrix(i,k);
					}
				}
				KRATOS_WATCH(perturbation_vector)
				
				// Compute perturbed RHS and LHS
				Vector error = ZeroVector(9);
				Vector RHS_obtained = ZeroVector(9);
				Vector RHS_perturbed = ZeroVector(9);
				Vector solution_increment = ZeroVector(9);
				pElement->CalculateRightHandSide(RHS_perturbed, modelPart.GetProcessInfo());

				KRATOS_WATCH(RHS_perturbed)

				solution_increment = prod(LHS, perturbation_vector);
				noalias(RHS_obtained) = RHS - solution_increment;
				noalias(error) = RHS_perturbed - RHS_obtained;
				error_norms.push_back(norm_2(error));
			}
				
			// // Check quadratic convergence (if FullNR has been selected when generating the element)
			// for(unsigned int i=1; i<error_norms.size(); ++i)
			// 	KRATOS_CHECK_NEAR(error_norms[i-1]/error_norms[i], 4.0, 1e-1);

			//Check quadratic convergence (if Picard has been selected when generating the element)
			// for(unsigned int i=1; i<error_norms.size(); ++i)
			// 	KRATOS_CHECK_NEAR(error_norms[i-1]/error_norms[i], 2.0, 2.5e-1);

			// for(unsigned int i=0;i<error_norms.size();++i){
			// 	std::cout << "Error norm "<< i << " " << error_norms[i] << std::endl;
			// }


	    }
	}

}