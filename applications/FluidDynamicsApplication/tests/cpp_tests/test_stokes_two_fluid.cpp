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
#include "containers/model.h"
#include "includes/properties.h"
#include "includes/model_part.h"
#include "custom_elements/stokes_3D_twofluid.h"
#include "custom_constitutive/newtonian_two_fluid_3d_law.h"

namespace Kratos {
	namespace Testing {

		typedef ModelPart::IndexType									 IndexType;
		typedef ModelPart::NodeIterator					          NodeIteratorType;

        void PrepareModelPart(ModelPart& modelPart)
        {
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
            bdf_coefs[0] = 3.0 / (2.0*delta_time);
            bdf_coefs[1] = -2.0 / delta_time;
            bdf_coefs[2] = 0.5*delta_time;
            modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
			pElemProp->SetValue(DENSITY_AIR, 1.0);
            NewtonianTwoFluid3DLaw::Pointer pConsLaw(new NewtonianTwoFluid3DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
            std::vector<ModelPart::IndexType> elemNodes{ 1, 2, 3, 4 };
            modelPart.CreateNewElement("StokesTwoFluid3D4N", 1, elemNodes, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);

            // Define the nodal values
            Matrix vel_original(4, 3);
            vel_original(0, 0) = 0.0; vel_original(0, 1) = 0.1; vel_original(0, 2) = 0.2;
            vel_original(1, 0) = 0.1; vel_original(1, 1) = 0.2; vel_original(1, 2) = 0.3;
            vel_original(2, 0) = 0.2; vel_original(2, 1) = 0.3; vel_original(2, 2) = 0.4;
            vel_original(3, 0) = 0.3; vel_original(3, 1) = 0.4; vel_original(3, 2) = 0.5;

            // Set the nodal BODY_FORCE, DENSITY and DYNAMIC_VISCOSITY values
            for (NodeIteratorType it_node = modelPart.NodesBegin(); it_node < modelPart.NodesEnd(); ++it_node) {
                it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = -9.81;
            }

            for (unsigned int i = 0; i < 4; i++) {
                pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
                for (unsigned int k = 0; k < 3; k++) {
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k] = vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
                }
            }
        }

	    // /** Checks the StokesTwoFluid3D4N element
	    //  * Checks the LHS and RHS for a cut element
	    //  */
	    KRATOS_TEST_CASE_IN_SUITE(ElementStokesTwoFluidCut3D4N, FluidDynamicsApplicationFastSuite)
		{
			Model current_model;
			ModelPart& modelPart = current_model.CreateModelPart("Main");
            PrepareModelPart(modelPart);
            
            auto pElement = modelPart.pGetElement(1);

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
            KRATOS_CHECK_NEAR(RHS(0), -8.86905, 1e-2);
            KRATOS_CHECK_NEAR(RHS(1), -0.37064, 1e-2);
            KRATOS_CHECK_NEAR(RHS(2), 7.71903, 1e-2);
            KRATOS_CHECK_NEAR(RHS(3), -318569, 1.0);
            KRATOS_CHECK_NEAR(RHS(4), 22.639, 1e-2);
            KRATOS_CHECK_NEAR(RHS(5), 8.65338, 1e-2);
            KRATOS_CHECK_NEAR(RHS(6), 12.1078, 1e-2);
            KRATOS_CHECK_NEAR(RHS(7), 2866.56, 1e-2);
            KRATOS_CHECK_NEAR(RHS(8), 12.9595, 1e-2);
            KRATOS_CHECK_NEAR(RHS(9), 39.2958, 1e-2);
            KRATOS_CHECK_NEAR(RHS(10), 29.5256, 1e-2);
            KRATOS_CHECK_NEAR(RHS(11), 316699, 1.0);
            KRATOS_CHECK_NEAR(RHS(12), 5.71452, 1e-2);
            KRATOS_CHECK_NEAR(RHS(13), 9.57769, 1e-2);
            KRATOS_CHECK_NEAR(RHS(14), 30.8809, 1e-2);
            KRATOS_CHECK_NEAR(RHS(15), -996.765, 1e-2);
		}

        // /** Checks the StokesTwoFluid3D4N element
        //  * Checks the LHS and RHS for a negative element (distance <= 0.0)
        //  */
        KRATOS_TEST_CASE_IN_SUITE(ElementStokesTwoFluidNegativeSide3D4N, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            PrepareModelPart(modelPart);
            
            auto pElement = modelPart.pGetElement(1);

            pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) = -1.0;

            // Compute RHS and LHS
            Vector RHS = ZeroVector(16);
            Matrix LHS = ZeroMatrix(16, 16);

			const auto& r_process_info = modelPart.GetProcessInfo();
			pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
			pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

            // Check the RHS values (the RHS is computed as the LHS x previous_solution,
            // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
            KRATOS_CHECK_NEAR(RHS(0), 14.8262, 1e-2);
            KRATOS_CHECK_NEAR(RHS(1), 27.1782, 1e-2);
            KRATOS_CHECK_NEAR(RHS(2), 39.1214, 1e-2);
            KRATOS_CHECK_NEAR(RHS(3), -2.04821e+06, 10.0);
            KRATOS_CHECK_NEAR(RHS(4), 17.2867, 1e-2);
            KRATOS_CHECK_NEAR(RHS(5), 29.6168, 1e-2);
            KRATOS_CHECK_NEAR(RHS(6), 41.549, 1e-2);
            KRATOS_CHECK_NEAR(RHS(7), 411458, 1.0);
            KRATOS_CHECK_NEAR(RHS(8), 19.7418, 1e-2);
            KRATOS_CHECK_NEAR(RHS(9), 32.0938, 1e-2);
            KRATOS_CHECK_NEAR(RHS(10), 44.015, 1e-2);
            KRATOS_CHECK_NEAR(RHS(11), 685764, 1.0);
            KRATOS_CHECK_NEAR(RHS(12), 22.2078, 1e-2);
            KRATOS_CHECK_NEAR(RHS(13), 34.5488, 1e-2);
            KRATOS_CHECK_NEAR(RHS(14), 46.492, 1e-2);
            KRATOS_CHECK_NEAR(RHS(15), 950986, 1.0);
        }

        // /** Checks the StokesTwoFluid3D4N element
        //  * Checks the LHS and RHS for a positive element (distance > 0.0)
        //  */
        KRATOS_TEST_CASE_IN_SUITE(ElementStokesTwoFluidPositiveSide3D4N, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            
            PrepareModelPart(modelPart);
            
            auto pElement = modelPart.pGetElement(1);

            pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = 1.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = 1.0;
            pElement->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) = 1.0;

            // Compute RHS and LHS
            Vector RHS = ZeroVector(16);
            Matrix LHS = ZeroMatrix(16, 16);

			const auto& r_process_info = modelPart.GetProcessInfo();
			pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
			pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

            // Check the RHS values (the RHS is computed as the LHS x previous_solution,
            // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
            KRATOS_CHECK_NEAR(RHS(0), 0.039501, 1e-2);
            KRATOS_CHECK_NEAR(RHS(1), 0.0600744, 1e-2);
            KRATOS_CHECK_NEAR(RHS(2), -0.328102, 1e-2);
            KRATOS_CHECK_NEAR(RHS(3), 7026.02, 1e-2);
            KRATOS_CHECK_NEAR(RHS(4), 0.0117953, 1e-2);
            KRATOS_CHECK_NEAR(RHS(5), 0.0213953, 1e-2);
            KRATOS_CHECK_NEAR(RHS(6), -0.377754, 1e-2);
            KRATOS_CHECK_NEAR(RHS(7), 411.433, 1e-2);
            KRATOS_CHECK_NEAR(RHS(8), 0.0115203, 1e-2);
            KRATOS_CHECK_NEAR(RHS(9), 0.0211214, 1e-2);
            KRATOS_CHECK_NEAR(RHS(10), -0.378029, 1e-2);
            KRATOS_CHECK_NEAR(RHS(11), 685.739, 1e-2);
            KRATOS_CHECK_NEAR(RHS(12), 0.0112459, 1e-2);
            KRATOS_CHECK_NEAR(RHS(13), 0.0208464, 1e-2);
            KRATOS_CHECK_NEAR(RHS(14), -0.378303, 1e-2);
            KRATOS_CHECK_NEAR(RHS(15), -8123.29, 1e-2);
        }

	} // namespace Testing
}  // namespace Kratos.