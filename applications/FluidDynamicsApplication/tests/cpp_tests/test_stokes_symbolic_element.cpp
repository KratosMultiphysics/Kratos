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
#include "custom_elements/symbolic_stokes.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"

namespace Kratos {
	namespace Testing {

        static constexpr double TOLERANCE = 1.0e-12;

		typedef ModelPart::IndexType									 IndexType;
		typedef ModelPart::NodeIterator					          NodeIteratorType;

	    /** Checks the SymbolicStokes2D3N element.
	     * Checks the LHS and RHS computation using a small perturbation.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(ElementSymbolicStokes2D3N, FluidDynamicsApplicationFastSuite)
		{
			Model model;
			ModelPart& r_model_part = model.CreateModelPart("Main", 3);

			// Variables addition
			r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
			r_model_part.AddNodalSolutionStepVariable(PRESSURE);
			r_model_part.AddNodalSolutionStepVariable(VELOCITY);

			// Process info creation
			double delta_time = 0.1;
			r_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.0);
			r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
			Vector bdf_coefs(3);
			bdf_coefs[0] = 3.0/(2.0*delta_time);
			bdf_coefs[1] = -2.0/delta_time;
			bdf_coefs[2] = 0.5*delta_time;
			r_model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

			// Set the element properties
			auto p_elem_prop = r_model_part.CreateNewProperties(1);
			p_elem_prop->SetValue(DENSITY, 1000.0);
			p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
			Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
			p_elem_prop->SetValue(CONSTITUTIVE_LAW, pConsLaw);

			// Geometry creation
			r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
			r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
			r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
			std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3};
			r_model_part.CreateNewElement("SymbolicStokes2D3N", 1, elem_nodes, p_elem_prop);

			auto p_element = r_model_part.pGetElement(1);

			// Define the nodal values
			Matrix vel_original(3,2);
			vel_original(0,0) = 0.0; vel_original(0,1) = 0.1;
			vel_original(1,0) = 0.1; vel_original(1,1) = 0.2;
			vel_original(2,0) = 0.2; vel_original(2,1) = 0.3;

			for(unsigned int i=0; i<3; i++){
				p_element->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE)    = 0.0;
				for(unsigned int k=0; k<2; k++){
					p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel_original(i,k);
					p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i,k);
					p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*vel_original(i,k);
				}
			}

			// Compute RHS and LHS
			Vector RHS = ZeroVector(9);
			Matrix LHS = ZeroMatrix(9,9);

			p_element->Initialize(); // Initialize the element to initialize the constitutive law
			p_element->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

			// Check values
			std::vector<double> RHS_expected({224.531,273.906,-0.0944375,-138.125,98.75,-0.0351875,61.7187,-76.4063,-0.020375});
			std::vector<double> LHS_row_0_expected({1875,625,0.166667,-1.16667e-05,-5e-06,0.166667,625,-625,0.166667});
			KRATOS_CHECK_VECTOR_NEAR(RHS, RHS_expected, 1.0e-2);
			KRATOS_CHECK_VECTOR_NEAR(row(LHS,0), LHS_row_0_expected, 1.0e-2);
	    }

		/** Checks the SymbolicStokes2D3N element.
		 * Checks the LHS and RHS stationary solid rigid movements.
		 */
	    KRATOS_TEST_CASE_IN_SUITE(ElementSymbolicStokes2D3NStationary, FluidDynamicsApplicationFastSuite)
		{
			Model model;
			ModelPart& r_model_part = model.CreateModelPart("Main", 3);

			// Variables addition
			r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
			r_model_part.AddNodalSolutionStepVariable(PRESSURE);
			r_model_part.AddNodalSolutionStepVariable(VELOCITY);

			// Process info creation
			double delta_time = 0.1;
			r_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.0);
			r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
			Vector bdf_coefs(3);
			bdf_coefs[0] = 0.0;
			bdf_coefs[1] = 0.0;
			bdf_coefs[2] = 0.0;
			r_model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

			// Set the element properties
			Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(1);
			p_elem_prop->SetValue(DENSITY, 1.0);
			p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.0);
			Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
			p_elem_prop->SetValue(CONSTITUTIVE_LAW, pConsLaw);

			// Geometry creation
			r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
			r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
			r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
			std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3};
			r_model_part.CreateNewElement("SymbolicStokes2D3N", 1, elem_nodes, p_elem_prop);

			auto p_element = r_model_part.pGetElement(1);

			// Define the nodal values
			array_1d<double, 3> velocity_values;
			velocity_values[0] = 0.0;
			velocity_values[1] = 0.0;
			velocity_values[2] = 0.0;
			double pressure_value = 0.0;

			// Set the nodal values
			for (NodeIteratorType it_node=r_model_part.NodesBegin(); it_node<r_model_part.NodesEnd(); ++it_node){
				it_node->FastGetSolutionStepValue(PRESSURE) = pressure_value;
				it_node->FastGetSolutionStepValue(VELOCITY) = velocity_values;
			}

			// Compute RHS and LHS
			Vector RHS = ZeroVector(9);
			Matrix LHS = ZeroMatrix(9,9);

			p_element->Initialize(); // Initialize the element to initialize the constitutive law
			p_element->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

			// Check obtained RHS
			double sum_RHS = 0.0;
			for (unsigned int i=0; i<RHS.size(); ++i) {
				sum_RHS += RHS[i];
            }
			KRATOS_CHECK_NEAR(sum_RHS, 0.0, TOLERANCE);

			// Check rigid movement modes
			Vector a(9);
			Vector rhs(9);
			double sum_rhs;
			// Mode 1 check
			a[0] = 1.0; a[1] = 0.0; a[2] = 0.0;	a[3] = 1.0;	a[4] = 0.0;	a[5] = 0.0;	a[6] = 1.0;	a[7] = 0.0;	a[8] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i) {
				sum_rhs += rhs[i];
            }
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
			// Mode 2 check
			a[0] = 0.0; a[1] = 1.0; a[2] = 0.0;	a[3] = 0.0;	a[4] = 1.0;	a[5] = 0.0;	a[6] = 0.0;	a[7] = 1.0;	a[8] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i) {
				sum_rhs += rhs[i];
            }
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
			// Mode 3 check
			a[0] = 1.0; a[1] = 1.0; a[2] = 0.0;	a[3] = 1.0;	a[4] = 1.0;	a[5] = 0.0;	a[6] = 1.0;	a[7] = 1.0;	a[8] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i) {
				sum_rhs += rhs[i];
            }
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);

	    }

        /** Checks the SymbolicStokes2D4N element.
         * Checks the LHS and RHS computation using a small perturbation.
         */
        KRATOS_TEST_CASE_IN_SUITE(ElementSymbolicStokes2D4N, FluidDynamicsApplicationFastSuite)
        {
            Model model;
            ModelPart &r_model_part = model.CreateModelPart("Main", 3);

            // Variables addition
            r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
            r_model_part.AddNodalSolutionStepVariable(PRESSURE);
            r_model_part.AddNodalSolutionStepVariable(VELOCITY);

            // Process info creation
            double delta_time = 0.1;
            r_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.0);
            r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0 / (2.0 * delta_time);
            bdf_coefs[1] = -2.0 / delta_time;
            bdf_coefs[2] = 0.5 * delta_time;
            r_model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

            // Set the element properties
            auto p_elem_prop = r_model_part.CreateNewProperties(1);
            p_elem_prop->SetValue(DENSITY, 1000.0);
            p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
            r_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
            r_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
            std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3, 4};
            r_model_part.CreateNewElement("SymbolicStokes2D4N", 1, elem_nodes, p_elem_prop);

            auto p_element = r_model_part.pGetElement(1);

            // Define the nodal values
            Matrix vel_original(4,2);
            vel_original(0,0) = 0.0; vel_original(0,1) = 0.1;
            vel_original(1,0) = 0.1; vel_original(1,1) = 0.2;
            vel_original(2,0) = 0.2; vel_original(2,1) = 0.3;
            vel_original(3,0) = 0.4; vel_original(3,1) = 0.4;

            for (unsigned int i = 0; i < 4; i++)
            {
                p_element->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
                for (unsigned int k = 0; k < 2; k++)
                {
                    p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k] = vel_original(i, k);
                    p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9 * vel_original(i, k);
                    p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75 * vel_original(i, k);
                }
            }

            // Compute RHS and LHS
            Vector RHS = ZeroVector(12);
            Matrix LHS = ZeroMatrix(12, 12);

            p_element->Initialize(); // Initialize the element to initialize the constitutive law
            p_element->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

            // Check values
            std::vector<double> RHS_expected({161.2500018, 213.6458374, -0.1151145824, 36.25000067, 201.0416671, -0.05894791639, 116.8749998, 169.2708304, 0.05245833223, 204.0624977, 156.6666651, -0.02839583341});
            std::vector<double> LHS_row_0_expected({1875.000011, 156.2500033, 0.1666666667, 624.9999939, 156.2499983, 0.1666666667, 312.4999944, -156.2500033, 0.08333333333, 937.5000006, -156.2499983, 0.08333333333});
            KRATOS_CHECK_VECTOR_NEAR(RHS, RHS_expected, 1.0e-2);
        }

        /** Checks the SymbolicStokes2D4N element.
         * Checks the LHS and RHS stationary solid rigid movements.
         */
        KRATOS_TEST_CASE_IN_SUITE(ElementSymbolicStokes2D4NStationary, FluidDynamicsApplicationFastSuite)
        {
            Model model;
            ModelPart &r_model_part = model.CreateModelPart("Main", 3);

            // Variables addition
            r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
            r_model_part.AddNodalSolutionStepVariable(PRESSURE);
            r_model_part.AddNodalSolutionStepVariable(VELOCITY);

            // Process info creation
            double delta_time = 0.1;
            r_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.0);
            r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 0.0;
            bdf_coefs[1] = 0.0;
            bdf_coefs[2] = 0.0;
            r_model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

            // Set the element properties
            Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(1);
            p_elem_prop->SetValue(DENSITY, 1.0);
            p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.0);
            Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
            r_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
            r_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
            std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3, 4};
            r_model_part.CreateNewElement("SymbolicStokes2D4N", 1, elem_nodes, p_elem_prop);

            auto p_element = r_model_part.pGetElement(1);

            // Define the nodal values
            array_1d<double, 3> velocity_values;
            velocity_values[0] = 0.0;
            velocity_values[1] = 0.0;
            velocity_values[2] = 0.0;
            double pressure_value = 0.0;

            // Set the nodal values
            for (auto &node : r_model_part.Nodes())
            {
                node.FastGetSolutionStepValue(PRESSURE) = pressure_value;
                node.FastGetSolutionStepValue(VELOCITY) = velocity_values;
            }

            // Compute RHS and LHS
            Vector RHS = ZeroVector(12);
            Matrix LHS = ZeroMatrix(12, 12);

            p_element->Initialize(); // Initialize the element to initialize the constitutive law
            p_element->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

            // Check obtained RHS
            double sum_RHS = 0.0;
            for (const auto &i_RHS : RHS)
            {
                sum_RHS += i_RHS;
            }
            KRATOS_CHECK_NEAR(sum_RHS, 0.0, TOLERANCE);

            // Check rigid movement modes
            Vector a(12);
            Vector rhs(12);
            double sum_rhs;
            // Mode 1 check
            a[0] = 1.0; a[1] = 0.0; a[2] = 0.0; a[3] = 1.0; a[4] = 0.0; a[5] = 0.0; a[6] = 1.0; a[7] = 0.0; a[8] = 0.0; a[9] = 1.0; a[10] = 0.0; a[11] = 0.0;
            sum_rhs = 0.0;
            rhs = prod(LHS, a);
            for (const auto &i_rhs : rhs)
            {
                sum_rhs += i_rhs;
            }
            KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);

            // Mode 2 check
            a[0] = 0.0; a[1] = 1.0; a[2] = 0.0; a[3] = 0.0; a[4] = 1.0; a[5] = 0.0; a[6] = 0.0; a[7] = 1.0; a[8] = 0.0; a[9] = 0.0; a[10] = 1.0; a[11] = 0.0;
            sum_rhs = 0.0;
            rhs = prod(LHS, a);
            for (const auto &i_rhs : rhs)
            {
                sum_rhs += i_rhs;
            }
            KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);

            // Mode 3 check
            a[0] = 1.0; a[1] = 1.0; a[2] = 0.0; a[3] = 1.0; a[4] = 1.0; a[5] = 0.0; a[6] = 1.0; a[7] = 1.0; a[8] = 0.0; a[9] = 1.0; a[10] = 1.0; a[11] = 0.0;
            sum_rhs = 0.0;
            rhs = prod(LHS, a);
            for (const auto &i_rhs : rhs)
            {
                sum_rhs += i_rhs;
            }
            KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
        }

	    // /** Checks the SymbolicStokes3D4N element.
	    //  * Checks the LHS and RHS computation using a small perturbation.
	    //  */
	    KRATOS_TEST_CASE_IN_SUITE(ElementSymbolicStokes3D4N, FluidDynamicsApplicationFastSuite)
		{
			Model model;
			ModelPart& r_model_part = model.CreateModelPart("Main", 3);

			// Variables addition
			r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
			r_model_part.AddNodalSolutionStepVariable(PRESSURE);
			r_model_part.AddNodalSolutionStepVariable(VELOCITY);

			// Process info creation
			double delta_time = 0.1;
			r_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.0);
			r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
			Vector bdf_coefs(3);
			bdf_coefs[0] = 3.0/(2.0*delta_time);
			bdf_coefs[1] = -2.0/delta_time;
			bdf_coefs[2] = 0.5*delta_time;
			r_model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

			// Set the element properties
			auto p_elem_prop = r_model_part.CreateNewProperties(1);
			p_elem_prop->SetValue(DENSITY, 1000.0);
			p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
			Newtonian3DLaw::Pointer pConsLaw(new Newtonian3DLaw());
			p_elem_prop->SetValue(CONSTITUTIVE_LAW, pConsLaw);

			// Geometry creation
			r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
			r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
			r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
			r_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
			std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3, 4};
			r_model_part.CreateNewElement("SymbolicStokes3D4N", 1, elem_nodes, p_elem_prop);

			auto p_element = r_model_part.pGetElement(1);

			// Define the nodal values
			Matrix vel_original(4,3);
			vel_original(0,0) = 0.0; vel_original(0,1) = 0.1; vel_original(0,2) = 0.2;
			vel_original(1,0) = 0.1; vel_original(1,1) = 0.2; vel_original(1,2) = 0.3;
			vel_original(2,0) = 0.2; vel_original(2,1) = 0.3; vel_original(2,2) = 0.4;
			vel_original(3,0) = 0.3; vel_original(3,1) = 0.4; vel_original(3,2) = 0.5;

			for(unsigned int i=0; i<4; i++){
				p_element->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE)    = 0.0;
				for(unsigned int k=0; k<3; k++){
					p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel_original(i,k);
					p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i,k);
					p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*vel_original(i,k);
				}
			}

			// Compute RHS and LHS
			Vector RHS = ZeroVector(16);
			Matrix LHS = ZeroMatrix(16,16);

			p_element->Initialize(); // Initialize the element to initialize the constitutive law
			p_element->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

			// Check values
			std::vector<double> RHS_expected({98.1458,110.49,122.833,-0.0620312,-66.0521,29.625,41.9687,-0.0175938,19.75,-51.2396,44.4375,-0.0126563,22.2187,34.5625,-36.4271,-0.00771875});
			std::vector<double> LHS_row_0_expected({388.889,138.889,138.889,0.0416667,-13.8889,-1.66667e-06,-1.66667e-06,0.0416667,125,-138.889,0,0.0416667,125,0,-138.889,0.0416667});
			KRATOS_CHECK_VECTOR_NEAR(RHS, RHS_expected, 1.0e-2);
			KRATOS_CHECK_VECTOR_NEAR(row(LHS,0), LHS_row_0_expected, 1.0e-2);
		}

		KRATOS_TEST_CASE_IN_SUITE(ElementSymbolicStokes3D4NStationary, FluidDynamicsApplicationFastSuite)
		{
			Model model;
			ModelPart& r_model_part = model.CreateModelPart("Main", 3);

			// Variables addition
			r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
			r_model_part.AddNodalSolutionStepVariable(PRESSURE);
			r_model_part.AddNodalSolutionStepVariable(VELOCITY);

			// Process info creation
			double delta_time = 0.1;
			r_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.0);
			r_model_part.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+03);
			r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
			Vector bdf_coefs(3);
			bdf_coefs[0] = 0.0;
			bdf_coefs[1] = 0.0;
			bdf_coefs[2] = 0.0;
			r_model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

			// Set the element properties
			auto p_elem_prop = r_model_part.CreateNewProperties(1);
			p_elem_prop->SetValue(DENSITY, 1.0);
			p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.0);
			Newtonian3DLaw::Pointer pConsLaw(new Newtonian3DLaw());
			p_elem_prop->SetValue(CONSTITUTIVE_LAW, pConsLaw);

			// Geometry creation
			r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
			r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
			r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
			r_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
			std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3, 4};
			r_model_part.CreateNewElement("SymbolicStokes3D4N", 1, elem_nodes, p_elem_prop);

			auto p_element = r_model_part.pGetElement(1);

			// Define the nodal values
			array_1d<double, 3> velocity_values;
			velocity_values[0] = 0.0;
			velocity_values[1] = 0.0;
			velocity_values[2] = 0.0;
			double pressure_value = 0.0;

			// Set the nodal values
			for (NodeIteratorType it_node=r_model_part.NodesBegin(); it_node<r_model_part.NodesEnd(); ++it_node){
				it_node->FastGetSolutionStepValue(PRESSURE) = pressure_value;
				it_node->FastGetSolutionStepValue(VELOCITY) = velocity_values;
			}

			// Compute RHS and LHS
			Vector RHS = ZeroVector(16);
			Matrix LHS = ZeroMatrix(16,16);

			p_element->Initialize(); // Initialize the element to initialize the constitutive law
			p_element->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

			// Check obtained RHS
			double sum_RHS = 0.0;
			for (unsigned int i=0; i<RHS.size(); ++i) {
				sum_RHS += RHS[i];
            }
			KRATOS_CHECK_NEAR(sum_RHS, 0.0, TOLERANCE);

			// Check modes
			Vector a(16);
			Vector rhs(16);
			double sum_rhs;
			// Mode 1 check
			a[0] = 1.0; a[1] = 0.0; a[2] = 0.0;	a[3] = 0.0;	a[4] = 1.0;	a[5] = 0.0;	a[6] = 0.0;	a[7] = 0.0;	a[8] = 1.0; a[9] = 0.0; a[10] = 0.0; a[11] = 0.0; a[12] = 1.0; a[13] = 0.0; a[14] = 0.0; a[15] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i) {
				sum_rhs += rhs[i];
            }
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
			// Mode 2 check
			a[0] = 0.0; a[1] = 1.0; a[2] = 0.0;	a[3] = 0.0;	a[4] = 0.0;	a[5] = 1.0;	a[6] = 0.0;	a[7] = 0.0;	a[8] = 0.0; a[9] = 1.0; a[10] = 0.0; a[11] = 0.0; a[12] = 0.0; a[13] = 1.0; a[14] = 0.0; a[15] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i) {
				sum_rhs += rhs[i];
            }
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
			// Mode 3 check
			a[0] = 0.0; a[1] = 0.0; a[2] = 1.0;	a[3] = 0.0;	a[4] = 0.0;	a[5] = 0.0;	a[6] = 1.0;	a[7] = 0.0;	a[8] = 0.0; a[9] = 0.0; a[10] = 1.0; a[11] = 0.0; a[12] = 0.0; a[13] = 0.0; a[14] = 1.0; a[15] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i) {
				sum_rhs += rhs[i];
            }
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
			// Mode 4 check
			a[0] = 1.0; a[1] = 0.0; a[2] = 1.0;	a[3] = 0.0;	a[4] = 1.0;	a[5] = 0.0;	a[6] = 1.0;	a[7] = 0.0;	a[8] = 1.0; a[9] = 0.0; a[10] = 1.0; a[11] = 0.0; a[12] = 1.0; a[13] = 0.0; a[14] = 1.0; a[15] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i) {
				sum_rhs += rhs[i];
            }
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
			// Mode 5 check
			a[0] = 0.0; a[1] = 1.0; a[2] = 1.0;	a[3] = 0.0;	a[4] = 0.0;	a[5] = 1.0;	a[6] = 1.0;	a[7] = 0.0;	a[8] = 0.0; a[9] = 1.0; a[10] = 1.0; a[11] = 0.0; a[12] = 0.0; a[13] = 1.0; a[14] = 1.0; a[15] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i) {
				sum_rhs += rhs[i];
            }
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
			// Mode 6 check
			a[0] = 1.0; a[1] = 1.0; a[2] = 0.0;	a[3] = 0.0;	a[4] = 1.0;	a[5] = 1.0;	a[6] = 0.0;	a[7] = 0.0;	a[8] = 1.0; a[9] = 1.0; a[10] = 0.0; a[11] = 0.0; a[12] = 1.0; a[13] = 1.0; a[14] = 0.0; a[15] = 0.0;
			sum_rhs = 0.0;
			rhs = prod(LHS,a);
			for (unsigned int i=0; i<rhs.size(); ++i) {
				sum_rhs += rhs[i];
            }
			KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);

		}

        // /** Checks the SymbolicStokes3D6N element.
        //  * Checks the LHS and RHS computation using a small perturbation.
        //  */
        KRATOS_TEST_CASE_IN_SUITE(ElementSymbolicStokes3D6N, FluidDynamicsApplicationFastSuite)
        {
            Model model;
            ModelPart &r_model_part = model.CreateModelPart("Main", 3);

            // Variables addition
            r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
            r_model_part.AddNodalSolutionStepVariable(PRESSURE);
            r_model_part.AddNodalSolutionStepVariable(VELOCITY);

            // Process info creation
            double delta_time = 0.1;
            r_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.0);
            r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0 / (2.0 * delta_time);
            bdf_coefs[1] = -2.0 / delta_time;
            bdf_coefs[2] = 0.5 * delta_time;
            r_model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

            // Set the element properties
            auto p_elem_prop = r_model_part.CreateNewProperties(1);
            p_elem_prop->SetValue(DENSITY, 1000.0);
            p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            Newtonian3DLaw::Pointer pConsLaw(new Newtonian3DLaw());
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
            r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
            r_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
            r_model_part.CreateNewNode(5, 1.0, 0.0, 1.0);
            r_model_part.CreateNewNode(6, 0.0, 1.0, 1.0);
            std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3, 4, 5, 6};
            r_model_part.CreateNewElement("SymbolicStokes3D6N", 1, elem_nodes, p_elem_prop);

            auto p_element = r_model_part.pGetElement(1);

            // Define the nodal values
            Matrix vel_original(6,3);
            vel_original(0,0) = 0.0; vel_original(0,1) = 0.1; vel_original(0,2) = 0.2;
            vel_original(1,0) = 0.1; vel_original(1,1) = 0.2; vel_original(1,2) = 0.3;
            vel_original(2,0) = 0.2; vel_original(2,1) = 0.3; vel_original(2,2) = 0.4;
            vel_original(3,0) = 0.3; vel_original(3,1) = 0.4; vel_original(3,2) = 0.5;
            vel_original(4,0) = 0.4; vel_original(4,1) = 0.5; vel_original(4,2) = 0.6;
            vel_original(5,0) = 0.5; vel_original(5,1) = 0.6; vel_original(5,2) = 0.7;

            for (unsigned int i = 0; i < 6; i++)
            {
                p_element->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
                for (unsigned int k = 0; k < 3; k++)
                {
                    p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k] = vel_original(i, k);
                    p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9 * vel_original(i, k);
                    p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75 * vel_original(i, k);
                }
            }

            // Compute RHS and LHS
            Vector RHS = ZeroVector(24);
            Matrix LHS = ZeroMatrix(24, 24);

            p_element->Initialize(); // Initialize the element to initialize the constitutive law
            p_element->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

            // Check values
            std::vector<double> RHS_expected({230.7031274, 255.3906281, 217.5781286, -0.1080156245, -138.1250003, 74.06250008, 223.7500003, -0.05740624994, 55.54687492, -107.2656257, 229.9218751, -0.05123437499, 255.3906261, 280.0781264, -7.734374083, -0.08085937475, -113.4375017, 98.74999842, -1.562502333, -0.005562500355, 80.23437358, -82.57812733, 4.609372417, 0.003078124575});
            std::vector<double> LHS_row_0_expected({625.0000064, 208.3333356, 104.1666678, 0.05555555556, -3.472222214e-06, -1.666666667e-06, 104.1666661, 0.05555555556, 208.3333321, -208.3333339, 104.1666669, 0.05555555556, 312.5000019, 104.1666678, -104.1666661, 0.02777777778, -2.361111063e-06, -8.333333333e-07, -104.1666678, 0.02777777778, 104.1666654, -104.1666669, -104.1666669, 0.02777777778});
            KRATOS_CHECK_VECTOR_NEAR(RHS, RHS_expected, 1.0e-2);
            KRATOS_CHECK_VECTOR_NEAR(row(LHS, 0), LHS_row_0_expected, 1.0e-2);
        }

        KRATOS_TEST_CASE_IN_SUITE(ElementSymbolicStokes3D6NStationary, FluidDynamicsApplicationFastSuite)
        {
            Model model;
            ModelPart &r_model_part = model.CreateModelPart("Main", 3);

            // Variables addition
            r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
            r_model_part.AddNodalSolutionStepVariable(PRESSURE);
            r_model_part.AddNodalSolutionStepVariable(VELOCITY);

            // Process info creation
            double delta_time = 0.1;
            r_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.0);
            r_model_part.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+03);
            r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 0.0;
            bdf_coefs[1] = 0.0;
            bdf_coefs[2] = 0.0;
            r_model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

            // Set the element properties
            auto p_elem_prop = r_model_part.CreateNewProperties(1);
            p_elem_prop->SetValue(DENSITY, 1.0);
            p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.0);
            Newtonian3DLaw::Pointer pConsLaw(new Newtonian3DLaw());
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
            r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
            r_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
            r_model_part.CreateNewNode(5, 1.0, 0.0, 1.0);
            r_model_part.CreateNewNode(6, 0.0, 1.0, 1.0);
            std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3, 4, 5, 6};
            r_model_part.CreateNewElement("SymbolicStokes3D6N", 1, elem_nodes, p_elem_prop);

            auto p_element = r_model_part.pGetElement(1);

            // Define the nodal values
            array_1d<double, 3> velocity_values;
            velocity_values[0] = 0.0;
            velocity_values[1] = 0.0;
            velocity_values[2] = 0.0;
            double pressure_value = 0.0;

            // Set the nodal values
            for (auto &node : r_model_part.Nodes())
            {
                node.FastGetSolutionStepValue(PRESSURE) = pressure_value;
                node.FastGetSolutionStepValue(VELOCITY) = velocity_values;
            }

            // Compute RHS and LHS
            Vector RHS = ZeroVector(24);
            Matrix LHS = ZeroMatrix(24, 24);

            p_element->Initialize(); // Initialize the element to initialize the constitutive law
            p_element->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

            // Check obtained RHS
            double sum_RHS = 0.0;
            for (auto &i_RHS : RHS)
            {
                sum_RHS += i_RHS;
            }
            KRATOS_CHECK_NEAR(sum_RHS, 0.0, TOLERANCE);
            // Check modes
            Vector a(24);
            Vector rhs(24);
            double sum_rhs;
            // Mode 1 check
            for (unsigned int i = 0; i < a.size(); ++i)
            {
                a[i] = (i % 4 == 0) ? 1.0 : 0.0;
            }
            sum_rhs = 0.0;
            rhs = prod(LHS, a);
            for (const auto &i_rhs : rhs)
            {
                sum_rhs += i_rhs;
            }
            KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
            // Mode 2 check
            for (unsigned int i = 0; i < a.size(); ++i)
            {
                a[i] = (i % 4 == 1) ? 1.0 : 0.0;
            }
            sum_rhs = 0.0;
            rhs = prod(LHS, a);
            for (const auto &i_rhs : rhs)
            {
                sum_rhs += i_rhs;
            }
            KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
            // Mode 3 check
            for (unsigned int i = 0; i < a.size(); ++i)
            {
                a[i] = (i % 4 == 2) ? 1.0 : 0.0;
            }
            sum_rhs = 0.0;
            rhs = prod(LHS, a);
            for (const auto &i_rhs : rhs)
            {
                sum_rhs += i_rhs;
            }
            KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
            // Mode 4 check
            for (unsigned int i = 0; i < a.size(); ++i)
            {
                a[i] = (i % 4 == 0 || i % 4 == 2) ? 1.0 : 0.0;
            }
            sum_rhs = 0.0;
            rhs = prod(LHS, a);
            for (const auto &i_rhs : rhs)
            {
                sum_rhs += i_rhs;
            }
            KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
            // Mode 5 check
            for (unsigned int i = 0; i < a.size(); ++i)
            {
                a[i] = (i % 4 == 1 || i % 4 == 2) ? 1.0 : 0.0;
            }
            sum_rhs = 0.0;
            rhs = prod(LHS, a);
            for (const auto &i_rhs : rhs)
            {
                sum_rhs += i_rhs;
            }
            KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
            // Mode 6 check
            for (unsigned int i = 0; i < a.size(); ++i)
            {
                a[i] = (i % 4 == 0 || i % 4 == 1) ? 1.0 : 0.0;
            }
            sum_rhs = 0.0;
            rhs = prod(LHS, a);
            for (const auto &i_rhs : rhs)
            {
                sum_rhs += i_rhs;
            }
            KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
        }

        // /** Checks the SymbolicStokes3D8N element.
        //  * Checks the LHS and RHS computation using a small perturbation.
        //  */
        KRATOS_TEST_CASE_IN_SUITE(ElementSymbolicStokes3D8N, FluidDynamicsApplicationFastSuite)
        {
            Model model;
            ModelPart &r_model_part = model.CreateModelPart("Main", 3);

            // Variables addition
            r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
            r_model_part.AddNodalSolutionStepVariable(PRESSURE);
            r_model_part.AddNodalSolutionStepVariable(VELOCITY);

            // Process info creation
            double delta_time = 0.1;
            r_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.0);
            r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 3.0 / (2.0 * delta_time);
            bdf_coefs[1] = -2.0 / delta_time;
            bdf_coefs[2] = 0.5 * delta_time;
            r_model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

            // Set the element properties
            auto p_elem_prop = r_model_part.CreateNewProperties(1);
            p_elem_prop->SetValue(DENSITY, 1000.0);
            p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            Newtonian3DLaw::Pointer pConsLaw(new Newtonian3DLaw());
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
            r_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
            r_model_part.CreateNewNode(4, 0.0, 1.0, 1.0);
            r_model_part.CreateNewNode(5, 0.0, 0.0, 1.0);
            r_model_part.CreateNewNode(6, 1.0, 0.0, 1.0);
            r_model_part.CreateNewNode(7, 1.0, 1.0, 1.0);
            r_model_part.CreateNewNode(8, 0.0, 1.0, 1.0);
            std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3, 4, 5, 6, 7, 8};
            r_model_part.CreateNewElement("SymbolicStokes3D8N", 1, elem_nodes, p_elem_prop);

            auto p_element = r_model_part.pGetElement(1);

            // Define the nodal values
            Matrix vel_original(8, 3);
            vel_original(0, 0) = 0.0; vel_original(0, 1) = 0.1; vel_original(0, 2) = 0.2;
            vel_original(1, 0) = 0.1; vel_original(1, 1) = 0.2; vel_original(1, 2) = 0.3;
            vel_original(2, 0) = 0.2; vel_original(2, 1) = 0.3; vel_original(2, 2) = 0.4;
            vel_original(3, 0) = 0.3; vel_original(3, 1) = 0.4; vel_original(3, 2) = 0.5;
            vel_original(4, 0) = 0.4; vel_original(4, 1) = 0.5; vel_original(4, 2) = 0.6;
            vel_original(5, 0) = 0.5; vel_original(5, 1) = 0.6; vel_original(5, 2) = 0.7;
            vel_original(6, 0) = 0.6; vel_original(6, 1) = 0.7; vel_original(6, 2) = 0.8;
            vel_original(7, 0) = 0.7; vel_original(7, 1) = 0.8; vel_original(7, 2) = 0.9;

            for (unsigned int i = 0; i < 8; i++)
            {
                p_element->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
                for (unsigned int k = 0; k < 3; k++)
                {
                    p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k] = vel_original(i, k);
                    p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9 * vel_original(i, k);
                    p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75 * vel_original(i, k);
                }
            }

            // Compute RHS and LHS
            Vector RHS = ZeroVector(32);
            Matrix LHS = ZeroMatrix(32, 32);

            p_element->Initialize(); // Initialize the element to initialize the constitutive law
            p_element->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

            // Check values
            std::vector<double> RHS_expected({312.2168837, 190.8386773, 389.472495, -0.1311111106, -67.20285686, 280.1368874, 383.5463436, -0.1183055552, -7.053952332, -142.4946575, 419.4204074, -0.08101736106, 386.3828838, -226.9384326, 468.0936225, -0.1025659719, 255.1268689, 410.6717412, -101.4516571, -0.08564583323, -117.4846442, 397.8839468, -42.03592709, -0.04259722243, -157.8939673, 11.7134052, -90.25374564, 0.03069444374, 144.7629509, 49.23009891, -233.562372, -0.0194513892});
            std::vector<double> LHS_row_0_expected({635.2831243, 50.48076971, 146.0336552, 0.05555555556, 109.6420927, 86.53846153, 129.8076919, 0.05555555556, 53.68589666, -137.0192312, 75.72115373, 0.02777777778, 292.6014983, -165.8653845, 102.7644241, 0.02777777778, 276.1752134, 125.3004811, -146.0336536, 0.02777777778, 20.56623719, 89.24278847, -129.8076935, 0.02777777778, -14.62339902, -38.7620196, -75.72115456, 0.01388888889, 85.00266994, -9.915865479, -102.7644232, 0.01388888889});
            KRATOS_CHECK_VECTOR_NEAR(RHS, RHS_expected, 1.0e-2);
            KRATOS_CHECK_VECTOR_NEAR(row(LHS, 0), LHS_row_0_expected, 1.0e-2);
        }

        KRATOS_TEST_CASE_IN_SUITE(ElementSymbolicStokes3D8NStationary, FluidDynamicsApplicationFastSuite)
        {
            Model model;
            ModelPart &r_model_part = model.CreateModelPart("Main", 3);

            // Variables addition
            r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
            r_model_part.AddNodalSolutionStepVariable(PRESSURE);
            r_model_part.AddNodalSolutionStepVariable(VELOCITY);

            // Process info creation
            double delta_time = 0.1;
            r_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.0);
            r_model_part.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+03);
            r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            Vector bdf_coefs(3);
            bdf_coefs[0] = 0.0;
            bdf_coefs[1] = 0.0;
            bdf_coefs[2] = 0.0;
            r_model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

            // Set the element properties
            auto p_elem_prop = r_model_part.CreateNewProperties(1);
            p_elem_prop->SetValue(DENSITY, 1.0);
            p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.0);
            Newtonian3DLaw::Pointer pConsLaw(new Newtonian3DLaw());
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
            r_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
            r_model_part.CreateNewNode(4, 0.0, 1.0, 1.0);
            r_model_part.CreateNewNode(5, 0.0, 0.0, 1.0);
            r_model_part.CreateNewNode(6, 1.0, 0.0, 1.0);
            r_model_part.CreateNewNode(7, 1.0, 1.0, 1.0);
            r_model_part.CreateNewNode(8, 0.0, 1.0, 1.0);
            std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3, 4, 5, 6, 7, 8};
            r_model_part.CreateNewElement("SymbolicStokes3D8N", 1, elem_nodes, p_elem_prop);

            auto p_element = r_model_part.pGetElement(1);

            // Define the nodal values
            array_1d<double, 3> velocity_values;
            velocity_values[0] = 0.0;
            velocity_values[1] = 0.0;
            velocity_values[2] = 0.0;
            double pressure_value = 0.0;

            // Set the nodal values
            for (auto &node : r_model_part.Nodes())
            {
                node.FastGetSolutionStepValue(PRESSURE) = pressure_value;
                node.FastGetSolutionStepValue(VELOCITY) = velocity_values;
            }

            // Compute RHS and LHS
            Vector RHS = ZeroVector(32);
            Matrix LHS = ZeroMatrix(32, 32);

            p_element->Initialize(); // Initialize the element to initialize the constitutive law
            p_element->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

            // Check obtained RHS
            double sum_RHS = 0.0;
            for (const auto &i_RHS : RHS)
            {
                sum_RHS += i_RHS;
            }
            KRATOS_CHECK_NEAR(sum_RHS, 0.0, TOLERANCE);
            // Check modes
            Vector a(32);
            Vector rhs(32);
            double sum_rhs;
            // Mode 1 check
            for (unsigned int i = 0; i < a.size(); ++i)
            {
                a[i] = (i % 4 == 0) ? 1.0 : 0.0;
            }
            sum_rhs = 0.0;
            rhs = prod(LHS, a);
            for (const auto &i_rhs : rhs)
            {
                sum_rhs += i_rhs;
            }
            KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
            // Mode 2 check
            for (unsigned int i = 0; i < a.size(); ++i)
            {
                a[i] = (i % 4 == 1) ? 1.0 : 0.0;
            }
            sum_rhs = 0.0;
            rhs = prod(LHS, a);
            for (const auto &i_rhs : rhs)
            {
                sum_rhs += i_rhs;
            }
            KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
            // Mode 3 check
            for (unsigned int i = 0; i < a.size(); ++i)
            {
                a[i] = (i % 4 == 2) ? 1.0 : 0.0;
            }
            sum_rhs = 0.0;
            rhs = prod(LHS, a);
            for (const auto &i_rhs : rhs)
            {
                sum_rhs += i_rhs;
            }
            KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
            // Mode 4 check
            for (unsigned int i = 0; i < a.size(); ++i)
            {
                a[i] = (i % 4 == 0 || i % 4 == 2) ? 1.0 : 0.0;
            }
            sum_rhs = 0.0;
            rhs = prod(LHS, a);
            for (const auto &i_rhs : rhs)
            {
                sum_rhs += i_rhs;
            }
            KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
            // Mode 5 check
            for (unsigned int i = 0; i < a.size(); ++i)
            {
                a[i] = (i % 4 == 1 || i % 4 == 2) ? 1.0 : 0.0;
            }
            sum_rhs = 0.0;
            rhs = prod(LHS, a);
            for (const auto &i_rhs : rhs)
            {
                sum_rhs += i_rhs;
            }
            KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
            // Mode 6 check
            for (unsigned int i = 0; i < a.size(); ++i)
            {
                a[i] = (i % 4 == 0 || i % 4 == 1) ? 1.0 : 0.0;
            }
            sum_rhs = 0.0;
            rhs = prod(LHS, a);
            for (const auto &i_rhs : rhs)
            {
                sum_rhs += i_rhs;
            }
            KRATOS_CHECK_NEAR(sum_rhs, 0.0, TOLERANCE);
        }

	} // namespace Testing
}  // namespace Kratos.
