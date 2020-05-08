//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Mohit Tyagi
//
//

// System includes

// External includes
// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "spaces/ublas_space.h"
#include "includes/properties.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"
#include "includes/global_pointer_variables.h"
#include "custom_elements/symbolic_stokes.h"
#include "custom_constitutive/newtonian_3d_law.h"

#include "processes/find_nodal_neighbours_process.h"
#include "utilities/normal_calculation_utils.h"
#include "utilities/time_discretization.h"

namespace Kratos
{
    namespace Testing
    {
        typedef ModelPart::IndexType IndexType;
        typedef ModelPart::NodeIterator NodeIteratorType;

        /** Checks the SymbolicStokes3D4N element.
         * Checks the LHS and RHS computation
         */
        KRATOS_TEST_CASE_IN_SUITE(SymbolicStokes3D4N, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart &model_part = current_model.CreateModelPart("model_part");
            model_part.SetBufferSize(3);
            // variable addition
            model_part.AddNodalSolutionStepVariable(BODY_FORCE);
            model_part.AddNodalSolutionStepVariable(DENSITY);
            model_part.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            model_part.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            model_part.AddNodalSolutionStepVariable(PRESSURE);
            model_part.AddNodalSolutionStepVariable(VELOCITY);

            // Process info creation
            double delta_time = 0.1;
            ProcessInfo &process_info = model_part.GetProcessInfo();
            process_info.SetValue(STEP, 0);
            process_info.SetValue(DYNAMIC_TAU, 0.0);
            process_info.SetValue(DELTA_TIME, delta_time);
            process_info.SetValue(START_TIME, 0.0);
            model_part.CloneTimeStep(delta_time);
            process_info[STEP] = 1;
            TimeDiscretization::BDF time_disc_BDF2(2);
            time_disc_BDF2.ComputeAndSaveBDFCoefficients(process_info);

            // Set element properties
            Properties::Pointer p_elem_prop = model_part.CreateNewProperties(0);
            p_elem_prop->SetValue(DENSITY, 1000.0);
            p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1000.0);
            Newtonian3DLaw::Pointer p_const_law(new Newtonian3DLaw());
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_const_law);

            // Geometry creation
            model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
            model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
            model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
            model_part.CreateNewElement("SymbolicStokes3D4N", 1, {1, 2, 3, 4}, p_elem_prop);

            // Set the nodal DENSITY and DYNAMICS_VISCOSITY
            Element::Pointer p_element = model_part.pGetElement(1);
            for (NodeIteratorType it_node = model_part.NodesBegin(); it_node != model_part.NodesEnd(); ++it_node)
            {
                it_node->FastGetSolutionStepValue(DENSITY) = p_elem_prop->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = p_elem_prop->GetValue(DYNAMIC_VISCOSITY);
                if (it_node->Z() == 0.0)
                {
                    it_node->Fix(VELOCITY_X);
                    it_node->Fix(VELOCITY_Y);
                    it_node->Fix(VELOCITY_Z);
                }
            }
            GeometricalObject::NodeType &node4 = p_element->GetGeometry()[3];
            node4.FastGetSolutionStepValue(VELOCITY_X) = 0.1;
            node4.FastGetSolutionStepValue(VELOCITY_Y) = 0.1;
            node4.FastGetSolutionStepValue(VELOCITY_Z) = 0.1;

            // Compute RHS and LHS
            Vector RHS = ZeroVector(16);
            Matrix LHS = ZeroMatrix(16, 16);
            p_element->Initialize();
            p_element->CalculateLocalSystem(LHS, RHS, process_info);
            // Check the RHS values
            KRATOS_CHECK_NEAR(RHS(0), 9.722222285, 1e-07);
            KRATOS_CHECK_NEAR(RHS(1), 9.722222285, 1e-07);
            KRATOS_CHECK_NEAR(RHS(2), 59.72222229, 1e-07);
            KRATOS_CHECK_NEAR(RHS(3), 0.01145833333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(4), -18.05555549, 1e-07);
            KRATOS_CHECK_NEAR(RHS(5), -12.49999994, 1e-07);
            KRATOS_CHECK_NEAR(RHS(6), -29.1666666, 1e-07);
            KRATOS_CHECK_NEAR(RHS(7), -0.009375, 1e-07);
            KRATOS_CHECK_NEAR(RHS(8), -12.49999994, 1e-07);
            KRATOS_CHECK_NEAR(RHS(9), -18.05555549, 1e-07);
            KRATOS_CHECK_NEAR(RHS(10), -29.1666666, 1e-07);
            KRATOS_CHECK_NEAR(RHS(11), -0.009375, 1e-07);
            KRATOS_CHECK_NEAR(RHS(12), -41.66666686, 1e-07);
            KRATOS_CHECK_NEAR(RHS(13), -41.66666686, 1e-07);
            KRATOS_CHECK_NEAR(RHS(14), -63.88888908, 1e-07);
            KRATOS_CHECK_NEAR(RHS(15), -0.009375, 1e-07);
        }

        /** Checks the SymbolicStokes3D6N element.
         * Checks the LHS and RHS computation
         */
        KRATOS_TEST_CASE_IN_SUITE(SymbolicStokes3D6N, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart &model_part = current_model.CreateModelPart("model_part");
            model_part.SetBufferSize(3);
            // variable addition
            model_part.AddNodalSolutionStepVariable(BODY_FORCE);
            model_part.AddNodalSolutionStepVariable(DENSITY);
            model_part.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            model_part.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            model_part.AddNodalSolutionStepVariable(PRESSURE);
            model_part.AddNodalSolutionStepVariable(VELOCITY);

            // Process info creation
            double delta_time = 0.1;
            ProcessInfo &process_info = model_part.GetProcessInfo();
            process_info.SetValue(STEP, 0);
            process_info.SetValue(DYNAMIC_TAU, 0.0);
            process_info.SetValue(DELTA_TIME, delta_time);
            process_info.SetValue(START_TIME, 0.0);
            model_part.CloneTimeStep(delta_time);
            process_info[STEP] = 1;
            TimeDiscretization::BDF time_disc_BDF2(2);
            time_disc_BDF2.ComputeAndSaveBDFCoefficients(process_info);

            // Set element properties
            Properties::Pointer p_elem_prop = model_part.CreateNewProperties(0);
            p_elem_prop->SetValue(DENSITY, 1000.0);
            p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1000.0);
            Newtonian3DLaw::Pointer p_const_law(new Newtonian3DLaw());
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_const_law);

            // Geometry creation
            model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
            model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
            model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
            model_part.CreateNewNode(5, 1.0, 0.0, 1.0);
            model_part.CreateNewNode(6, 0.0, 1.0, 1.0);
            model_part.CreateNewElement("SymbolicStokes3D6N", 1, {1, 2, 3, 4, 5, 6}, p_elem_prop);

            // Set the nodal DENSITY and DYNAMICS_VISCOSITY
            Element::Pointer p_element = model_part.pGetElement(1);
            for (NodeIteratorType it_node = model_part.NodesBegin(); it_node != model_part.NodesEnd(); ++it_node)
            {
                it_node->FastGetSolutionStepValue(DENSITY) = p_elem_prop->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = p_elem_prop->GetValue(DYNAMIC_VISCOSITY);
                if (it_node->Z() == 0.0)
                {
                    it_node->Fix(VELOCITY_X);
                    it_node->Fix(VELOCITY_Y);
                    it_node->Fix(VELOCITY_Z);
                }
                if (it_node->Z() == 1.0)
                {
                    it_node->FastGetSolutionStepValue(VELOCITY_X) = 0.1;
                    it_node->FastGetSolutionStepValue(VELOCITY_Y) = 0.1;
                    it_node->FastGetSolutionStepValue(VELOCITY_Z) = 0.1;
                }
            }

            // Compute RHS and LHS
            Vector RHS = ZeroVector(24);
            Matrix LHS = ZeroMatrix(24, 24);
            p_element->Initialize();
            p_element->CalculateLocalSystem(LHS, RHS, process_info);

            // Check the RHS values
            KRATOS_CHECK_NEAR(RHS(0), -16.66666667, 1e-07);
            KRATOS_CHECK_NEAR(RHS(1), -16.66666667, 1e-07);
            KRATOS_CHECK_NEAR(RHS(2), 47.22222222, 1e-07);
            KRATOS_CHECK_NEAR(RHS(3), 0.03854166667, 1e-07);
            KRATOS_CHECK_NEAR(RHS(4), -33.33333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(5), -25, 1e-07);
            KRATOS_CHECK_NEAR(RHS(6), -27.77777778, 1e-07);
            KRATOS_CHECK_NEAR(RHS(7), -0.008333333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(8), -25, 1e-07);
            KRATOS_CHECK_NEAR(RHS(9), -33.33333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(10), -27.77777778, 1e-07);
            KRATOS_CHECK_NEAR(RHS(11), -0.008333333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(12), -91.66666667, 1e-07);
            KRATOS_CHECK_NEAR(RHS(13), -91.66666667, 1e-07);
            KRATOS_CHECK_NEAR(RHS(14), -72.22222222, 1e-07);
            KRATOS_CHECK_NEAR(RHS(15), 0.03854166667, 1e-07);
            KRATOS_CHECK_NEAR(RHS(16), -108.3333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(17), -100, 1e-07);
            KRATOS_CHECK_NEAR(RHS(18), -147.2222222, 1e-07);
            KRATOS_CHECK_NEAR(RHS(19), -0.05520833333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(20), -100, 1e-07);
            KRATOS_CHECK_NEAR(RHS(21), -108.3333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(22), -147.2222222, 1e-07);
            KRATOS_CHECK_NEAR(RHS(23), -0.05520833333, 1e-07);
        }

        /** Checks the SymbolicStokes3D8N element.
         * Checks the LHS and RHS computation
         */
        KRATOS_TEST_CASE_IN_SUITE(SymbolicStokes3D8N, FluidDynamicsApplicationFastSuite)
        {
                        Model current_model;
            ModelPart &model_part = current_model.CreateModelPart("model_part");
            model_part.SetBufferSize(3);
            // variable addition
            model_part.AddNodalSolutionStepVariable(BODY_FORCE);
            model_part.AddNodalSolutionStepVariable(DENSITY);
            model_part.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            model_part.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            model_part.AddNodalSolutionStepVariable(PRESSURE);
            model_part.AddNodalSolutionStepVariable(VELOCITY);

            // Process info creation
            double delta_time = 0.1;
            ProcessInfo &process_info = model_part.GetProcessInfo();
            process_info.SetValue(STEP, 0);
            process_info.SetValue(DYNAMIC_TAU, 0.0);
            process_info.SetValue(DELTA_TIME, delta_time);
            process_info.SetValue(START_TIME, 0.0);
            model_part.CloneTimeStep(delta_time);
            process_info[STEP] = 1;
            TimeDiscretization::BDF time_disc_BDF2(2);
            time_disc_BDF2.ComputeAndSaveBDFCoefficients(process_info);

            // Set element properties
            Properties::Pointer p_elem_prop = model_part.CreateNewProperties(0);
            p_elem_prop->SetValue(DENSITY, 1000.0);
            p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1000.0);
            Newtonian3DLaw::Pointer p_const_law(new Newtonian3DLaw());
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_const_law);

            // Geometry creation
            model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
            model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
            model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
            model_part.CreateNewNode(5, 0.0, 0.0, 1.0);
            model_part.CreateNewNode(6, 1.0, 0.0, 1.0);
            model_part.CreateNewNode(7, 1.0, 1.0, 1.0);
            model_part.CreateNewNode(8, 0.0, 1.0, 1.0);
            model_part.CreateNewElement("SymbolicStokes3D8N", 1, {1, 2, 3, 4, 5, 6, 7, 8}, p_elem_prop);

            // Set the nodal DENSITY and DYNAMICS_VISCOSITY
            Element::Pointer p_element = model_part.pGetElement(1);
            for (NodeIteratorType it_node = model_part.NodesBegin(); it_node != model_part.NodesEnd(); ++it_node)
            {
                it_node->FastGetSolutionStepValue(DENSITY) = p_elem_prop->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = p_elem_prop->GetValue(DYNAMIC_VISCOSITY);
                if (it_node->Z() == 0.0)
                {
                    it_node->Fix(VELOCITY_X);
                    it_node->Fix(VELOCITY_Y);
                    it_node->Fix(VELOCITY_Z);
                }
                if (it_node->Z() == 1.0)
                {
                    it_node->FastGetSolutionStepValue(VELOCITY_X) = 0.1;
                    it_node->FastGetSolutionStepValue(VELOCITY_Y) = 0.1;
                    it_node->FastGetSolutionStepValue(VELOCITY_Z) = 0.1;
                }
            }

            // Compute RHS and LHS
            Vector RHS = ZeroVector(32);
            Matrix LHS = ZeroMatrix(32, 32);
            p_element->Initialize();
            p_element->CalculateLocalSystem(LHS, RHS, process_info);
            // Check the RHS values
            KRATOS_CHECK_NEAR(RHS(0), -29.16666667, 1e-07);
            KRATOS_CHECK_NEAR(RHS(1), -29.16666667, 1e-07);
            KRATOS_CHECK_NEAR(RHS(2), 45.83333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(3), 0.096875, 1e-07);
            KRATOS_CHECK_NEAR(RHS(4), -45.83333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(5), -29.16666667, 1e-07);
            KRATOS_CHECK_NEAR(RHS(6), -4.166666667, 1e-07);
            KRATOS_CHECK_NEAR(RHS(7), 0.034375, 1e-07);
            KRATOS_CHECK_NEAR(RHS(8), -45.83333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(9), -45.83333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(10), -54.16666667, 1e-07);
            KRATOS_CHECK_NEAR(RHS(11), -0.028125, 1e-07);
            KRATOS_CHECK_NEAR(RHS(12), -29.16666667, 1e-07);
            KRATOS_CHECK_NEAR(RHS(13), -45.83333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(14), -4.166666667, 1e-07);
            KRATOS_CHECK_NEAR(RHS(15), 0.034375, 1e-07);
            KRATOS_CHECK_NEAR(RHS(16), -141.6666667, 1e-07);
            KRATOS_CHECK_NEAR(RHS(17), -141.6666667, 1e-07);
            KRATOS_CHECK_NEAR(RHS(18), -133.3333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(19), 0.065625, 1e-07);
            KRATOS_CHECK_NEAR(RHS(20), -158.3333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(21), -141.6666667, 1e-07);
            KRATOS_CHECK_NEAR(RHS(22), -183.3333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(23), -0.059375, 1e-07);
            KRATOS_CHECK_NEAR(RHS(24), -158.3333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(25), -158.3333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(26), -233.3333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(27), -0.184375, 1e-07);
            KRATOS_CHECK_NEAR(RHS(28), -141.6666667, 1e-07);
            KRATOS_CHECK_NEAR(RHS(29), -158.3333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(30), -183.3333333, 1e-07);
            KRATOS_CHECK_NEAR(RHS(31), -0.059375, 1e-07);
        }


    } // namespace Testing
} // namespace Kratos.