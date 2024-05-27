//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Uxue Chasco
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "spaces/ublas_space.h"
#include "includes/properties.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"
#include "includes/global_pointer_variables.h"
#include "custom_elements/two_fluid_navier_stokes_alpha_method.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"
#include "custom_constitutive/newtonian_two_fluid_2d_law.h"
#include "custom_constitutive/newtonian_two_fluid_3d_law.h"

#include "processes/find_nodal_neighbours_process.h"
#include "utilities/normal_calculation_utils.h"

namespace Kratos {
    namespace Testing {

        typedef ModelPart::IndexType									 IndexType;
        typedef ModelPart::NodeIterator					          NodeIteratorType;

        /** Checks the TwoFluidNavierStokesAlphaMethod2D3N element.
         * Checks the LHS and RHS computation
         */
        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesAlphaMethod2D3N, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            //FIXME: SET TO 2 WHEN A NEW DATA CONTAINER IS CREATED
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
            modelPart.GetProcessInfo().SetValue(SPECTRAL_RADIUS_LIMIT, 0.0);
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.0);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-03);
            auto p_cons_law = Kratos::make_shared<NewtonianTwoFluid2DLaw>();
            pElemProp->SetValue(CONSTITUTIVE_LAW, p_cons_law);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};
            modelPart.CreateNewElement("TwoFluidNavierStokesAlphaMethod2D3N", 1, elemNodes, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);
             // Fake time advance to set the previous ProcessInfo container
            modelPart.CloneSolutionStep();

            // Define the nodal values
            Matrix vel_original(3, 3);
            vel_original(0, 0) = 0.0;
            vel_original(0, 1) = 0.1;
            vel_original(0, 2) = 0.2;
            vel_original(1, 0) = 0.1;
            vel_original(1, 1) = 0.2;
            vel_original(1, 2) = 0.3;
            vel_original(2, 0) = 0.2;
            vel_original(2, 1) = 0.3;
            vel_original(2, 2) = 0.4;

            // Set the nodal BODY_FORCE, DENSITY and DYNAMIC_VISCOSITY values
            for (NodeIteratorType it_node = modelPart.NodesBegin(); it_node < modelPart.NodesEnd(); ++it_node)
            {
                it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = -9.81;
            }

            for (unsigned int i = 0; i < 3; i++)
            {
                pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
                for (unsigned int k = 0; k < 3; k++)
                {
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k] = vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] =0.9*vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                }
            }
            pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -1.0;

            // Compute RHS and LHS
            Vector RHS = ZeroVector(9);
            Matrix LHS = ZeroMatrix(9, 9);

            const auto &r_process_info = modelPart.GetProcessInfo();
            pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
            pElement->CalculateLocalSystem(LHS, RHS, r_process_info);
            // Check the RHS values (the RHS is computed as the LHS x previous_solution,
            // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
            Vector reference_RHS = ZeroVector(9);
            reference_RHS[0] = -14.67064121;
            reference_RHS[1] = -38.30659006;
            reference_RHS[2] =-0.02477803762;
            reference_RHS[3] = -44.16968252;
            reference_RHS[4] = -61.67597243;
            reference_RHS[5] = -0.05915402772;
            reference_RHS[6] =-43.90588458;
            reference_RHS[7] = -79.83781791;
            reference_RHS[8] = -0.06606793466;

            KRATOS_EXPECT_VECTOR_NEAR(reference_RHS, RHS, 1e-2);
        }
        // /** Checks the TwoFluidNavierStokesAlphaMethod3D4N element
        //  * Checks the LHS and RHS for a cut element
        //  */
        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesCutAlphaMethod3D4N, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
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

            modelPart.GetProcessInfo().SetValue(SPECTRAL_RADIUS_LIMIT, 0.0);
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
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
            modelPart.CreateNewElement("TwoFluidNavierStokesAlphaMethod3D4N", 1, elemNodes, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);
            // Fake time advance to set the previous ProcessInfo container
            modelPart.CloneSolutionStep();

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
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.9*vel_original(i,k);
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

            KRATOS_EXPECT_NEAR(RHS(0), -119.5094945, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(1), 13.86753717, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(2), -264.4471963, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(3), -0.008525833469, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(4), 31.81898305, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(5), -23.48989075, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(6), -515.5487839, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(7), 0.1753941423, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(8), 52.86326183, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(9), -58.7570097, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(10), -510.5819917, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(11), -0.03898369627, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(12), 89.23040167, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(13), -30.84418781, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(14), -581.986129, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(15), -0.2278846126, 1e-7);

        }
                // /** Checks the TwoFluidNavierStokesAlphaMethod3D4N element
        //  * Checks the LHS and RHS for a cut element with surface tension
        //  */
        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesSurfaceTensionAlphaMethod3D4N, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(2);

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

            modelPart.GetProcessInfo().SetValue(SPECTRAL_RADIUS_LIMIT, 0.0);
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            modelPart.GetProcessInfo().SetValue(SURFACE_TENSION, true);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            pElemProp->SetValue(SURFACE_TENSION_COEFFICIENT, 1.0);
            NewtonianTwoFluid3DLaw::Pointer pConsLaw(new NewtonianTwoFluid3DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
            std::vector<ModelPart::IndexType> elemNodes {1, 2, 3, 4};
            modelPart.CreateNewElement("TwoFluidNavierStokesAlphaMethod3D4N", 1, elemNodes, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);
            // Fake time advance to set the previous ProcessInfo container
            modelPart.CloneSolutionStep();

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
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
            }
            }

            pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) =  1.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) =  1.0;

            pElement->GetGeometry()[0].SetValue(CURVATURE, 1.0);
            pElement->GetGeometry()[1].SetValue(CURVATURE, 1.0);
            pElement->GetGeometry()[2].SetValue(CURVATURE, 1.0);
            pElement->GetGeometry()[3].SetValue(CURVATURE, 1.0);

            Vector inner_pressure_grad = ZeroVector(3);
            Vector outer_pressure_grad = ZeroVector(3);

            for(unsigned int k=0; k<3; k++){
                inner_pressure_grad[k] = 0.5*k;
                outer_pressure_grad[k] = 1.5*k + 1.0;
            }

            pElement->GetGeometry()[0].SetValue(PRESSURE_GRADIENT, inner_pressure_grad);
            pElement->GetGeometry()[1].SetValue(PRESSURE_GRADIENT, outer_pressure_grad);
            pElement->GetGeometry()[2].SetValue(PRESSURE_GRADIENT, inner_pressure_grad);
            pElement->GetGeometry()[3].SetValue(PRESSURE_GRADIENT, outer_pressure_grad);

            // Compute RHS and LHS
            Vector RHS = ZeroVector(16);
            Vector reference_RHS = ZeroVector(16);
            Matrix LHS = ZeroMatrix(16,16);

            const auto& r_process_info = modelPart.GetProcessInfo();

            pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
            pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

            reference_RHS[0] = -241.33267395784;
            reference_RHS[1] = 17.632097705363;
            reference_RHS[2]= -182.29720750899;
            reference_RHS[3]= 0.020140666295797;
            reference_RHS[4]= 101.35024256784;
            reference_RHS[5]= -35.407695645367;
            reference_RHS[6]= -564.88597999491;
            reference_RHS[7]= 0.40838413082497;
            reference_RHS[8]= 148.25438916454;
            reference_RHS[9]= -52.14562036006;
            reference_RHS[10]= -605.8561165501;
            reference_RHS[11]= -0.070140666295799;
            reference_RHS[12]= 185.67518738712;
            reference_RHS[13]= -47.99651481075;
            reference_RHS[14]= -658.07450777431;
            reference_RHS[15]= -0.45838413082497;

            // Check the RHS values (the RHS is computed as the LHS x previous_solution,
            // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
            KRATOS_EXPECT_VECTOR_NEAR(RHS, reference_RHS, 1e-7);

        }
       // /** Checks the TwoFluidNavierStokesAlphaMethod3D4N element
        //  * Checks the LHS and RHS for a negative element (distance <= 0.0)
        //  */
        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesNegativeSideAlphaMethod3D4N, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(2);

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

            modelPart.GetProcessInfo().SetValue(SPECTRAL_RADIUS_LIMIT, 0.0);
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            NewtonianTwoFluid3DLaw::Pointer pConsLaw(new NewtonianTwoFluid3DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
            std::vector<ModelPart::IndexType> elemNodes{ 1, 2, 3, 4 };
            modelPart.CreateNewElement("TwoFluidNavierStokesAlphaMethod3D4N", 1, elemNodes, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);
            // Fake time advance to set the previous ProcessInfo container
            modelPart.CloneSolutionStep();

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
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                }
            }
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
            KRATOS_EXPECT_NEAR(RHS(0), 16.6700148724, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(1), 17.5504263651, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(2), 76.0097167815, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(3),0.951210720699, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(4), -35.762027016, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(5), -35.3837331297, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(6), -701.815365218, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(7), -0.0572209011027, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(8), -30.4898366163, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(9), -52.0920660863, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(10), -784.514819936, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(11), -0.0700797902459, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(12), -35.584025235, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(13), -47.95583463, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(14), -879.858882572, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(15), -0.923910029351, 1e-7);

        }

        // /** Checks the TwoFluidNavierStokesAlphaMethod3D4N element
        //  * Checks the LHS and RHS for a positive element (distance > 0.0)
        //  */
        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesPositiveSideAlphaMethod3D4N, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(2);

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
            modelPart.GetProcessInfo().SetValue(SPECTRAL_RADIUS_LIMIT, 0.0);
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            NewtonianTwoFluid3DLaw::Pointer pConsLaw(new NewtonianTwoFluid3DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
            std::vector<ModelPart::IndexType> elemNodes{ 1, 2, 3, 4 };
            modelPart.CreateNewElement("TwoFluidNavierStokesAlphaMethod3D4N", 1, elemNodes, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);

            // Fake time advance to set the previous ProcessInfo container
            modelPart.CloneSolutionStep();

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
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z, 1) = -9.81;
            }

            for (unsigned int i = 0; i < 4; i++) {
                pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
                for (unsigned int k = 0; k < 3; k++) {
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k] = vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;

                }
            }
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
            KRATOS_EXPECT_NEAR(RHS(0), 16.6700148724, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(1), 17.5504263651, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(2), 76.0097167815, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(3), 0.951210720699, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(4), -35.762027016, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(5), -35.3837331297, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(6),-701.815365218, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(7), -0.0572209011027, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(8), -30.4898366163, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(9), -52.0920660863, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(10), -784.514819936, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(11),-0.0700797902459, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(12), -35.584025235, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(13), -47.95583463, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(14), -879.858882572, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(15), -0.923910029351, 1e-7);

        }

        /** Checks the TwoFluidNavierStokesAlphaMethod2D3N element in a hydrostatic case.
         *  Checks the computation of the RHS
         */
        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesAlphaMethod2D3NHydrostatic, FluidDynamicsApplicationFastSuite)
        {

            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(2);

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

            modelPart.GetProcessInfo().SetValue(SPECTRAL_RADIUS_LIMIT, 0.0);
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-03);
            Newtonian2DLaw::Pointer pConsLaw(new Newtonian2DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 2.0, 0.0, 0.0);  // 0 = node 1
            modelPart.CreateNewNode(2, 2.0, 2.0, 0.0);	// 1 = node 2
            modelPart.CreateNewNode(3, 0.0, 2.0, 0.0);	// 2 = node 3

            std::vector<ModelPart::IndexType> elemNodes1 {1, 2, 3};

            modelPart.CreateNewElement("TwoFluidNavierStokesAlphaMethod2D3N", 1, elemNodes1, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);
            // Fake time advance to set the previous ProcessInfo container
            modelPart.CloneSolutionStep();

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
                    // pElement->GetGeometry()[i].Fix(VELOCITY);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
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

            const auto& r_process_info = modelPart.GetProcessInfo();
            pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
            pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

            double det;
            MathUtils<double>::InvertMatrix(LHS, LHS, det);

            const Vector solVec = prod(LHS, RHS);

            // The remaining residuals in the velocities have the size of the boundary integrals over the enriched pressure.
            // If the "standard" pressure shape functions are used, the results do not hold.

            KRATOS_EXPECT_NEAR(RHS(0), 0.0, 1e-7);		// U_x at node 1
            KRATOS_EXPECT_NEAR(RHS(1), -17.5, 1e-7); 	// U_y at node 1
            KRATOS_EXPECT_NEAR(RHS(2), 0.0, 1e-7);		// P   at node 1

            KRATOS_EXPECT_NEAR(RHS(3), 7.5, 1e-7);		// U_x at node 2
            KRATOS_EXPECT_NEAR(RHS(4), 0.0, 1e-7);		// U_y at node 2
            KRATOS_EXPECT_NEAR(RHS(5), 0.0, 1e-7);		// P   at node 2

            KRATOS_EXPECT_NEAR(RHS(6), -7.5, 1e-7);		// U_x at node 3
            KRATOS_EXPECT_NEAR(RHS(7), -7.5, 1e-7);		// U_y at node 3
            KRATOS_EXPECT_NEAR(RHS(8), 0.0, 1e-7);		// P   at node 3
        }

        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesDarcyAlphaMethod3D4N, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(2);

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


            modelPart.GetProcessInfo().SetValue(SPECTRAL_RADIUS_LIMIT,0.0);
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
            pElemProp->SetValue(DENSITY, 1000.0);
            pElemProp->SetValue(LIN_DARCY_COEF, 1.0 / 4.339E-08);
            pElemProp->SetValue(NONLIN_DARCY_COEF, 1.0 / 5.086E-04);
            pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            NewtonianTwoFluid3DLaw::Pointer pConsLaw(new NewtonianTwoFluid3DLaw());
            pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

            // Geometry creation
            modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
            modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
            std::vector<ModelPart::IndexType> elemNodes{ 1, 2, 3, 4 };
            modelPart.CreateNewElement("TwoFluidNavierStokesAlphaMethod3D4N", 1, elemNodes, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);
            // Fake time advance to set the previous ProcessInfo container
            modelPart.CloneSolutionStep();

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
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                }
            }
            pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) = 1.0;

            // Compute RHS and LHS
            Vector RHS = ZeroVector(16);
            Matrix LHS = ZeroMatrix(16, 16);

            const auto& r_process_info = modelPart.GetProcessInfo();
            pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
            pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

            // Check the RHS values (the RHS is computed as the LHS x previous_solution,
            // hence, it is assumed that if the RHS is correct, the LHS is correct as well)


            KRATOS_EXPECT_NEAR(RHS(0), -119.5094945, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(1), 13.86753717, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(2), -264.4471963, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(3), -0.008525833469, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(4), 31.81898305, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(5), -23.48989075, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(6), -515.5487839, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(7), 0.1753941423, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(8), 52.86326183, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(9), -58.7570097, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(10), -510.5819917, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(11), -0.03898369627, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(12), 89.23040167, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(13), -30.84418781, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(14), -581.986129, 1e-7);
            KRATOS_EXPECT_NEAR(RHS(15), -0.2278846126, 1e-7);

        }

        // Giving a value different to zero to source term in order to test it.

        // /** Checks the TwoFluidNavierStokesAlphaMethod2D3N element with a source term in mass conservation equation
        //  * Checks the LHS and RHS for a cut element
        //  */

        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesAlphaMethod2D3NError, FluidDynamicsApplicationFastSuite)
        {
             Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(2);

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
            modelPart.GetProcessInfo().SetValue(SPECTRAL_RADIUS_LIMIT, 0.0);
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR,10.0);


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
            modelPart.CreateNewElement("TwoFluidNavierStokesAlphaMethod2D3N", 1, elemNodes, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);
            // Fake time advance to set the previous ProcessInfo container
            modelPart.CloneSolutionStep();


            // Define the nodal values
            Matrix vel_original(3, 2);
            vel_original(0, 0) = 0.0;
            vel_original(0, 1) = 0.1;
            vel_original(1, 0) = 0.1;
            vel_original(1, 1) = 0.2;
            vel_original(2, 0) = 0.2;
            vel_original(2, 1) = 0.3;

            // Set the nodal DENSITY and DYNAMIC_VISCOSITY values
            for (NodeIteratorType it_node = modelPart.NodesBegin(); it_node < modelPart.NodesEnd(); ++it_node)
            {
                it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = -9.81;
            }

            for (unsigned int i = 0; i < 3; i++)
            {
                pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
                for (unsigned int k = 0; k < 2; k++)
                {
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k] = vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9 * vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                }
            }
            pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = 0.5;

            // Compute RHS and LHS
            Vector RHS = ZeroVector(9);
            Matrix LHS = ZeroMatrix(9, 9);

            const auto &r_process_info = modelPart.GetProcessInfo();
            pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
            pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

            // Check the RHS values (the RHS is computed as the LHS x previous_solution,
            // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
            Vector reference_RHS = ZeroVector(9);
            reference_RHS[0] = 5836.450756;
            reference_RHS[1] = 4764.408449;
            reference_RHS[2] = -28.37919168;
            reference_RHS[3] = -6019.528438;
            reference_RHS[4] = 1996.025722;
            reference_RHS[5] = -16.19446094;
            reference_RHS[6] = 239.7394;
            reference_RHS[7] = -3593.338385;
            reference_RHS[8] = -5.576347377;

            KRATOS_EXPECT_VECTOR_NEAR(reference_RHS, RHS, 1e-2);
        }

        // /** Checks the TwoFluidNavierStokesAlphaMethod3D4N element with a source term in mass conservation equation
        //  * Checks the LHS and RHS for a cut element
        //  */

        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesAlphaMethodCut3D4NError, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart& modelPart = current_model.CreateModelPart("Main");
            modelPart.SetBufferSize(2);

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

            modelPart.GetProcessInfo().SetValue(SPECTRAL_RADIUS_LIMIT, 0.0);
            modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            modelPart.GetProcessInfo().SetValue(VOLUME_ERROR, 10.0);

            // Set the element properties
            Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
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
            modelPart.CreateNewElement("TwoFluidNavierStokesAlphaMethod3D4N", 1, elemNodes, pElemProp);

            Element::Pointer pElement = modelPart.pGetElement(1);
             // Fake time advance to set the previous ProcessInfo container
            modelPart.CloneSolutionStep();

            // Define the nodal values
            Matrix vel_original(4, 3);
            vel_original(0, 0) = 0.0;
            vel_original(0, 1) = 0.1;
            vel_original(0, 2) = 0.2;
            vel_original(1, 0) = 0.1;
            vel_original(1, 1) = 0.2;
            vel_original(1, 2) = 0.3;
            vel_original(2, 0) = 0.2;
            vel_original(2, 1) = 0.3;
            vel_original(2, 2) = 0.4;
            vel_original(3, 0) = 0.3;
            vel_original(3, 1) = 0.4;
            vel_original(3, 2) = 0.5;

            // Set the nodal BODY_FORCE, DENSITY and DYNAMIC_VISCOSITY values
            for (NodeIteratorType it_node = modelPart.NodesBegin(); it_node < modelPart.NodesEnd(); ++it_node)
            {
                it_node->FastGetSolutionStepValue(DENSITY) = pElemProp->GetValue(DENSITY);
                it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = pElemProp->GetValue(DYNAMIC_VISCOSITY);
                it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = -9.81;
            }

            for (unsigned int i = 0; i < 4; i++)
            {
                pElement->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
                for (unsigned int k = 0; k < 3; k++)
                {
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k] = vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9 * vel_original(i, k);
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k] = 0.0;
                    pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
                }
            }
            pElement->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
            pElement->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -1.0;
            pElement->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) = 1.0;

            // Compute RHS and LHS
            Vector RHS = ZeroVector(16);
            Matrix LHS = ZeroMatrix(16, 16);

            const auto &r_process_info = modelPart.GetProcessInfo();
            pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
            pElement->CalculateLocalSystem(LHS, RHS, r_process_info);
            // std::cout << pElement->Info() << std::setprecision(10) << std::endl;

            // Check the RHS values (the RHS is computed as the LHS x previous_solution,
            // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
            Vector reference_RHS = ZeroVector(16);
            reference_RHS[0] =3948.998936;
            reference_RHS[1] = 2873.057834;
            reference_RHS[2] = 3804.061234;
            reference_RHS[3] = -5.479993813;
            reference_RHS[4] = -3931.287229;
            reference_RHS[5] = 473.7242468;
            reference_RHS[6] = -728.3575273;
            reference_RHS[7] = -3.988810599;
            reference_RHS[8] = 130.4431454;
            reference_RHS[9] = -3350.834823;
            reference_RHS[10] = -433.0021081;
            reference_RHS[11] = -2.9057729;
            reference_RHS[12] = -92.27454488;
            reference_RHS[13] = 684.7556692;
            reference_RHS[14] = -4513.788544;
            reference_RHS[15] = -4.392089354;
            KRATOS_EXPECT_VECTOR_NEAR(reference_RHS, RHS, 1e-2);
        }

        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesAlphaMethodArtificialDynamicViscosity3D4N, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart &r_model_part = current_model.CreateModelPart("Main");
            r_model_part.SetBufferSize(3);

            // Variables addition
            r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
            r_model_part.AddNodalSolutionStepVariable(DENSITY);
            r_model_part.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            r_model_part.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            r_model_part.AddNodalSolutionStepVariable(PRESSURE);
            r_model_part.AddNodalSolutionStepVariable(VELOCITY);
            r_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
            r_model_part.AddNodalSolutionStepVariable(DISTANCE);

            // Process info creation
            double delta_time = 0.1;

            r_model_part.GetProcessInfo().SetValue(SPECTRAL_RADIUS_LIMIT, 0.0);
            r_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            r_model_part.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            auto p_elem_prop = r_model_part.CreateNewProperties(0);
            p_elem_prop->SetValue(DENSITY, 1000.0);
            p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            auto p_cons_law = Kratos::make_shared<NewtonianTwoFluid3DLaw>();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_cons_law);

            // Geometry creation
            r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
            r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
            r_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
            std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3, 4};
            auto p_element = r_model_part.CreateNewElement("TwoFluidNavierStokesAlphaMethod3D4N", 1, elem_nodes, p_elem_prop);

            // Fake time advance to set the previous ProcessInfo container
            r_model_part.CloneSolutionStep();

            array_1d<double,3> aux_v;
            for (auto& r_node : r_model_part.Nodes()) {
                aux_v[0] = r_node.Id();
                aux_v[1] = 2.0*r_node.Id();
                aux_v[2] = std::pow(r_node.Id(),2);
                r_node.GetSolutionStepValue(VELOCITY) = aux_v;
                r_node.GetSolutionStepValue(PRESSURE) = std::pow(r_node.Id(), 3.0);
                r_node.GetSolutionStepValue(DISTANCE) = 1.0;
                r_node.GetSolutionStepValue(DENSITY) = p_elem_prop->GetValue(DENSITY);
                r_node.GetSolutionStepValue(DYNAMIC_VISCOSITY) = p_elem_prop->GetValue(DYNAMIC_VISCOSITY);
            }

            // The integration points
            double art_dyn_visc;
            const auto &r_process_info = r_model_part.GetProcessInfo();
            p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law

            // Obtain the artificial dynamic viscosity in the element.
            p_element->Calculate(ARTIFICIAL_DYNAMIC_VISCOSITY, art_dyn_visc, r_process_info);

            // std::cout << std::setprecision(12) << art_dyn_visc << std::endl;
            const double tolerance = 1.0e-8;
            double exact_art_dyn_visc =  3019.8895531;

            KRATOS_EXPECT_NEAR(art_dyn_visc, exact_art_dyn_visc, tolerance);
        }

        KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesAlphaMethodArtificialDynamicViscosity2D3N, FluidDynamicsApplicationFastSuite)
        {
            Model current_model;
            ModelPart &r_model_part = current_model.CreateModelPart("Main");
            r_model_part.SetBufferSize(3);

            // Variables addition
            r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
            r_model_part.AddNodalSolutionStepVariable(DENSITY);
            r_model_part.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
            r_model_part.AddNodalSolutionStepVariable(DYNAMIC_TAU);
            r_model_part.AddNodalSolutionStepVariable(PRESSURE);
            r_model_part.AddNodalSolutionStepVariable(VELOCITY);
            r_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
            r_model_part.AddNodalSolutionStepVariable(DISTANCE);

            // Process info creation
            double delta_time = 0.1;

            r_model_part.GetProcessInfo().SetValue(SPECTRAL_RADIUS_LIMIT, 0.0);
            r_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
            r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
            r_model_part.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

            // Set the element properties
            auto p_elem_prop = r_model_part.CreateNewProperties(0);
            p_elem_prop->SetValue(DENSITY, 1000.0);
            p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
            auto p_cons_law = Kratos::make_shared<NewtonianTwoFluid3DLaw>();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_cons_law);

            // Geometry creation
            r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
            r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
            r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
            std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3};
            auto p_element = r_model_part.CreateNewElement("TwoFluidNavierStokesAlphaMethod2D3N", 1, elem_nodes, p_elem_prop);

            // Fake time advance to set the previous ProcessInfo container
            r_model_part.CloneSolutionStep();

            array_1d<double, 3> aux_v;
            for (auto &r_node : r_model_part.Nodes())
            {
                aux_v[0] = r_node.Id();
                aux_v[1] = 2.0 * r_node.Id();
                aux_v[2] = std::pow(r_node.Id(), 2);
                r_node.GetSolutionStepValue(VELOCITY) = aux_v;
                r_node.GetSolutionStepValue(PRESSURE) = std::pow(r_node.Id(), 3.0);
                r_node.GetSolutionStepValue(DISTANCE) = 1.0;
                r_node.GetSolutionStepValue(DENSITY) = p_elem_prop->GetValue(DENSITY);
                r_node.GetSolutionStepValue(DYNAMIC_VISCOSITY) = p_elem_prop->GetValue(DYNAMIC_VISCOSITY);
            }

            // The integration points
            double art_dyn_visc;
            const auto &r_process_info = r_model_part.GetProcessInfo();
            p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law

            // Obtain the artificial dynamic viscosity in the element.
            p_element->Calculate(ARTIFICIAL_DYNAMIC_VISCOSITY, art_dyn_visc, r_process_info);

            const double tolerance = 1.0e-8;
            double exact_art_dyn_visc = 3772.34868808;
            KRATOS_EXPECT_NEAR(art_dyn_visc, exact_art_dyn_visc, tolerance);
        }
    } // namespace Testing
}  // namespace Kratos.