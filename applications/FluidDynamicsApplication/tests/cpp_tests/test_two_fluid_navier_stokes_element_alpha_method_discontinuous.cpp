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
#include "spaces/ublas_space.h"
#include "includes/properties.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"
#include "includes/global_pointer_variables.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"
#include "custom_constitutive/newtonian_two_fluid_2d_law.h"
#include "custom_constitutive/newtonian_two_fluid_3d_law.h"
#include "processes/find_nodal_neighbours_process.h"
#include "utilities/normal_calculation_utils.h"

// Application includes
#include "custom_elements/two_fluid_navier_stokes_alpha_method_discontinuous.h"

namespace Kratos::Testing {

    // KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesAlphaMethod2D3N, FluidDynamicsApplicationFastSuite)
    // {
    //     Model current_model;
    //     r_model_part& r_model_part = current_model.Creater_model_part("Main");
    //     //FIXME: SET TO 2 WHEN A NEW DATA CONTAINER IS CREATED
    //     r_model_part.SetBufferSize(3);

    //     // Variables addition
    //     r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    //     r_model_part.AddNodalSolutionStepVariable(DENSITY);
    //     r_model_part.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
    //     r_model_part.AddNodalSolutionStepVariable(DYNAMIC_TAU);
    //     r_model_part.AddNodalSolutionStepVariable(PRESSURE);
    //     r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    //     r_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    //     r_model_part.AddNodalSolutionStepVariable(DISTANCE);

    //     // Process info creation
    //     double delta_time = 0.1;
    //     r_model_part.GetProcessInfo().SetValue(SPECTRAL_RADIUS_LIMIT, 0.0);
    //     r_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.0);
    //     r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    //     r_model_part.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

    //     // Set the element properties
    //     Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);
    //     p_elem_prop->SetValue(DENSITY, 1000.0);
    //     p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.0e-03);
    //     auto p_cons_law = Kratos::make_shared<NewtonianTwoFluid2DLaw>();
    //     p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_cons_law);

    //     // Geometry creation
    //     r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    //     r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    //     r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    //     std::vector<r_model_part::IndexType> elemNodes {1, 2, 3};
    //     r_model_part.CreateNewElement("TwoFluidNavierStokesAlphaMethod2D3N", 1, elemNodes, p_elem_prop);

    //     Element::Pointer p_element = r_model_part.pGetElement(1);
    //         // Fake time advance to set the previous ProcessInfo container
    //     r_model_part.CloneSolutionStep();

    //     // Define the nodal values
    //     Matrix vel_original(3, 3);
    //     vel_original(0, 0) = 0.0;
    //     vel_original(0, 1) = 0.1;
    //     vel_original(0, 2) = 0.2;
    //     vel_original(1, 0) = 0.1;
    //     vel_original(1, 1) = 0.2;
    //     vel_original(1, 2) = 0.3;
    //     vel_original(2, 0) = 0.2;
    //     vel_original(2, 1) = 0.3;
    //     vel_original(2, 2) = 0.4;

    //     // Set the nodal BODY_FORCE, DENSITY and DYNAMIC_VISCOSITY values
    //     for (NodeIteratorType it_node = r_model_part.NodesBegin(); it_node < r_model_part.NodesEnd(); ++it_node)
    //     {
    //         it_node->FastGetSolutionStepValue(DENSITY) = p_elem_prop->GetValue(DENSITY);
    //         it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = p_elem_prop->GetValue(DYNAMIC_VISCOSITY);
    //         it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = -9.81;
    //     }

    //     for (unsigned int i = 0; i < 3; i++)
    //     {
    //         p_element->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
    //         for (unsigned int k = 0; k < 3; k++)
    //         {
    //             p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k] = vel_original(i, k);
    //             p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] =0.9*vel_original(i, k);
    //             p_element->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k] = 0.0;
    //             p_element->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
    //         }
    //     }
    //     p_element->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
    //     p_element->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
    //     p_element->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -1.0;

    //     // Compute RHS and LHS
    //     Vector RHS = ZeroVector(9);
    //     Matrix LHS = ZeroMatrix(9, 9);

    //     const auto &r_process_info = r_model_part.GetProcessInfo();
    //     p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
    //     p_element->CalculateLocalSystem(LHS, RHS, r_process_info);
    //     // Check the RHS values (the RHS is computed as the LHS x previous_solution,
    //     // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
    //     Vector reference_RHS = ZeroVector(9);
    //     reference_RHS[0] = -14.67064121;
    //     reference_RHS[1] = -38.30659006;
    //     reference_RHS[2] =-0.02477803762;
    //     reference_RHS[3] = -44.16968252;
    //     reference_RHS[4] = -61.67597243;
    //     reference_RHS[5] = -0.05915402772;
    //     reference_RHS[6] =-43.90588458;
    //     reference_RHS[7] = -79.83781791;
    //     reference_RHS[8] = -0.06606793466;

    //     KRATOS_EXPECT_VECTOR_NEAR(reference_RHS, RHS, 1e-2);
    // }

    // KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesCutAlphaMethod3D4N, FluidDynamicsApplicationFastSuite)
    // {
    //     Model current_model;
    //     r_model_part& r_model_part = current_model.Creater_model_part("Main");
    //     r_model_part.SetBufferSize(3);

    //     // Variables addition
    //     r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    //     r_model_part.AddNodalSolutionStepVariable(DENSITY);
    //     r_model_part.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
    //     r_model_part.AddNodalSolutionStepVariable(DYNAMIC_TAU);
    //     r_model_part.AddNodalSolutionStepVariable(PRESSURE);
    //     r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    //     r_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    //     r_model_part.AddNodalSolutionStepVariable(DISTANCE);

    //     // Process info creation
    //     double delta_time = 0.1;

    //     r_model_part.GetProcessInfo().SetValue(SPECTRAL_RADIUS_LIMIT, 0.0);
    //     r_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
    //     r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    //     r_model_part.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

    //     // Set the element properties
    //     Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);
    //     p_elem_prop->SetValue(DENSITY, 1000.0);
    //     p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.0e-05);
    //     NewtonianTwoFluid3DLaw::Pointer p_cons_law(new NewtonianTwoFluid3DLaw());
    //     p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_cons_law);

    //     // Geometry creation
    //     r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    //     r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    //     r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    //     r_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
    //     std::vector<r_model_part::IndexType> elemNodes {1, 2, 3, 4};
    //     r_model_part.CreateNewElement("TwoFluidNavierStokesAlphaMethod3D4N", 1, elemNodes, p_elem_prop);

    //     Element::Pointer p_element = r_model_part.pGetElement(1);
    //     // Fake time advance to set the previous ProcessInfo container
    //     r_model_part.CloneSolutionStep();

    //     // Define the nodal values
    //     Matrix vel_original(4,3);
    //     vel_original(0,0) = 0.0; vel_original(0,1) = 0.1; vel_original(0,2) = 0.2;
    //     vel_original(1,0) = 0.1; vel_original(1,1) = 0.2; vel_original(1,2) = 0.3;
    //     vel_original(2,0) = 0.2; vel_original(2,1) = 0.3; vel_original(2,2) = 0.4;
    //     vel_original(3,0) = 0.3; vel_original(3,1) = 0.4; vel_original(3,2) = 0.5;

    //     // Set the nodal BODY_FORCE, DENSITY and DYNAMIC_VISCOSITY values
    //     for (NodeIteratorType it_node=r_model_part.NodesBegin(); it_node<r_model_part.NodesEnd(); ++it_node){
    //         it_node->FastGetSolutionStepValue(DENSITY) = p_elem_prop->GetValue(DENSITY);
    //         it_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = p_elem_prop->GetValue(DYNAMIC_VISCOSITY);
    //         it_node->FastGetSolutionStepValue(BODY_FORCE_Z) = -9.81;
    //     }

    //     for(unsigned int i=0; i<4; i++){
    //         p_element->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE)    = 0.0;
    //         for(unsigned int k=0; k<3; k++){
    //             p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel_original(i,k);
    //             p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i,k);
    //             p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.9*vel_original(i,k);
    //             p_element->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
    //             p_element->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
    //             p_element->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
    //         }
    //     }
    //     p_element->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
    //     p_element->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) =  1.0;
    //     p_element->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = -1.0;
    //     p_element->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) =  1.0;

    //     // Compute RHS and LHS
    //     Vector RHS = ZeroVector(16);
    //     Matrix LHS = ZeroMatrix(16,16);

    //     const auto& r_process_info = r_model_part.GetProcessInfo();
    //     p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
    //     p_element->CalculateLocalSystem(LHS, RHS, r_process_info);

    //     // Check the RHS values (the RHS is computed as the LHS x previous_solution,
    //     // hence, it is assumed that if the RHS is correct, the LHS is correct as well)

    //     KRATOS_EXPECT_NEAR(RHS(0), -119.5094945, 1e-7);
    //     KRATOS_EXPECT_NEAR(RHS(1), 13.86753717, 1e-7);
    //     KRATOS_EXPECT_NEAR(RHS(2), -264.4471963, 1e-7);
    //     KRATOS_EXPECT_NEAR(RHS(3), -0.008525833469, 1e-7);
    //     KRATOS_EXPECT_NEAR(RHS(4), 31.81898305, 1e-7);
    //     KRATOS_EXPECT_NEAR(RHS(5), -23.48989075, 1e-7);
    //     KRATOS_EXPECT_NEAR(RHS(6), -515.5487839, 1e-7);
    //     KRATOS_EXPECT_NEAR(RHS(7), 0.1753941423, 1e-7);
    //     KRATOS_EXPECT_NEAR(RHS(8), 52.86326183, 1e-7);
    //     KRATOS_EXPECT_NEAR(RHS(9), -58.7570097, 1e-7);
    //     KRATOS_EXPECT_NEAR(RHS(10), -510.5819917, 1e-7);
    //     KRATOS_EXPECT_NEAR(RHS(11), -0.03898369627, 1e-7);
    //     KRATOS_EXPECT_NEAR(RHS(12), 89.23040167, 1e-7);
    //     KRATOS_EXPECT_NEAR(RHS(13), -30.84418781, 1e-7);
    //     KRATOS_EXPECT_NEAR(RHS(14), -581.986129, 1e-7);
    //     KRATOS_EXPECT_NEAR(RHS(15), -0.2278846126, 1e-7);
    // }

    KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesAlphaMethodDiscontinuous2D3NHydrostatic, FluidDynamicsApplicationFastSuite)
    {

        Model current_model;
        auto& r_model_part = current_model.CreateModelPart("Main");
        r_model_part.SetBufferSize(2);

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
        r_model_part.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
        r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
        r_model_part.GetProcessInfo().SetValue(VOLUME_ERROR, 0.0);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(DENSITY, 1000.0);
        p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.0e-03);
        auto p_cons_law = Kratos::make_shared<Newtonian2DLaw>();
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_cons_law);

        // Geometry creation
        r_model_part.CreateNewNode(1, 2.0, 0.0, 0.0);  // 0 = node 1
        r_model_part.CreateNewNode(2, 2.0, 2.0, 0.0);	// 1 = node 2
        r_model_part.CreateNewNode(3, 0.0, 2.0, 0.0);	// 2 = node 3

        std::vector<std::size_t> elemNodes1 {1, 2, 3};
        auto p_element = r_model_part.CreateNewElement("TwoFluidNavierStokesAlphaMethodDiscontinuous2D3N", 1, elemNodes1, p_elem_prop);

        // Fake time advance to set the previous ProcessInfo container
        r_model_part.CloneSolutionStep();

        // Define the nodal values as 0 for hydrostatic case
        Matrix vel_original(3,2);
        vel_original(0,0) = 0.0; vel_original(0,1) = 0.0;
        vel_original(1,0) = 0.0; vel_original(1,1) = 0.0;
        vel_original(2,0) = 0.0; vel_original(2,1) = 0.0;

        // Setting equal nodal values for DENSITY, DYNAMIC_VISCOSITY, BODY_FORCE
        for (auto& r_node : r_model_part.Nodes()){
            r_node.FastGetSolutionStepValue(DENSITY) = p_elem_prop->GetValue(DENSITY);
            r_node.FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = p_elem_prop->GetValue(DYNAMIC_VISCOSITY);
            r_node.FastGetSolutionStepValue(BODY_FORCE_X) = 0.0;
            r_node.FastGetSolutionStepValue(BODY_FORCE_Y) = -10.0;
            r_node.FastGetSolutionStepValue(BODY_FORCE_Z) = 0.0;
        }

        // Setting the density (different for nodes since element cut by surface)
        p_element->GetGeometry()[0].FastGetSolutionStepValue(DENSITY) = 2.0;
        p_element->GetGeometry()[1].FastGetSolutionStepValue(DENSITY) = 1.0;
        p_element->GetGeometry()[2].FastGetSolutionStepValue(DENSITY) = 1.0;

        for(unsigned int i=0; i<3; i++){
            for(unsigned int k=0; k<2; k++){
                p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel_original(i,k);
                p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = vel_original(i,k);
                // p_element->GetGeometry()[i].Fix(VELOCITY);
                p_element->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
                p_element->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
            }
        }

        p_element->GetGeometry()[0].Fix(VELOCITY_X);
        p_element->GetGeometry()[0].Fix(VELOCITY_Y);

        // Setting the density (different for nodes to define the position of the surface)
        p_element->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
        p_element->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = 1.0;
        p_element->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) = 1.0;

        // Simon : Setting the pressure
        p_element->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE) = 30.0;
        p_element->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE) = 0.0;
        p_element->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE) = 0.0;

        p_element->GetGeometry()[0].Fix(PRESSURE);

        // Compute RHS and LHS
        Vector RHS = ZeroVector(9);
        Matrix LHS = ZeroMatrix(9,9);

        const auto& r_process_info = r_model_part.GetProcessInfo();
        p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
        p_element->CalculateLocalSystem(LHS, RHS, r_process_info);

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

}  // namespace Kratos::Testing