//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Uxue Chasco
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "containers/model.h"
#include "spaces/ublas_space.h"
#include "includes/properties.h"
#include "includes/model_part.h"
#include "includes/global_pointer_variables.h"
#include "utilities/math_utils.h"
#include "utilities/normal_calculation_utils.h"

// Application includes
#include "custom_elements/two_fluid_navier_stokes_fractional.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"
#include "custom_constitutive/newtonian_two_fluid_3d_law.h"
#include "processes/find_nodal_neighbours_process.h"
#include "tests/cpp_tests/fluid_dynamics_fast_suite.h"

namespace Kratos::Testing
{

typedef ModelPart::IndexType                              IndexType;
typedef ModelPart::NodeIterator                           NodeIteratorType;

/** Checks the TwoFluidNavierStokesFractional2D3N element.
 * Checks the LHS and RHS computation
 */

KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesFractional2D3N, FluidDynamicsApplicationFastSuite)
{

    Model current_model;
    ModelPart& modelPart = current_model.CreateModelPart("Main");
    // Variables addition
    modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
    modelPart.AddNodalSolutionStepVariable(DENSITY);
    modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
    modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
    modelPart.AddNodalSolutionStepVariable(PRESSURE);
    modelPart.AddNodalSolutionStepVariable(VELOCITY);
    modelPart.AddNodalSolutionStepVariable(FRACTIONAL_VELOCITY);
    modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
    modelPart.AddNodalSolutionStepVariable(DISTANCE);

    // Process info creation
    const int buffer_size = 4;
    const double delta_time = 0.1;
    modelPart.SetBufferSize(buffer_size);
    modelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    modelPart.GetProcessInfo().SetValue(TIME, 0.6);
    modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    modelPart.CloneTimeStep(modelPart.GetProcessInfo().GetValue(TIME) + delta_time);
    modelPart.GetProcessInfo().GetPreviousTimeStepInfo(1).SetValue(TIME, 0.0);
    modelPart.GetProcessInfo().GetPreviousTimeStepInfo(1).SetValue(DELTA_TIME, delta_time);
    modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
    Vector bdf_coefs(3);
    bdf_coefs[0] = 3.0/(2.0*delta_time);
    bdf_coefs[1] = -2.0/delta_time;
    bdf_coefs[2] = 0.5/delta_time;
    modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
    modelPart.GetProcessInfo().SetValue(WATER_VOLUME_ERROR,0.0);
    modelPart.GetProcessInfo().SetValue(AIR_VOLUME_ERROR, 0.0);

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
    modelPart.CreateNewElement("TwoFluidNavierStokesFractional2D3N", 1, elemNodes, pElemProp);

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
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY)[k] = 0.1 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 1)[k] = 0.2 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 2)[k] = 0.3 * vel_original(i, k);
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

    KRATOS_EXPECT_NEAR(RHS(0), 285.327812743, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(1),-75.7286708218, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(2), 0.678487936315, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(3),-758.586246712, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(4), -529.046820664, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(5), -0.808596445975, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(6), -549.320499823, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(7), -931.191955616, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(8), -0.0198914903401, 1e-7);
}

// /** Checks the TwoFluidNavierStokesFractional3D4N element
//  * Checks the LHS and RHS for a cut element
//  */
KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesFractionalCut3D4N, FluidDynamicsApplicationFastSuite)
{
    Model current_model;
    ModelPart& modelPart = current_model.CreateModelPart("Main");

    // Process info creation
    const int buffer_size = 4;
    const double delta_time = 0.1;
    modelPart.SetBufferSize(buffer_size);
    modelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
    modelPart.GetProcessInfo().SetValue(TIME, 0.6);
    modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    modelPart.CloneTimeStep(modelPart.GetProcessInfo().GetValue(TIME) + delta_time);
    modelPart.GetProcessInfo().GetPreviousTimeStepInfo(1).SetValue(TIME, 0.0);
    modelPart.GetProcessInfo().GetPreviousTimeStepInfo(1).SetValue(DELTA_TIME, delta_time);
    modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
    modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    Vector bdf_coefs(3);
    bdf_coefs[0] = 3.0/(2.0*delta_time);
    bdf_coefs[1] = -2.0/delta_time;
    bdf_coefs[2] = 0.5/delta_time;
    modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
    modelPart.GetProcessInfo().SetValue(AIR_VOLUME_ERROR, 0.0);
    modelPart.GetProcessInfo().SetValue(WATER_VOLUME_ERROR, 0.0);

    // Variables addition
    modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
    modelPart.AddNodalSolutionStepVariable(DENSITY);
    modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
    modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
    modelPart.AddNodalSolutionStepVariable(PRESSURE);
    modelPart.AddNodalSolutionStepVariable(VELOCITY);
    modelPart.AddNodalSolutionStepVariable(FRACTIONAL_VELOCITY);
    modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
    modelPart.AddNodalSolutionStepVariable(DISTANCE);

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
    modelPart.CreateNewElement("TwoFluidNavierStokesFractional3D4N", 1, elemNodes, pElemProp);

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
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY)[k] = 0.1 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 1)[k] = 0.2 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 2)[k] = 0.3 * vel_original(i, k);
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
    std::cout << std::setprecision(12) << RHS << std::endl;
    // Check the RHS values (the RHS is computed as the LHS x previous_solution,
    // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
    KRATOS_EXPECT_NEAR(RHS(0), -182.780712656, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(1), -5.97373013, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(2), -366.425812874, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(3), 0.136180451508, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(4), -65.5996318562, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(5),-195.157657105, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(6), -732.8993499, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(7), 0.228653068676, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(8), -30.1350973711, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(9), -309.557107994, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(10), -744.055344314, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(11), -0.182507276303, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(12), 4.70728308942, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(13), -259.815870553, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(14), -872.015239241, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(15), -0.282326243881, 1e-7);
}


// /** Checks the TwoFluidNavierStokesFractional3D4N element
//  * Checks the LHS and RHS for a negative element (distance <= 0.0)
//  */
KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesFractionalNegativeSide3D4N, FluidDynamicsApplicationFastSuite)
{
    Model current_model;
    ModelPart& modelPart = current_model.CreateModelPart("Main");

    // Process info creation
    const int buffer_size = 4;
    const double delta_time = 0.1;
    modelPart.SetBufferSize(buffer_size);
    modelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
    modelPart.GetProcessInfo().SetValue(TIME, 0.6);
    modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    modelPart.CloneTimeStep(modelPart.GetProcessInfo().GetValue(TIME) + delta_time);
    modelPart.GetProcessInfo().GetPreviousTimeStepInfo(1).SetValue(TIME, 0.0);
    modelPart.GetProcessInfo().GetPreviousTimeStepInfo(1).SetValue(DELTA_TIME, delta_time);
    modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
    modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    Vector bdf_coefs(3);
    bdf_coefs[0] = 3.0/(2.0*delta_time);
    bdf_coefs[1] = -2.0/delta_time;
    bdf_coefs[2] = 0.5/delta_time;
    modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
    modelPart.GetProcessInfo().SetValue(AIR_VOLUME_ERROR, 0.0);
    modelPart.GetProcessInfo().SetValue(WATER_VOLUME_ERROR, 0.0);

    // Variables addition
    modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
    modelPart.AddNodalSolutionStepVariable(DENSITY);
    modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
    modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
    modelPart.AddNodalSolutionStepVariable(PRESSURE);
    modelPart.AddNodalSolutionStepVariable(VELOCITY);
    modelPart.AddNodalSolutionStepVariable(FRACTIONAL_VELOCITY);
    modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
    modelPart.AddNodalSolutionStepVariable(DISTANCE);

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
    modelPart.CreateNewElement("TwoFluidNavierStokesFractional3D4N", 1, elemNodes, pElemProp);

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
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY)[k] = 0.1*vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 1)[k] = 0.2 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 2)[k] = 0.3 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k] = 0.0;
            pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
            pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
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

    KRATOS_EXPECT_NEAR(RHS(0),55.5608619747, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(1),64.3649724019, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(2),130.747961753, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(3),1.82139199133, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(4),-173.148878249, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(5),-263.100173027, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(6),-1019.86139378, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(7),-0.231551322358, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(8),-198.517992416, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(9),-320.806050475, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(10),-1153.91463291, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(11),-0.36014021379, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(12),-233.824784725, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(13),-357.542877174, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(14),-1300.79220817, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(15),-1.32970045518 , 1e-7);
}

// /** Checks the TwoFluidNavierStokesFractional3D4N element
//  * Checks the LHS and RHS for a positive element (distance > 0.0)
//  */
KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesFractionalPositiveSide3D4N, FluidDynamicsApplicationFastSuite)
{
    Model current_model;
    ModelPart& modelPart = current_model.CreateModelPart("Main");

    // Process info creation
    const int buffer_size = 4;
    const double delta_time = 0.1;
    modelPart.SetBufferSize(buffer_size);
    modelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
    modelPart.GetProcessInfo().SetValue(TIME, 0.6);
    modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    modelPart.CloneTimeStep(modelPart.GetProcessInfo().GetValue(TIME) + delta_time);
    modelPart.GetProcessInfo().GetPreviousTimeStepInfo(1).SetValue(TIME, 0.0);
    modelPart.GetProcessInfo().GetPreviousTimeStepInfo(1).SetValue(DELTA_TIME, delta_time);
    modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
    modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    Vector bdf_coefs(3);
    bdf_coefs[0] = 3.0/(2.0*delta_time);
    bdf_coefs[1] = -2.0/delta_time;
    bdf_coefs[2] = 0.5/delta_time;
    modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
    modelPart.GetProcessInfo().SetValue(AIR_VOLUME_ERROR, 0.0);
    modelPart.GetProcessInfo().SetValue(WATER_VOLUME_ERROR, 0.0);

    // Variables addition
    modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
    modelPart.AddNodalSolutionStepVariable(DENSITY);
    modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
    modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
    modelPart.AddNodalSolutionStepVariable(PRESSURE);
    modelPart.AddNodalSolutionStepVariable(VELOCITY);
    modelPart.AddNodalSolutionStepVariable(FRACTIONAL_VELOCITY);
    modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
    modelPart.AddNodalSolutionStepVariable(DISTANCE);

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
    modelPart.CreateNewElement("TwoFluidNavierStokesFractional3D4N", 1, elemNodes, pElemProp);

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
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY)[k] = 0.1 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 1)[k] = 0.2 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 2)[k] = 0.3 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k] = 0.0;
            pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
            pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
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
    // std::cout<< std::setprecision(12) << RHS << std::endl;
    KRATOS_EXPECT_NEAR(RHS(0),55.5608619747, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(1), 64.3649724019, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(2), 130.747961753, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(3), 1.82139199133, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(4), -173.148878249, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(5), -263.100173027, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(6), -1019.86139378, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(7), -0.231551322358, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(8), -198.517992416, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(9), -320.806050475, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(10), -1153.91463291, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(11), -0.36014021379, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(12), -233.824784725, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(13), -357.542877174, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(14), -1300.79220817, 1e-7);
    KRATOS_EXPECT_NEAR(RHS(15), -1.32970045518, 1e-7);
}

/** Checks the TwoFluidNavierStokesFractional2D3N element in a hydrostatic case.
    *  Checks the computation of the RHS
    */
KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesFractional2D3NHydrostatic, FluidDynamicsApplicationFastSuite)
{

    Model current_model;
    ModelPart& modelPart = current_model.CreateModelPart("Main");

    // Process info creation
    const int buffer_size = 4;
    const double delta_time = 0.1;
    modelPart.SetBufferSize(buffer_size);
    modelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    modelPart.GetProcessInfo().SetValue(TIME, 0.6);
    modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    modelPart.CloneTimeStep(modelPart.GetProcessInfo().GetValue(TIME) + delta_time);
    modelPart.GetProcessInfo().GetPreviousTimeStepInfo(1).SetValue(TIME, 0.0);
    modelPart.GetProcessInfo().GetPreviousTimeStepInfo(1).SetValue(DELTA_TIME, delta_time);
    modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
    modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
    modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);

    Vector bdf_coefs(3);
    bdf_coefs[0] = 3.0/(2.0*delta_time);
    bdf_coefs[1] = -2.0/delta_time;
    bdf_coefs[2] = 0.5/delta_time;
    modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
    modelPart.GetProcessInfo().SetValue(WATER_VOLUME_ERROR, 0.0);
    modelPart.GetProcessInfo().SetValue(AIR_VOLUME_ERROR, 0.0);

    // Variables addition
    modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
    modelPart.AddNodalSolutionStepVariable(DENSITY);
    modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
    modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
    modelPart.AddNodalSolutionStepVariable(PRESSURE);
    modelPart.AddNodalSolutionStepVariable(VELOCITY);
    modelPart.AddNodalSolutionStepVariable(FRACTIONAL_VELOCITY);
    modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
    modelPart.AddNodalSolutionStepVariable(DISTANCE);

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

    modelPart.CreateNewElement("TwoFluidNavierStokesFractional2D3N", 1, elemNodes1, pElemProp);

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
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY)[k] = 0.1 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 1)[k] = 0.2 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 2)[k] = 0.3 * vel_original(i, k);
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

    const auto& r_process_info = modelPart.GetProcessInfo();
    pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
    pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

    double det;
    MathUtils<double>::InvertMatrix(LHS, LHS, det);

    const Vector solVec = prod(LHS, RHS);

    // The remaining residuals in the velocities have the size of the boundary integrals over the enriched pressure.
    // If the "standard" pressure shape functions are used, the results do not hold.
    // std::cout<< std::setprecision(12) << RHS << std::endl;
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


// Giving a value different to zero to source term in order to test it.

// /** Checks the TwoFluidNavierStokeFractional2D3N element with a source term in mass conservation equation
//  * Checks the LHS and RHS for a cut element
//  */

KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesFractional2D3NError, FluidDynamicsApplicationFastSuite)
{
    Model current_model;
    ModelPart& modelPart = current_model.CreateModelPart("Main");

    // Process info creation
    const int buffer_size = 4;
    const double delta_time = 0.1;
    modelPart.SetBufferSize(buffer_size);
    modelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    modelPart.GetProcessInfo().SetValue(TIME, 0.6);
    modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    modelPart.CloneTimeStep(modelPart.GetProcessInfo().GetValue(TIME) + delta_time);
    modelPart.GetProcessInfo().GetPreviousTimeStepInfo(1).SetValue(TIME, 0.0);
    modelPart.GetProcessInfo().GetPreviousTimeStepInfo(1).SetValue(DELTA_TIME, delta_time);

    // Variables addition
    modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
    modelPart.AddNodalSolutionStepVariable(DENSITY);
    modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
    modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
    modelPart.AddNodalSolutionStepVariable(PRESSURE);
    modelPart.AddNodalSolutionStepVariable(VELOCITY);
    modelPart.AddNodalSolutionStepVariable(FRACTIONAL_VELOCITY);
    modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
    modelPart.AddNodalSolutionStepVariable(DISTANCE);

    modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
    modelPart.GetProcessInfo().SetValue(SOUND_VELOCITY, 1.0e+3);
    modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    Vector bdf_coefs(3);
    bdf_coefs[0] = 3.0/(2.0*delta_time);
    bdf_coefs[1] = -2.0/delta_time;
    bdf_coefs[2] = 0.5/delta_time;
    modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
    modelPart.GetProcessInfo().SetValue(WATER_VOLUME_ERROR,10.0);
    modelPart.GetProcessInfo().SetValue(AIR_VOLUME_ERROR, 10.0);

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
    modelPart.CreateNewElement("TwoFluidNavierStokesFractional2D3N", 1, elemNodes, pElemProp);

    Element::Pointer pElement = modelPart.pGetElement(1);

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
            pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY)[k] = 0.1 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 1)[k] = 0.2 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 2)[k] = 0.3 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k] = 0.0;
            pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
            pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
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
    reference_RHS[0] = 6075.49;
    reference_RHS[1] = 4684.29;
    reference_RHS[2] = -27.713;
    reference_RHS[3] = -6671.02;
    reference_RHS[4] = 1521.5;
    reference_RHS[5] = -16.8578;
    reference_RHS[6] = -241.804;
    reference_RHS[7] = -4409.55;
    reference_RHS[8] = -5.5792;

    KRATOS_EXPECT_VECTOR_NEAR(reference_RHS, RHS, 1e-2);
}

// /** Checks the TwoFluidNavierStokesFractional3D4N element with a source term in mass conservation equation
//  * Checks the LHS and RHS for a cut element
//  */

KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesFractionalCut3D4NError, FluidDynamicsApplicationFastSuite)
{
    Model current_model;
    ModelPart& modelPart = current_model.CreateModelPart("Main");

    // Process info creation
    const int buffer_size = 4;
    const double delta_time = 0.1;
    modelPart.SetBufferSize(buffer_size);
    modelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
    modelPart.GetProcessInfo().SetValue(TIME, 0.6);
    modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    modelPart.CloneTimeStep(modelPart.GetProcessInfo().GetValue(TIME) + delta_time);
    modelPart.GetProcessInfo().GetPreviousTimeStepInfo(1).SetValue(TIME, 0.0);
    modelPart.GetProcessInfo().GetPreviousTimeStepInfo(1).SetValue(DELTA_TIME, delta_time);
    // Variables addition
    modelPart.AddNodalSolutionStepVariable(BODY_FORCE);
    modelPart.AddNodalSolutionStepVariable(DENSITY);
    modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
    modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
    modelPart.AddNodalSolutionStepVariable(PRESSURE);
    modelPart.AddNodalSolutionStepVariable(VELOCITY);
    modelPart.AddNodalSolutionStepVariable(FRACTIONAL_VELOCITY);
    modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
    modelPart.AddNodalSolutionStepVariable(DISTANCE);

    // Process info creation

    modelPart.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
    modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    Vector bdf_coefs(3);
    bdf_coefs[0] = 3.0/(2.0*delta_time);
    bdf_coefs[1] = -2.0/delta_time;
    bdf_coefs[2] = 0.5/delta_time;
    modelPart.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);
    modelPart.GetProcessInfo().SetValue(AIR_VOLUME_ERROR, 10.0);
    modelPart.GetProcessInfo().SetValue(WATER_VOLUME_ERROR, 10.0);

    // Set the element properties
    Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);
    pElemProp->SetValue(DENSITY, 1000.0);
    pElemProp->SetValue(DYNAMIC_VISCOSITY, 1.0e-03);
    NewtonianTwoFluid3DLaw::Pointer pConsLaw(new NewtonianTwoFluid3DLaw());
    pElemProp->SetValue(CONSTITUTIVE_LAW, pConsLaw);

    // Geometry creation
    modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
    std::vector<ModelPart::IndexType> elemNodes {1, 2, 3, 4};
    modelPart.CreateNewElement("TwoFluidNavierStokesFractional3D4N", 1, elemNodes, pElemProp);

    Element::Pointer pElement = modelPart.pGetElement(1);
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
            pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY)[k] = 0.1 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 1)[k] = 0.2 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 2)[k] = 0.3 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY)[k] = 0.0;
            pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
            pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
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
    // Check the RHS values (the RHS is computed as the LHS x previous_solution,
    // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
    Vector reference_RHS = ZeroVector(16);
    reference_RHS[0] = 3885.76546946;
    reference_RHS[1] = 2853.252194541;
    reference_RHS[2] = 3702.11733421;
    reference_RHS[3] = -5.33528927901;
    reference_RHS[4] = -4028.74316648;
    reference_RHS[5] = 302.056712497;
    reference_RHS[6] = -945.707036998;
    reference_RHS[7] = -3.93555403796;
    reference_RHS[8] = 47.4442472421;
    reference_RHS[9] = -3601.67069179;
    reference_RHS[10] = -666.474404348;
    reference_RHS[11] = -3.04929475341;
    reference_RHS[12] = -176.798973004;
    reference_RHS[13] = 455.784932702;
    reference_RHS[14] = -4803.85305063;
    reference_RHS[15] = -4.44652859628;

    KRATOS_EXPECT_VECTOR_NEAR(reference_RHS, RHS, 1e-7);
}

}  // namespace Kratos::Testing.
