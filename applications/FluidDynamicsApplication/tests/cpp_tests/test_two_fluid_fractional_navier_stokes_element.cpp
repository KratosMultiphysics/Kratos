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
    const int buffer_size = 3;
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
    const std::vector<double> rhs_ref = {285.327812743,-75.7286708218,0.678487936315,-758.586246712,-529.046820664,-0.808596445975,-549.320499823,-931.191955616,-0.0198914903401};
    KRATOS_EXPECT_VECTOR_NEAR(RHS, rhs_ref, 1.0e-7);

}

// /** Checks the TwoFluidNavierStokesFractional3D4N element
//  * Checks the LHS and RHS for a cut element
//  */
KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesFractionalCut3D4N, FluidDynamicsApplicationFastSuite)
{
    Model current_model;
    ModelPart& modelPart = current_model.CreateModelPart("Main");

    // Process info creation
    const int buffer_size = 3;
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
    // Check the RHS values (the RHS is computed as the LHS x previous_solution,
    // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
    const std::vector<double> rhs_ref = {-182.780712656374,-5.97373013000447,-366.425812873773,0.13618045150788,-65.5996318561992,-195.157657104678,-732.899349900163,0.228653068676312,-30.1350973710576,-309.557107994407,-744.055344313566,-0.182507276303036,4.70728308941767,-259.815870552732,-872.015239241193,-0.282326243881156};
    KRATOS_EXPECT_VECTOR_NEAR(RHS, rhs_ref, 1.0e-7);
}


// /** Checks the TwoFluidNavierStokesFractional3D4N element
//  * Checks the LHS and RHS for a negative element (distance <= 0.0)
//  */
KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesFractionalNegativeSide3D4N, FluidDynamicsApplicationFastSuite)
{
    Model current_model;
    ModelPart& modelPart = current_model.CreateModelPart("Main");

    // Process info creation
    const int buffer_size = 3;
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
    const std::vector<double> rhs_ref={55.5608619747481,64.3649724019167,130.747961752768,1.82139199133129,-173.148878249225,-263.100173026906,-1019.86139377883,-0.231551322358482,-198.517992416071,-320.806050474686,-1153.91463291277,-0.360140213789913,-233.824784724542,-357.542877174272,-1300.79220817092,-1.3297004551829};
    KRATOS_EXPECT_VECTOR_NEAR(RHS, rhs_ref, 1.0e-7);
}

// /** Checks the TwoFluidNavierStokesFractional3D4N element
//  * Checks the LHS and RHS for a positive element (distance > 0.0)
//  */
KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesFractionalPositiveSide3D4N, FluidDynamicsApplicationFastSuite)
{
    Model current_model;
    ModelPart& modelPart = current_model.CreateModelPart("Main");

    // Process info creation
    const int buffer_size = 3;
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
    const std::vector<double> rhs_ref={55.5608619747481,64.3649724019167,130.747961752768,1.82139199133129,-173.148878249225,-263.100173026906,-1019.86139377883,-0.231551322358482,-198.517992416071,-320.806050474686,-1153.91463291277,-0.360140213789913,-233.824784724542,-357.542877174272,-1300.79220817092,-1.3297004551829};
    KRATOS_EXPECT_VECTOR_NEAR(RHS, rhs_ref, 1.0e-7);
}

/** Checks the TwoFluidNavierStokesFractional2D3N element in a hydrostatic case.
    *  Checks the computation of the RHS
    */
KRATOS_TEST_CASE_IN_SUITE(ElementTwoFluidNavierStokesFractional2D3NHydrostatic, FluidDynamicsApplicationFastSuite)
{

    Model current_model;
    ModelPart& modelPart = current_model.CreateModelPart("Main");

    // Process info creation
    const int buffer_size = 3;
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
    const std::vector<double> rhs_ref={0,-17.5,-5.6843418860808e-14,7.5,4.44089209850063e-16,2.8421709430404e-14,-7.5,-7.5,7.105427357601e-15};
    KRATOS_EXPECT_VECTOR_NEAR(RHS, rhs_ref, 1.0e-7);
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
    const std::vector<double> rhs_ref = {6075.49466849675,4684.29040473282,-27.7130224974724,-6671.02384528602,1521.49790451337,-16.8577767803168,-241.803857305756,-4409.54852590944,-5.57920072221086};
    KRATOS_EXPECT_VECTOR_NEAR(RHS, rhs_ref, 1.0e-10);
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
    // hence, it is assumd that if the RHS is correct, the LHS is correct as well)
    const std::vector<double> rhs_ref = {3885.76546946416,2853.25219454485,3702.11733420833,-5.33528927901302,-4028.74316647894,302.056712496537,-945.70703699794,-3.93555403796184,47.44424724214,-3601.67069178779,-666.474404348206,-3.04929475340806,-176.798973004457,455.784932701568,-4803.85305063243,-4.44652859628375};
    KRATOS_EXPECT_VECTOR_NEAR(RHS, rhs_ref, 1.0e-10);
}

}  // namespace Kratos::Testing.
