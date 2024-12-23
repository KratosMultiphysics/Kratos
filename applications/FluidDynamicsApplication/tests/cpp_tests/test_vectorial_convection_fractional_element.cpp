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
#include "containers/model.h"
#include "spaces/ublas_space.h"
#include "includes/properties.h"
#include "includes/model_part.h"
#include "includes/global_pointer_variables.h"
#include "utilities/math_utils.h"
#include "utilities/normal_calculation_utils.h"

// Application includes
#include "custom_elements/vectorial_convection_fractional_element.h"
#include "processes/find_nodal_neighbours_process.h"
#include "tests/cpp_tests/fluid_dynamics_fast_suite.h"

namespace Kratos::Testing {

typedef ModelPart::IndexType                             IndexType;
typedef ModelPart::NodeIterator                          NodeIteratorType;

/** Checks the VectorialConvectionFractional2D3N element.
 * Checks the LHS and RHS computation
 */

KRATOS_TEST_CASE_IN_SUITE(ElementVectorialConvectionFractional2D3N, FluidDynamicsApplicationFastSuite)
{

    Model current_model;
    ModelPart& modelPart = current_model.CreateModelPart("Main");
    // Variables addition
    modelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
    modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
    modelPart.AddNodalSolutionStepVariable(VELOCITY);
    modelPart.AddNodalSolutionStepVariable(FRACTIONAL_VELOCITY);
    modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);

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

    // Set the element properties
    Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);

    // Geometry creation
    modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};
    modelPart.CreateNewElement("VectorialConvectionFractionalElement2D3N", 1, elemNodes, pElemProp);

    Element::Pointer pElement = modelPart.pGetElement(1);

    // Define the nodal values
    Matrix vel_original(3,2);
    vel_original(0,0) = 0.0; vel_original(0,1) = 0.1;
    vel_original(1,0) = 0.1; vel_original(1,1) = 0.2;
    vel_original(2,0) = 0.2; vel_original(2,1) = 0.3;

    for(unsigned int i=0; i<3; i++){
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

    // Compute RHS and LHS
    Vector RHS = ZeroVector(6);
    Matrix LHS = ZeroMatrix(6,6);

    const auto& r_process_info = modelPart.GetProcessInfo();
    pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
    pElement->CalculateLocalSystem(LHS, RHS, r_process_info);
    // Check the RHS values (the RHS is computed as the LHS x previous_solution,
    // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
    const std::vector<double> rhs_ref={0.0897168467616033,0.176984086785688,0.39861648657173,0.786349246547645,0.475841396524262,0.938690536488135};
    KRATOS_EXPECT_VECTOR_NEAR(RHS, rhs_ref, 1.0e-10)
}

// /** Checks the VectorialConvectionFractional3D4N element
//  * Checks the LHS and RHS for a cut element
//  */
KRATOS_TEST_CASE_IN_SUITE(ElementVectorialConvectionFractional3D4N, FluidDynamicsApplicationFastSuite)
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

    // Variables addition

    modelPart.AddNodalSolutionStepVariable(DYNAMIC_TAU);
    modelPart.AddNodalSolutionStepVariable(VELOCITY);
    modelPart.AddNodalSolutionStepVariable(FRACTIONAL_VELOCITY);
    modelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);

    // Set the element properties
    Properties::Pointer pElemProp = modelPart.CreateNewProperties(0);

    // Geometry creation
    modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    modelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
    std::vector<ModelPart::IndexType> elemNodes {1, 2, 3, 4};
    modelPart.CreateNewElement("VectorialConvectionFractionalElement3D4N", 1, elemNodes, pElemProp);

    Element::Pointer pElement = modelPart.pGetElement(1);

    // Define the nodal values
    Matrix vel_original(4,3);
    vel_original(0,0) = 0.0; vel_original(0,1) = 0.1; vel_original(0,2) = 0.2;
    vel_original(1,0) = 0.1; vel_original(1,1) = 0.2; vel_original(1,2) = 0.3;
    vel_original(2,0) = 0.2; vel_original(2,1) = 0.3; vel_original(2,2) = 0.4;
    vel_original(3,0) = 0.3; vel_original(3,1) = 0.4; vel_original(3,2) = 0.5;

    for(unsigned int i=0; i<4; i++){
        for(unsigned int k=0; k<3; k++){
            pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel_original(i,k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i,k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*vel_original(i,k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY)[k] = 0.1 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 1)[k] = 0.2 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 2)[k] = 0.3 * vel_original(i, k);
            pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY,0)[k]    = 0.0;
            pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
            pElement->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY, 2)[k] = 0.0;
        }
    }

    // Compute RHS and LHS
    Vector RHS = ZeroVector(12);
    Matrix LHS = ZeroMatrix(12,12);

    const auto& r_process_info = modelPart.GetProcessInfo();
    pElement->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
    pElement->CalculateLocalSystem(LHS, RHS, r_process_info);

    // Check the RHS values (the RHS is computed as the LHS x previous_solution,
    // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
    const std::vector<double> rhs_ref={-0.0137395522001275,-0.0223513234273902,-0.0309630946546529,0.148963526100064,0.242331911713695,0.335700297327326,0.167041645911196,0.271741160062705,0.376440674214213,0.185119765722328,0.301150408411714,0.4171810511011};
    KRATOS_EXPECT_VECTOR_NEAR(RHS, rhs_ref, 1.0e-10)
    
}

} // namespace Kratos::Testing