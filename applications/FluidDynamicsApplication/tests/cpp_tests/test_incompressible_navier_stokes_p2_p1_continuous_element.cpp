//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes
#include <iomanip> // for std::setprecision

// External includes


// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

// Application includes
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"
#include "custom_elements/incompressible_navier_stokes_p2_p1_continuous.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(IncompressibleNavierStokesP2P1Continuous2D6N, FluidDynamicsApplicationFastSuite)
{
    Model model;
    unsigned int buffer_size = 3;
    auto& r_model_part = model.CreateModelPart("TestModelPart", buffer_size);

    // Variables addition
    r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    r_model_part.AddNodalSolutionStepVariable(PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    r_model_part.AddNodalSolutionStepVariable(REACTION);
    r_model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);

    // ProcessInfo container fill
    double delta_time = 0.1;
    r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    Vector bdf_coefs(3);
    bdf_coefs[0] = 3.0 / (2.0 * delta_time);
    bdf_coefs[1] = -2.0 / delta_time;
    bdf_coefs[2] = 0.5 * delta_time;
    r_model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

    // Set the element properties
    auto p_properties = r_model_part.CreateNewProperties(0);
    p_properties->SetValue(DENSITY, 1.0);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 1.0e-3);
    p_properties->SetValue(CONSTITUTIVE_LAW, Kratos::make_shared<Newtonian2DLaw>());

    // Element creation
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.1, 0.0);
    r_model_part.CreateNewNode(3, 1.0, 0.9, 0.0);
    r_model_part.CreateNewNode(4, 0.5, 0.05, 0.0);
    r_model_part.CreateNewNode(5, 1.0, 0.5, 0.0);
    r_model_part.CreateNewNode(6, 0.5, 0.45, 0.0);

    for (auto it_node = r_model_part.NodesBegin(); it_node < r_model_part.NodesEnd(); ++it_node){
        it_node->AddDof(VELOCITY_X,REACTION_X);
        it_node->AddDof(VELOCITY_Y,REACTION_Y);
        it_node->AddDof(PRESSURE,REACTION_WATER_PRESSURE);
    }

    std::vector<ModelPart::IndexType> element_nodes {1, 2, 3, 4, 5, 6};
    auto p_elem = r_model_part.CreateNewElement("IncompressibleNavierStokesP2P1Continuous2D6N", 1, element_nodes, p_properties);

    // Define and set the nodal values
    Matrix reference_velocity(6,2);
    reference_velocity(0,0) = 0.0; reference_velocity(0,1) = 0.1;
    reference_velocity(1,0) = 0.1; reference_velocity(1,1) = 0.2;
    reference_velocity(2,0) = 0.2; reference_velocity(2,1) = 0.3;
    reference_velocity(3,0) = 0.4; reference_velocity(3,1) = 0.4;
    reference_velocity(4,0) = 0.5; reference_velocity(4,1) = 0.5;
    reference_velocity(5,0) = 0.6; reference_velocity(5,1) = 0.6;

    auto& r_geometry = r_model_part.ElementsBegin()->GetGeometry();
    for(std::size_t i = 0; i < 6; ++i){
        r_geometry[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
        for(std::size_t k = 0; k < 2; ++k){
            r_geometry[i].FastGetSolutionStepValue(VELOCITY)[k] = reference_velocity(i,k);
            r_geometry[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*reference_velocity(i,k);
            r_geometry[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.7*reference_velocity(i,k);
        }
    }

    // Calculate RHS and LHS
    Vector RHS = ZeroVector(15);
    Matrix LHS = ZeroMatrix(15,15);
    const auto& r_process_info = r_model_part.GetProcessInfo();
    p_elem->Initialize(r_process_info);
    p_elem->CalculateLocalSystem(LHS, RHS, r_process_info);

    // std::cout << p_elem->Info() << std::setprecision(12) << std::endl;
    // KRATOS_WATCH(RHS)
    // KRATOS_WATCH(row(LHS,0))

    // Check values
    const std::vector<double> rhs_ref = {-0.00906903569079,-0.0216587001242,-0.0178476996397,-0.0172703814363,0.0834031923607,0.0915008945538,-0.00823181605452,0.0498746884469,0.351823061946,0.322803675622,0.158034060439,0.142938341794,-0.3377346826,-0.0527051752743,0.295439857875};
    const std::vector<double> lhs_0_ref = {0.275040358619,2.13488594544e-06,0.029996360203,-0.0100000622721,0.0293109266242,0.0095065886921,-0.18426639514,0.041493102637,-0.0301319829049,0.000495608465964,-0.127955089696,-0.0414973724089,0.13386705482,-0.000600436672153,6.67151857946e-05};
    KRATOS_EXPECT_VECTOR_NEAR(RHS, rhs_ref, 1.0e-10)
    KRATOS_EXPECT_VECTOR_NEAR(row(LHS,0), lhs_0_ref, 1.0e-10)
}

KRATOS_TEST_CASE_IN_SUITE(IncompressibleNavierStokesP2P1Continuous3D10N, FluidDynamicsApplicationFastSuite)
{
    Model model;
    unsigned int buffer_size = 3;
    auto& r_model_part = model.CreateModelPart("TestModelPart", buffer_size);

    // Variables addition
    r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    r_model_part.AddNodalSolutionStepVariable(PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    r_model_part.AddNodalSolutionStepVariable(REACTION);
    r_model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);

    // ProcessInfo container fill
    double delta_time = 0.1;
    r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    Vector bdf_coefs(3);
    bdf_coefs[0] = 3.0 / (2.0 * delta_time);
    bdf_coefs[1] = -2.0 / delta_time;
    bdf_coefs[2] = 0.5 * delta_time;
    r_model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

    // Set the element properties
    auto p_properties = r_model_part.CreateNewProperties(0);
    p_properties->SetValue(DENSITY, 1.0);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 1.0e-3);
    p_properties->SetValue(CONSTITUTIVE_LAW, Kratos::make_shared<Newtonian3DLaw>());

    // Element creation
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.1, 0.0);
    r_model_part.CreateNewNode(3, 1.0, 0.9, 0.0);
    r_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
    r_model_part.CreateNewNode(5, 0.5, 0.05, 0.0);
    r_model_part.CreateNewNode(6, 1.0, 0.5, 0.0);
    r_model_part.CreateNewNode(7, 0.5, 0.45, 0.0);
    r_model_part.CreateNewNode(8, 0.0, 0.0, 0.5);
    r_model_part.CreateNewNode(9, 0.5, 0.05, 0.5);
    r_model_part.CreateNewNode(10, 0.5, 0.45, 0.5);

    for (auto it_node = r_model_part.NodesBegin(); it_node < r_model_part.NodesEnd(); ++it_node){
        it_node->AddDof(VELOCITY_X,REACTION_X);
        it_node->AddDof(VELOCITY_Y,REACTION_Y);
        it_node->AddDof(PRESSURE,REACTION_WATER_PRESSURE);
    }

    std::vector<ModelPart::IndexType> element_nodes {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    auto p_elem = r_model_part.CreateNewElement("IncompressibleNavierStokesP2P1Continuous3D10N", 1, element_nodes, p_properties);

    // Define and set the nodal values
    Matrix reference_velocity(10,3);
    reference_velocity(0,0) = 0.0; reference_velocity(0,1) = 0.1; reference_velocity(0,2) = 0.2;
    reference_velocity(1,0) = 0.1; reference_velocity(1,1) = 0.2; reference_velocity(1,2) = 0.4;
    reference_velocity(2,0) = 0.2; reference_velocity(2,1) = 0.3; reference_velocity(2,2) = 0.6;
    reference_velocity(3,0) = 0.4; reference_velocity(3,1) = 0.4; reference_velocity(3,2) = 0.8;
    reference_velocity(4,0) = 0.5; reference_velocity(4,1) = 0.5; reference_velocity(4,2) = 1.0;
    reference_velocity(5,0) = 0.6; reference_velocity(5,1) = 0.6; reference_velocity(5,2) = 1.2;
    reference_velocity(6,0) = 0.0; reference_velocity(6,1) = 0.1; reference_velocity(6,2) = 0.2;
    reference_velocity(7,0) = 0.1; reference_velocity(7,1) = 0.2; reference_velocity(7,2) = 0.4;
    reference_velocity(8,0) = 0.2; reference_velocity(8,1) = 0.3; reference_velocity(8,2) = 0.6;
    reference_velocity(9,0) = 0.4; reference_velocity(9,1) = 0.4; reference_velocity(9,2) = 0.8;

    auto& r_geometry = r_model_part.ElementsBegin()->GetGeometry();
    for(std::size_t i = 0; i < 10; ++i){
        r_geometry[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
        for(std::size_t k = 0; k < 3; ++k){
            r_geometry[i].FastGetSolutionStepValue(VELOCITY)[k] = reference_velocity(i,k);
            r_geometry[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*reference_velocity(i,k);
            r_geometry[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.7*reference_velocity(i,k);
        }
    }

    // Calculate RHS and LHS
    Vector RHS = ZeroVector(34);
    Matrix LHS = ZeroMatrix(34,34);
    const auto& r_process_info = r_model_part.GetProcessInfo();
    p_elem->Initialize(r_process_info);
    p_elem->CalculateLocalSystem(LHS, RHS, r_process_info);

    // std::cout << p_elem->Info() << std::setprecision(12) << std::endl;
    // KRATOS_WATCH(RHS)
    // KRATOS_WATCH(row(LHS,0))

    // Check values
    const std::vector<double> rhs_ref = {-0.00657701934375,-0.00620550235972,-0.0136913620863,-0.000170114993796,-0.0143517156581,-0.0157786851078,-0.00860833486305,-0.0150132458839,-0.0189512865026,-0.0130363146785,-0.011378009438,-0.0250112171768,-0.0242266633507,-0.0170353120096,-0.0450749598297,0.0346452605859,0.0609125107186,0.098362051396,0.0160967983649,0.0035688376663,0.0250792755284,0.0237285184092,0.0190965941479,0.0416218137988,0.0449979595942,0.0685522075529,0.120269383925,0.0474728594943,0.0547584464252,0.120148974574,-0.119401045188,-0.0123259565984,0.00846909310257,0.0482579086839};
    const std::vector<double> lhs_0_ref = {0.0257086073875,3.0423641244e-06,0.00522857129039,0.0178570580167,-0.00186127586489,3.38504920949e-05,0.0158010322097,0.00253394525089,-3.28579639389e-06,0.018188732577,0,0.00162630966951,-0.0629346872215,0.00488268545444,-0.00284951737883,0.0131204630345,-0.000669627021879,2.58110017567e-05,-0.0763856921432,-0.00488877018269,-0.00327692855409,-0.0802590504228,1.90236453582e-20,-0.00679850526245,0.0139363248065,-0.00357526747101,0.00275600539288,0.00355834628902,0.00357526747101,0.00325768914512,0.0207605910311,-0.00752233157666,-0.00657159278778,-0.00666666666667};
    KRATOS_EXPECT_VECTOR_NEAR(RHS, rhs_ref, 1.0e-10)
    KRATOS_EXPECT_VECTOR_NEAR(row(LHS,0), lhs_0_ref, 1.0e-10)
}

}  // namespace Kratos::Testing