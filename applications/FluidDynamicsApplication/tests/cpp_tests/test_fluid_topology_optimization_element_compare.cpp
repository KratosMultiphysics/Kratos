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
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "tests/cpp_tests/fluid_dynamics_fast_suite.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(FluidTopologyOptimizationElement2D3N_compare, FluidDynamicsApplicationFastSuite)
{
    Model model;
    std::size_t buffer_size = 3;
    auto& r_model_part = model.CreateModelPart("TestModelPart", buffer_size);

    // Variables addition
    r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    r_model_part.AddNodalSolutionStepVariable(DENSITY);
    r_model_part.AddNodalSolutionStepVariable(PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(REACTION);
    r_model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);

    // ProcessInfo container fill
    double delta_time = 1e15;
    r_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 1.0);
    r_model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    Vector bdf_coefs(3);
    bdf_coefs[0] = 3.0 / (2.0 * delta_time);
    bdf_coefs[1] = -2.0 / delta_time;
    bdf_coefs[2] = 0.5 / delta_time;
    for (int i = 0; i< 3; i++)
    {
        bdf_coefs[i] = 0.0;
    }
    r_model_part.GetProcessInfo().SetValue(BDF_COEFFICIENTS, bdf_coefs);

    // Set the element properties
    auto p_properties = r_model_part.CreateNewProperties(0);
    p_properties->SetValue(DENSITY, 1.0);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 1.0);
    auto p_cons_law = Kratos::make_shared<Newtonian2DLaw>();
    p_properties->SetValue(CONSTITUTIVE_LAW, p_cons_law);

    // Element creation
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 0.5, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 0.0, 0.5, 0.0);

    for (auto it_node = r_model_part.NodesBegin(); it_node < r_model_part.NodesEnd(); ++it_node){
        it_node->AddDof(VELOCITY_X);
        it_node->AddDof(VELOCITY_Y);
        it_node->AddDof(VELOCITY_Z);
        it_node->AddDof(PRESSURE);
    }

    std::vector<ModelPart::IndexType> element_nodes {1, 2, 3};
    auto p_elem = r_model_part.CreateNewElement("WeaklyCompressibleNavierStokes2D3N", 1, element_nodes, p_properties);
    p_elem->SetValue(RESISTANCE, 0.0); // IT IS A NODAL VALUE

    // Define and set the nodal values
    Matrix reference_velocity(3,2);
    reference_velocity(0,0) = 0.1; reference_velocity(0,1) = 0.1;
    reference_velocity(1,0) = 0.2; reference_velocity(1,1) = 0.2;
    reference_velocity(2,0) = 0.3; reference_velocity(2,1) = 0.3;

    auto& r_geometry = r_model_part.ElementsBegin()->GetGeometry();
    for (std::size_t i=0; i<3; i++){
        r_geometry[i].SetValue(SOUND_VELOCITY, 1e15);
        r_geometry[i].FastGetSolutionStepValue(DENSITY) = 1.0;
        r_geometry[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
        r_geometry[i].FastGetSolutionStepValue(PRESSURE, 1) = 0.0;
        r_geometry[i].FastGetSolutionStepValue(PRESSURE, 2) = 0.0;
        for (std::size_t k=0; k<2; k++){
            r_geometry[i].FastGetSolutionStepValue(VELOCITY)[k] = 1.0*reference_velocity(i,k);
            r_geometry[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.0*reference_velocity(i,k);
            r_geometry[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.0*reference_velocity(i,k);
        }
    }


    // Calculate RHS and LHS
    Vector RHS = ZeroVector(9);
    Matrix LHS = ZeroMatrix(9,9);
    const auto& r_process_info = r_model_part.GetProcessInfo();
    p_elem->Initialize(r_process_info);
    p_elem->CalculateLocalSystem(LHS, RHS, r_process_info);

    std::cout << "\n";
    std::cout << "\n-------------------\n";
    std::cout << p_elem->Info() << std::setprecision(12) << std::endl;
    KRATOS_WATCH(RHS)
    for (int i = 0; i < 9; i++)
    {
        KRATOS_WATCH(row(LHS,i))
    }
    std::cout << "-------------------\n";
    std::cout << "\n";

    // Check values
    const std::vector<double> rhs_ref = {0.303417198205,0.403417198205,-0.0232176605582,-0.16277452246,-0.15527452246,-0.0258911697209,-0.155910026661,-0.263410026661,-0.0258911697209}; // WeaklyCompressibleNavierStokes2D3N_RHS
    const std::vector<double> lhs_0_ref = {1.66444798803,0.691666666667,0.0877773348295,-1.17805732735,-0.5,0.0788893318372,-0.486390660682,-0.191666666667,0.0833333333333}; // WeaklyCompressibleNavierStokes2D3N_LHS_0
    const std::vector<double> lhs_1_ref = {0.691666666667,1.66444798803,0.0877773348295,-0.191666666667,-0.486390660682,0.0833333333333,-0.5,-1.17805732735,0.0788893318372}; // WeaklyCompressibleNavierStokes2D3N_LHS_1
    const std::vector<double> lhs_2_ref = {-0.0773922018607,-0.0773922018607,0.0297647171318,0.080362767597,-0.00297056573633,-0.0148823585659,-0.00297056573633,0.080362767597,-0.0148823585659}; // WeaklyCompressibleNavierStokes2D3N_LHS_2
    const std::vector<double> lhs_3_ref = {-1.22683014973,-0.191666666667,-0.087792064617,1.2092484082,0,-0.0788746020497,0.017581741532,0.191666666667,-0.0833333333333}; // WeaklyCompressibleNavierStokes2D3N_LHS_3
    const std::vector<double> lhs_4_ref = {-0.5,-0.535163483064,-0.00445873128368,0,0.517581741532,0,0.5,0.017581741532,0.00445873128368}; // WeaklyCompressibleNavierStokes2D3N_LHS_4
    const std::vector<double> lhs_5_ref = {-0.089274464806,-0.0833333333333,-0.0148823585659,0.0863038990697,0,0.0148823585659,0.00297056573633,0.0833333333333,0}; // WeaklyCompressibleNavierStokes2D3N_LHS_5
    const std::vector<double> lhs_6_ref = {-0.5394001777417794,-0.5,-0.00444997778227757,0.01970008887088974,0.5,0.00444997778227757,0.5197000888708898,0,0}; // WeaklyCompressibleNavierStokes2D3N_LHS_6
    const std::vector<double> lhs_7_ref = {-0.191666666667,-1.23106684441,-0.0877833111156,0.191666666667,0.0197000888709,-0.0833333333333,0,1.21136675554,-0.0788833555511}; // WeaklyCompressibleNavierStokes2D3N_LHS_7
    const std::vector<double> lhs_8_ref = {-0.0833333333333,-0.089274464806,-0.0148823585659,0.0833333333333,0.00297056573633,0,0,0.0863038990697,0.0148823585659}; // WeaklyCompressibleNavierStokes2D3N_LHS_8

    KRATOS_EXPECT_VECTOR_NEAR(RHS, rhs_ref, 1.0e-8)
    KRATOS_EXPECT_VECTOR_NEAR(row(LHS,0), lhs_0_ref, 1.0e-8)
    KRATOS_EXPECT_VECTOR_NEAR(row(LHS,1), lhs_1_ref, 1.0e-8)
    KRATOS_EXPECT_VECTOR_NEAR(row(LHS,2), lhs_2_ref, 1.0e-8)
    KRATOS_EXPECT_VECTOR_NEAR(row(LHS,3), lhs_3_ref, 1.0e-8)
    KRATOS_EXPECT_VECTOR_NEAR(row(LHS,4), lhs_4_ref, 1.0e-8)
    KRATOS_EXPECT_VECTOR_NEAR(row(LHS,5), lhs_5_ref, 1.0e-8)
    KRATOS_EXPECT_VECTOR_NEAR(row(LHS,6), lhs_6_ref, 1.0e-8)
    KRATOS_EXPECT_VECTOR_NEAR(row(LHS,7), lhs_7_ref, 1.0e-8)
    KRATOS_EXPECT_VECTOR_NEAR(row(LHS,8), lhs_8_ref, 1.0e-8)
}

}  // namespace Kratos::Testing