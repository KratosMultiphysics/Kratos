//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Gianmarco Boscolo
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

KRATOS_TEST_CASE_IN_SUITE(FluidTopologyOptimizationElement2D3N, FluidDynamicsApplicationFastSuite)
{
    Model model;
    std::size_t buffer_size = 3;
    auto& r_model_part = model.CreateModelPart("TestModelPart", buffer_size);
    // Variables addition
    // Navier-Stokes
    r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    r_model_part.AddNodalSolutionStepVariable(DENSITY);
    r_model_part.AddNodalSolutionStepVariable(PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(REACTION);
    r_model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);
    // ADJOINT Navier-Stokes
    r_model_part.AddNodalSolutionStepVariable(BODY_FORCE_ADJ);
    r_model_part.AddNodalSolutionStepVariable(PRESSURE_ADJ);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY_ADJ);
    r_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY_ADJ);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION_ADJ);

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

    std::vector<ModelPart::IndexType> element_nodes {1, 2, 3};
    auto p_elem = r_model_part.CreateNewElement("FluidTopologyOptimizationElement2D3N", 1, element_nodes, p_properties);

    auto& r_geometry = r_model_part.ElementsBegin()->GetGeometry();

    const auto& r_process_info = r_model_part.GetProcessInfo();

    for (auto it_node = r_model_part.NodesBegin(); it_node < r_model_part.NodesEnd(); ++it_node){
        it_node->AddDof(VELOCITY_X);
        it_node->AddDof(VELOCITY_Y);
        it_node->AddDof(VELOCITY_Z);
        it_node->AddDof(PRESSURE);
    }

    // Define and set the nodal values
    // NS
    Matrix reference_velocity(3,2);
    reference_velocity(0,0) = 0.1; reference_velocity(0,1) = 0.1;
    reference_velocity(1,0) = 0.2; reference_velocity(1,1) = 0.2;
    reference_velocity(2,0) = 0.3; reference_velocity(2,1) = 0.3;

    for (std::size_t i=0; i<3; i++){
        r_geometry[i].SetValue(SOUND_VELOCITY, 1e15);
        r_geometry[i].FastGetSolutionStepValue(DENSITY) = 1.0;
        r_geometry[i].FastGetSolutionStepValue(PRESSURE) = 0.0;
        r_geometry[i].SetValue(RESISTANCE, 0.0); // RESISTANCE IS A NODAL VARIABLE
        r_geometry[i].FastGetSolutionStepValue(PRESSURE, 1) = 0.0;
        r_geometry[i].FastGetSolutionStepValue(PRESSURE, 2) = 0.0;
        for (std::size_t k=0; k<2; k++){
            r_geometry[i].FastGetSolutionStepValue(VELOCITY)[k] = 1.0*reference_velocity(i,k);
            r_geometry[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.0*reference_velocity(i,k);
            r_geometry[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.0*reference_velocity(i,k);
        }
    }

    // Set the elemental values
    // p_elem->SetValue(RESISTANCE, 0.0); IT IS A NODAL VALUE

    // Calculate RHS and LHS
    Vector RHS = ZeroVector(9);
    Matrix LHS = ZeroMatrix(9,9);

    // SOLVE NAVIER-STOKES
    r_model_part.GetProcessInfo().SetValue(FLUID_TOP_OPT_PROBLEM_STAGE, 1);
    p_elem->Initialize(r_process_info);
    p_elem->CalculateLocalSystem(LHS, RHS, r_process_info);

    std::cout << "\n";
    std::cout << "\n-------------------\n";
    std::cout << p_elem->Info() << std::setprecision(12) << std::endl;
    std::cout << "NAVIER-STOKES\n";
    bool print_rhs_and_lhs = false;
        if (print_rhs_and_lhs)
        {
            KRATOS_WATCH(RHS)
            for (int i = 0; i < 9; i++)
            {
                KRATOS_WATCH(row(LHS,i))
            }
            std::cout << "-------------------\n";
            std::cout << "\n";
        }

    // Check values
    bool check_ns = false;
    if (check_ns)
    {
        const std::vector<double> rhs_ref = {0.303417198205,0.403417198205,-0.0232176605582,-0.16277452246,-0.15527452246,-0.0258911697209,-0.155910026661,-0.263410026661,-0.0258911697209}; // FluidTopologyOptimization2D3N-NS_RHS
        const std::vector<double> lhs_0_ref = {1.66444798803,0.691666666667,0.0877773348295,-1.17805732735,-0.5,0.0788893318372,-0.486390660682,-0.191666666667,0.0833333333333}; // FluidTopologyOptimization2D3N-NS_LHS_0
        const std::vector<double> lhs_1_ref = {0.691666666667,1.66444798803,0.0877773348295,-0.191666666667,-0.486390660682,0.0833333333333,-0.5,-1.17805732735,0.0788893318372}; // FluidTopologyOptimization2D3N-NS_LHS_1
        const std::vector<double> lhs_2_ref = {-0.0773922018607,-0.0773922018607,0.0297647171318,0.080362767597,-0.00297056573633,-0.0148823585659,-0.00297056573633,0.080362767597,-0.0148823585659}; // FluidTopologyOptimization2D3N-NS_LHS_2
        const std::vector<double> lhs_3_ref = {-1.22683014973,-0.191666666667,-0.087792064617,1.2092484082,0,-0.0788746020497,0.017581741532,0.191666666667,-0.0833333333333}; // FluidTopologyOptimization2D3N-NS_LHS_3
        const std::vector<double> lhs_4_ref = {-0.5,-0.535163483064,-0.00445873128368,0,0.517581741532,0,0.5,0.017581741532,0.00445873128368}; // FluidTopologyOptimization2D3N-NS_LHS_4
        const std::vector<double> lhs_5_ref = {-0.089274464806,-0.0833333333333,-0.0148823585659,0.0863038990697,0,0.0148823585659,0.00297056573633,0.0833333333333,0}; // FluidTopologyOptimization2D3N-NS_LHS_5
        const std::vector<double> lhs_6_ref = {-0.5394001777417794,-0.5,-0.00444997778227757,0.01970008887088974,0.5,0.00444997778227757,0.5197000888708898,0,0}; // FluidTopologyOptimization2D3N-NS_LHS_6
        const std::vector<double> lhs_7_ref = {-0.191666666667,-1.23106684441,-0.0877833111156,0.191666666667,0.0197000888709,-0.0833333333333,0,1.21136675554,-0.0788833555511}; // FluidTopologyOptimization2D3N-NS_LHS_7
        const std::vector<double> lhs_8_ref = {-0.0833333333333,-0.089274464806,-0.0148823585659,0.0833333333333,0.00297056573633,0,0,0.0863038990697,0.0148823585659}; // FluidTopologyOptimization2D3N-NS_LHS_8

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
    
    bool test_adjoint = true;
    if (test_adjoint)
    {

        for (auto it_node = r_model_part.NodesBegin(); it_node < r_model_part.NodesEnd(); ++it_node){
            it_node->AddDof(VELOCITY_ADJ_X);
            it_node->AddDof(VELOCITY_ADJ_Y);
            it_node->AddDof(VELOCITY_ADJ_Z);
            it_node->AddDof(PRESSURE_ADJ);
        }

        // Define and set the nodal values
        // ADJ NS
        Matrix reference_velocity_adj(3,2);
        reference_velocity_adj(0,0) = 0.1; reference_velocity_adj(0,1) = 0.1;
        reference_velocity_adj(1,0) = 0.2; reference_velocity_adj(1,1) = 0.2;
        reference_velocity_adj(2,0) = 0.3; reference_velocity_adj(2,1) = 0.3;

        for (std::size_t i=0; i<3; i++){
            r_geometry[i].SetValue(SOUND_VELOCITY, 1e15);
            r_geometry[i].FastGetSolutionStepValue(DENSITY) = 1.0;
            r_geometry[i].FastGetSolutionStepValue(PRESSURE_ADJ) = 0.0;
            r_geometry[i].SetValue(RESISTANCE, 0.0); // RESISTANCE IS A NODAL VARIABLE
            r_geometry[i].FastGetSolutionStepValue(PRESSURE_ADJ, 1) = 0.0;
            r_geometry[i].FastGetSolutionStepValue(PRESSURE_ADJ, 2) = 0.0;
            for (std::size_t k=0; k<2; k++){
                r_geometry[i].FastGetSolutionStepValue(VELOCITY_ADJ)[k] = 1.0*reference_velocity_adj(i,k);
                r_geometry[i].FastGetSolutionStepValue(VELOCITY_ADJ, 1)[k] = 0.0*reference_velocity_adj(i,k);
                r_geometry[i].FastGetSolutionStepValue(VELOCITY_ADJ, 2)[k] = 0.0*reference_velocity_adj(i,k);
            }
        }

        // Calculate RHS_ADJ and LHS_ADJ
        Vector RHS_ADJ = ZeroVector(9);
        Matrix LHS_ADJ = ZeroMatrix(9,9);

        // SOLVE ADJ NAVIER-STOKES
        r_model_part.GetProcessInfo().SetValue(FLUID_TOP_OPT_PROBLEM_STAGE, 2);
        p_elem->Initialize(r_process_info);
        p_elem->CalculateLocalSystem(LHS_ADJ, RHS_ADJ, r_process_info);

        std::cout << "\n";
        std::cout << "\n-------------------\n";
        std::cout << p_elem->Info() << std::setprecision(12) << std::endl;
        std::cout << "ADJOINT NAVIER-STOKES\n";
        bool print_rhs_and_lhs_adj = false;
        if (print_rhs_and_lhs_adj)
        {
            KRATOS_WATCH(RHS_ADJ)
            for (int i = 0; i < 9; i++)
            {
                KRATOS_WATCH(row(LHS_ADJ,i))
            }
            std::cout << "-------------------\n";
            std::cout << "\n";
        }
        

        // Check values
        bool check_adj = true;
        if (check_adj)
        {
            // const std::vector<double> rhs_adj_ref = {0.303417198205,0.403417198205,-0.0232176605582,-0.16277452246,-0.15527452246,-0.0258911697209,-0.155910026661,-0.263410026661,-0.0258911697209}; // FluidTopologyOptimization2D3N-NS_RHS
            // const std::vector<double> lhs_adj_0_ref = {1.66444798803,0.691666666667,0.0877773348295,-1.17805732735,-0.5,0.0788893318372,-0.486390660682,-0.191666666667,0.0833333333333}; // FluidTopologyOptimization2D3N-NS_LHS_0
            // const std::vector<double> lhs_adj_1_ref = {0.691666666667,1.66444798803,0.0877773348295,-0.191666666667,-0.486390660682,0.0833333333333,-0.5,-1.17805732735,0.0788893318372}; // FluidTopologyOptimization2D3N-NS_LHS_1
            // const std::vector<double> lhs_adj_2_ref = {-0.0773922018607,-0.0773922018607,0.0297647171318,0.080362767597,-0.00297056573633,-0.0148823585659,-0.00297056573633,0.080362767597,-0.0148823585659}; // FluidTopologyOptimization2D3N-NS_LHS_2
            // const std::vector<double> lhs_adj_3_ref = {-1.22683014973,-0.191666666667,-0.087792064617,1.2092484082,0,-0.0788746020497,0.017581741532,0.191666666667,-0.0833333333333}; // FluidTopologyOptimization2D3N-NS_LHS_3
            // const std::vector<double> lhs_adj_4_ref = {-0.5,-0.535163483064,-0.00445873128368,0,0.517581741532,0,0.5,0.017581741532,0.00445873128368}; // FluidTopologyOptimization2D3N-NS_LHS_4
            // const std::vector<double> lhs_adj_5_ref = {-0.089274464806,-0.0833333333333,-0.0148823585659,0.0863038990697,0,0.0148823585659,0.00297056573633,0.0833333333333,0}; // FluidTopologyOptimization2D3N-NS_LHS_5
            // const std::vector<double> lhs_adj_6_ref = {-0.5394001777417794,-0.5,-0.00444997778227757,0.01970008887088974,0.5,0.00444997778227757,0.5197000888708898,0,0}; // FluidTopologyOptimization2D3N-NS_LHS_6
            // const std::vector<double> lhs_adj_7_ref = {-0.191666666667,-1.23106684441,-0.0877833111156,0.191666666667,0.0197000888709,-0.0833333333333,0,1.21136675554,-0.0788833555511}; // FluidTopologyOptimization2D3N-NS_LHS_7
            // const std::vector<double> lhs_adj_8_ref = {-0.0833333333333,-0.089274464806,-0.0148823585659,0.0833333333333,0.00297056573633,0,0,0.0863038990697,0.0148823585659}; // FluidTopologyOptimization2D3N-NS_LHS_8

            KRATOS_EXPECT_VECTOR_NEAR(RHS_ADJ, RHS, 1.0e-8)
            KRATOS_EXPECT_VECTOR_NEAR(row(LHS_ADJ,0), row(LHS,0), 1.0e-8)
            KRATOS_EXPECT_VECTOR_NEAR(row(LHS_ADJ,1), row(LHS,1), 1.0e-8)
            KRATOS_EXPECT_VECTOR_NEAR(row(LHS_ADJ,2), row(LHS,2), 1.0e-8)
            KRATOS_EXPECT_VECTOR_NEAR(row(LHS_ADJ,3), row(LHS,3), 1.0e-8)
            KRATOS_EXPECT_VECTOR_NEAR(row(LHS_ADJ,4), row(LHS,4), 1.0e-8)
            KRATOS_EXPECT_VECTOR_NEAR(row(LHS_ADJ,5), row(LHS,5), 1.0e-8)
            KRATOS_EXPECT_VECTOR_NEAR(row(LHS_ADJ,6), row(LHS,6), 1.0e-8)
            KRATOS_EXPECT_VECTOR_NEAR(row(LHS_ADJ,7), row(LHS,7), 1.0e-8)
            KRATOS_EXPECT_VECTOR_NEAR(row(LHS_ADJ,8), row(LHS,8), 1.0e-8)
        }
    }
}

}  // namespace Kratos::Testing