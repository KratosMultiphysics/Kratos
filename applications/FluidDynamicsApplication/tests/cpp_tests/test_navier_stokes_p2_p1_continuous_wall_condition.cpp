//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//


// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

// Application includes
#include "fluid_dynamics_application.h"

namespace Kratos::Testing {

namespace
{
    auto CreateTestingNavierStokesWallCondition2D3N(
        const std::string ConditionName,
        ModelPart& rModelPart)
    {
        // Add required nodal variables
        rModelPart.AddNodalSolutionStepVariable(NORMAL);
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(PRESSURE);
        rModelPart.AddNodalSolutionStepVariable(BODY_FORCE);
        rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);

        // Set the required properties
        // Note that the condition DENSITY is retrieved from the element one
        auto p_properties_0 = rModelPart.CreateNewProperties(0);
        auto p_properties_1 = rModelPart.CreateNewProperties(1);
        p_properties_1->SetValue(DENSITY, 1000.0);
        p_properties_1->SetValue(DYNAMIC_VISCOSITY, 2.0);
        auto p_cons_law = Kratos::make_shared<Newtonian2DLaw>();
        p_properties_1->SetValue(CONSTITUTIVE_LAW, p_cons_law);

        // Create a fake element to serve as parent of current testing condition
        rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
        rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
        rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
        rModelPart.CreateNewNode(4, 0.5, 0.0, 0.0);
        rModelPart.CreateNewNode(5, 0.5, 0.5, 0.0);
        rModelPart.CreateNewNode(6, 0.0, 0.5, 0.0);
        auto p_element = rModelPart.CreateNewElement("IncompressibleNavierStokesP2P1Continuous2D6N", 1, {1,2,3,4,5,6}, p_properties_1);
        p_element->Initialize(rModelPart.GetProcessInfo()); // Initialize constitutive law

        // Create the testing condition
        auto p_test_condition = rModelPart.CreateNewCondition(ConditionName, 1, {{3,1,2}}, p_properties_0);

        // Add DOFs
        for (auto& r_node : rModelPart.Nodes()){
            r_node.AddDof(VELOCITY_X);
            r_node.AddDof(VELOCITY_Y);
            r_node.AddDof(PRESSURE);
        }

        // Manually set the NORMALS
        array_1d<double,3> aux_normal = ZeroVector(3);
        aux_normal[0] = -1.0;
        rModelPart.pGetNode(1)->FastGetSolutionStepValue(NORMAL) = aux_normal;
        rModelPart.pGetNode(3)->FastGetSolutionStepValue(NORMAL) = aux_normal;

        // Manually set the element as condition neighbour
        GlobalPointersVector<Element> neigh_vect;
        neigh_vect.resize(1);
        neigh_vect(0) = p_element;
        p_test_condition->SetValue(NEIGHBOUR_ELEMENTS, neigh_vect);

        return p_test_condition;
    }
}

KRATOS_TEST_CASE_IN_SUITE(NavierStokesP2P1ContinuousWallCondition2D3NZero, FluidDynamicsApplicationFastSuite)
{
    // Create the test model part
    Model model;
    std::size_t buffer_size = 2;
    auto& r_model_part = model.CreateModelPart("TestModelPart",buffer_size);

    // Create the testing condition
    auto p_test_condition = CreateTestingNavierStokesWallCondition2D3N("NavierStokesP2P1ContinuousWallCondition2D3N", r_model_part);

    // Set the testing nodal values
    array_1d<double,3> aux_v = ZeroVector(3);
    for (auto& r_node: r_model_part.Nodes()) {
        aux_v[0] = r_node.Id();
        aux_v[1] = 2.0*r_node.Id();
        r_node.FastGetSolutionStepValue(VELOCITY) = aux_v;
    }

    // Calculate the RHS and LHS
    // Note that in this case it must have zero contribution
    Vector RHS;
    Matrix LHS;
    p_test_condition->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

    // Check results
    KRATOS_EXPECT_VECTOR_NEAR(RHS, ZeroVector(8), 1.0e-12)
    KRATOS_EXPECT_MATRIX_NEAR(LHS, ZeroMatrix(8,8), 1.0e-12)
}

KRATOS_TEST_CASE_IN_SUITE(NavierStokesP2P1ContinuousWallCondition2D3NOutletInflow, FluidDynamicsApplicationFastSuite)
{
    // Create the test model part
    Model model;
    std::size_t buffer_size = 2;
    auto& r_model_part = model.CreateModelPart("TestModelPart",buffer_size);

    // Create the testing condition
    auto p_test_condition = CreateTestingNavierStokesWallCondition2D3N("NavierStokesP2P1ContinuousWallCondition2D3N", r_model_part);

    // Set the testing nodal values
    array_1d<double,3> aux_v = ZeroVector(3);
    for (auto& r_node: r_model_part.Nodes()) {
        aux_v[0] = r_node.Id();
        aux_v[1] = 2.0*r_node.Id();
        r_node.FastGetSolutionStepValue(VELOCITY) = aux_v;
    }

    // Activate the outlet inflow contribution and set required values
    p_test_condition->Set(OUTLET, true);
    r_model_part.GetProcessInfo().SetValue(CHARACTERISTIC_VELOCITY, 1.0);
    r_model_part.GetProcessInfo().SetValue(OUTLET_INFLOW_CONTRIBUTION_SWITCH, true);

    // Calculate the RHS and LHS
    // Note that in this case the LHS must have zero contribution
    Vector RHS;
    Matrix LHS;
    p_test_condition->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

    // Check results
    std::vector<double> rhs_out = {-14456.1331739, 0, -382.691792867, 0, -14336.9944223, 0, 0, 0};
    KRATOS_EXPECT_VECTOR_NEAR(RHS, rhs_out, 1.0e-8)
    KRATOS_EXPECT_MATRIX_NEAR(LHS, ZeroMatrix(8,8), 1.0e-12)
}

}  // namespace Kratos::Testing