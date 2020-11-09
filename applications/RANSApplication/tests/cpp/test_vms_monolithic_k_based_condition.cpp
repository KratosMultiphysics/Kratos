//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "containers/global_pointers_vector.h"
#include "containers/model.h"
#include "testing/testing.h"

// Application includes
#include "custom_utilities/test_utilities.h"
#include "includes/cfd_variables.h"
#include "rans_application_variables.h"

namespace Kratos
{
namespace Testing
{
namespace
{
ModelPart& RansVMSMonolithicKBasedWall2D2NSetUp(Model& rModel)
{
    const auto add_variables_function = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(PRESSURE);
        rModelPart.AddNodalSolutionStepVariable(DENSITY);
        rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
        rModelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
    };

    using namespace RansApplicationTestUtilities;

    auto& r_model_part = CreateTestModelPart(
        rModel, "Element2D3N", "RansVMSMonolithicKBasedWall2D2N", add_variables_function,
        [](ModelPart::NodeType& rNode) {
            rNode.AddDof(VELOCITY_X).SetEquationId(rNode.Id() * 4);
            rNode.AddDof(VELOCITY_Y).SetEquationId(rNode.Id() * 4 + 1);
            rNode.AddDof(VELOCITY_Z).SetEquationId(rNode.Id() * 4 + 2);
            rNode.AddDof(PRESSURE).SetEquationId(rNode.Id() * 4 + 3);
        },
        1);

    // set nodal historical variables
    RandomFillNodalHistoricalVariable(r_model_part, VELOCITY, -10.0, 10.0);
    RandomFillNodalHistoricalVariable(r_model_part, DENSITY, 2.0, 2.0);
    RandomFillNodalHistoricalVariable(r_model_part, KINEMATIC_VISCOSITY, 1e-5, 1e-3);
    RandomFillNodalHistoricalVariable(r_model_part, PRESSURE, 10.0, 100.0);
    RandomFillNodalHistoricalVariable(r_model_part, MESH_VELOCITY, 1e-3, 1e-1);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_KINETIC_ENERGY, 10.0, 40.0);
    RandomFillNodalHistoricalVariable(r_model_part, ACCELERATION, 1.0, 1000.0);
    RandomFillNodalHistoricalVariable(r_model_part, EXTERNAL_PRESSURE, 1.0, 1000.0);

    RandomFillContainerVariable<ModelPart::ConditionsContainerType>(
        r_model_part, NORMAL, -2.0, -1.0);

    auto& r_condition = r_model_part.Conditions().front();
    auto& r_element = r_model_part.Elements().front();
    r_condition.SetValue(NEIGHBOUR_ELEMENTS, GlobalPointersVector<Element>{&r_element});

    // set process info variables
    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(WALL_VON_KARMAN, 3.1);
    r_process_info.SetValue(WALL_SMOOTHNESS_BETA, 4.2);
    r_process_info.SetValue(RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT, 12.0);
    r_process_info.SetValue(TURBULENCE_RANS_C_MU, 0.09);

    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansVMSMonolithicKBasedWall2D2N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansVMSMonolithicKBasedWall2D2NSetUp(model);

    // Test:
    auto eqn_ids = std::vector<IndexType>{};
    for (const auto& r_condition : r_model_part.Conditions())
    {
        r_condition.EquationIdVector(eqn_ids, r_model_part.GetProcessInfo());
        KRATOS_CHECK_EQUAL(eqn_ids.size(), r_condition.GetGeometry().PointsNumber() * 3);
        const auto& r_geometry = r_condition.GetGeometry();
        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
        {
            KRATOS_CHECK_EQUAL(eqn_ids[local_index++], r_geometry[i_node].Id() * 4);
            KRATOS_CHECK_EQUAL(eqn_ids[local_index++], r_geometry[i_node].Id() * 4 + 1);
            KRATOS_CHECK_EQUAL(eqn_ids[local_index++], r_geometry[i_node].Id() * 4 + 3);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansVMSMonolithicKBasedWall2D2N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansVMSMonolithicKBasedWall2D2NSetUp(model);

    // Test:
    auto dofs = Element::DofsVectorType{};
    for (const auto& r_condition : r_model_part.Conditions())
    {
        r_condition.GetDofList(dofs, r_model_part.GetProcessInfo());
        KRATOS_CHECK_EQUAL(dofs.size(), r_condition.GetGeometry().PointsNumber() * 3);
        const auto& r_geometry = r_condition.GetGeometry();
        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
        {
            KRATOS_CHECK_EQUAL(dofs[local_index]->GetVariable(), VELOCITY_X);
            KRATOS_CHECK_EQUAL(dofs[local_index++]->EquationId(),
                               r_geometry[i_node].Id() * 4);
            KRATOS_CHECK_EQUAL(dofs[local_index]->GetVariable(), VELOCITY_Y);
            KRATOS_CHECK_EQUAL(dofs[local_index++]->EquationId(),
                               r_geometry[i_node].Id() * 4 + 1);
            KRATOS_CHECK_EQUAL(dofs[local_index]->GetVariable(), PRESSURE);
            KRATOS_CHECK_EQUAL(dofs[local_index++]->EquationId(),
                               r_geometry[i_node].Id() * 4 + 3);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansVMSMonolithicKBasedWall2D2N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansVMSMonolithicKBasedWall2D2NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS;
    auto& r_condition = r_model_part.Conditions().front();

    const auto& r_process_info = r_model_part.GetProcessInfo();

    // checking for no-wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 0);
    r_condition.Initialize();
    r_condition.CalculateLocalSystem(LHS, RHS, r_process_info);

    // setting reference values
    ref_RHS = ZeroVector(6);
    ref_LHS = ZeroMatrix(6, 6);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);

    // checking for wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 1);
    r_condition.Initialize();
    r_condition.CalculateLocalSystem(LHS, RHS, r_process_info);

    // setting reference values
    ref_RHS = ZeroVector(6);
    ref_LHS = ZeroMatrix(6, 6);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansVMSMonolithicKBasedWall2D2N_CalculateLeftHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansVMSMonolithicKBasedWall2D2NSetUp(model);

    const auto& r_process_info = r_model_part.GetProcessInfo();

    // Test:
    Matrix LHS, ref_LHS;
    auto& r_condition = r_model_part.Conditions().front();

    // checking for no-wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 0);
    r_condition.Initialize();
    r_condition.CalculateLeftHandSide(LHS, r_process_info);

    // setting reference values
    ref_LHS = ZeroMatrix(6, 6);

    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);

    // checking for wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 1);
    r_condition.Initialize();
    r_condition.CalculateLeftHandSide(LHS, r_process_info);

    // setting reference values
    ref_LHS = ZeroMatrix(6, 6);

    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansVMSMonolithicKBasedWall2D2N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansVMSMonolithicKBasedWall2D2NSetUp(model);

    const auto& r_process_info = r_model_part.GetProcessInfo();

    // Test:
    Vector RHS, ref_RHS;
    auto& r_condition = r_model_part.Conditions().front();

    // checking for no-wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 0);
    r_condition.Initialize();
    r_condition.CalculateRightHandSide(RHS, r_process_info);

    // setting reference values
    ref_RHS = ZeroVector(6);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);

    // checking for wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 1);
    r_condition.Initialize();
    r_condition.CalculateRightHandSide(RHS, r_process_info);

    // setting reference values
    ref_RHS = ZeroVector(6);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansVMSMonolithicKBasedWall2D2N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansVMSMonolithicKBasedWall2D2NSetUp(model);

    const auto& r_process_info = r_model_part.GetProcessInfo();

    // Test:
    Matrix LHS, ref_LHS;
    auto& r_condition = r_model_part.Conditions().front();

    // checking for no-wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 0);
    r_condition.Initialize();
    r_condition.CalculateDampingMatrix(LHS, r_process_info);

    // setting reference values
    ref_LHS = ZeroMatrix(6, 6);

    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);

    // checking for wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 1);
    r_condition.Initialize();
    r_condition.CalculateDampingMatrix(LHS, r_process_info);

    // setting reference values
    ref_LHS(0, 0) = 7.5197560293399790e-01;
    ref_LHS(0, 3) = 7.5197560293399790e-01;
    ref_LHS(1, 1) = 7.5197560293399790e-01;
    ref_LHS(1, 4) = 7.5197560293399790e-01;
    ref_LHS(3, 0) = 7.5197560293399790e-01;
    ref_LHS(3, 3) = 7.5197560293399790e-01;
    ref_LHS(4, 1) = 7.5197560293399790e-01;
    ref_LHS(4, 4) = 7.5197560293399790e-01;

    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansVMSMonolithicKBasedWall2D2N_CalculateLocalVelocityContribution,
                          KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansVMSMonolithicKBasedWall2D2NSetUp(model);

    const auto& r_process_info = r_model_part.GetProcessInfo();

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS;
    auto& r_condition = r_model_part.Conditions().front();

    // checking for no-wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 0);
    r_condition.Initialize();
    r_condition.CalculateLocalVelocityContribution(LHS, RHS, r_process_info);

    // setting reference values
    ref_RHS = ZeroVector(6);
    ref_LHS = ZeroMatrix(6, 6);
    ref_RHS[0] = -1.3614337272157539e+02;
    ref_RHS[3] = -2.0529174563896703e+02;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);

    // checking for wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 1);
    r_condition.Initialize();
    r_condition.CalculateLocalVelocityContribution(LHS, RHS, r_process_info);

    // setting reference values
    ref_RHS[0] = -1.3560206485323209e+02;
    ref_RHS[1] = 1.4876293019322391e+00;
    ref_RHS[3] = -2.0475043777062373e+02;
    ref_RHS[4] = 1.4876293019322391e+00;

    ref_LHS(0, 0) = 7.5197560293399790e-01;
    ref_LHS(0, 3) = 7.5197560293399790e-01;
    ref_LHS(1, 1) = 7.5197560293399790e-01;
    ref_LHS(1, 4) = 7.5197560293399790e-01;
    ref_LHS(3, 0) = 7.5197560293399790e-01;
    ref_LHS(3, 3) = 7.5197560293399790e-01;
    ref_LHS(4, 1) = 7.5197560293399790e-01;
    ref_LHS(4, 4) = 7.5197560293399790e-01;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);

    const array_1d<double, 3>& friction_velocity = r_condition.GetValue(FRICTION_VELOCITY);
    array_1d<double, 3> ref_friction_velocity;
    ref_friction_velocity[0] = -1.9591460562647886e-01;
    ref_friction_velocity[1] = -5.3841505924981037e-01;
    ref_friction_velocity[2] = -2.7029201972010961e+00;
    KRATOS_CHECK_VECTOR_NEAR(friction_velocity, ref_friction_velocity, 1e-12);

    KRATOS_CHECK_NEAR(r_condition.GetValue(RANS_Y_PLUS), 4.12772678625194032520e+02, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
