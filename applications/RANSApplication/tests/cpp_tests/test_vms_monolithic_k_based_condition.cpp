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
#include <iomanip>

// External includes

// Project includes
#include "containers/global_pointers_vector.h"
#include "containers/model.h"
#include "testing/testing.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/fluid_test_utilities.h"
#include "custom_utilities/rans_variable_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
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
        rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
    };

    const auto& set_element_properties = [](Properties& rProperties) {
        rProperties.SetValue(DENSITY, 2.0);
        rProperties.SetValue(DYNAMIC_VISCOSITY, 2e-2);
        rProperties.SetValue(CONSTITUTIVE_LAW, KratosComponents<ConstitutiveLaw>::Get("RansKEpsilonNewtonian2DLaw").Clone());
    };

    const auto& set_condition_properties = [](Properties& rProperties) {
        rProperties.SetValue(WALL_SMOOTHNESS_BETA, 4.2);
        rProperties.SetValue(RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT, 12.0);
    };

    auto& r_model_part = FluidTestUtilities::CreateTestModelPart(
        rModel, "test", "Element2D3N", "RansVMSMonolithicKBasedWall2D2N", set_element_properties, set_condition_properties, add_variables_function,
        [](ModelPart::NodeType& rNode) {
            rNode.AddDof(VELOCITY_X).SetEquationId(rNode.Id() * 4);
            rNode.AddDof(VELOCITY_Y).SetEquationId(rNode.Id() * 4 + 1);
            rNode.AddDof(VELOCITY_Z).SetEquationId(rNode.Id() * 4 + 2);
            rNode.AddDof(PRESSURE).SetEquationId(rNode.Id() * 4 + 3);
        },
        1);

    // set nodal historical variables
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, VELOCITY, -10.0, 10.0);
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, PRESSURE, 10.0, 100.0);
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, MESH_VELOCITY, 1e-3, 1e-1);
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, TURBULENT_KINETIC_ENERGY, 10.0, 40.0);
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, ACCELERATION, 1.0, 1000.0);

    FluidTestUtilities::RandomFillNonHistoricalVariable(r_model_part.Conditions(), NORMAL, 2, -2.0, -1.0);

    auto& r_condition = r_model_part.Conditions().front();
    auto& r_element = r_model_part.Elements().front();
    r_condition.SetValue(NEIGHBOUR_ELEMENTS, GlobalPointersVector<Element>{&r_element});

    // set process info variables
    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(TURBULENCE_RANS_C_MU, 0.09);
    r_process_info.SetValue(VON_KARMAN, 3.1);

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansVMSMonolithicKBasedWall2D2N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansVMSMonolithicKBasedWall2D2NSetUp(model);

    // Test:
    FluidTestUtilities::RunEntityEquationIdVectorTest(r_model_part.Conditions(), r_model_part.GetProcessInfo(), {&VELOCITY_X, &VELOCITY_Y, &PRESSURE});
}

KRATOS_TEST_CASE_IN_SUITE(RansVMSMonolithicKBasedWall2D2N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansVMSMonolithicKBasedWall2D2NSetUp(model);

    // Test:
    FluidTestUtilities::RunEntityGetDofListTest(r_model_part.Conditions(), r_model_part.GetProcessInfo(), {&VELOCITY_X, &VELOCITY_Y, &PRESSURE});
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
    r_condition.Initialize(r_process_info);
    r_condition.CalculateLocalSystem(LHS, RHS, r_process_info);

    // setting reference values
    ref_RHS = ZeroVector(6);
    ref_LHS = ZeroMatrix(6, 6);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);

    // checking for wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 1);
    r_condition.Initialize(r_process_info);
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
    r_condition.Initialize(r_process_info);
    r_condition.CalculateLeftHandSide(LHS, r_process_info);

    // setting reference values
    ref_LHS = ZeroMatrix(6, 6);

    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);

    // checking for wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 1);
    r_condition.Initialize(r_process_info);
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
    r_condition.Initialize(r_process_info);
    r_condition.CalculateRightHandSide(RHS, r_process_info);

    // setting reference values
    ref_RHS = ZeroVector(6);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);

    // checking for wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 1);
    r_condition.Initialize(r_process_info);
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
    r_condition.Initialize(r_process_info);
    r_condition.CalculateDampingMatrix(LHS, r_process_info);

    // setting reference values
    ref_LHS = ZeroMatrix(6, 6);

    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);

    // checking for wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 1);
    r_condition.Initialize(r_process_info);
    r_condition.CalculateDampingMatrix(LHS, r_process_info);

    // setting reference values
    ref_LHS(0, 0) =  1.2878652580682977e+00;
    ref_LHS(0, 3) =  6.4718969471920806e-01;
    ref_LHS(1, 1) =  1.2878652580682977e+00;
    ref_LHS(1, 4) =  6.4718969471920806e-01;
    ref_LHS(3, 0) =  6.4718969471920806e-01;
    ref_LHS(3, 3) =  1.3008935208085346e+00;
    ref_LHS(4, 1) =  6.4718969471920806e-01;
    ref_LHS(4, 4) =  1.3008935208085346e+00;

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

    // set wall distance
    r_condition.SetValue(DISTANCE, RansCalculationUtilities::CalculateWallHeight(r_condition, r_condition.GetValue(NORMAL)));

    // checking for no-wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 0);
    r_condition.Initialize(r_process_info);
    r_condition.CalculateLocalVelocityContribution(LHS, RHS, r_process_info);

    // setting reference values
    ref_RHS = ZeroVector(6);
    ref_LHS = ZeroMatrix(6, 6);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);

    // checking for wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 1);
    r_condition.Initialize(r_process_info);
    r_condition.CalculateLocalVelocityContribution(LHS, RHS, r_process_info);

    // setting reference values
    ref_RHS[0] =  1.0664335127604463e-01;
    ref_RHS[1] =  8.2557769254555886e+00;
    ref_RHS[3] =  3.3530730354550062e+00;
    ref_RHS[4] =  9.2724325134262724e+00;

    ref_LHS(0, 0) =  1.2878652580682977e+00;
    ref_LHS(0, 3) =  6.4718969471920806e-01;
    ref_LHS(1, 1) =  1.2878652580682977e+00;
    ref_LHS(1, 4) =  6.4718969471920806e-01;
    ref_LHS(3, 0) =  6.4718969471920806e-01;
    ref_LHS(3, 3) =  1.3008935208085346e+00;
    ref_LHS(4, 1) =  6.4718969471920806e-01;
    ref_LHS(4, 4) =  1.3008935208085346e+00;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);

    const Vector& gauss_rans_y_plus = r_condition.GetValue(GAUSS_RANS_Y_PLUS);

    Vector ref_gauss_rans_y_plus(2);

    ref_gauss_rans_y_plus[0] =  2.8522374483513367e+01;
    ref_gauss_rans_y_plus[1] =  3.7167376010071798e+01;

    KRATOS_CHECK_VECTOR_NEAR(gauss_rans_y_plus, ref_gauss_rans_y_plus, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
