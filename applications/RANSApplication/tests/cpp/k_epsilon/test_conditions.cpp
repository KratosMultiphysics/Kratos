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
#include "containers/model.h"
#include "testing/testing.h"

// Application includes
#include "custom_utilities/test_utilities.h"
#include "rans_application_variables.h"
#include "test_utilities.h"

namespace Kratos
{
namespace Testing
{
namespace
{
ModelPart& RansKEpsilonEpsilonKBasedWall2D2NSetUp(
    Model& rModel)
{
    auto& r_model_part = KEpsilonTestUtilities::RansKEpsilonEpsilon2D2NSetUp(
        rModel, "RansKEpsilonEpsilonKBasedWall2D2N");

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

ModelPart& RansKEpsilonEpsilonUBasedWall2D2NSetUp(
    Model& rModel)
{
    auto& r_model_part = KEpsilonTestUtilities::RansKEpsilonEpsilon2D2NSetUp(
        rModel, "RansKEpsilonEpsilonUBasedWall2D2N");

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonKBasedWall2D2N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonKBasedWall2D2NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestEquationIdVector<ModelPart::ConditionsContainerType>(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonKBasedWall2D2N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonKBasedWall2D2NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestGetDofList<ModelPart::ConditionsContainerType>(
        r_model_part, TURBULENT_ENERGY_DISSIPATION_RATE);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonKBasedWall2D2N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonKBasedWall2D2NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS;
    auto& r_condition = r_model_part.Conditions().front();

    const auto& r_process_info = r_model_part.GetProcessInfo();

    // checking for no-wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 0);
    r_condition.CalculateLocalSystem(LHS, RHS, r_process_info);
    // setting reference values
    ref_RHS = ZeroVector(2);
    ref_LHS = ZeroMatrix(2, 2);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);

    // checking for wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 1);
    r_condition.CalculateLocalSystem(LHS, RHS, r_process_info);
    // setting reference values
    ref_RHS[0] = 6.5112436747195190e+01;
    ref_RHS[1] = 6.5112436747195190e+01;
    ref_LHS = ZeroMatrix(2, 2);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonKBasedWall2D2N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonKBasedWall2D2NSetUp(model);

    const auto& r_process_info = r_model_part.GetProcessInfo();

    // Test:
    Vector RHS, ref_RHS;
    auto& r_condition = r_model_part.Conditions().front();

    // checking for no-wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 0);
    r_condition.CalculateRightHandSide(RHS, r_process_info);
    // setting reference values
    ref_RHS = ZeroVector(2);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);

    // checking for wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 1);
    r_condition.CalculateRightHandSide(RHS, r_process_info);
    // setting reference values
    ref_RHS[0] = 6.5112436747195190e+01;
    ref_RHS[1] = 6.5112436747195190e+01;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonUBasedWall2D2N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonUBasedWall2D2NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestEquationIdVector<ModelPart::ConditionsContainerType>(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonUBasedWall2D2N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonUBasedWall2D2NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestGetDofList<ModelPart::ConditionsContainerType>(
        r_model_part, TURBULENT_ENERGY_DISSIPATION_RATE);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonUBasedWall2D2N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonUBasedWall2D2NSetUp(model);

    const auto& r_process_info = r_model_part.GetProcessInfo();

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS;
    auto& r_condition = r_model_part.Conditions().front();

    // checking for no-wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 0);
    r_condition.CalculateLocalSystem(LHS, RHS, r_process_info);
    // setting reference values
    ref_RHS = ZeroVector(2);
    ref_LHS = ZeroMatrix(2, 2);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);

    // checking for wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 1);
    r_condition.CalculateLocalSystem(LHS, RHS, r_process_info);
    // setting reference values
    ref_RHS[0] = 3.8869408853394147e+04;
    ref_RHS[1] = 3.8869408853394147e+04;
    ref_LHS = ZeroMatrix(2, 2);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonUBasedWall2D2N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonUBasedWall2D2NSetUp(model);

    const auto& r_process_info = r_model_part.GetProcessInfo();

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_condition = r_model_part.Conditions().front();

    // checking for no-wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 0);
    r_condition.CalculateRightHandSide(RHS, r_process_info);
    // setting reference values
    ref_RHS = ZeroVector(2);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);

    // checking for wall function
    r_condition.SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 1);
    r_condition.CalculateRightHandSide(RHS, r_process_info);
    // setting reference values
    ref_RHS[0] = 3.8869408853394147e+04;
    ref_RHS[1] = 3.8869408853394147e+04;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
