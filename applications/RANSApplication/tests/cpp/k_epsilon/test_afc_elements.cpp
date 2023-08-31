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
ModelPart& RansKEpsilonKAFC2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KEpsilonTestUtilities::RansKEpsilonK2D3NSetUp(
        rModel, "RansKEpsilonKAFC2D3N");

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

ModelPart& RansKEpsilonEpsilonAFC2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KEpsilonTestUtilities::RansKEpsilonEpsilon2D3NSetUp(
        rModel, "RansKEpsilonEpsilonAFC2D3N");

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKAFC2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKAFC2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestEquationIdVector<ModelPart::ElementsContainerType>(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKAFC2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKAFC2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestGetDofList<ModelPart::ElementsContainerType>(
        r_model_part, TURBULENT_KINETIC_ENERGY);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKAFC2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKAFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = 1.2327968933223705e+00;
    ref_RHS[1] = 7.8878546835975738e-01;
    ref_RHS[2] = 8.1631178311485209e-01;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_EXPECT_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKAFC2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKAFC2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = 1.2327968933223705e+00;
    ref_RHS[1] = 7.8878546835975738e-01;
    ref_RHS[2] = 8.1631178311485209e-01;

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKAFC2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKAFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS, ref_RHS;
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) = 1.4131020272484841e+02;
    ref_LHS(0, 1) = 8.4468027013988021e+01;
    ref_LHS(0, 2) = 1.2522456766544926e+02;
    ref_LHS(1, 0) = 8.5505614426575448e+01;
    ref_LHS(1, 1) = 1.7170851645741465e+02;
    ref_LHS(1, 2) = 1.3124961554377217e+02;
    ref_LHS(2, 0) = 1.2442629649217810e+02;
    ref_LHS(2, 1) = 1.3057061398977064e+02;
    ref_LHS(2, 2) = 3.6990083768146155e+02;

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_EXPECT_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKAFC2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKAFC2D3NSetUp(model);

    // Test:
    Matrix M, ref_M;
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(
        M, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_M = ZeroMatrix(3, 3);
    ref_M(0, 0) = 1.66666666666666657415e-01;
    ref_M(1, 1) = 1.66666666666666657415e-01;
    ref_M(2, 2) = 1.66666666666666657415e-01;

    KRATOS_EXPECT_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKAFC2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKAFC2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) = 1.4131020272484841e+02;
    ref_D(0, 1) = 8.4468027013988021e+01;
    ref_D(0, 2) = 1.2522456766544926e+02;
    ref_D(1, 0) = 8.5505614426575448e+01;
    ref_D(1, 1) = 1.7170851645741465e+02;
    ref_D(1, 2) = 1.3124961554377217e+02;
    ref_D(2, 0) = 1.2442629649217810e+02;
    ref_D(2, 1) = 1.3057061398977064e+02;
    ref_D(2, 2) = 3.6990083768146155e+02;

    KRATOS_EXPECT_MATRIX_NEAR(D, ref_D, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonAFC2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonAFC2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestEquationIdVector<ModelPart::ElementsContainerType>(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonAFC2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonAFC2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestGetDofList<ModelPart::ElementsContainerType>(
        r_model_part, TURBULENT_ENERGY_DISSIPATION_RATE);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonAFC2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonAFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = 3.6382647412713577e+03;
    ref_RHS[1] = 2.8060674875212635e+03;
    ref_RHS[2] = 4.4737912354588825e+03;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_EXPECT_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonAFC2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonAFC2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = 3.6382647412713577e+03;
    ref_RHS[1] = 2.8060674875212635e+03;
    ref_RHS[2] = 4.4737912354588825e+03;

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonAFC2D3N_CalculateLocalVelocityContribution,
                          KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonAFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS, ref_RHS;
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) = 3.1107947789649791e+02;
    ref_LHS(0, 1) = 1.8725080934712418e+02;
    ref_LHS(0, 2) = 2.7527586704580659e+02;
    ref_LHS(1, 0) = 1.8828839675971159e+02;
    ref_LHS(1, 1) = 3.7892045735263218e+02;
    ref_LHS(1, 2) = 2.8881138802873329e+02;
    ref_LHS(2, 0) = 2.7447759587253546e+02;
    ref_LHS(2, 1) = 2.8813238647473179e+02;
    ref_LHS(2, 2) = 8.1356506361223580e+02;

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_EXPECT_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonAFC2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonAFC2D3NSetUp(model);

    // Test:
    Matrix M, ref_M;
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(
        M, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_M = ZeroMatrix(3, 3);
    ref_M(0, 0) = 1.66666666666666657415e-01;
    ref_M(1, 1) = 1.66666666666666657415e-01;
    ref_M(2, 2) = 1.66666666666666657415e-01;

    KRATOS_EXPECT_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonAFC2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonAFC2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) = 3.1107947789649791e+02;
    ref_D(0, 1) = 1.8725080934712418e+02;
    ref_D(0, 2) = 2.7527586704580659e+02;
    ref_D(1, 0) = 1.8828839675971159e+02;
    ref_D(1, 1) = 3.7892045735263218e+02;
    ref_D(1, 2) = 2.8881138802873329e+02;
    ref_D(2, 0) = 2.7447759587253546e+02;
    ref_D(2, 1) = 2.8813238647473179e+02;
    ref_D(2, 2) = 8.1356506361223580e+02;


    KRATOS_EXPECT_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
