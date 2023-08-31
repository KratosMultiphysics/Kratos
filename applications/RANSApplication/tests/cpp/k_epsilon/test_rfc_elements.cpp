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
#include "../stabilization_method_test_utilities.h"
#include "custom_utilities/test_utilities.h"
#include "rans_application_variables.h"
#include "test_utilities.h"

namespace Kratos
{
namespace Testing
{
namespace
{
ModelPart& RansKEpsilonKRFC2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KEpsilonTestUtilities::RansKEpsilonK2D3NSetUp(
        rModel, "RansKEpsilonKRFC2D3N");

    StabilizationMethodTestUtilities::InitializeResidualBasedFluxCorrectedConstants(
        r_model_part.GetProcessInfo());

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

ModelPart& RansKEpsilonEpsilonRFC2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KEpsilonTestUtilities::RansKEpsilonEpsilon2D3NSetUp(
        rModel, "RansKEpsilonEpsilonRFC2D3N");

    StabilizationMethodTestUtilities::InitializeResidualBasedFluxCorrectedConstants(
        r_model_part.GetProcessInfo());

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKRFC2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKRFC2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestEquationIdVector<ModelPart::ElementsContainerType>(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKRFC2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKRFC2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestGetDofList<ModelPart::ElementsContainerType>(
        r_model_part, TURBULENT_KINETIC_ENERGY);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKRFC2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKRFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = 1.8924955742800251e+00;
    ref_RHS[1] = 1.8878284694563217e+00;
    ref_RHS[2] = 1.8954495378356369e+00;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_EXPECT_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKRFC2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKRFC2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = 1.8924955742800251e+00;
    ref_RHS[1] = 1.8878284694563217e+00;
    ref_RHS[2] = 1.8954495378356369e+00;

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKRFC2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKRFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS(3, 0.0), ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) = 1.2194714293253678e+03;
    ref_LHS(0, 1) = -2.1289725996884789e+02;
    ref_LHS(0, 2) = -2.1346518602756640e+02;
    ref_LHS(1, 0) = -2.1224529156687947e+02;
    ref_LHS(1, 1) = 1.2163061499430414e+03;
    ref_LHS(1, 2) = -2.1290777011619701e+02;
    ref_LHS(2, 0) = -2.1387783819021848e+02;
    ref_LHS(2, 1) = -2.1397239068081763e+02;
    ref_LHS(2, 2) = 1.2221971622719593e+03;

    ref_RHS[0] = -2.9907675429503579e+04;
    ref_RHS[1] = 1.9583263472729057e+04;
    ref_RHS[2] = -7.9850125960094752e+04;

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_EXPECT_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKRFC2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKRFC2D3NSetUp(model);

    // Test:
    Matrix M, ref_M(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(
        M, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_M(0, 0) = 2.2225547131388934e-01;
    ref_M(0, 1) = 5.5588804647222656e-02;
    ref_M(0, 2) = 5.5588804647222656e-02;
    ref_M(1, 0) = 5.5314710278559671e-02;
    ref_M(1, 1) = 2.2198137694522629e-01;
    ref_M(1, 2) = 5.5314710278559651e-02;
    ref_M(2, 0) = 5.5762287953557865e-02;
    ref_M(2, 1) = 5.5762287953557844e-02;
    ref_M(2, 2) = 2.2242895462022449e-01;

    KRATOS_EXPECT_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKRFC2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKRFC2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) = 1.2194714293253678e+03;
    ref_D(0, 1) = -2.1289725996884789e+02;
    ref_D(0, 2) = -2.1346518602756640e+02;
    ref_D(1, 0) = -2.1224529156687947e+02;
    ref_D(1, 1) = 1.2163061499430414e+03;
    ref_D(1, 2) = -2.1290777011619701e+02;
    ref_D(2, 0) = -2.1387783819021848e+02;
    ref_D(2, 1) = -2.1397239068081763e+02;
    ref_D(2, 2) = 1.2221971622719593e+03;

    KRATOS_EXPECT_MATRIX_NEAR(D, ref_D, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonRFC2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonRFC2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestEquationIdVector<ModelPart::ElementsContainerType>(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonRFC2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonRFC2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestGetDofList<ModelPart::ElementsContainerType>(
        r_model_part, TURBULENT_ENERGY_DISSIPATION_RATE);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonRFC2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonRFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = 7.2797422181906841e+03;
    ref_RHS[1] = 7.2715936386889534e+03;
    ref_RHS[2] = 7.2848997211368023e+03;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_EXPECT_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonRFC2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonRFC2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = 7.2797422181906841e+03;
    ref_RHS[1] = 7.2715936386889534e+03;
    ref_RHS[2] = 7.2848997211368023e+03;

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonRFC2D3N_CalculateLocalVelocityContribution,
                          KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonRFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS(3, 0.0), ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) = 2.6659507621547718e+03;
    ref_LHS(0, 1) = -4.5903337346355045e+02;
    ref_LHS(0, 2) = -4.5956119036827181e+02;
    ref_LHS(1, 0) = -4.5838140506158209e+02;
    ref_LHS(1, 1) = 2.6628125412181089e+03;
    ref_LHS(1, 2) = -4.5903084100612341e+02;
    ref_LHS(2, 0) = -4.5997384253092400e+02;
    ref_LHS(2, 1) = -4.6009546157074396e+02;
    ref_LHS(2, 2) = 2.6686634576255929e+03;

    ref_RHS[0] = 3.2092070925074752e+05;
    ref_RHS[1] = -1.3816269382337187e+06;
    ref_RHS[2] = -2.0722915797429322e+06;

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_EXPECT_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonRFC2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonRFC2D3NSetUp(model);

    // Test:
    Matrix M, ref_M(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(
        M, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_M(0, 0) = 2.2223738419925163e-01;
    ref_M(0, 1) = 5.5570717532584954e-02;
    ref_M(0, 2) = 5.5570717532584954e-02;
    ref_M(1, 0) = 5.5446328338043610e-02;
    ref_M(1, 1) = 2.2211299500471024e-01;
    ref_M(1, 2) = 5.5446328338043589e-02;
    ref_M(2, 0) = 5.5649447529288888e-02;
    ref_M(2, 1) = 5.5649447529288867e-02;
    ref_M(2, 2) = 2.2231611419595554e-01;

    KRATOS_EXPECT_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonRFC2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonRFC2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) = 2.6659507621547718e+03;
    ref_D(0, 1) = -4.5903337346355045e+02;
    ref_D(0, 2) = -4.5956119036827181e+02;
    ref_D(1, 0) = -4.5838140506158209e+02;
    ref_D(1, 1) = 2.6628125412181089e+03;
    ref_D(1, 2) = -4.5903084100612341e+02;
    ref_D(2, 0) = -4.5997384253092400e+02;
    ref_D(2, 1) = -4.6009546157074396e+02;
    ref_D(2, 2) = 2.6686634576255929e+03;

    KRATOS_EXPECT_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
