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
ModelPart& RansKOmegaSSTKCWD2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KOmegaSSTTestUtilities::RansKOmegaSSTK2D3NSetUp(
        rModel, "RansKOmegaSSTKCWD2D3N");

    StabilizationMethodTestUtilities::InitializeCrossWindStabilizationConstants(
        r_model_part.GetProcessInfo());

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

ModelPart& RansKOmegaSSTOmegaCWD2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KOmegaSSTTestUtilities::RansKOmegaSSTOmega2D3NSetUp(
        rModel, "RansKOmegaSSTOmegaCWD2D3N");

    StabilizationMethodTestUtilities::InitializeCrossWindStabilizationConstants(
        r_model_part.GetProcessInfo());

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestEquationIdVector<ModelPart::ElementsContainerType>(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestGetDofList<ModelPart::ElementsContainerType>(
        r_model_part, TURBULENT_KINETIC_ENERGY);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = 2.4695761769301994e+00;
    ref_RHS[1] = 1.5697917787021267e+00;
    ref_RHS[2] = 1.6364014138735621e+00;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_EXPECT_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = 2.4695761769301994e+00;
    ref_RHS[1] = 1.5697917787021267e+00;
    ref_RHS[2] = 1.6364014138735621e+00;

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS(3, 0.0), ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) = 1.4194205606480323e+03;
    ref_LHS(0, 1) = -9.6750384967553202e+02;
    ref_LHS(0, 2) = 2.5032611229520009e+02;
    ref_LHS(1, 0) = -9.6646626226294461e+02;
    ref_LHS(1, 1) = 2.6162092923767191e+03;
    ref_LHS(1, 2) = -8.7453288473790849e+02;
    ref_LHS(2, 0) = 2.4952784112192893e+02;
    ref_LHS(2, 1) = -8.7521188629191010e+02;
    ref_LHS(2, 2) = 1.8769558951312886e+03;

    ref_RHS[0] = -6.8684326732130939e+04;
    ref_RHS[1] = 9.1585664588081621e+04;
    ref_RHS[2] = -1.4272006927454341e+05;

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_EXPECT_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3NSetUp(model);

    // Test:
    Matrix M, ref_M(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(
        M, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_M(0, 0) = 2.5029361539588230e-01;
    ref_M(0, 1) = 4.1823659192556498e-02;
    ref_M(0, 2) = 4.1665035727527236e-02;
    ref_M(1, 0) = 4.1127052768545523e-02;
    ref_M(1, 1) = 2.4972648471388814e-01;
    ref_M(1, 2) = 4.1479037216497421e-02;
    ref_M(2, 0) = 4.1911314394537547e-02;
    ref_M(2, 1) = 4.1782743351944499e-02;
    ref_M(2, 2) = 2.5018882921263130e-01;

    KRATOS_EXPECT_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) = 1.4194205606480323e+03;
    ref_D(0, 1) = -9.6750384967553202e+02;
    ref_D(0, 2) = 2.5032611229520009e+02;
    ref_D(1, 0) = -9.6646626226294461e+02;
    ref_D(1, 1) = 2.6162092923767191e+03;
    ref_D(1, 2) = -8.7453288473790849e+02;
    ref_D(2, 0) = 2.4952784112192893e+02;
    ref_D(2, 1) = -8.7521188629191010e+02;
    ref_D(2, 2) = 1.8769558951312886e+03;

    KRATOS_EXPECT_MATRIX_NEAR(D, ref_D, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestEquationIdVector<ModelPart::ElementsContainerType>(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestGetDofList<ModelPart::ElementsContainerType>(
        r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = -4.4653425825643089e+03;
    ref_RHS[1] = -4.4594309354935785e+03;
    ref_RHS[2] = -4.5000424077124380e+03;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_EXPECT_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = -4.4653425825643089e+03;
    ref_RHS[1] = -4.4594309354935785e+03;
    ref_RHS[2] = -4.5000424077124380e+03;

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS(3, 0.0), ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) = 1.1018055571251605e+03;
    ref_LHS(0, 1) = -6.4850066305868836e+02;
    ref_LHS(0, 2) = 1.1062713457913765e+02;
    ref_LHS(1, 0) = -6.4746307564610106e+02;
    ref_LHS(1, 1) = 1.8220246384015632e+03;
    ref_LHS(1, 2) = -6.8661079910453986e+02;
    ref_LHS(2, 0) = 1.0982886340586641e+02;
    ref_LHS(2, 1) = -6.8728980065854125e+02;
    ref_LHS(2, 2) = 9.0674940624315116e+02;

    ref_RHS[0] = -5.1385198942209681e+05;
    ref_RHS[1] = -4.6800235702979070e+05;
    ref_RHS[2] = 2.0254898170882917e+05;

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_EXPECT_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3NSetUp(model);

    // Test:
    Matrix M, ref_M(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(
        M, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_M(0, 0) = 2.5004265721772956e-01;
    ref_M(0, 1) = 4.1649022583236869e-02;
    ref_M(0, 2) = 4.0923771964438999e-02;
    ref_M(1, 0) = 4.1197290445312833e-02;
    ref_M(1, 1) = 2.4967444841500339e-01;
    ref_M(1, 2) = 4.1303363776237487e-02;
    ref_M(2, 0) = 4.2091296598204285e-02;
    ref_M(2, 1) = 4.2008202013356553e-02;
    ref_M(2, 2) = 2.5110068355697002e-01;

    KRATOS_EXPECT_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) = 1.1018055571251605e+03;
    ref_D(0, 1) = -6.4850066305868836e+02;
    ref_D(0, 2) = 1.1062713457913765e+02;
    ref_D(1, 0) = -6.4746307564610106e+02;
    ref_D(1, 1) = 1.8220246384015632e+03;
    ref_D(1, 2) = -6.8661079910453986e+02;
    ref_D(2, 0) = 1.0982886340586641e+02;
    ref_D(2, 1) = -6.8728980065854125e+02;
    ref_D(2, 2) = 9.0674940624315116e+02;

    KRATOS_EXPECT_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
