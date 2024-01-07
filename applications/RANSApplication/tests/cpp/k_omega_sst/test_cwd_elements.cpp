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
#include "custom_utilities/fluid_test_utilities.h"
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
    FluidTestUtilities::RunEntityEquationIdVectorTest(r_model_part.Elements(), r_model_part.GetProcessInfo(), {&TURBULENT_KINETIC_ENERGY});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3NSetUp(model);

    // Test:
    FluidTestUtilities::RunEntityGetDofListTest(r_model_part.Elements(), r_model_part.GetProcessInfo(), {&TURBULENT_KINETIC_ENERGY});
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
    ref_RHS[0] =  3.3339662429003218e+00;
    ref_RHS[1] =  5.5303468180366711e+00;
    ref_RHS[2] =  5.6949045987938263e+00;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
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
    ref_RHS[0] =  3.3339662429003218e+00;
    ref_RHS[1] =  5.5303468180366711e+00;
    ref_RHS[2] =  5.6949045987938263e+00;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
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
    ref_LHS(0, 0) =  7.6926959655028395e+02;
    ref_LHS(0, 1) =  -6.0872470624058030e+02;
    ref_LHS(0, 2) =  8.3983173248588443e+01;
    ref_LHS(1, 0) =  -6.1472090087943457e+02;
    ref_LHS(1, 1) =  1.4026652752830455e+03;
    ref_LHS(1, 2) =  -5.9017200673821287e+02;
    ref_LHS(2, 0) =  7.2766421748272307e+01;
    ref_LHS(2, 1) =  -5.9727836840027533e+02;
    ref_LHS(2, 2) =  7.9656488500325793e+02;

    ref_RHS[0] =  3.9844756180300928e+03;
    ref_RHS[1] =  6.1953394588344891e+03;
    ref_RHS[2] =  -5.1347221688072110e+04;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
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
    ref_M(0, 0) =  2.4197000392273543e-01;
    ref_M(0, 1) =  3.0895205202373426e-02;
    ref_M(0, 2) =  3.4262732263680040e-02;
    ref_M(1, 0) =  4.1585032797009272e-02;
    ref_M(1, 1) =  2.4900513297086099e-01;
    ref_M(1, 2) =  4.0553557299337337e-02;
    ref_M(2, 0) =  4.9566808896052195e-02;
    ref_M(2, 1) =  5.3009110300893988e-02;
    ref_M(2, 2) =  2.5828448595606734e-01;

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
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
    ref_D(0, 0) =  7.6926959655028395e+02;
    ref_D(0, 1) =  -6.0872470624058030e+02;
    ref_D(0, 2) =  8.3983173248588443e+01;
    ref_D(1, 0) =  -6.1472090087943457e+02;
    ref_D(1, 1) =  1.4026652752830455e+03;
    ref_D(1, 2) =  -5.9017200673821287e+02;
    ref_D(2, 0) =  7.2766421748272307e+01;
    ref_D(2, 1) =  -5.9727836840027533e+02;
    ref_D(2, 2) =  7.9656488500325793e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3NSetUp(model);

    // Test:
    FluidTestUtilities::RunEntityEquationIdVectorTest(r_model_part.Elements(), r_model_part.GetProcessInfo(), {&TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3NSetUp(model);

    // Test:
    FluidTestUtilities::RunEntityGetDofListTest(r_model_part.Elements(), r_model_part.GetProcessInfo(), {&TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE});
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
    ref_RHS[0] =  -2.1587022794893951e+03;
    ref_RHS[1] =  -2.2874097051165322e+03;
    ref_RHS[2] =  -2.4426877707931167e+03;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
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
    ref_RHS[0] =  -2.1587022794893951e+03;
    ref_RHS[1] =  -2.2874097051165322e+03;
    ref_RHS[2] =  -2.4426877707931167e+03;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
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
    ref_LHS(0, 0) =  9.1346382617329527e+02;
    ref_LHS(0, 1) =  -6.5879829404579345e+02;
    ref_LHS(0, 2) =  1.1149262641313474e+02;
    ref_LHS(1, 0) =  -6.6478033272885523e+02;
    ref_LHS(1, 1) =  1.5647649234727155e+03;
    ref_LHS(1, 2) =  -6.4338357045053033e+02;
    ref_LHS(2, 0) =  1.0042804207222419e+02;
    ref_LHS(2, 1) =  -6.5029943423125962e+02;
    ref_LHS(2, 2) =  9.3513492794875947e+02;

    ref_RHS[0] =  -4.0422844584480737e+05;
    ref_RHS[1] =  4.1067389045164717e+05;
    ref_RHS[2] =  -3.7286942267926800e+05;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
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
    ref_M(0, 0) =  2.4424017835675169e-01;
    ref_M(0, 1) =  3.2670874847167249e-02;
    ref_M(0, 2) =  3.6277727132438965e-02;
    ref_M(1, 0) =  4.1629895492917776e-02;
    ref_M(1, 1) =  2.4929391600980519e-01;
    ref_M(1, 2) =  4.0939163141299029e-02;
    ref_M(2, 0) =  4.7341722770322454e-02;
    ref_M(2, 1) =  5.1064053997629880e-02;
    ref_M(2, 2) =  2.5598453014493261e-01;

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
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
    ref_D(0, 0) =  9.1346382617329527e+02;
    ref_D(0, 1) =  -6.5879829404579345e+02;
    ref_D(0, 2) =  1.1149262641313474e+02;
    ref_D(1, 0) =  -6.6478033272885523e+02;
    ref_D(1, 1) =  1.5647649234727155e+03;
    ref_D(1, 2) =  -6.4338357045053033e+02;
    ref_D(2, 0) =  1.0042804207222419e+02;
    ref_D(2, 1) =  -6.5029943423125962e+02;
    ref_D(2, 2) =  9.3513492794875947e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
