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
ModelPart& RansKEpsilonKCWD2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KEpsilonTestUtilities::RansKEpsilonK2D3NSetUp(
        rModel, "RansKEpsilonKCWD2D3N");

    StabilizationMethodTestUtilities::InitializeCrossWindStabilizationConstants(
        r_model_part.GetProcessInfo());

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

ModelPart& RansKEpsilonEpsilonCWD2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KEpsilonTestUtilities::RansKEpsilonEpsilon2D3NSetUp(
        rModel, "RansKEpsilonEpsilonCWD2D3N");

    StabilizationMethodTestUtilities::InitializeCrossWindStabilizationConstants(
        r_model_part.GetProcessInfo());

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKCWD2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKCWD2D3NSetUp(model);

    // Test:
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityEquationIdVectorTest(r_model_part, {&TURBULENT_KINETIC_ENERGY});
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKCWD2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKCWD2D3NSetUp(model);

    // Test:
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityGetDofListTest(r_model_part, {&TURBULENT_KINETIC_ENERGY});
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKCWD2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKCWD2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  1.5392770768777893e+02;
    ref_RHS[1] =  2.0085159669087230e+02;
    ref_RHS[2] =  2.9241615576967308e+02;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKCWD2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKCWD2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  1.5392770768777893e+02;
    ref_RHS[1] =  2.0085159669087230e+02;
    ref_RHS[2] =  2.9241615576967308e+02;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKCWD2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKCWD2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS(3, 0.0), ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) =  1.8303674024026321e+01;
    ref_LHS(0, 1) =  -1.4831345683544091e+01;
    ref_LHS(0, 2) =  1.7157693836332222e+00;
    ref_LHS(1, 0) =  -1.3362222776831324e+01;
    ref_LHS(1, 1) =  2.8680968667196144e+01;
    ref_LHS(1, 2) =  -1.3193592114789842e+01;
    ref_LHS(2, 0) =  1.2825521977500876e+00;
    ref_LHS(2, 1) =  -1.4611254949402884e+01;
    ref_LHS(2, 2) =  1.6423737968131007e+01;

    ref_RHS[0] =  1.3748588321166096e+02;
    ref_RHS[1] =  2.5052632272961273e+02;
    ref_RHS[2] =  -9.5186959391440428e+02;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKCWD2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKCWD2D3NSetUp(model);

    // Test:
    Matrix M, ref_M(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(
        M, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_M(0, 0) =  2.1097776352072176e-01;
    ref_M(0, 1) =  1.2982186516959494e-02;
    ref_M(0, 2) =  1.2058657173943594e-02;
    ref_M(1, 0) =  -7.3186368869048975e-03;
    ref_M(1, 1) =  1.6584887405294227e-01;
    ref_M(1, 2) =  -2.4073245010721067e-03;
    ref_M(2, 0) =  1.6475739438823456e-02;
    ref_M(2, 1) =  7.8870451934899520e-03;
    ref_M(2, 2) =  1.7701446350611955e-01;

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonKCWD2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonKCWD2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) =  1.8303674024026321e+01;
    ref_D(0, 1) =  -1.4831345683544091e+01;
    ref_D(0, 2) =  1.7157693836332222e+00;
    ref_D(1, 0) =  -1.3362222776831324e+01;
    ref_D(1, 1) =  2.8680968667196144e+01;
    ref_D(1, 2) =  -1.3193592114789842e+01;
    ref_D(2, 0) =  1.2825521977500876e+00;
    ref_D(2, 1) =  -1.4611254949402884e+01;
    ref_D(2, 2) =  1.6423737968131007e+01;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonCWD2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonCWD2D3NSetUp(model);

    // Test:
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityEquationIdVectorTest(r_model_part, {&TURBULENT_ENERGY_DISSIPATION_RATE});
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonCWD2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonCWD2D3NSetUp(model);

    // Test:
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityGetDofListTest(r_model_part, {&TURBULENT_ENERGY_DISSIPATION_RATE});
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonCWD2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonCWD2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  4.1870674572650069e+03;
    ref_RHS[1] =  3.7306016874891466e+03;
    ref_RHS[2] =  5.7060562773292304e+03;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonCWD2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonCWD2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  4.1870674572650069e+03;
    ref_RHS[1] =  3.7306016874891466e+03;
    ref_RHS[2] =  5.7060562773292304e+03;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonCWD2D3N_CalculateLocalVelocityContribution,
                          KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonCWD2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS(3, 0.0), ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) =  3.2919202184477911e+01;
    ref_LHS(0, 1) =  -1.1866604480782106e+01;
    ref_LHS(0, 2) =  6.2559402951690091e+00;
    ref_LHS(1, 0) =  -1.0397481574069340e+01;
    ref_LHS(1, 1) =  3.5625576592780980e+01;
    ref_LHS(1, 2) =  -1.1592707967032670e+01;
    ref_LHS(2, 0) =  5.8227231092858744e+00;
    ref_LHS(2, 1) =  -1.3010370801645713e+01;
    ref_LHS(2, 2) =  2.4422535365010507e+01;

    ref_RHS[0] =  -3.2102016201773553e+04;
    ref_RHS[1] =  9.9970721225592497e+03;
    ref_RHS[2] =  -1.5827717525576181e+04;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonCWD2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonCWD2D3NSetUp(model);

    // Test:
    Matrix M, ref_M(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(
        M, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_M(0, 0) =  2.5185423167416821e-01;
    ref_M(0, 1) =  3.7162365341543359e-02;
    ref_M(0, 2) =  3.4172408168318721e-02;
    ref_M(1, 0) =  2.0807409333467024e-02;
    ref_M(1, 1) =  2.0930631655676540e-01;
    ref_M(1, 2) =  1.7944000637448443e-02;
    ref_M(2, 0) =  3.9477874974911217e-02;
    ref_M(2, 1) =  3.6579646506471400e-02;
    ref_M(2, 2) =  2.3059927825111598e-01;

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonEpsilonCWD2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKEpsilonEpsilonCWD2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) =  3.2919202184477911e+01;
    ref_D(0, 1) =  -1.1866604480782106e+01;
    ref_D(0, 2) =  6.2559402951690091e+00;
    ref_D(1, 0) =  -1.0397481574069340e+01;
    ref_D(1, 1) =  3.5625576592780980e+01;
    ref_D(1, 2) =  -1.1592707967032670e+01;
    ref_D(2, 0) =  5.8227231092858744e+00;
    ref_D(2, 1) =  -1.3010370801645713e+01;
    ref_D(2, 2) =  2.4422535365010507e+01;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
