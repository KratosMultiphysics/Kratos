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
ModelPart& RansKOmegaKCWD2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KOmegaTestUtilities::RansKOmegaK2D3NSetUp(
        rModel, "RansKOmegaKCWD2D3N");

    StabilizationMethodTestUtilities::InitializeCrossWindStabilizationConstants(
        r_model_part.GetProcessInfo());

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

ModelPart& RansKOmegaOmegaCWD2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KOmegaTestUtilities::RansKOmegaOmega2D3NSetUp(
        rModel, "RansKOmegaOmegaCWD2D3N");

    StabilizationMethodTestUtilities::InitializeCrossWindStabilizationConstants(
        r_model_part.GetProcessInfo());

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaKCWD2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaKCWD2D3NSetUp(model);

    // Test:
    FluidTestUtilities::RunEntityEquationIdVectorTest(r_model_part.Elements(), r_model_part.GetProcessInfo(), {&TURBULENT_KINETIC_ENERGY});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaKCWD2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaKCWD2D3NSetUp(model);

    // Test:
    FluidTestUtilities::RunEntityGetDofListTest(r_model_part.Elements(), r_model_part.GetProcessInfo(), {&TURBULENT_KINETIC_ENERGY});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaKCWD2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaKCWD2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  3.3339662209298471e+00;
    ref_RHS[1] =  5.5303467118285479e+00;
    ref_RHS[2] =  5.6949045237285443e+00;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaKCWD2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaKCWD2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  3.3339662209298471e+00;
    ref_RHS[1] =  5.5303467118285479e+00;
    ref_RHS[2] =  5.6949045237285443e+00;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaKCWD2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaKCWD2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS(3, 0.0), ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) =  7.6927353725327703e+02;
    ref_LHS(0, 1) =  -6.0872864713943716e+02;
    ref_LHS(0, 2) =  8.3983172584421411e+01;
    ref_LHS(1, 0) =  -6.1472484177829142e+02;
    ref_LHS(1, 1) =  1.4026731548263288e+03;
    ref_LHS(1, 2) =  -5.9017594830839494e+02;
    ref_LHS(2, 0) =  7.2766421084105261e+01;
    ref_LHS(2, 1) =  -5.9728230997045739e+02;
    ref_LHS(2, 2) =  7.9656882444117241e+02;

    ref_RHS[0] =  3.9845741862450031e+03;
    ref_RHS[1] =  6.1954440579272632e+03;
    ref_RHS[2] =  -5.1347424438975184e+04;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaKCWD2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaKCWD2D3NSetUp(model);

    // Test:
    Matrix M, ref_M(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(
        M, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_M(0, 0) =  2.4197000358750095e-01;
    ref_M(0, 1) =  3.0895204628414958e-02;
    ref_M(0, 2) =  3.4262731836817599e-02;
    ref_M(1, 0) =  4.1585031814360865e-02;
    ref_M(1, 1) =  2.4900512944120617e-01;
    ref_M(1, 2) =  4.0553556100190888e-02;
    ref_M(2, 0) =  4.9566808151758013e-02;
    ref_M(2, 1) =  5.3009108625078398e-02;
    ref_M(2, 2) =  2.5828448413427296e-01;

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaKCWD2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaKCWD2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) =  7.6927353725327703e+02;
    ref_D(0, 1) =  -6.0872864713943716e+02;
    ref_D(0, 2) =  8.3983172584421411e+01;
    ref_D(1, 0) =  -6.1472484177829142e+02;
    ref_D(1, 1) =  1.4026731548263288e+03;
    ref_D(1, 2) =  -5.9017594830839494e+02;
    ref_D(2, 0) =  7.2766421084105261e+01;
    ref_D(2, 1) =  -5.9728230997045739e+02;
    ref_D(2, 2) =  7.9656882444117241e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaCWD2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaCWD2D3NSetUp(model);

    // Test:
    FluidTestUtilities::RunEntityEquationIdVectorTest(r_model_part.Elements(), r_model_part.GetProcessInfo(), {&TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaCWD2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaCWD2D3NSetUp(model);

    // Test:
    FluidTestUtilities::RunEntityGetDofListTest(r_model_part.Elements(), r_model_part.GetProcessInfo(), {&TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaCWD2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaCWD2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  1.1533905140853946e+02;
    ref_RHS[1] =  1.2121288782511556e+02;
    ref_RHS[2] =  1.2844272535636901e+02;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaCWD2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaCWD2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  1.1533905140853946e+02;
    ref_RHS[1] =  1.2121288782511556e+02;
    ref_RHS[2] =  1.2844272535636901e+02;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaCWD2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaCWD2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS(3, 0.0), ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) =  9.2212938564768797e+02;
    ref_LHS(0, 1) =  -6.6146914654211207e+02;
    ref_LHS(0, 2) =  1.1329425410908260e+02;
    ref_LHS(1, 0) =  -6.6746534118096622e+02;
    ref_LHS(1, 1) =  1.6096468937074219e+03;
    ref_LHS(1, 2) =  -6.4482370288684615e+02;
    ref_LHS(2, 0) =  1.0207750260876654e+02;
    ref_LHS(2, 1) =  -6.5193006454890860e+02;
    ref_LHS(2, 2) =  9.4759383210742897e+02;

    ref_RHS[0] =  -4.0878900992645137e+05;
    ref_RHS[1] =  4.0744143687354401e+05;
    ref_RHS[2] =  -3.7878621679629531e+05;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaCWD2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaCWD2D3NSetUp(model);

    // Test:
    Matrix M, ref_M(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(
        M, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_M(0, 0) =  2.4464627024944774e-01;
    ref_M(0, 1) =  3.4506802377775231e-02;
    ref_M(0, 2) =  3.6726973383662591e-02;
    ref_M(1, 0) =  4.1631248550572944e-02;
    ref_M(1, 1) =  2.4939457549163468e-01;
    ref_M(1, 2) =  4.0940947972314019e-02;
    ref_M(2, 0) =  4.6961466139465793e-02;
    ref_M(2, 1) =  4.9244084328260308e-02;
    ref_M(2, 2) =  2.5556163939429732e-01;

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaCWD2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaCWD2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) =  9.2212938564768797e+02;
    ref_D(0, 1) =  -6.6146914654211207e+02;
    ref_D(0, 2) =  1.1329425410908260e+02;
    ref_D(1, 0) =  -6.6746534118096622e+02;
    ref_D(1, 1) =  1.6096468937074219e+03;
    ref_D(1, 2) =  -6.4482370288684615e+02;
    ref_D(2, 0) =  1.0207750260876654e+02;
    ref_D(2, 1) =  -6.5193006454890860e+02;
    ref_D(2, 2) =  9.4759383210742897e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
