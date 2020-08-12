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

    return r_model_part;
}

ModelPart& RansKEpsilonEpsilonRFC2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KEpsilonTestUtilities::RansKEpsilonEpsilon2D3NSetUp(
        rModel, "RansKEpsilonEpsilonRFC2D3N");

    StabilizationMethodTestUtilities::InitializeResidualBasedFluxCorrectedConstants(
        r_model_part.GetProcessInfo());

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
    ref_RHS[0] = 4.50994145217781028379e+00;
    ref_RHS[1] = 4.49055195644842175540e+00;
    ref_RHS[2] = 4.49880538668604224029e+00;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
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
    ref_RHS[0] = 4.50994145217781028379e+00;
    ref_RHS[1] = 4.49055195644842175540e+00;
    ref_RHS[2] = 4.49880538668604224029e+00;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
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
    ref_LHS(0, 0) = 9.75406602650585227821e+02;
    ref_LHS(0, 1) = -1.70080032921509769039e+02;
    ref_LHS(0, 2) = -1.70048687938027256905e+02;
    ref_LHS(1, 0) = -1.69169605947782457633e+02;
    ref_LHS(1, 1) = 9.70052072207196033560e+02;
    ref_LHS(1, 2) = -1.68335821205377698107e+02;
    ref_LHS(2, 0) = -1.69525797887331805214e+02;
    ref_LHS(2, 1) = -1.68723358128409557821e+02;
    ref_LHS(2, 2) = 9.71958393031216814961e+02;

    ref_RHS[0] = -2.98311916208666880266e+04;
    ref_RHS[1] = -5.37915542593372520059e+04;
    ref_RHS[2] = -1.30869844126473908545e+04;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
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
    ref_M(0, 0) = 2.22472573700204456770e-01;
    ref_M(0, 1) = 5.58059070335377715999e-02;
    ref_M(0, 2) = 5.58059070335377715999e-02;
    ref_M(1, 0) = 5.53271330569383990028e-02;
    ref_M(1, 1) = 2.21993799723605028662e-01;
    ref_M(1, 2) = 5.53271330569383781861e-02;
    ref_M(2, 0) = 5.55309303949713742310e-02;
    ref_M(2, 1) = 5.55309303949713534143e-02;
    ref_M(2, 2) = 2.22197597061637996951e-01;

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
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
    ref_D(0, 0) = 9.75406602650585227821e+02;
    ref_D(0, 1) = -1.70080032921509769039e+02;
    ref_D(0, 2) = -1.70048687938027256905e+02;
    ref_D(1, 0) = -1.69169605947782457633e+02;
    ref_D(1, 1) = 9.70052072207196033560e+02;
    ref_D(1, 2) = -1.68335821205377698107e+02;
    ref_D(2, 0) = -1.69525797887331805214e+02;
    ref_D(2, 1) = -1.68723358128409557821e+02;
    ref_D(2, 2) = 9.71958393031216814961e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
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
    ref_RHS[0] = 1.38886658671788154606e+04;
    ref_RHS[1] = 1.38615066269478447794e+04;
    ref_RHS[2] = 1.38730673656943035894e+04;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
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
    ref_RHS[0] = 1.38886658671788154606e+04;
    ref_RHS[1] = 1.38615066269478447794e+04;
    ref_RHS[2] = 1.38730673656943035894e+04;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
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
    ref_LHS(0, 0) = 2.08324515220861621856e+03;
    ref_LHS(0, 1) = -3.43306415149386680241e+02;
    ref_LHS(0, 2) = -3.43223474424340679434e+02;
    ref_LHS(1, 0) = -3.42395988175659340413e+02;
    ref_LHS(1, 1) = 2.07800761461857746326e+03;
    ref_LHS(1, 2) = -3.41627635890958970322e+02;
    ref_LHS(2, 0) = -3.42700584373645199321e+02;
    ref_LHS(2, 1) = -3.42015172813990830036e+02;
    ref_LHS(2, 2) = 2.07986235474655586586e+03;

    ref_RHS[0] = 2.74716727342334634159e+05;
    ref_RHS[1] = -1.40242927389094606042e+06;
    ref_RHS[2] = -8.13739162544821156189e+04;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
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
    ref_M(0, 0) = 2.22336181869094179353e-01;
    ref_M(0, 1) = 5.56695152024275080604e-02;
    ref_M(0, 2) = 5.56695152024275080604e-02;
    ref_M(1, 0) = 5.54520149439209908415e-02;
    ref_M(1, 1) = 2.22118681610587620501e-01;
    ref_M(1, 2) = 5.54520149439209700248e-02;
    ref_M(2, 0) = 5.55445971974928079384e-02;
    ref_M(2, 1) = 5.55445971974927871218e-02;
    ref_M(2, 2) = 2.22211263864159458414e-01;

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
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
    ref_D(0, 0) = 2.08324515220861621856e+03;
    ref_D(0, 1) = -3.43306415149386680241e+02;
    ref_D(0, 2) = -3.43223474424340679434e+02;
    ref_D(1, 0) = -3.42395988175659340413e+02;
    ref_D(1, 1) = 2.07800761461857746326e+03;
    ref_D(1, 2) = -3.41627635890958970322e+02;
    ref_D(2, 0) = -3.42700584373645199321e+02;
    ref_D(2, 1) = -3.42015172813990830036e+02;
    ref_D(2, 2) = 2.07986235474655586586e+03;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
