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
    RansApplicationTestUtilities::TestEquationIdVector<ModelPart::ElementsContainerType>(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaKCWD2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaKCWD2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestGetDofList<ModelPart::ElementsContainerType>(
        r_model_part, TURBULENT_KINETIC_ENERGY);
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
    ref_RHS[0] = 2.4695761696119565e+00;
    ref_RHS[1] = 1.5697917767132448e+00;
    ref_RHS[2] = 1.6364014119815309e+00;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_EXPECT_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
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
    ref_RHS[0] = 2.4695761696119565e+00;
    ref_RHS[1] = 1.5697917767132448e+00;
    ref_RHS[2] = 1.6364014119815309e+00;

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
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
    ref_LHS(0, 0) = 1.4194218759467108e+03;
    ref_LHS(0, 1) = -9.6750516598254501e+02;
    ref_LHS(0, 2) = 2.5032611209123323e+02;
    ref_LHS(1, 0) = -9.6646757856995748e+02;
    ref_LHS(1, 1) = 2.6162119244428650e+03;
    ref_LHS(1, 2) = -8.7453420090158500e+02;
    ref_LHS(2, 0) = 2.4952784091796218e+02;
    ref_LHS(2, 1) = -8.7521320245558672e+02;
    ref_LHS(2, 2) = 1.8769572111640177e+03;

    ref_RHS[0] = -6.8684372210569287e+04;
    ref_RHS[1] = 9.1585801355519128e+04;
    ref_RHS[2] = -1.4272016049218399e+05;

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_EXPECT_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
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
    ref_M(0, 0) = 2.5029361486252505e-01;
    ref_M(0, 1) = 4.1823659049689606e-02;
    ref_M(0, 2) = 4.1665035592019277e-02;
    ref_M(1, 0) = 4.1127052629825391e-02;
    ref_M(1, 1) = 2.4972648464155722e-01;
    ref_M(1, 2) = 4.1479037174023022e-02;
    ref_M(2, 0) = 4.1911314258101878e-02;
    ref_M(2, 1) = 4.1782743308181054e-02;
    ref_M(2, 2) = 2.5018882917553120e-01;

    KRATOS_EXPECT_MATRIX_NEAR(M, ref_M, 1e-12);
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
    ref_D(0, 0) = 1.4194218759467108e+03;
    ref_D(0, 1) = -9.6750516598254501e+02;
    ref_D(0, 2) = 2.5032611209123323e+02;
    ref_D(1, 0) = -9.6646757856995748e+02;
    ref_D(1, 1) = 2.6162119244428650e+03;
    ref_D(1, 2) = -8.7453420090158500e+02;
    ref_D(2, 0) = 2.4952784091796218e+02;
    ref_D(2, 1) = -8.7521320245558672e+02;
    ref_D(2, 2) = 1.8769572111640177e+03;

    KRATOS_EXPECT_MATRIX_NEAR(D, ref_D, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaCWD2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaCWD2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestEquationIdVector<ModelPart::ElementsContainerType>(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaCWD2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaCWD2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestGetDofList<ModelPart::ElementsContainerType>(
        r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
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
    ref_RHS[0] = 2.3727362021772407e+02;
    ref_RHS[1] = 2.3658299713380416e+02;
    ref_RHS[2] = 2.3732254636130006e+02;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_EXPECT_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
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
    ref_RHS[0] = 2.3727362021772407e+02;
    ref_RHS[1] = 2.3658299713380416e+02;
    ref_RHS[2] = 2.3732254636130006e+02;

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
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
    ref_LHS(0, 0) = 2.5872332277389887e+03;
    ref_LHS(0, 1) = -1.9130604217267542e+03;
    ref_LHS(0, 2) = 3.7083315566389865e+02;
    ref_LHS(1, 0) = -1.9120228343141666e+03;
    ref_LHS(1, 1) = 4.8427540304619652e+03;
    ref_LHS(1, 2) = -1.7770812587735570e+03;
    ref_LHS(2, 0) = 3.7003488449062763e+02;
    ref_LHS(2, 1) = -1.7777602603275584e+03;
    ref_LHS(2, 2) = 3.2626123139209249e+03;

    ref_RHS[0] = -9.9086526170477737e+05;
    ref_RHS[1] = -1.0969821411709106e+06;
    ref_RHS[2] = 3.3969010895097250e+05;

    KRATOS_EXPECT_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_EXPECT_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
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
    ref_M(0, 0) = 2.5019672508174201e-01;
    ref_M(0, 1) = 4.1771908671012055e-02;
    ref_M(0, 2) = 4.1665131342759218e-02;
    ref_M(1, 0) = 4.1305266280965751e-02;
    ref_M(1, 1) = 2.4981657650204664e-01;
    ref_M(1, 2) = 4.1540827508372923e-02;
    ref_M(2, 0) = 4.1830733853424423e-02;
    ref_M(2, 1) = 4.1744647059419077e-02;
    ref_M(2, 2) = 2.5012717991932987e-01;

    KRATOS_EXPECT_MATRIX_NEAR(M, ref_M, 1e-12);
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
    ref_D(0, 0) = 2.5872332277389887e+03;
    ref_D(0, 1) = -1.9130604217267542e+03;
    ref_D(0, 2) = 3.7083315566389865e+02;
    ref_D(1, 0) = -1.9120228343141666e+03;
    ref_D(1, 1) = 4.8427540304619652e+03;
    ref_D(1, 2) = -1.7770812587735570e+03;
    ref_D(2, 0) = 3.7003488449062763e+02;
    ref_D(2, 1) = -1.7777602603275584e+03;
    ref_D(2, 2) = 3.2626123139209249e+03;

    KRATOS_EXPECT_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
