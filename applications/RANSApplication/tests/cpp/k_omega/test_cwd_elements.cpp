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

    return r_model_part;
}

ModelPart& RansKOmegaOmegaCWD2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KOmegaTestUtilities::RansKOmegaOmega2D3NSetUp(
        rModel, "RansKOmegaOmegaCWD2D3N");

    StabilizationMethodTestUtilities::InitializeCrossWindStabilizationConstants(
        r_model_part.GetProcessInfo());

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
    ref_RHS[0] = 2.4695761705848365e+00;
    ref_RHS[1] = 1.5697917770959631e+00;
    ref_RHS[2] = 1.6364014122654156e+00;
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
    ref_RHS[0] = 2.4695761705848365e+00;
    ref_RHS[1] = 1.5697917770959631e+00;
    ref_RHS[2] = 1.6364014122654156e+00;

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
    ref_LHS(0, 0) = 1.6511120905670732e+03;
    ref_LHS(0, 1) = -1.1991217462257839e+03;
    ref_LHS(0, 2) = 2.5025247788852712e+02;
    ref_LHS(1, 0) = -1.1980841588131966e+03;
    ref_LHS(1, 1) = 3.0794757605509021e+03;
    ref_LHS(1, 2) = -1.1061814566478565e+03;
    ref_LHS(2, 0) = 2.4945420671525602e+02;
    ref_LHS(2, 1) = -1.1068604582018577e+03;
    ref_LHS(2, 2) = 2.1086781011832550e+03;

    ref_RHS[0] = -7.6693476963781155e+04;
    ref_RHS[1] = 1.1565496460867976e+05;
    ref_RHS[2] = -1.5878021900419163e+05;

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
    ref_M(0, 0) = 2.5029361493243962e-01;
    ref_M(0, 1) = 4.1823659073663109e-02;
    ref_M(0, 2) = 4.1665035611182719e-02;
    ref_M(1, 0) = 4.1127052653224049e-02;
    ref_M(1, 1) = 2.4972648467294564e-01;
    ref_M(1, 2) = 4.1479037185367829e-02;
    ref_M(2, 0) = 4.1911314277381484e-02;
    ref_M(2, 1) = 4.1782743319712837e-02;
    ref_M(2, 2) = 2.5018882918324814e-01;

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
    ref_D(0, 0) = 1.6511120905670732e+03;
    ref_D(0, 1) = -1.1991217462257839e+03;
    ref_D(0, 2) = 2.5025247788852712e+02;
    ref_D(1, 0) = -1.1980841588131966e+03;
    ref_D(1, 1) = 3.0794757605509021e+03;
    ref_D(1, 2) = -1.1061814566478565e+03;
    ref_D(2, 0) = 2.4945420671525602e+02;
    ref_D(2, 1) = -1.1068604582018577e+03;
    ref_D(2, 2) = 2.1086781011832550e+03;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
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
    ref_RHS[0] = 2.3727362034143948e+02;
    ref_RHS[1] = 2.3658299719726426e+02;
    ref_RHS[2] = 2.3732254640107863e+02;
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
    ref_RHS[0] = 2.3727362034143948e+02;
    ref_RHS[1] = 2.3658299719726426e+02;
    ref_RHS[2] = 2.3732254640107863e+02;

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
    ref_LHS(0, 0) = 2.9223258912741758e+03;
    ref_LHS(0, 1) = -2.2480998381747859e+03;
    ref_LHS(0, 2) = 3.7077990897416737e+02;
    ref_LHS(1, 0) = -2.2470622507621983e+03;
    ref_LHS(1, 1) = 5.5128389912631137e+03;
    ref_LHS(1, 2) = -2.1121268028925451e+03;
    ref_LHS(2, 0) = 3.6998163780089618e+02;
    ref_LHS(2, 1) = -2.1128058044465465e+03;
    ref_LHS(2, 2) = 3.5977111048774441e+03;

    ref_RHS[0] = -1.0617797415564791e+06;
    ref_RHS[1] = -1.1811532836006940e+06;
    ref_RHS[2] = 4.9477573075801518e+05;

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
    ref_M(0, 0) = 2.5019672519114761e-01;
    ref_M(0, 1) = 4.1771908706167697e-02;
    ref_M(0, 2) = 4.1665131372155939e-02;
    ref_M(1, 0) = 4.1305266315527639e-02;
    ref_M(1, 1) = 2.4981657654144007e-01;
    ref_M(1, 2) = 4.1540827523649765e-02;
    ref_M(2, 0) = 4.1830733882944858e-02;
    ref_M(2, 1) = 4.1744647074886808e-02;
    ref_M(2, 2) = 2.5012717993027489e-01;

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
    ref_D(0, 0) = 2.9223258912741758e+03;
    ref_D(0, 1) = -2.2480998381747859e+03;
    ref_D(0, 2) = 3.7077990897416737e+02;
    ref_D(1, 0) = -2.2470622507621983e+03;
    ref_D(1, 1) = 5.5128389912631137e+03;
    ref_D(1, 2) = -2.1121268028925451e+03;
    ref_D(2, 0) = 3.6998163780089618e+02;
    ref_D(2, 1) = -2.1128058044465465e+03;
    ref_D(2, 2) = 3.5977111048774441e+03;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
