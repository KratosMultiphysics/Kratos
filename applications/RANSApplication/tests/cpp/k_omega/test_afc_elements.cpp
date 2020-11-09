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
ModelPart& RansKOmegaKAFC2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KOmegaTestUtilities::RansKOmegaK2D3NSetUp(
        rModel, "RansKOmegaKAFC2D3N");

    return r_model_part;
}

ModelPart& RansKOmegaOmegaAFC2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KOmegaTestUtilities::RansKOmegaOmega2D3NSetUp(
        rModel, "RansKOmegaOmegaAFC2D3N");

    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaKAFC2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaKAFC2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestEquationIdVector<ModelPart::ElementsContainerType>(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaKAFC2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaKAFC2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestGetDofList<ModelPart::ElementsContainerType>(
        r_model_part, TURBULENT_KINETIC_ENERGY);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaKAFC2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaKAFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = 2.67829732581686741355e+00;
    ref_RHS[1] = 2.21576728463429173388e+00;
    ref_RHS[2] = 1.85563938248039850265e+00;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaKAFC2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaKAFC2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = 2.67829732581686741355e+00;
    ref_RHS[1] = 2.21576728463429173388e+00;
    ref_RHS[2] = 1.85563938248039850265e+00;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaKAFC2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaKAFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS, ref_RHS;
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) = 1.22727546832704533131e+02;
    ref_LHS(0, 1) = 7.87491041488959808703e+01;
    ref_LHS(0, 2) = 8.08415536207580629480e+01;
    ref_LHS(1, 0) = 7.99118115689366987908e+01;
    ref_LHS(1, 1) = 1.89309264135145298269e+02;
    ref_LHS(1, 2) = 9.45686406772083216765e+01;
    ref_LHS(2, 0) = 8.11121632251401223357e+01;
    ref_LHS(2, 1) = 9.44333842004898542655e+01;
    ref_LHS(2, 2) = 1.98034136968214227181e+02;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaKAFC2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaKAFC2D3NSetUp(model);

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

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaKAFC2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaKAFC2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) = 1.22727546832704533131e+02;
    ref_D(0, 1) = 7.87491041488959808703e+01;
    ref_D(0, 2) = 8.08415536207580629480e+01;
    ref_D(1, 0) = 7.99118115689366987908e+01;
    ref_D(1, 1) = 1.89309264135145298269e+02;
    ref_D(1, 2) = 9.45686406772083216765e+01;
    ref_D(2, 0) = 8.11121632251401223357e+01;
    ref_D(2, 1) = 9.44333842004898542655e+01;
    ref_D(2, 2) = 1.98034136968214227181e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaAFC2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaAFC2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestEquationIdVector<ModelPart::ElementsContainerType>(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaAFC2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaAFC2D3NSetUp(model);

    // Test:
    RansApplicationTestUtilities::TestGetDofList<ModelPart::ElementsContainerType>(
        r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaAFC2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaAFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = 1.67913396150496367909e+02;
    ref_RHS[1] = 1.67913396150496339487e+02;
    ref_RHS[2] = 1.67913396150496339487e+02;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaAFC2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaAFC2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = 1.67913396150496367909e+02;
    ref_RHS[1] = 1.67913396150496339487e+02;
    ref_RHS[2] = 1.67913396150496339487e+02;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaAFC2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaAFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS, ref_RHS;
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) = 1.81440253624390123832e+02;
    ref_LHS(0, 1) = 1.16830954029190266397e+02;
    ref_LHS(0, 2) = 1.19717268052484996588e+02;
    ref_LHS(1, 0) = 1.17993661449230984317e+02;
    ref_LHS(1, 1) = 2.80368818881605193383e+02;
    ref_LHS(1, 2) = 1.39893465144128413158e+02;
    ref_LHS(2, 0) = 1.19987877656867055975e+02;
    ref_LHS(2, 1) = 1.39758208667409974169e+02;
    ref_LHS(2, 2) = 2.92961716693029018188e+02;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaAFC2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaAFC2D3NSetUp(model);

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

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaAFC2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaAFC2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) = 1.81440253624390123832e+02;
    ref_D(0, 1) = 1.16830954029190266397e+02;
    ref_D(0, 2) = 1.19717268052484996588e+02;
    ref_D(1, 0) = 1.17993661449230984317e+02;
    ref_D(1, 1) = 2.80368818881605193383e+02;
    ref_D(1, 2) = 1.39893465144128413158e+02;
    ref_D(2, 0) = 1.19987877656867055975e+02;
    ref_D(2, 1) = 1.39758208667409974169e+02;
    ref_D(2, 2) = 2.92961716693029018188e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
