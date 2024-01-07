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
ModelPart& RansKOmegaSSTKAFC2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KOmegaSSTTestUtilities::RansKOmegaSSTK2D3NSetUp(
        rModel, "RansKOmegaSSTKAFC2D3N");

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

ModelPart& RansKOmegaSSTOmegaAFC2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KOmegaSSTTestUtilities::RansKOmegaSSTOmega2D3NSetUp(
        rModel, "RansKOmegaSSTOmegaAFC2D3N");

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKAFC2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKAFC2D3NSetUp(model);

    // Test:
    FluidTestUtilities::RunEntityEquationIdVectorTest(r_model_part.Elements(), r_model_part.GetProcessInfo(), {&TURBULENT_KINETIC_ENERGY});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKAFC2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKAFC2D3NSetUp(model);

    // Test:
    FluidTestUtilities::RunEntityGetDofListTest(r_model_part.Elements(), r_model_part.GetProcessInfo(), {&TURBULENT_KINETIC_ENERGY});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKAFC2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKAFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  1.8682173192276441e+00;
    ref_RHS[1] =  2.7875608478310872e+00;
    ref_RHS[2] =  2.6311202701767287e+00;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKAFC2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKAFC2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  1.8682173192276441e+00;
    ref_RHS[1] =  2.7875608478310872e+00;
    ref_RHS[2] =  2.6311202701767287e+00;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKAFC2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKAFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS, ref_RHS;
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) =  6.3451838654161087e+01;
    ref_LHS(0, 1) =  2.8521338254976225e+01;
    ref_LHS(0, 2) =  3.8954331113741716e+01;
    ref_LHS(1, 0) =  2.2525143616121976e+01;
    ref_LHS(1, 1) =  4.2905226978547020e+01;
    ref_LHS(1, 2) =  3.4120445312612631e+01;
    ref_LHS(2, 0) =  2.7737579613425574e+01;
    ref_LHS(2, 1) =  2.7014083650550191e+01;
    ref_LHS(2, 2) =  7.2208961159504938e+01;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKAFC2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKAFC2D3NSetUp(model);

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

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKAFC2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKAFC2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) =  6.3451838654161087e+01;
    ref_D(0, 1) =  2.8521338254976225e+01;
    ref_D(0, 2) =  3.8954331113741716e+01;
    ref_D(1, 0) =  2.2525143616121976e+01;
    ref_D(1, 1) =  4.2905226978547020e+01;
    ref_D(1, 2) =  3.4120445312612631e+01;
    ref_D(2, 0) =  2.7737579613425574e+01;
    ref_D(2, 1) =  2.7014083650550191e+01;
    ref_D(2, 2) =  7.2208961159504938e+01;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaAFC2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaAFC2D3NSetUp(model);

    // Test:
    FluidTestUtilities::RunEntityEquationIdVectorTest(r_model_part.Elements(), r_model_part.GetProcessInfo(), {&TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaAFC2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaAFC2D3NSetUp(model);

    // Test:
    FluidTestUtilities::RunEntityGetDofListTest(r_model_part.Elements(), r_model_part.GetProcessInfo(), {&TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaAFC2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaAFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  -1.1487756872879527e+03;
    ref_RHS[1] =  -1.1487721715993939e+03;
    ref_RHS[2] =  -1.1487748491231873e+03;

    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaAFC2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaAFC2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  -1.1487756872879527e+03;
    ref_RHS[1] =  -1.1487721715993939e+03;
    ref_RHS[2] =  -1.1487748491231873e+03;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaAFC2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaAFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS, ref_RHS;
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) =  9.7941220320961847e+01;
    ref_LHS(0, 1) =  3.9349486203577364e+01;
    ref_LHS(0, 2) =  5.4353123658674932e+01;
    ref_LHS(1, 0) =  3.3367447520515491e+01;
    ref_LHS(1, 1) =  5.1391121704115264e+01;
    ref_LHS(1, 2) =  4.4097984859249166e+01;
    ref_LHS(2, 0) =  4.3288539317764361e+01;
    ref_LHS(2, 1) =  3.7182121078519771e+01;
    ref_LHS(2, 2) =  1.0323748808230647e+02;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaAFC2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaAFC2D3NSetUp(model);

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

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaAFC2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaAFC2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) =  9.7941220320961847e+01;
    ref_D(0, 1) =  3.9349486203577364e+01;
    ref_D(0, 2) =  5.4353123658674932e+01;
    ref_D(1, 0) =  3.3367447520515491e+01;
    ref_D(1, 1) =  5.1391121704115264e+01;
    ref_D(1, 2) =  4.4097984859249166e+01;
    ref_D(2, 0) =  4.3288539317764361e+01;
    ref_D(2, 1) =  3.7182121078519771e+01;
    ref_D(2, 2) =  1.0323748808230647e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
