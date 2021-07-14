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
ModelPart& RansKOmegaKAFC2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KOmegaTestUtilities::RansKOmegaK2D3NSetUp(
        rModel, "RansKOmegaKAFC2D3N");

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

ModelPart& RansKOmegaOmegaAFC2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KOmegaTestUtilities::RansKOmegaOmega2D3NSetUp(
        rModel, "RansKOmegaOmegaAFC2D3N");

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaKAFC2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaKAFC2D3NSetUp(model);

    // Test:
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityEquationIdVectorTest(r_model_part, {&TURBULENT_KINETIC_ENERGY});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaKAFC2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaKAFC2D3NSetUp(model);

    // Test:
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityGetDofListTest(r_model_part, {&TURBULENT_KINETIC_ENERGY});
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
    ref_RHS[0] =  1.8682173192276443e+00;
    ref_RHS[1] =  2.7875608478310872e+00;
    ref_RHS[2] =  2.6311202701767296e+00;
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
    ref_RHS[0] =  1.8682173192276443e+00;
    ref_RHS[1] =  2.7875608478310872e+00;
    ref_RHS[2] =  2.6311202701767296e+00;

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
    ref_LHS(0, 0) =  7.0027238672003577e+01;
    ref_LHS(0, 1) =  2.7295946055994204e+01;
    ref_LHS(0, 2) =  3.3604323294881240e+01;
    ref_LHS(1, 0) =  2.8765068962706966e+01;
    ref_LHS(1, 1) =  4.2391079988133619e+01;
    ref_LHS(1, 2) =  2.8394666956441021e+01;
    ref_LHS(2, 0) =  3.3171106108998103e+01;
    ref_LHS(2, 1) =  2.6977004121827981e+01;
    ref_LHS(2, 2) =  6.6812514192654589e+01;

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
    ref_D(0, 0) =  7.0027238672003577e+01;
    ref_D(0, 1) =  2.7295946055994204e+01;
    ref_D(0, 2) =  3.3604323294881240e+01;
    ref_D(1, 0) =  2.8765068962706966e+01;
    ref_D(1, 1) =  4.2391079988133619e+01;
    ref_D(1, 2) =  2.8394666956441021e+01;
    ref_D(2, 0) =  3.3171106108998103e+01;
    ref_D(2, 1) =  2.6977004121827981e+01;
    ref_D(2, 2) =  6.6812514192654589e+01;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaAFC2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaAFC2D3NSetUp(model);

    // Test:
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityEquationIdVectorTest(r_model_part, {&TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaAFC2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaAFC2D3NSetUp(model);

    // Test:
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityGetDofListTest(r_model_part, {&TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE});
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
    ref_RHS[0] = 6.0855934008867820e+01;
    ref_RHS[1] = 6.0855934008867806e+01;
    ref_RHS[2] = 6.0855934008867806e+01;
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
    ref_RHS[0] = 6.0855934008867820e+01;
    ref_RHS[1] = 6.0855934008867806e+01;
    ref_RHS[2] = 6.0855934008867806e+01;

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
    ref_LHS(0, 0) =  1.0429645399079686e+02;
    ref_LHS(0, 1) =  4.1397396148019652e+01;
    ref_LHS(0, 2) =  4.9931445876593408e+01;
    ref_LHS(1, 0) =  4.2866519054732422e+01;
    ref_LHS(1, 1) =  6.4259517868503622e+01;
    ref_LHS(1, 2) =  4.2181285016768193e+01;
    ref_LHS(2, 0) =  4.9498228690710278e+01;
    ref_LHS(2, 1) =  4.0763622182155153e+01;
    ref_LHS(2, 2) =  9.9507569352956480e+01;

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
    ref_D(0, 0) =  1.0429645399079686e+02;
    ref_D(0, 1) =  4.1397396148019652e+01;
    ref_D(0, 2) =  4.9931445876593408e+01;
    ref_D(1, 0) =  4.2866519054732422e+01;
    ref_D(1, 1) =  6.4259517868503622e+01;
    ref_D(1, 2) =  4.2181285016768193e+01;
    ref_D(2, 0) =  4.9498228690710278e+01;
    ref_D(2, 1) =  4.0763622182155153e+01;
    ref_D(2, 2) =  9.9507569352956480e+01;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
