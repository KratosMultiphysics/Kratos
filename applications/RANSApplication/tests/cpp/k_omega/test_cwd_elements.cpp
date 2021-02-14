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
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityEquationIdVectorTest(r_model_part, {&TURBULENT_KINETIC_ENERGY});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaKCWD2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaKCWD2D3NSetUp(model);

    // Test:
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityGetDofListTest(r_model_part, {&TURBULENT_KINETIC_ENERGY});
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
    ref_RHS[0] =  3.7539002230212812e+00;
    ref_RHS[1] =  5.5128164584891115e+00;
    ref_RHS[2] =  5.3066252612017948e+00;
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
    ref_RHS[0] =  3.7539002230212812e+00;
    ref_RHS[1] =  5.5128164584891115e+00;
    ref_RHS[2] =  5.3066252612017948e+00;

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
    ref_LHS(0, 0) =  4.6304773467248617e+02;
    ref_LHS(0, 1) =  -2.6712642334829155e+02;
    ref_LHS(0, 2) =  6.6962955362638638e+01;
    ref_LHS(1, 0) =  -2.6565730044157880e+02;
    ref_LHS(1, 1) =  7.2692462653433211e+02;
    ref_LHS(1, 2) =  -2.6505940838625685e+02;
    ref_LHS(2, 0) =  6.6529738176755487e+01;
    ref_LHS(2, 1) =  -2.6647707122086979e+02;
    ref_LHS(2, 2) =  4.5571445146833622e+02;

    ref_RHS[0] =  -3.6125320879998380e+03;
    ref_RHS[1] =  -1.7516548125893678e+03;
    ref_RHS[2] =  -3.3224029802140183e+04;

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
    ref_M(0, 0) =  2.5073038992536423e-01;
    ref_M(0, 1) =  4.2234902458487322e-02;
    ref_M(0, 2) =  4.1845753571334500e-02;
    ref_M(1, 0) =  4.0181502727893316e-02;
    ref_M(1, 1) =  2.4831138382970699e-01;
    ref_M(1, 2) =  4.0559574541282617e-02;
    ref_M(2, 0) =  4.2412484039581116e-02;
    ref_M(2, 1) =  4.2773816381014559e-02;
    ref_M(2, 2) =  2.5092088828537007e-01;

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
    ref_D(0, 0) =  4.6304773467248617e+02;
    ref_D(0, 1) =  -2.6712642334829155e+02;
    ref_D(0, 2) =  6.6962955362638638e+01;
    ref_D(1, 0) =  -2.6565730044157880e+02;
    ref_D(1, 1) =  7.2692462653433211e+02;
    ref_D(1, 2) =  -2.6505940838625685e+02;
    ref_D(2, 0) =  6.6529738176755487e+01;
    ref_D(2, 1) =  -2.6647707122086979e+02;
    ref_D(2, 2) =  4.5571445146833622e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaCWD2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaCWD2D3NSetUp(model);

    // Test:
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityEquationIdVectorTest(r_model_part, {&TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaOmegaCWD2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaOmegaCWD2D3NSetUp(model);

    // Test:
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityGetDofListTest(r_model_part, {&TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE});
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
    ref_RHS[0] =  1.2207271695330226e+02;
    ref_RHS[1] =  1.2066982567876937e+02;
    ref_RHS[2] =  1.2238810882376549e+02;
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
    ref_RHS[0] =  1.2207271695330226e+02;
    ref_RHS[1] =  1.2066982567876937e+02;
    ref_RHS[2] =  1.2238810882376549e+02;

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
    ref_LHS(0, 0) =  5.8976102030279617e+02;
    ref_LHS(0, 1) =  -2.9705653005623549e+02;
    ref_LHS(0, 2) =  9.9577487952062597e+01;
    ref_LHS(1, 0) =  -2.9558740714952279e+02;
    ref_LHS(1, 1) =  8.8755535648540069e+02;
    ref_LHS(1, 2) =  -2.9624497807878333e+02;
    ref_LHS(2, 0) =  9.9144270766179488e+01;
    ref_LHS(2, 1) =  -2.9766264091339644e+02;
    ref_LHS(2, 2) =  5.7990465798288551e+02;

    ref_RHS[0] =  -2.8739732645200915e+05;
    ref_RHS[1] =  1.6508181176165611e+05;
    ref_RHS[2] =  -2.6068941413931243e+05;

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
    ref_M(0, 0) =  2.5048963472432273e-01;
    ref_M(0, 1) =  4.2045494975105534e-02;
    ref_M(0, 2) =  4.1786463698278500e-02;
    ref_M(1, 0) =  4.0674175483377553e-02;
    ref_M(1, 1) =  2.4887773572381777e-01;
    ref_M(1, 2) =  4.0927571865032361e-02;
    ref_M(2, 0) =  4.2165425334750165e-02;
    ref_M(2, 1) =  4.2403927052567070e-02;
    ref_M(2, 2) =  2.5061600742199719e-01;

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
    ref_D(0, 0) =  5.8976102030279617e+02;
    ref_D(0, 1) =  -2.9705653005623549e+02;
    ref_D(0, 2) =  9.9577487952062597e+01;
    ref_D(1, 0) =  -2.9558740714952279e+02;
    ref_D(1, 1) =  8.8755535648540069e+02;
    ref_D(1, 2) =  -2.9624497807878333e+02;
    ref_D(2, 0) =  9.9144270766179488e+01;
    ref_D(2, 1) =  -2.9766264091339644e+02;
    ref_D(2, 2) =  5.7990465798288551e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
