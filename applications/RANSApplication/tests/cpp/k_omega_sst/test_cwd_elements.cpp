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
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityEquationIdVectorTest(r_model_part, {&TURBULENT_KINETIC_ENERGY});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3NSetUp(model);

    // Test:
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityGetDofListTest(r_model_part, {&TURBULENT_KINETIC_ENERGY});
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
    ref_RHS[0] =  3.7539009241285344e+00;
    ref_RHS[1] =  5.5128188214882536e+00;
    ref_RHS[2] =  5.3066261412215487e+00;
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
    ref_RHS[0] =  3.7539009241285344e+00;
    ref_RHS[1] =  5.5128188214882536e+00;
    ref_RHS[2] =  5.3066261412215487e+00;

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
    ref_LHS(0, 0) =  4.6305475568786983e+02;
    ref_LHS(0, 1) =  -2.6713342691836431e+02;
    ref_LHS(0, 2) =  6.6962960549350072e+01;
    ref_LHS(1, 0) =  -2.6566430401165161e+02;
    ref_LHS(1, 1) =  7.2693869418204736e+02;
    ref_LHS(1, 2) =  -2.6506641023053049e+02;
    ref_LHS(2, 0) =  6.6529743363466949e+01;
    ref_LHS(2, 1) =  -2.6648407306514360e+02;
    ref_LHS(2, 2) =  4.5572147710329466e+02;

    ref_RHS[0] =  -3.6123578793624638e+03;
    ref_RHS[1] =  -1.7514720854860825e+03;
    ref_RHS[2] =  -3.3224392747523350e+04;

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
    ref_M(0, 0) =  2.5073039926503526e-01;
    ref_M(0, 1) =  4.2234925013615872e-02;
    ref_M(0, 2) =  4.1845762212368076e-02;
    ref_M(1, 0) =  4.0181524374683836e-02;
    ref_M(1, 1) =  2.4831146555249350e-01;
    ref_M(1, 2) =  4.0559597423850263e-02;
    ref_M(2, 0) =  4.2412492843330041e-02;
    ref_M(2, 1) =  4.2773840689482853e-02;
    ref_M(2, 2) =  2.5092090407515100e-01;

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
    ref_D(0, 0) =  4.6305475568786983e+02;
    ref_D(0, 1) =  -2.6713342691836431e+02;
    ref_D(0, 2) =  6.6962960549350072e+01;
    ref_D(1, 0) =  -2.6566430401165161e+02;
    ref_D(1, 1) =  7.2693869418204736e+02;
    ref_D(1, 2) =  -2.6506641023053049e+02;
    ref_D(2, 0) =  6.6529743363466949e+01;
    ref_D(2, 1) =  -2.6648407306514360e+02;
    ref_D(2, 2) =  4.5572147710329466e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3NSetUp(model);

    // Test:
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityEquationIdVectorTest(r_model_part, {&TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3NSetUp(model);

    // Test:
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityGetDofListTest(r_model_part, {&TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE});
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
    ref_RHS[0] =  -2.3068106857354496e+03;
    ref_RHS[1] =  -2.2764275815564715e+03;
    ref_RHS[2] =  -2.3092960441756359e+03;
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
    ref_RHS[0] =  -2.3068106857354496e+03;
    ref_RHS[1] =  -2.2764275815564715e+03;
    ref_RHS[2] =  -2.3092960441756359e+03;

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
    ref_LHS(0, 0) =  5.8880588430706837e+02;
    ref_LHS(0, 1) =  -3.0208441226714029e+02;
    ref_LHS(0, 2) =  9.7764022847602305e+01;
    ref_LHS(1, 0) =  -3.0060113294723993e+02;
    ref_LHS(1, 1) =  8.5829931210466884e+02;
    ref_LHS(1, 2) =  -3.0270118236808219e+02;
    ref_LHS(2, 0) =  9.7482973442346577e+01;
    ref_LHS(2, 1) =  -3.0392834572037737e+02;
    ref_LHS(2, 2) =  5.7536653923988752e+02;

    ref_RHS[0] =  -2.8559792396780843e+05;
    ref_RHS[1] =  1.7354500941188674e+05;
    ref_RHS[2] =  -2.5725103871187189e+05;

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
    ref_M(0, 0) =  2.5056225334414944e-01;
    ref_M(0, 1) =  4.2249734458007122e-02;
    ref_M(0, 2) =  4.1864260500512757e-02;
    ref_M(1, 0) =  4.0662003298607093e-02;
    ref_M(1, 1) =  2.4869616565612798e-01;
    ref_M(1, 2) =  4.0910054414649392e-02;
    ref_M(2, 0) =  4.2104322993142894e-02;
    ref_M(2, 1) =  4.2377596573536985e-02;
    ref_M(2, 2) =  2.5055494546784829e-01;

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
    ref_D(0, 0) =  5.8880588430706837e+02;
    ref_D(0, 1) =  -3.0208441226714029e+02;
    ref_D(0, 2) =  9.7764022847602305e+01;
    ref_D(1, 0) =  -3.0060113294723993e+02;
    ref_D(1, 1) =  8.5829931210466884e+02;
    ref_D(1, 2) =  -3.0270118236808219e+02;
    ref_D(2, 0) =  9.7482973442346577e+01;
    ref_D(2, 1) =  -3.0392834572037737e+02;
    ref_D(2, 2) =  5.7536653923988752e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
