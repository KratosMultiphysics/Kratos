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
ModelPart& RansKOmegaSSTKRFC2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KOmegaSSTTestUtilities::RansKOmegaSSTK2D3NSetUp(
        rModel, "RansKOmegaSSTKRFC2D3N");

    StabilizationMethodTestUtilities::InitializeResidualBasedFluxCorrectedConstants(
        r_model_part.GetProcessInfo());

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

ModelPart& RansKOmegaSSTOmegaRFC2D3NSetUp(
    Model& rModel)
{
    auto& r_model_part = KOmegaSSTTestUtilities::RansKOmegaSSTOmega2D3NSetUp(
        rModel, "RansKOmegaSSTOmegaRFC2D3N");

    StabilizationMethodTestUtilities::InitializeResidualBasedFluxCorrectedConstants(
        r_model_part.GetProcessInfo());

    RansApplicationTestUtilities::CheckElementsAndConditions(r_model_part);

    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKRFC2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKRFC2D3NSetUp(model);

    // Test:
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityEquationIdVectorTest(r_model_part, {&TURBULENT_KINETIC_ENERGY});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKRFC2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKRFC2D3NSetUp(model);

    // Test:
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityGetDofListTest(r_model_part, {&TURBULENT_KINETIC_ENERGY});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKRFC2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKRFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  3.7537067121785825e+00;
    ref_RHS[1] =  5.5124731718459481e+00;
    ref_RHS[2] =  5.3063887570340711e+00;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKRFC2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKRFC2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  3.7537067121785825e+00;
    ref_RHS[1] =  5.5124731718459481e+00;
    ref_RHS[2] =  5.3063887570340711e+00;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKRFC2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKRFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS(3, 0.0), ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) =  3.6406859649465594e+02;
    ref_LHS(0, 1) =  -4.6945346280565772e+01;
    ref_LHS(0, 2) =  -5.4251676963416585e+01;
    ref_LHS(1, 0) =  -4.5476223373853003e+01;
    ref_LHS(1, 1) =  2.8682490489108380e+02;
    ref_LHS(1, 2) =  -4.5151867891760084e+01;
    ref_LHS(2, 0) =  -5.4684894149299723e+01;
    ref_LHS(2, 1) =  -4.6569530726373131e+01;
    ref_LHS(2, 2) =  3.5701070055836755e+02;

    ref_RHS[0] =  1.5598793172390560e+02;
    ref_RHS[1] =  -7.5718756793591456e+03;
    ref_RHS[2] =  -3.1170541173020167e+04;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKRFC2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKRFC2D3NSetUp(model);

    // Test:
    Matrix M, ref_M(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(
        M, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_M(0, 0) =  2.5072226261875297e-01;
    ref_M(0, 1) =  4.2230116748414516e-02;
    ref_M(0, 2) =  4.1841969415818095e-02;
    ref_M(1, 0) =  4.0176944727397114e-02;
    ref_M(1, 1) =  2.4830038434544086e-01;
    ref_M(1, 2) =  4.0555289868167484e-02;
    ref_M(2, 0) =  4.2408646489386220e-02;
    ref_M(2, 1) =  4.2769282082160680e-02;
    ref_M(2, 2) =  2.5091434447457339e-01;

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKRFC2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKRFC2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) =  3.6406859649465594e+02;
    ref_D(0, 1) =  -4.6945346280565772e+01;
    ref_D(0, 2) =  -5.4251676963416585e+01;
    ref_D(1, 0) =  -4.5476223373853003e+01;
    ref_D(1, 1) =  2.8682490489108380e+02;
    ref_D(1, 2) =  -4.5151867891760084e+01;
    ref_D(2, 0) =  -5.4684894149299723e+01;
    ref_D(2, 1) =  -4.6569530726373131e+01;
    ref_D(2, 2) =  3.5701070055836755e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaRFC2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaRFC2D3NSetUp(model);

    // Test:
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityEquationIdVectorTest(r_model_part, {&TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaRFC2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaRFC2D3NSetUp(model);

    // Test:
    FluidTestUtilities::Testing<ModelPart::ElementsContainerType>::RunEntityGetDofListTest(r_model_part, {&TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaRFC2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaRFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  -2.3067274247976061e+03;
    ref_RHS[1] =  -2.2762635523464137e+03;
    ref_RHS[2] =  -2.3092144219327747e+03;

    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaRFC2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaRFC2D3NSetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] =  -2.3067274247976061e+03;
    ref_RHS[1] =  -2.2762635523464137e+03;
    ref_RHS[2] =  -2.3092144219327747e+03;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaRFC2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaRFC2D3NSetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS(3, 0.0), ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) =  5.3838341225546185e+02;
    ref_LHS(0, 1) =  -6.9125819472129692e+01;
    ref_LHS(0, 2) =  -8.4783400409987962e+01;
    ref_LHS(1, 0) =  -6.7642540152229387e+01;
    ref_LHS(1, 1) =  3.8843111836835084e+02;
    ref_LHS(1, 2) =  -6.5806133441652463e+01;
    ref_LHS(2, 0) =  -8.5064449815243719e+01;
    ref_LHS(2, 1) =  -6.7033296793947585e+01;
    ref_LHS(2, 2) =  5.2100829120014055e+02;

    ref_RHS[0] =  -2.1033171775697754e+05;
    ref_RHS[1] =  1.5801320141793793e+04;
    ref_RHS[2] =  -1.7476201668847815e+05;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaRFC2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaRFC2D3NSetUp(model);

    // Test:
    Matrix M, ref_M(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(
        M, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_M(0, 0) =  2.5055758647313475e-01;
    ref_M(0, 1) =  4.2244947702503644e-02;
    ref_M(0, 2) =  4.1861634480182014e-02;
    ref_M(1, 0) =  4.0657401084665416e-02;
    ref_M(1, 1) =  2.4868153953415034e-01;
    ref_M(1, 2) =  4.0905485096895247e-02;
    ref_M(2, 0) =  4.2101682661747554e-02;
    ref_M(2, 1) =  4.2372836252294606e-02;
    ref_M(2, 2) =  2.5055050421763458e-01;

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaRFC2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaRFC2D3NSetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) =  5.3838341225546185e+02;
    ref_D(0, 1) =  -6.9125819472129692e+01;
    ref_D(0, 2) =  -8.4783400409987962e+01;
    ref_D(1, 0) =  -6.7642540152229387e+01;
    ref_D(1, 1) =  3.8843111836835084e+02;
    ref_D(1, 2) =  -6.5806133441652463e+01;
    ref_D(2, 0) =  -8.5064449815243719e+01;
    ref_D(2, 1) =  -6.7033296793947585e+01;
    ref_D(2, 2) =  5.2100829120014055e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
